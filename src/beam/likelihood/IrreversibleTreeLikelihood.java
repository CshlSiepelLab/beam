package beam.likelihood;

import java.util.Arrays;
import java.util.List;
import java.util.HashSet;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;

import beam.substitutionmodel.CrisprSubstitutionModel;


/*
 * Calculates the tree likelihood while considering the origin of the cell division process
 * and using a simplified pruning algorithm to save on computations given the irreversible
 * assumptions of the substitution model.
 *
 * @author Stephen Staklinski
 */
@Description("Calculates the tree likelihood while considering the origin of the cell division process " +
            "and using a simplified pruning algorithm to save on computations given the irreversible " +
            "assumptions of the substitution model.")
public class IrreversibleTreeLikelihood extends GenericTreeLikelihood {

    /* Input for the origin time of the cell division process */
    public Input<RealParameter> originInput = new Input<>("origin",
            "Start of the cell division process, usually start of the experiment.",
            Input.Validate.REQUIRED);

    /* Alignment data */
    protected Alignment data;
    /* Height of the origin node */
    protected Double originHeight;
    /* Substitution model for mutations */
    protected CrisprSubstitutionModel substitutionModel;
    /* Number of sites in the alignment */
    protected int nrOfSites;
    /* Transition probability matrix */
    protected double[][] probabilities;
    /* Whether the substitution model is global */
    protected boolean isGlobalModel;
    /*
     * Flag indicating the state of the tree:
     * CLEAN=0: nothing needs to be recalculated for the node</li>
     * DIRTY=1: node partial needs to be recalculated</li>
     * FILTHY=2: indices for the node need to be recalculated</li>
     */
    protected int hasDirt;
    /* Number of nodes in the tree */
    protected int nrOfNodes;
    /* Number of nodes excluding the origin node */
    protected int numNodesNoOrigin;
    /* Size of the partials array per site */
    protected int[] partialsSizes;
    /* Size of the transition matrix per site */
    protected int[] matrixSizes;
    /* Partial likelihoods for each node */
    protected double[][][][] partials;
    /* Transition matrices for each node */
    protected double[][][][] matrices;
    /* Current matrix indices for each node */
    protected int[] currentMatrixIndex;
    /* Stored matrix indices for each node */
    protected int[] storedMatrixIndex;
    /* Current partials indices for each node */
    protected int[] currentPartialsIndex;
    /* Stored partials indices for each node */
    protected int[] storedPartialsIndex;
    /* Ancestral states for each node and site */
    protected int[][] ancestralStates;
    /* Stored ancestral states for each node and site */
    protected int[][] storedAncestralStates;
    /* All possible states per site*/
    private int[][] allStatesPerSite;
    /* State representing unedited sequence */
    private static final int[] UNEDITED_STATE = new int[]{0};
    /* Whether to use scaling for numerical stability */
    protected boolean useScaling = false;
    /* Log scaling factors sum for each site */
    protected double[][] logScalingFactorsSum;
    // Threshold for scaling partials
    private static final double SCALING_THRESHOLD = 1.0E-100;
    // State representing missing data
    private int[] missingDataStatePerSite;
    // Number of state per site
    private int[] nrOfStatesPerSite;

    public IrreversibleTreeLikelihood() {
        branchRateModelInput.setRule(Validate.REQUIRED);
        siteModelInput.setType(SiteModel.Base.class);
    }

    @Override
    public void initAndValidate() {
        substitutionModel = (CrisprSubstitutionModel) ((SiteModel.Base) siteModelInput.get()).substModelInput.get();
        originHeight = originInput.get().getValue();
        nrOfSites = substitutionModel.getSiteCount();
        data = substitutionModel.getData(); // We use the data from the substitution model since the missing data state (-1) is replaced there during the rate matrix setup
        nrOfStatesPerSite = substitutionModel.getNrOfStatesPerSite();    // Reference, not copy
        isGlobalModel = substitutionModel.isGlobalModel();

        // Initialize transition probability matrices per site
        probabilities = new double[nrOfSites][];
        for (int i = 0; i < nrOfSites; i++) {
            probabilities[i] = new double[nrOfStatesPerSite[i] * nrOfStatesPerSite[i]];
        }

        // Initialize data structures for likelihood calculations
        numNodesNoOrigin = treeInput.get().getNodeCount();
        nrOfNodes = numNodesNoOrigin + 1; // +1 for origin node
        missingDataStatePerSite = substitutionModel.getMissingStates();

        matrixSizes = new int[nrOfSites];
        partialsSizes = new int[nrOfSites];
        for (int i = 0; i < nrOfSites; i++) {
            matrixSizes[i] = nrOfStatesPerSite[i] * nrOfStatesPerSite[i];
            partialsSizes[i] = nrOfStatesPerSite[i];
        }
        partials = new double[2][nrOfNodes][nrOfSites][];
        matrices = new double[2][nrOfNodes][nrOfSites][];
        logScalingFactorsSum = new double[2][nrOfSites];
        // Initialize all arrays to 0
        for (int i = 0; i < nrOfNodes; i++) {
            for (int j = 0; j < nrOfSites; j++) {
                partials[0][i][j] = new double[partialsSizes[j]];
                partials[1][i][j] = new double[partialsSizes[j]];
                matrices[0][i][j] = new double[matrixSizes[j]];
                matrices[1][i][j] = new double[matrixSizes[j]];
                Arrays.fill(partials[0][i][j], 0.0);
                Arrays.fill(partials[1][i][j], 0.0);
                Arrays.fill(matrices[0][i][j], 0.0);
                Arrays.fill(matrices[1][i][j], 0.0);
                logScalingFactorsSum[0][j] = 0.0;
                logScalingFactorsSum[1][j] = 0.0;
            }
        }

        // Initialize arrays
        allStatesPerSite = new int[nrOfSites][];
        for (int i = 0; i < nrOfSites; i++) {
            allStatesPerSite[i] = new int[nrOfStatesPerSite[i]];
            for (int j = 0; j < nrOfStatesPerSite[i]; j++) {
                allStatesPerSite[i][j] = j;
            }
        }

        currentMatrixIndex = new int[nrOfNodes];
        storedMatrixIndex = new int[nrOfNodes];
        currentPartialsIndex = new int[nrOfNodes];
        storedPartialsIndex = new int[nrOfNodes];

        ancestralStates = new int[numNodesNoOrigin][nrOfSites];
        storedAncestralStates = new int[numNodesNoOrigin][nrOfSites];
        for (int i = 0; i < numNodesNoOrigin; i++) {
            for (int j = 0; j < nrOfSites; j++) {
                ancestralStates[i][j] = -1;
                storedAncestralStates[i][j] = -1;
            }
        }

        // Set initial states
        setStates(treeInput.get().getRoot());
    }

    // Sets the observed data in likelihood core partials for leaf nodes across all alignment sites
    protected void setStates(Node node) {
        if (node.isLeaf()) {
            int taxonIndex = data.getTaxonIndex(node.getID());
            if (taxonIndex == -1) throw new RuntimeException("Could not find sequence " + node.getID() + " in the alignment");
            List<Integer> states = data.getCounts().get(taxonIndex);
            int leafIndex = node.getNr();
            for (int i = 0; i < nrOfSites; i++) {
                // Set the known state to 1.0 for tips
                partials[currentPartialsIndex[leafIndex]][leafIndex][i][states.get(i)] = 1.0;
                // Set the possible ancestral state for tips
                ancestralStates[leafIndex][i] = (states.get(i) == missingDataStatePerSite[i]) ? -1 : states.get(i);
            }
        } else {
            setStates(node.getLeft());
            setStates(node.getRight());
        }
    }

    @Override
    public double calculateLogP() {
        Node root = treeInput.get().getRoot();
        // Return -infinity if tree exceeds origin time
        if (root.getHeight() >= originHeight) {
            return Double.NEGATIVE_INFINITY;
        }
        // Calculate likelihood with optional scaling
        traverse(root);
        if (logP == Double.NEGATIVE_INFINITY) {
            setUseScaling(true);
            traverse(root);
            if (logP == Double.NEGATIVE_INFINITY) {
                throw new RuntimeException("Likelihood is still negative infinity after scaling");
            }
        }
        return logP;
    }

    /*
     * Computes the transition probabilities pre-order, then does a post-order traversal
     * calculating the partial likelihoods by a modified pruning algorithm
     * taking advantage of the irreversibility of the substitution model
     * to reduce the number of ancestral states to propagate partial
     * likelihoods for.
     */
    protected int traverse(Node node) {
        final int nodeIndex = node.getNr();
        int update = (node.isDirty() | hasDirt);

        // Update transition probabilities if needed
        if (update != Tree.IS_CLEAN) {
            double parentHeight = node.isRoot() ? originHeight : node.getParent().getHeight();
            substitutionModel.getTransitionProbabilitiesAllSites(node, parentHeight, node.getHeight(), 
                    branchRateModelInput.get().getRateForBranch(node), probabilities);
            // Set node transition probability matrix for all sites
            currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
            for (int siteNum = 0; siteNum < nrOfSites; siteNum++) {
                System.arraycopy(probabilities[siteNum], 0, matrices[currentMatrixIndex[nodeIndex]][nodeIndex][siteNum],
                    0, matrixSizes[siteNum]);
            }
        }

        // Process internal nodes
        if (!node.isLeaf()) {
            Node child1 = node.getLeft();
            Node child2 = node.getRight();
            int update1 = traverse(child1);
            int update2 = traverse(child2);

            // Calculate partials if either child is dirty
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {
                int childIndex1 = child1.getNr();
                int childIndex2 = child2.getNr();
                
                setPossibleAncestralStates(childIndex1, childIndex2, nodeIndex);
                calculatePartials(childIndex1, childIndex2, nodeIndex);

                // Calculate origin partials at root
                if (node.isRoot()) {
                    logP = calculateLogLikelihoods(nodeIndex, nodeIndex + 1);
                }

                update |= (update1 | update2);
            }
        }

        return update;
    }

    public void setPossibleAncestralStates(int childIndex1, int childIndex2, int parentIndex) {
        for (int i = 0; i < nrOfSites; i++) {
            final int child1State = ancestralStates[childIndex1][i];
            final int child2State = ancestralStates[childIndex2][i];

            if (child1State == 0 || child2State == 0 || 
                (child1State > 0 && child2State > 0 && child1State != child2State)) {
                ancestralStates[parentIndex][i] = 0;
                continue;
            } else if (child1State > 0 || child2State > 0) {
                ancestralStates[parentIndex][i] = child1State >= 0 ? child1State : child2State;
            }
            // The -1 initialization is the fallback here, which indicates consideration of all possible states for that site
        }
    }

    /*
     * Convert the possible ancestral state integer to an array of states to calculate partials for in the pruning algorithm
     */
    private int[] getPossibleStates(int state, int siteNum) {
        return state == 0 ? UNEDITED_STATE :
               state > 0 ? new int[]{0, state} :
               allStatesPerSite[siteNum];
    }

    /**
     * Calculates partial likelihoods at a node while first checking which partials
     * need to be calculated in the irreversible model based on the possible ancestral 
     * states given the states of the children.
     */
    public void calculatePartials(int childIndex1, int childIndex2, int parentIndex) {
        currentPartialsIndex[parentIndex] = 1 - currentPartialsIndex[parentIndex];
        
        final double[][] partials1 = partials[currentPartialsIndex[childIndex1]][childIndex1];
        final double[][] matrices1 = matrices[currentMatrixIndex[childIndex1]][childIndex1];
        final double[][] partials2 = partials[currentPartialsIndex[childIndex2]][childIndex2];
        final double[][] matrices2 = matrices[currentMatrixIndex[childIndex2]][childIndex2];
        final double[][] partials3 = partials[currentPartialsIndex[parentIndex]][parentIndex];

        for (int k = 0; k < nrOfSites; k++) {
            final int[] possibleStates = getPossibleStates(ancestralStates[parentIndex][k], k);
            final int[] child1States = getPossibleStates(ancestralStates[childIndex1][k], k);
            final int[] child2States = getPossibleStates(ancestralStates[childIndex2][k], k);

            // Calculate partials for all states
            for (int i : possibleStates) {
                double sum1 = 0.0;
                double sum2 = 0.0;
                // Row is the starting parent state, columns are the target child states for the matrix
                // Partials come from the children at the target child states
                final int rowOffset = i * nrOfStatesPerSite[k];
                for (int j : child1States) {
                    sum1 += matrices1[k][rowOffset + j] * partials1[k][j];
                }
                for (int j : child2States) {
                    sum2 += matrices2[k][rowOffset + j] * partials2[k][j];
                }
                partials3[k][i] = sum1 * sum2;
            }

            if (useScaling) {
                scalePartials(parentIndex, k, possibleStates);
            }
        }
    }

    /*
     * Calculates site partial likelihoods and total log likelihood at the origin node.
     * The origin is known to be in the unedited state, so we can only calculate the partials
     * for transitions from that state only.
     */
    public double calculateLogLikelihoods(int rootIndex, int originIndex) {
        double logP = 0.0;
        currentPartialsIndex[originIndex] = 1 - currentPartialsIndex[originIndex];

        final double[][] partials1 = partials[currentPartialsIndex[rootIndex]][rootIndex];
        final double[][] originBranchTransitionProbs = matrices[currentMatrixIndex[rootIndex]][rootIndex];
        final double[][] originPartials = partials[currentPartialsIndex[originIndex]][originIndex];

        for (int k = 0; k < nrOfSites; k++) {
            final int[] possibleStates = getPossibleStates(ancestralStates[rootIndex][k], k);
            double sum1 = 0.0;

            for (int j : possibleStates) {
                sum1 += originBranchTransitionProbs[k][j] * partials1[k][j];
            }
            originPartials[k][0] = sum1;

            if (useScaling) {
                scalePartials(originIndex, k, new int[]{0});
                logP += Math.log(originPartials[k][0]) + logScalingFactorsSum[currentPartialsIndex[originIndex]][k];
            } else {
                logP += Math.log(originPartials[k][0]);
            }
        }

        return logP;
    }

    protected void scalePartials(int nodeIndex, int siteNum, int[] possibleStates) {
        // Find the maximum partial value for scaling
        double scaleFactor = 0.0;
        for (int j : possibleStates) {
            scaleFactor = Math.max(scaleFactor, partials[currentPartialsIndex[nodeIndex]][nodeIndex][siteNum][j]);
        }

        // Scale partials if needed and set scaling factor
        if (scaleFactor < SCALING_THRESHOLD) {
            for (int j : possibleStates) {
                partials[currentPartialsIndex[nodeIndex]][nodeIndex][siteNum][j] /= scaleFactor;
            }
        }
        logScalingFactorsSum[currentPartialsIndex[nodeIndex]][siteNum] += scaleFactor < SCALING_THRESHOLD ? Math.log(scaleFactor) : 0.0;
    }

    public void setUseScaling(boolean status) {
        useScaling = status;
    }



    @Override
    public void store() {
        super.store();

        // Store likelihood core components
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, nrOfNodes);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nrOfNodes);

        for (int i = 0; i < numNodesNoOrigin; i++) {
            for (int j = 0; j < nrOfSites; j++) {
                storedAncestralStates[i][j] = ancestralStates[i][j];
            }
        }
    }

    @Override
    public void restore() {
        super.restore();
        
        // Restore likelihood core components
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);

        if (storedAncestralStates == null) {
            storedAncestralStates = new int[numNodesNoOrigin][nrOfSites];
            for (int i = 0; i < numNodesNoOrigin; i++) {
                for (int j = 0; j < nrOfSites; j++) {
                    storedAncestralStates[i][j] = ancestralStates[i][j];
                }
            }
        }

        for (int i = 0; i < numNodesNoOrigin; i++) {
            for (int j = 0; j < nrOfSites; j++) {
                ancestralStates[i][j] = storedAncestralStates[i][j];
            }
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        // Check site model first (most common case)
        if (((SiteModel.Base) siteModelInput.get()).isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }

        // Check branch rate model
        if (branchRateModelInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }

        // Check tree last (least common case)
        hasDirt = Tree.IS_CLEAN;
        return treeInput.get().somethingIsDirty();
    }
}