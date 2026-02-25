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
    /* Whether the substitution model is global */
    protected boolean isGlobalModel;
    // Recalculation flags
    private static final int IS_CLEAN = 0;
    private static final int IS_DIRTY = 1;
    protected int hasDirt;
    /* Number of nodes in the tree */
    protected int nrOfNodes;
    /* Number of nodes excluding the origin node */
    protected int numNodesNoOrigin;
    /* Size of the partials array per site */
    protected int[] partialsSizes;
    /* Partial likelihoods for each node */
    protected double[][][][] partials;
    /* Current matrix indices for each node */
    protected int[] currentMatrixIndex;
    /* Stored matrix indices for each node */
    protected int[] storedMatrixIndex;
    /* Current partials indices for each node */
    protected int[] currentPartialsIndex;
    /* Stored partials indices for each node */
    protected int[] storedPartialsIndex;
    /* Ancestral states for each node and site */
    protected int[][][] ancestralStates;
    /* All possible states per site*/
    private int[][] allStatesPerSite;
    /* State representing unedited sequence */
    private static final int[] UNEDITED_STATE = new int[]{0};
    /* Whether to use scaling for numerical stability */
    protected boolean useScaling = false;
    protected boolean storedUseScaling = false;
    /* Log scaling factors sum for each site */
    protected double[][][] logScalingFactors;
    // Track how many scaling attempt have been done
    protected int numScalingAttempts = 0;
    // After this many attempts, try to turn off scaling
    protected static final int MAX_SCALING_ATTEMPTS = 1000;
    // Threshold for scaling partials
    private static final double SCALING_THRESHOLD = 1.0E-100;
    // State representing missing data
    private int[] missingDataStatePerSite;
    // Number of state per site
    private int[] nrOfStatesPerSite;
    // Array to hold the core transition probability values that are consistent across sites for the crispr models, to avoid calculating the full sparse matrix upfront
    double[] coreTransitionValues = new double[4]; // 0: unedited to unedited (top left), 1: edited to same edited (diag), 2: any edit to missing data (last col), 3: unedited to any edited without editRate considered yet (first row)
    double[][][] nodeCoreTransitionValues;

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
        nrOfStatesPerSite = substitutionModel.getNrOfStatesPerSite();
        missingDataStatePerSite = substitutionModel.getMissingStates();

        // Initialize data structures for likelihood calculations
        numNodesNoOrigin = treeInput.get().getNodeCount();
        nrOfNodes = numNodesNoOrigin + 1; // +1 for origin node
        nodeCoreTransitionValues = new double[2][nrOfNodes][4]; // First index is for store/restore swapping
        partialsSizes = new int[nrOfSites];
        for (int i = 0; i < nrOfSites; i++) {
            partialsSizes[i] = nrOfStatesPerSite[i];
        }
        partials = new double[2][numNodesNoOrigin][nrOfSites][];
        // Initialize all arrays to 0
        for (int i = 0; i < numNodesNoOrigin; i++) {
            for (int j = 0; j < nrOfSites; j++) {
                partials[0][i][j] = new double[partialsSizes[j]];
                partials[1][i][j] = new double[partialsSizes[j]];
            }
        }

        logScalingFactors = new double[2][numNodesNoOrigin][nrOfSites]; // First index is for store/restore swapping
        numScalingAttempts = 0;

        // Initialize arrays
        allStatesPerSite = new int[nrOfSites][];
        for (int i = 0; i < nrOfSites; i++) {
            allStatesPerSite[i] = new int[nrOfStatesPerSite[i]];
            for (int j = 0; j < nrOfStatesPerSite[i]; j++) {
                allStatesPerSite[i][j] = j;
            }
        }

        ancestralStates = new int[2][numNodesNoOrigin][nrOfSites];

        // Initialize indices to 0
        currentMatrixIndex = new int[nrOfNodes];
        storedMatrixIndex = new int[nrOfNodes];
        currentPartialsIndex = new int[numNodesNoOrigin];
        storedPartialsIndex = new int[numNodesNoOrigin];

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
                // Set the observed states for tips that will never change, hence no need to initialize both partials index
                partials[currentPartialsIndex[leafIndex]][leafIndex][i][states.get(i)] = 1.0;
                ancestralStates[currentPartialsIndex[leafIndex]][leafIndex][i] = (states.get(i) == missingDataStatePerSite[i]) ? -1 : states.get(i);
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
        // If scaling is already on, check if we have tried scaling too many times without success and if so, turn it off and try again without scaling
        if (useScaling) {
            if (numScalingAttempts >= MAX_SCALING_ATTEMPTS) {
                setUseScaling(false);
                hasDirt = IS_DIRTY; // Set dirty to recalculate all partials without scaling
            } else {
                numScalingAttempts++;
            }
        }
        // Calculate likelihood with optional scaling
        traverse(root);
        if (logP == Double.NEGATIVE_INFINITY || Double.isNaN(logP)) {
            setUseScaling(true);
            hasDirt = IS_DIRTY; // Set dirty to recalculate all partials with scaling
            numScalingAttempts = 1;
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
        int updateTransitionProbs = (node.isDirty() | hasDirt);

        // Update transition probabilities if needed
        if (updateTransitionProbs != IS_CLEAN) {
            double parentHeight = node.isRoot() ? originHeight : node.getParent().getHeight();
            coreTransitionValues = substitutionModel.getCoreTransitionProbabilityValues(node, parentHeight, node.getHeight(), branchRateModelInput.get().getRateForBranch(node));
            currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
            nodeCoreTransitionValues[currentMatrixIndex[nodeIndex]][nodeIndex] = Arrays.copyOf(coreTransitionValues, coreTransitionValues.length);
        }

        // Process internal nodes
        if (!node.isLeaf()) {
            Node child1 = node.getLeft();
            Node child2 = node.getRight();
            int updateTransitionProbs1 = traverse(child1);
            int updateTransitionProbs2 = traverse(child2);

            // Calculate partials if either child is dirty or if scaling is on since we need to accumulate the log scaling factors
            if (updateTransitionProbs1 != IS_CLEAN || updateTransitionProbs2 != IS_CLEAN || useScaling) {
                int childIndex1 = child1.getNr();
                int childIndex2 = child2.getNr();
                setPossibleAncestralStates(childIndex1, childIndex2, nodeIndex);
                calculatePartials(childIndex1, childIndex2, nodeIndex);

                updateTransitionProbs |= (updateTransitionProbs1 | updateTransitionProbs2);
            }

            // Always calculate origin partials once at the root
            if (node.isRoot()) {
                logP = calculateLogLikelihoods(nodeIndex);
            }
        }
        return updateTransitionProbs;
    }

    public void setPossibleAncestralStates(int childIndex1, int childIndex2, int parentIndex) {
        currentPartialsIndex[parentIndex] = 1 - currentPartialsIndex[parentIndex];

        for (int i = 0; i < nrOfSites; i++) {
            final int child1State = ancestralStates[currentPartialsIndex[childIndex1]][childIndex1][i];
            final int child2State = ancestralStates[currentPartialsIndex[childIndex2]][childIndex2][i];

            if ((child1State == 0 || child2State == 0) || (child1State > 0 && child2State > 0 && child1State != child2State)) {
                ancestralStates[currentPartialsIndex[parentIndex]][parentIndex][i] = 0; // Must be unedited if either child is unedited OR two different edits, so parent must be unedited
            } else if (child1State > 0 && child2State > 0 && child1State == child2State) {
                ancestralStates[currentPartialsIndex[parentIndex]][parentIndex][i] = child1State; // Must be that edit if both children are edited with the same edit
            } else if (child1State > 0 && child2State <= 0) {
                ancestralStates[currentPartialsIndex[parentIndex]][parentIndex][i] = child1State;   // Must be 0 or that edit if one child is edited and the other is unedited or all possible states
            } else if (child2State > 0 && child1State <= 0) {
                ancestralStates[currentPartialsIndex[parentIndex]][parentIndex][i] = child2State;   // Must be 0 or that edit if one child is edited and the other is unedited or all possible states
            } else {
                ancestralStates[currentPartialsIndex[parentIndex]][parentIndex][i] = -1; // All states possible if both children are all states possible
            }
        }
    }

    /*
     * Convert the possible ancestral state integer to an array of states to calculate partials for in the pruning algorithm
     */
    private int[] getPossibleStates(int state, int[] allStates) {
        return state == 0 ? UNEDITED_STATE :
               state > 0 ? new int[]{0, state} :
               allStates;
    }

    /**
     * Calculates partial likelihoods at a node while first checking which partials
     * need to be calculated in the irreversible model based on the possible ancestral 
     * states given the states of the children.
     */
    public void calculatePartials(int childIndex1, int childIndex2, int parentIndex) {
        double[][] child1partials = partials[currentPartialsIndex[childIndex1]][childIndex1];
        double[][] child2partials = partials[currentPartialsIndex[childIndex2]][childIndex2];
        double[][] parentPartials = partials[currentPartialsIndex[parentIndex]][parentIndex];

        for (int k = 0; k < nrOfSites; k++) {
            int[] possibleParentStates = getPossibleStates(ancestralStates[currentPartialsIndex[parentIndex]][parentIndex][k], allStatesPerSite[k]);
            int[] child1States = getPossibleStates(ancestralStates[currentPartialsIndex[childIndex1]][childIndex1][k], allStatesPerSite[k]);
            int[] child2States = getPossibleStates(ancestralStates[currentPartialsIndex[childIndex2]][childIndex2][k], allStatesPerSite[k]);

            // Need to reset partials
            Arrays.fill(parentPartials[k], 0.0);

            // Calculate partials for all states
            for (int i : possibleParentStates) {
                double sum1 = 0.0;
                for (int j : child1States) {
                    sum1 += getTransitionProbFromCoreValues(i, j, nodeCoreTransitionValues[currentMatrixIndex[childIndex1]][childIndex1], k) * child1partials[k][j];
                }
                double sum2 = 0.0;
                for (int j : child2States) {
                    sum2 += getTransitionProbFromCoreValues(i, j, nodeCoreTransitionValues[currentMatrixIndex[childIndex2]][childIndex2], k) * child2partials[k][j];
                }
                parentPartials[k][i] = sum1 * sum2;
            }

            if (useScaling) {
                logScalingFactors[currentPartialsIndex[parentIndex]][parentIndex][k] = scalePartials(partials[currentPartialsIndex[parentIndex]][parentIndex][k], possibleParentStates);
            }
        }
    }

    public double getTransitionProbFromCoreValues(int parentState, int childState, double[] coreValues, int siteNum) {
        /* Reminder of core values orderering in the array:
        * 0: unedited to unedited (top left)
        * 1: edited to same edited (diag)
        * 2: any edit to missing data (last col)
        * 3: unedited to any edited without editRate considered yet (first row)
        */
        if (parentState == 0 && childState == 0) {
            return coreValues[0]; // unedited to unedited
        } else if (parentState == missingDataStatePerSite[siteNum] && childState == missingDataStatePerSite[siteNum]) {
            return 1.0; // Absorbing state once missing
        } else if (childState == missingDataStatePerSite[siteNum]) {
            return coreValues[2]; // Unedited or any edit to missing data
        } else if (parentState > 0 && childState != parentState) {
            return 0.0; // No transitions allowed when already edited, except to missing data state which is handled above
        } else if (childState == parentState) {
            return coreValues[1]; // Edit to the same edit, already accounting for the distinct [0][0] unedited to unedited case above
        } else if (parentState == 0 && childState > 0) {
            return coreValues[3] * substitutionModel.getEditRate(siteNum, childState); // unedited to any edited with editRate considered
        } else {
            throw new IllegalArgumentException("Invalid state combination: parentState=" + parentState + ", childState=" + childState);
        }
    }

    /*
     * Calculates site partial likelihoods and total log likelihood at the origin node.
     * The origin is known to be in the unedited state, so we can only calculate the partials
     * for transitions from that state only.
     */
    public double calculateLogLikelihoods(int rootIndex) {
        double lnl = 0.0;

        final double[][] rootPartials = partials[currentPartialsIndex[rootIndex]][rootIndex];

        for (int k = 0; k < nrOfSites; k++) {
            double sum = 0;    // Reset origin partials; Only the first position is needed, assuming the barcode starts unedited
            final int[] possibleRootStates = getPossibleStates(ancestralStates[currentPartialsIndex[rootIndex]][rootIndex][k], allStatesPerSite[k]);
            for (int rootState : possibleRootStates) {
                sum += getTransitionProbFromCoreValues(0, rootState, nodeCoreTransitionValues[currentMatrixIndex[rootIndex]][rootIndex], k) * rootPartials[k][rootState];
            }
            lnl += Math.log(sum);  // No scaling needed ever since its just one possible state at the origin
        }

        // Can just add the total log scaling factors at the end across all sites since we are in log space
        if (useScaling) lnl += accumulateLogScalingFactors();

        return lnl;
    }

    protected double scalePartials(double[] partials, int[] possibleStates) {
        // Find the maximum partial value for scaling
        double scaleFactor = 0.0;
        for (int j : possibleStates) {
            scaleFactor = Math.max(scaleFactor, partials[j]);
        }

        // Scale partials only if needed and set scaling factor
        if (scaleFactor < SCALING_THRESHOLD && scaleFactor > 0.0) {
            for (int j : possibleStates) {
                partials[j] /= scaleFactor;
            }
            return Math.log(scaleFactor);
        }
        return 0.0;
    }

    public double accumulateLogScalingFactors() {
        double totalLogScaling = 0.0;
        for (int i = 0; i < nrOfSites; i++) {
            for (int j = 0; j < numNodesNoOrigin; j++) {
                totalLogScaling += logScalingFactors[currentPartialsIndex[j]][j][i];
            }
        }
        return totalLogScaling;
    }

    public void setUseScaling(boolean status) {
        useScaling = status;
        if (useScaling) {
            logScalingFactors = new double[2][numNodesNoOrigin][nrOfSites]; // Reset log scaling factors when turning on scaling
        }
    }



    @Override
    public void store() {
        super.store();    // Store logP

        // Store likelihood core components
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, nrOfNodes);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, numNodesNoOrigin);
    }

    @Override
    public void restore() {
        super.restore();    // Restore logP

        // Restore likelihood core components
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, numNodesNoOrigin);
    }

    @Override
    protected boolean requiresRecalculation() {
        // Recalculate all transition probabilities and partials if any rates have changed
        if (((SiteModel.Base) siteModelInput.get()).isDirtyCalculation() ||
            branchRateModelInput.get().isDirtyCalculation()) {
            hasDirt = IS_DIRTY;
            return true;
        }

        // Check the tree last, and keep hasDirt as clean since recalculations will only need to be done node-wise
        hasDirt = IS_CLEAN;
        return treeInput.get().somethingIsDirty();
    }
}
