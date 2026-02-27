package beam.likelihood;

import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.datatype.DataType;
import beastclassic.evolution.tree.TreeTrait;
import beastclassic.evolution.tree.TreeTraitProvider;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;


/**
 * Calculate the tissue likelihood including the origin node.
 *
 * @author Stephen Staklinski
 */
@Description("Calculate the tissue likelihood including the origin node.")
public class TissueLikelihood extends GenericTreeLikelihood implements TreeTraitProvider {

    // Input for the origin parameter
    public Input<RealParameter> originInput = new Input<>("origin", "Start of the cell division process, usually start of the experiment.", Validate.REQUIRED);
    // Input for the label used for tissues written to posterior samples
    public Input<String> tagInput = new Input<>("tag", "label used to report trait in the output tree posterior samples. Default is location.", "location", Validate.OPTIONAL);

    // Substitution model
    protected SubstitutionModel.Base substitutionModel;
    // Site model
    SiteModel.Base siteModel;
    // Number of states in the model
    int stateCount;
    // Input frequencies
    private double[] rootFrequencies;
    /* Number of nodes in the tree */
    protected int nrOfNodes;
    /* Number of nodes excluding the origin node */
    protected int numNodesNoOrigin;
    /* Partial likelihoods for each node */
    protected double[][][] partials;
    // Transition matrix for each node
    protected double[][][] matrices;
    /* Current matrix indices for each node */
    protected int[] currentMatrixIndex;
    /* Stored matrix indices for each node */
    protected int[] storedMatrixIndex;
    /* Current partials indices for each node */
    protected int[] currentPartialsIndex;
    /* Stored partials indices for each node */
    protected int[] storedPartialsIndex;
    // Dirt flag for tree updates
    protected int hasDirt;
    protected int IS_CLEAN = 0;
    protected int IS_DIRTY = 1;
    // Data input
    protected Alignment data;
    /* Whether to use scaling for numerical stability */
    protected boolean useScaling = false;
    protected boolean storedUseScaling = false;
    /* Log scaling factors sum for each site */
    protected double[][] logScalingFactors;
    // Track how many scaling attempt have been done
    protected int numScalingAttempts = 0;
    protected int storedNumScalingAttempts = 0;
    // After this many attempts, try to turn off scaling
    protected static final int MAX_SCALING_ATTEMPTS = 1000;
    // Threshold for scaling partials
    private static final double SCALING_THRESHOLD = 1.0E-100;
    // Origin height
    protected Double originHeight;
    // Data type for the alignment
    protected DataType dataType;
    // Reconstructed states for each node
    private int[] reconstructedStates;
    // Stored reconstructed states for store/restore operations
    private int[] storedReconstructedStates;
    // Tag label for trait reporting
    private String tag;
    // Flag indicating if states have been redrawn
    private boolean areStatesRedrawn = false;
    // Stored flag for states redrawn status
    private boolean storedAreStatesRedrawn = false;
    // States at tip nodes
    private int[] tipStates;
    // Helper for tree traits
    protected TreeTraitProvider.Helper treeTraits = new Helper();
    // Origin node index (not physically in the tree, but we need to store partials for it)
    int originIndex;

    public TissueLikelihood() {
        branchRateModelInput.setRule(Validate.REQUIRED);
        siteModelInput.setType(SiteModel.Base.class);
    }

    @Override
    public void initAndValidate() {
        siteModel = (SiteModel.Base) siteModelInput.get();
        if (siteModel.getCategoryCount() > 1) throw new IllegalArgumentException("Site categories are not supported in the current implementation.");

        substitutionModel = (SubstitutionModel.Base) siteModel.substModelInput.get();
        if (!substitutionModel.canReturnComplexDiagonalization()) throw new IllegalArgumentException("Substitution model must be able to return transition probabilities.");

        data = dataInput.get();
        if (data.getSiteCount() > 1) throw new RuntimeException("Only one tissue per tip is allowed.");

        originHeight = originInput.get().getValue();
        rootFrequencies = substitutionModel.getFrequencies();
        numNodesNoOrigin = treeInput.get().getNodeCount();
        nrOfNodes = numNodesNoOrigin + 1; // +1 for origin node
        originIndex = numNodesNoOrigin; // Origin node is the last index in 0-based indexing
        stateCount = substitutionModel.getStateCount();
        
        partials = new double[2][nrOfNodes][];
        for (int i = 0; i < nrOfNodes; i++) {
            partials[0][i] = new double[stateCount];
            partials[1][i] = new double[stateCount];
        }
        currentPartialsIndex = new int[nrOfNodes];
        storedPartialsIndex = new int[nrOfNodes];
        logScalingFactors = new double[2][nrOfNodes];

        matrices = new double[2][numNodesNoOrigin][];
        for (int i = 0; i < numNodesNoOrigin; i++) {
            matrices[0][i] = new double[stateCount * stateCount];
            matrices[1][i] = new double[stateCount * stateCount];
        }
        currentMatrixIndex = new int[numNodesNoOrigin];
        storedMatrixIndex = new int[numNodesNoOrigin];

        // Set initial states
        setStates(treeInput.get().getRoot());

        // Initialize tissue input states and output tag
        tag = tagInput.get();
        tipStates = new int[treeInput.get().getLeafNodeCount()];
        dataType = data.getDataType();
        for (Node node : treeInput.get().getExternalNodes()) {
            int code = data.getPattern(data.getTaxonIndex(node.getID()), 0);
            tipStates[node.getNr()] = dataType.getStatesForCode(code)[0];
        }

        // Initialize reconstructed tissue states arrays
        reconstructedStates = new int[numNodesNoOrigin];
        storedReconstructedStates = new int[numNodesNoOrigin];
        
        // Initialize the input observed tissue state for each tip in the tree
        treeTraits.addTrait("states", new TreeTrait.IA() {
            @Override
            public String getTraitName() { return tag; }
            @Override
            public Intent getIntent() { return Intent.NODE; }
            @Override
            public int[] getTrait(TreeInterface tree, Node node) { return new int[] { getStateForNode(tree, node) }; }
            @Override
            public String getTraitString(TreeInterface tree, Node node) { return "\"" + dataType.getCode(getStateForNode(tree, node)) + "\""; }
        });
    }

    // Sets the observed data in likelihood core partials for leaf nodes across all alignment sites
    protected void setStates(Node node) {
        if (node.isLeaf()) {
            int taxonIndex = data.getTaxonIndex(node.getID());
            if (taxonIndex == -1) throw new RuntimeException("Could not find sequence " + node.getID() + " in the alignment");
            int leafIndex = node.getNr();
            partials[currentPartialsIndex[leafIndex]][leafIndex][data.getCounts().get(taxonIndex).get(0)] = 1.0;
        } else {
            setStates(node.getLeft());
            setStates(node.getRight());
        }
    }

    // Calculate the log likelihood of the current state.
    @Override
    public double calculateLogP() {
        final Node root = treeInput.get().getRoot();
        // Do not allow the tree root to be older than the origin 
        if (root.getHeight() >= originHeight) {
            areStatesRedrawn = false;   // Reset flag for ancestral state sampling
            return Double.NEGATIVE_INFINITY;
        }
        // Setup rescaling
        if (useScaling) {
            if (numScalingAttempts >= MAX_SCALING_ATTEMPTS) {
                setUseScaling(false);
                numScalingAttempts = 0;
                hasDirt = IS_DIRTY;
            } else {
                numScalingAttempts++;
            }
        }

        // Get new models
        siteModel = (SiteModel.Base) siteModelInput.get();
        substitutionModel = (SubstitutionModel.Base) siteModel.substModelInput.get();

        // Calculate likelihood with optional scaling
        traverse(root);
        if (logP == Double.NEGATIVE_INFINITY || Double.isNaN(logP)) {
            restore(); // Reset indicex to not double flip and overwrite restore state
            setUseScaling(true);
            hasDirt = IS_DIRTY;
            numScalingAttempts++;
            traverse(root);
            if (logP == Double.NEGATIVE_INFINITY || Double.isNaN(logP)) {
                throw new RuntimeException("Likelihood is still negative infinity after scaling");
            }
        }

        areStatesRedrawn = false;   // Reset flag for ancestral state sampling
        return logP;
    }


    // Traverse the tree to update transition probability matrices and partial likelihoods for each node
    private int traverse(Node node) {
        int nodeIndex = node.getNr();
        int updateTransitionProbs = (node.isDirty() | hasDirt);

        // Update transition probability matrix
        if (updateTransitionProbs != IS_CLEAN) {
            currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
            double parentHeight = node.isRoot() ? originHeight : node.getParent().getHeight();
            substitutionModel.getTransitionProbabilities(node, parentHeight, node.getHeight(), branchRateModelInput.get().getRateForBranch(node), matrices[currentMatrixIndex[nodeIndex]][nodeIndex]);
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
                currentPartialsIndex[nodeIndex] = 1 - currentPartialsIndex[nodeIndex];
                calculatePartials(childIndex1, childIndex2, nodeIndex);

                updateTransitionProbs |= (updateTransitionProbs1 | updateTransitionProbs2);
            }
        }

        // Always calculate origin partials once at the root
        if (node.isRoot()) {
            currentPartialsIndex[originIndex] = 1 - currentPartialsIndex[originIndex];
            logP = calculateLogLikelihood(nodeIndex, originIndex);
        }

        return updateTransitionProbs;
    }

    public void calculatePartials(int childIndex1, int childIndex2, int parentIndex) {
        double[] child1partials = partials[currentPartialsIndex[childIndex1]][childIndex1];
        double[] child2partials = partials[currentPartialsIndex[childIndex2]][childIndex2];
        double[] parentPartials = partials[currentPartialsIndex[parentIndex]][parentIndex];
        double[] child1Matrix = matrices[currentMatrixIndex[childIndex1]][childIndex1];
        double[] child2Matrix = matrices[currentMatrixIndex[childIndex2]][childIndex2];

        // Need to reset partials
        Arrays.fill(parentPartials, 0.0);

        // Calculate partials for all states
        for (int i = 0; i < stateCount; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            int rowOffset = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                sum1 += child1Matrix[rowOffset + j] * child1partials[j];
                sum2 += child2Matrix[rowOffset + j] * child2partials[j];
            }
            parentPartials[i] = sum1 * sum2;
        }

        if (useScaling) {
            logScalingFactors[currentPartialsIndex[parentIndex]][parentIndex] = scalePartials(parentPartials);
        }
    }

    // Calculate log likelihood after calculating partials for the origin
    public double calculateLogLikelihood(int rootIndex, int originIndex) {
        double lnl = 0.0;

        double[] rootPartials = partials[currentPartialsIndex[rootIndex]][rootIndex];
        double[] rootTransitionMatrix = matrices[currentMatrixIndex[rootIndex]][rootIndex];
        double[] originPartials = partials[currentPartialsIndex[originIndex]][originIndex];

        // Need to reset partials
        Arrays.fill(originPartials, 0.0);

        // Calculate partials for all states
        for (int i = 0; i < stateCount; i++) {
            int rowOffset = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                originPartials[i] += rootTransitionMatrix[rowOffset + j] * rootPartials[j];
            }
        }

        if (useScaling) {
            logScalingFactors[currentPartialsIndex[originIndex]][originIndex] = scalePartials(originPartials);
        }

        for (int i = 0; i < stateCount; i++) {
            lnl += rootFrequencies[i] * originPartials[i];
        }
        lnl = Math.log(lnl);

        if (useScaling) {
            lnl += accumulateLogScalingFactors();
        }

        return lnl;
    }

    protected double scalePartials(double[] partials) {
        // Find the maximum partial value for scaling
        double scaleFactor = 0.0;
        for (int j = 0; j < stateCount; j++) {
            scaleFactor = Math.max(scaleFactor, partials[j]);
        }

        // Scale partials only if needed and set scaling factor
        if (scaleFactor < SCALING_THRESHOLD && scaleFactor > 0.0) {
            for (int j = 0; j < stateCount; j++) {
                partials[j] /= scaleFactor;
            }
            return Math.log(scaleFactor);
        }
        return 0.0;
    }

    public double accumulateLogScalingFactors() {
        double totalLogScaling = 0.0;
        for (int j = 0; j < nrOfNodes; j++) {
            totalLogScaling += logScalingFactors[currentPartialsIndex[j]][j];
        }
        return totalLogScaling;
    }

    public void setUseScaling(boolean status) {
        useScaling = status;
        if (useScaling) {
            logScalingFactors = new double[2][nrOfNodes]; // Reset log scaling factors when turning on scaling
        }
    }

    @Override
    protected boolean requiresRecalculation() {
        siteModel = (SiteModel.Base) siteModelInput.get();
        substitutionModel = (SubstitutionModel.Base) siteModel.substModelInput.get();
        if (substitutionModel.isDirtyCalculation() || 
            siteModel.isDirtyCalculation() || 
            ((CalculationNode) branchRateModelInput.get()).isDirtyCalculation()) {
            hasDirt = IS_DIRTY;
            return true;
        }

        hasDirt = IS_CLEAN;
        return treeInput.get().somethingIsDirty();
    }

    // Stores the state
    @Override
    public void store() {
        super.store();    // Store logP

        // Store likelihood core components
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, currentMatrixIndex.length);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, currentPartialsIndex.length);
        storedUseScaling = useScaling;
        storedNumScalingAttempts = numScalingAttempts;

        // Store ancestral states
        System.arraycopy(reconstructedStates, 0, storedReconstructedStates, 0, reconstructedStates.length);
        storedAreStatesRedrawn = areStatesRedrawn;
    }

    // Restores the state
    @Override
    public void restore() {
        super.restore();    // Restore logP

        // Restore likelihood core components
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, storedMatrixIndex.length);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, storedPartialsIndex.length);
        useScaling = storedUseScaling;
        numScalingAttempts = storedNumScalingAttempts;

        // Restore ancestral states
        System.arraycopy(storedReconstructedStates, 0, reconstructedStates, 0, storedReconstructedStates.length);
        areStatesRedrawn = storedAreStatesRedrawn;
    }

    // Gets the state for a given node in the tree.
    public int getStateForNode(TreeInterface tree, Node node) {
        if (!areStatesRedrawn) {
            sampleAncestralStates(tree.getRoot(), -1);
            areStatesRedrawn = true;
        }
        return reconstructedStates[node.getNr()];
    }

    // Traverses the tree in pre-order to sample internal node states
    public void sampleAncestralStates(Node node, int parentState) {
        int nodeNum = node.getNr();
        double[] conditionalProbabilities = new double[stateCount];

        if (!node.isLeaf()) {
            // Handle an internal node
            double[] partialLikelihood = partials[currentPartialsIndex[nodeNum]][nodeNum];

            if (node.isRoot()) {
                // Handle root node (below origin) when origin is present since the origin is not physically in the tree
                double[] oPs = new double[stateCount];
                System.arraycopy(partials[currentPartialsIndex[originIndex]][originIndex], 0, oPs, 0, stateCount);
                for (int i = 0; i < stateCount; i++) {
                    oPs[i] *= rootFrequencies[i];
                }
                parentState = Randomizer.randomChoicePDF(oPs); // Looks like samples here, but really it will just pick the known origin tissue state since I have not implemented unknown origin state calculations yet
            }

            // Calculate conditional probabilities
            int rowOffset = parentState * stateCount;
            double[] probabilities = matrices[currentMatrixIndex[nodeNum]][nodeNum];
            // Calculate the conditional probability of being in a state at the recipient given a known sampled parent state and the partial likelihoods of the recipient state given the tree below which were precalculated already
            for (int i = 0; i < stateCount; i++) {
                conditionalProbabilities[i] = probabilities[rowOffset + i] * partialLikelihood[i];
            }

            // Sample state from conditional probabilities
            reconstructedStates[node.getNr()] = Randomizer.randomChoicePDF(conditionalProbabilities);

            // Traverse children
            Node child1 = node.getChild(0);
            sampleAncestralStates(child1, reconstructedStates[nodeNum]);
            Node child2 = node.getChild(1);
            sampleAncestralStates(child2, reconstructedStates[nodeNum]);
            
        } else {
            // Handle leaf node where tissue is known
            reconstructedStates[nodeNum] = tipStates[nodeNum];
        }
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }


    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }
}