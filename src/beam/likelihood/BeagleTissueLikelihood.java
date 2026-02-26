package beam.likelihood;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.datatype.DataType;
import beastclassic.evolution.tree.TreeTrait;
import beastclassic.evolution.tree.TreeTraitProvider;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.BeagleInfo;


/**
 * Uses the BEAGLE library to calculate the tissue likelihood including the origin node.
 *
 * @author Stephen Staklinski
 */
@Description("Uses the beagle library to calculate the tissue likelihood including the origin node.")
public class BeagleTissueLikelihood extends GenericTreeLikelihood implements TreeTraitProvider {

    // Input for the origin parameter
    public Input<RealParameter> originInput = new Input<>("origin", "Start of the cell division process, usually start of the experiment.", Input.Validate.REQUIRED);
    // Input for the label used for tissues written to posterior samples
    public Input<String> tagInput = new Input<>("tag", "label used to report trait in the output tree posterior samples. Default is location.", "location", Input.Validate.OPTIONAL);

    // Preferred flags list
    private static List<Integer> preferredOrder = null;
    // Threshold for scaling partials
    private double scalingThreshold = 1.0E-100;
    // Number of states in the model
    int stateCount;
    // Number of nodes in the tree
    int nodeCount;
    // Input frequencies
    private double[] rootFrequencies;
    // Operations array
    private int[] operations;
    // Operation counts
    private int operationCount;
    // Buffer index helper for partials
    protected BufferIndexHelper partialBufferHelper;
    // Buffer index helper for matrices
    protected BufferIndexHelper matrixBufferHelper;
    // Buffer index helper for scaling
    protected BufferIndexHelper scaleBufferHelper;
    // Number of tip nodes
    protected int tipCount;
    // Number of internal nodes
    protected int internalNodeCount;
    // BEAGLE instance
    protected Beagle beagle;
    // Partial likelihoods for origin node
    protected double[] originPartials;
    // Stored partial likelihoods for origin node
    protected double[] storedOriginPartials;
    // Transition matrix for root branch
    protected double[] rootTransitionMatrix;
    // Stored transition matrix for root branch
    protected double[] storedRootTransitionMatrix;
    // Dirt flag for tree updates
    protected int hasDirt;
    protected int IS_CLEAN = 0;
    protected int IS_DIRTY = 1;
    // Data input
    protected Alignment data;
    // Substitution model
    protected SubstitutionModel substitutionModel;
    // Branch rate model
    protected BranchRateModel.Base branchRateModel;
    // Probability array
    protected double[] probabilities;
    // Whether to use scale factors
    protected boolean useScaleFactors = false;
    // Track how many scaling attempt have been done
    protected int numScalingAttempts = 0;
    // After this many attempts, try to turn off scaling
    protected static final int MAX_SCALING_ATTEMPTS = 1000;
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
    // Array to store root partials
    double[] rootPartials;

    @Override
    public void initAndValidate() {
        // Validate and setup site model
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }

        // Validate site categories
        if (((SiteModel.Base) siteModelInput.get()).getCategoryCount() > 1) throw new IllegalArgumentException("Site categories are not supported in the current implementation.");

        // Validate and setup substitution model
        substitutionModel = (SubstitutionModel.Base) ((SiteModel.Base) siteModelInput.get()).substModelInput.get();
        if (!substitutionModel.canReturnComplexDiagonalization()) throw new IllegalArgumentException("Substitution model must be able to return transition probabilities.");

        // Validate and setup branch rate model
        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) throw new IllegalArgumentException("Branch rate model must be specified in the current implementation.");

        // Get state count and site count
        stateCount = substitutionModel.getStateCount();
        data = dataInput.get();
        if (data.getSiteCount() > 1) throw new RuntimeException("Only one tissue per tip is allowed.");

        // Setup origin
        originHeight = originInput.get().getValue();

        // Setup arrays
        // Initialize frequency array
        rootFrequencies = substitutionModel.getFrequencies();
        // Initialize probability arrays
        probabilities = new double[stateCount * stateCount];
        // Setup tree node counts
        nodeCount = treeInput.get().getNodeCount();
        tipCount = treeInput.get().getLeafNodeCount();
        internalNodeCount = nodeCount - tipCount;
        // Setup origin arrays
        int partialsSize = stateCount;
        originPartials = new double[partialsSize];
        storedOriginPartials = new double[partialsSize];
        rootTransitionMatrix = new double[stateCount * stateCount];
        storedRootTransitionMatrix = new double[stateCount * stateCount];
        rootPartials = new double[stateCount];

        // Initialize other arrays
        operations = new int[internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
        operationCount = 0;

        // Setup BEAGLE
        setupBeagle();

        // Initialize tissue input states and output tag
        tag = tagInput.get();
        tipStates = new int[tipCount];
        dataType = data.getDataType();
        for (Node node : treeInput.get().getExternalNodes()) {
            int code = data.getPattern(data.getTaxonIndex(node.getID()), 0);
            tipStates[node.getNr()] = dataType.getStatesForCode(code)[0];
        }

        // Initialize reconstructed tissue states arrays
        reconstructedStates = new int[nodeCount];
        storedReconstructedStates = new int[nodeCount];
        
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

    // Calculate the log likelihood of the current state.
    @Override
    public double calculateLogP() {
        final Node root = treeInput.get().getRoot();

        // Do not allow the tree root to be older than the origin 
        if (root.getHeight() >= originHeight) {
            return Double.NEGATIVE_INFINITY;
        }

        // Setup rescaling
        if (useScaleFactors) {
            if (numScalingAttempts <= MAX_SCALING_ATTEMPTS) {
                useScaleFactors = true;
                hasDirt = IS_DIRTY;
                numScalingAttempts++;
            } else {
                System.out.println("Turning off scaling...");
                // Turn off scaling after many attempts, since it may no longer be needed
                useScaleFactors = false;
                hasDirt = IS_DIRTY;
                numScalingAttempts = 0;
            }
        }

        // Reset counters
        operationCount = 0;

        // Setup updates for transition matrices and partial likelihoods
        traverse(root, true);

        if (logP == Double.NEGATIVE_INFINITY || Double.isNaN(logP)) {
            System.out.println("Likelihood is negative infinity or NaN, trying to scale partials...");
            useScaleFactors = true;
            hasDirt = IS_DIRTY;
            numScalingAttempts++;
            operationCount = 0;
            traverse(root, false);  // False because we do not need to swap the buffer indices, we just want to overwrite the underflowed partials with scaled partials for the same state
            if (logP == Double.NEGATIVE_INFINITY || Double.isNaN(logP)) {
                throw new RuntimeException("Likelihood is still negative infinity after scaling");
            }
        }

        areStatesRedrawn = false;   // Reset flag for ancestral state sampling
        return logP;
    }


    private int getScaleBufferIndex() {
        int cumulateScaleBufferIndex = Beagle.NONE;
        if (useScaleFactors) {
            scaleBufferHelper.flipOffset(internalNodeCount);
            cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
            beagle.resetScaleFactors(cumulateScaleBufferIndex);
            beagle.accumulateScaleFactors(scaleBufferHelper.getOffsetIndices(), internalNodeCount, cumulateScaleBufferIndex);
        }
        return cumulateScaleBufferIndex;
    }

    private double calculateLikelihoodWithOrigin(Node root, int rootIndex) {
        double lnl = 0.0;

        // Get root partials from BEAGLE
        Arrays.fill(rootPartials, 0.0);
        beagle.getPartials(rootIndex, Beagle.NONE, rootPartials);

        // Setup root transition matrix from origin to root
        double br = branchRateModel.getRateForBranch(root);
        substitutionModel.getTransitionProbabilities(root, originHeight, root.getHeight(), br, probabilities);
        System.arraycopy(probabilities, 0, rootTransitionMatrix, 0, probabilities.length);

        // Calculate origin partials
        calculateOriginPartials(rootPartials, rootTransitionMatrix, originPartials);

        // Scale origin partials, if needed
        double originScaleFactorsSum = scaleOriginPartials();

        // Replace root partials with origin partials in BEAGLE (just a trick to get BEAGLE to calculate the likelihood through the origin with scale factors, if needed)
        beagle.setPartials(partialBufferHelper.getOffsetIndex(root.getNr()), originPartials);

        // Calculate final likelihood
        double[] sumLogLikelihoods = new double[1];
        int cumulateScaleBufferIndex = getScaleBufferIndex();
        beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

        // Restore original root partials
        beagle.setPartials(partialBufferHelper.getOffsetIndex(root.getNr()), rootPartials);

        // Calculate final log likelihood
        lnl = sumLogLikelihoods[0] + originScaleFactorsSum;
        return lnl;
    }

    private double scaleOriginPartials() {
        double originScaleFactorsSum = 0.0;
        if (useScaleFactors) {
            double scaleFactor = 0.0;
            for (int j = 0; j < stateCount; j++) {
                scaleFactor = Math.max(scaleFactor, originPartials[0 + j]);
            }
            if (scaleFactor < scalingThreshold && scaleFactor > 0.0) {
                for (int j = 0; j < stateCount; j++) {
                    originPartials[0 + j] /= scaleFactor;
                }
                originScaleFactorsSum += Math.log(scaleFactor);
            }
        }
        return originScaleFactorsSum;
    }


    // Traverse the tree to update transition probability matrices and partial likelihoods for each node
    private int traverse(Node node, boolean flip) {
        int nodeNum = node.getNr();
        int updateTransitionProbs = (node.isDirty() | hasDirt);

        // Update transition probability matrix
        if (!node.isRoot() && updateTransitionProbs != IS_CLEAN) {
            // Get the new transition probability matrix and store it in beagle
            substitutionModel.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(), branchRateModel.getRateForBranch(node), probabilities);
            if (flip) matrixBufferHelper.flipOffset(nodeNum);
            int matrixIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
            beagle.setTransitionMatrix(matrixIndex, probabilities, 1);
            updateTransitionProbs |= IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {
            // Traverse down the two child nodes to enforce post-order traversal
            Node child1 = node.getLeft();
            final int updateTransitionProbs1 = traverse(child1, flip);
            Node child2 = node.getRight();
            final int updateTransitionProbs2 = traverse(child2, flip);

            // If either child node was dirty, then update the parent node
            if (updateTransitionProbs1 != IS_CLEAN || updateTransitionProbs2 != IS_CLEAN) {
                int x = operationCount * Beagle.OPERATION_TUPLE_SIZE;

                // Flip places a flag to calculate these values later in beagle updatePartials()
                if (flip) partialBufferHelper.flipOffset(nodeNum);
                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

                if (useScaleFactors) {
                    int n = nodeNum - tipCount; // Get the index of this scaling buffer, since only internal node have scaling buffers
                    scaleBufferHelper.flipOffset(n);
                    operations[x + 1] = scaleBufferHelper.getOffsetIndex(n); // Write new scaleFactor
                    operations[x + 2] = Beagle.NONE;
                } else {
                    operations[x + 1] = Beagle.NONE; // Not using scaleFactors
                    operations[x + 2] = Beagle.NONE;
                }

                // specify operations for beagle to perform later in updatePartials()
                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNr()); // source node 1
                operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNr()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNr()); // source node 2
                operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNr()); // source matrix 2
                operationCount++;
                updateTransitionProbs |= (updateTransitionProbs1 | updateTransitionProbs2);
            }
        }

        // Once at the root after post-order traversal, then complete the full likelihood calculation through the origin node
        if (node.isRoot()) {
            // Tell BEAGLE to perform all of the updates from traverse() to calculate the partials for all nodes in the tree
            beagle.updatePartials(operations, operationCount, Beagle.NONE);
            // Calculate likelihood based on origin
            logP = calculateLikelihoodWithOrigin(node, partialBufferHelper.getOffsetIndex(node.getNr()));
        }

        return updateTransitionProbs;
    }

    // Calculate partials for the origin node of degree 1 that goes from the start of the experiment to the root
    protected double[] calculateOriginPartials(double[] partials1, double[] matrices1, double[] partials3) {
        double sum1;
        int u = 0;
        int w = 0;
        for (int i = 0; i < stateCount; i++) {
            sum1 = 0.0;
            for (int j = 0; j < stateCount; j++) {
                sum1 += matrices1[w] * partials1[0 + j];
                w++;
            }
            partials3[u] = sum1;
            u++;
        }
        return partials3;
    }


    @Override
    protected boolean requiresRecalculation() {

        if (((CalculationNode) substitutionModel).isDirtyCalculation() || 
            ((SiteModel.Base) siteModelInput.get()).isDirtyCalculation() ||
            branchRateModel.isDirtyCalculation()) {
            hasDirt = IS_DIRTY;
            return true;
        }

        hasDirt = IS_CLEAN;
        return treeInput.get().somethingIsDirty();
    }

    // Stores the state
    @Override
    public void store() {
        partialBufferHelper.storeState();
        matrixBufferHelper.storeState();
        if (useScaleFactors) scaleBufferHelper.storeState();

        // Store origin partials
        System.arraycopy(originPartials, 0, storedOriginPartials, 0, originPartials.length);
        // Store root to origin branch transition matrix
        System.arraycopy(rootTransitionMatrix, 0, storedRootTransitionMatrix, 0, rootTransitionMatrix.length);
        // Store logP and reset isDirty to false
        super.store();
        // Store ancestral states
        System.arraycopy(reconstructedStates, 0, storedReconstructedStates, 0, reconstructedStates.length);
        storedAreStatesRedrawn = areStatesRedrawn;
    }

    // Restores the state
    @Override
    public void restore() {
        partialBufferHelper.restoreState();
        matrixBufferHelper.restoreState();
        if (useScaleFactors) scaleBufferHelper.restoreState();

        // restore origin partials
        System.arraycopy(storedOriginPartials, 0, originPartials, 0, originPartials.length);
        // Restore root to origin branch transition matrix
        System.arraycopy(storedRootTransitionMatrix, 0, rootTransitionMatrix, 0, rootTransitionMatrix.length);
        // Restore logP and reset isDirty to false
        logP = storedLogP;
        // Reset isDirty to false
        super.restore(); 
        // Restore ancestral states
        System.arraycopy(storedReconstructedStates, 0, reconstructedStates, 0, storedReconstructedStates.length);
        areStatesRedrawn = storedAreStatesRedrawn;
    }

    private static List<Integer> parseSystemPropertyIntegerArray(String propertyName) {
        List<Integer> order = new ArrayList<>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    int n = Integer.parseInt(part.trim());
                    order.add(n);
                } catch (NumberFormatException nfe) {
                	Log.warning.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    // Basic initialization of a beagle instance
    private void setupBeagle() {
        // Initialize buffer helpers
        partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount); // Partials for tips are already set
        matrixBufferHelper = new BufferIndexHelper(nodeCount, 0); // Need transition matrix for all nodes
        scaleBufferHelper = new BufferIndexHelper(internalNodeCount + 1, 0);    // +1 because we need one buffer for each internal node and one buffer for the cumulative scale factors

        // Load system properties
        preferredOrder = parseSystemPropertyIntegerArray("beagle.preferred.flags");

        // Initialize BEAGLE instance
        long preferenceFlags = 0;
        if (preferredOrder.size() > 0) preferenceFlags = preferredOrder.get(0);
        long requirementFlags = 0;
        requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
        try {
            beagle = BeagleFactory.loadBeagleInstance(
                tipCount,
                partialBufferHelper.getBufferCount(),
                tipCount,
                stateCount,
                1,
                1,
                matrixBufferHelper.getBufferCount(),
                1,
                scaleBufferHelper.getBufferCount(),
                null,
                preferenceFlags,
                requirementFlags
            );
        } catch (Exception e) {
            Log.warning.println("Error setting up BEAGLE. Check install.");
            System.exit(1);
        }
        Log.info.println("Using BEAGLE version: " + BeagleInfo.getVersion());

        // Setup tip states and pattern weights
        Node[] nodes = treeInput.get().getNodesAsArray();
        for (int i = 0; i < tipCount; i++) {
            int taxon = dataInput.get().getTaxonIndex(nodes[i].getID());
            beagle.setTipStates(i, data.getDataType().getStatesForCode(data.getPattern(taxon, 0)));
        }

        beagle.setPatternWeights(new double[]{1.0});
        beagle.setCategoryWeights(0, new double[]{1.0});
        beagle.setStateFrequencies(0, rootFrequencies); // Root tissue frequencies can not change during a run in the current implementation
    }

    /**
     * Helper class to manage buffer indices for BEAGLE operations.
     * Handles both single and mirrored buffer indices for store/restore operations.
     */
    public class BufferIndexHelper {
        private final int maxIndexValue;
        private final int minIndexValue;
        private final int offsetCount;
        private int[] indexOffsets;
        private int[] storedIndexOffsets;

        public BufferIndexHelper(int maxIndexValue, int minIndexValue) {
            this.maxIndexValue = maxIndexValue;
            this.minIndexValue = minIndexValue;
            this.offsetCount = maxIndexValue - minIndexValue;
            this.indexOffsets = new int[offsetCount];
            this.storedIndexOffsets = new int[offsetCount];
        }

        public int getBufferCount() {
            return 2 * offsetCount + minIndexValue;
        }

        public int getOffsetIndex(int i) {
            return i < minIndexValue ? i : indexOffsets[i - minIndexValue] + i;
        }

        public void flipOffset(int i) {
            if (i >= minIndexValue) {
                indexOffsets[i - minIndexValue] = offsetCount - indexOffsets[i - minIndexValue];
            }
        }

        public int[] getOffsetIndices() {
            return indexOffsets;
        }

        public void storeState() {
            System.arraycopy(indexOffsets, 0, storedIndexOffsets, 0, indexOffsets.length);
        }

        public void restoreState() {
            System.arraycopy(storedIndexOffsets, 0, indexOffsets, 0, storedIndexOffsets.length);
        }
    }

    // Gets the state for a given node in the tree.
    public int getStateForNode(TreeInterface tree, Node node) {
        if (!areStatesRedrawn) {
            sampleAncestralStates(tree, tree.getRoot(), -1);
            areStatesRedrawn = true;
        }
        return reconstructedStates[node.getNr()];
    }

    // Traverses the tree in pre-order to sample internal node states
    public void sampleAncestralStates(TreeInterface tree, Node node, int parentState) {
        int nodeNum = node.getNr();
        double[] conditionalProbabilities = new double[stateCount];

        if (!node.isLeaf()) {
            // Handle an internal node
            double[] partialLikelihood = new double[stateCount];
            beagle.getPartials(partialBufferHelper.getOffsetIndex(node.getNr()), Beagle.NONE, partialLikelihood);

            if (node.isRoot()) {
                // Handle root node (below origin) when origin is present since the origin is not physically in the tree
                double[] oPs = new double[stateCount];
                System.arraycopy(originPartials, 0, oPs, 0, originPartials.length);
                for (int i = 0; i < stateCount; i++) {
                    oPs[i] *= rootFrequencies[i];
                }

                parentState = Randomizer.randomChoicePDF(oPs);
                probabilities = rootTransitionMatrix;
            }
            else {
                // Handle any other internal node with a parent in the tree
                beagle.getTransitionMatrix(matrixBufferHelper.getOffsetIndex(node.getNr()), probabilities);
            }

            // Calculate conditional probabilities
            int parentIndex = parentState * stateCount;
            for (int i = 0; i < stateCount; i++) {
                conditionalProbabilities[i] = partialLikelihood[i] * probabilities[parentIndex + i];
            }

            // Sample state
            reconstructedStates[node.getNr()] = Randomizer.randomChoicePDF(conditionalProbabilities);

            // Traverse children
            Node child1 = node.getChild(0);
            sampleAncestralStates(treeInput.get(), child1, reconstructedStates[node.getNr()]);
            Node child2 = node.getChild(1);
            sampleAncestralStates(treeInput.get(), child2, reconstructedStates[node.getNr()]);
            
        } else {
            // Handle leaf node where tissue is known
            reconstructedStates[node.getNr()] = tipStates[node.getNr()];
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