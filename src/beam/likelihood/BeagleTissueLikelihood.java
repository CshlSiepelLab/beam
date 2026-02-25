package beam.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
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
    // Dimensions of the transition matrix
    int matrixDimensions;
    // Transition matrices
    private double[] matrices;
    // Current frequencies
    private double[] currentFreqs;
    // Matrix update indices
    private int[][] matrixUpdateIndices;
    // Branch lengths
    private double[][] branchLengths;
    // Branch update counts
    private int[] branchUpdateCount;
    // Operations array
    private int[][] operations;
    // Operation counts
    private int operationCount;
    // Buffer index helper for partials
    protected BufferIndexHelper partialBufferHelper;
    // Buffer index helper for eigen values
    private BufferIndexHelper eigenBufferHelper;
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
    // Branch lengths
    protected double[] m_branchLengths;
    // Stored branch lengths
    protected double[] storedBranchLengths;
    // Dirt flag for tree updates
    protected int hasDirt;
    // Substitution model
    protected SubstitutionModel substitutionModel;
    // Site model
    protected SiteModel.Base siteModel;
    // Branch rate model
    protected BranchRateModel.Base branchRateModel;
    // Pattern log likelihoods
    protected double[] patternLogLikelihoods;
    // Probability array
    protected double[] probabilities;
    // Current rescaling frequency
    private int rescalingFrequency = 10000;
    // Number of times to rescale
    private static final int RESCALE_TIMES = 1;
    // Whether to use scale factors
    protected boolean useScaleFactors = false;
    // Whether to use auto scaling
    private boolean useAutoScaling = false;
    // Whether to recompute scale factors
    private boolean recomputeScaleFactors = false;
    // Whether underflow has occurred
    private boolean everUnderflowed = false;
    // Rescaling counter
    private int rescalingCount = 0;
    // Inner rescaling counter
    private int rescalingCountInner = 0;
    // Scale buffer indices
    private int[] scaleBufferIndices;
    // Stored scale buffer indices
    private int[] storedScaleBufferIndices;
    // Origin height
    protected Double originHeight;
    // Key for states trait
    public static final String STATES_KEY = "states";
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

    @Override
    public void initAndValidate() {
        // Validate and setup site model
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        siteModel = (SiteModel.Base) siteModelInput.get();
        Alignment data = dataInput.get();
        dataType = data.getDataType();
        siteModel.setDataType(dataType);
        // Validate site categories
        if (siteModel.getCategoryCount() > 1) {
            throw new IllegalArgumentException("Site categories are not supported in the current implementation.");
        }
        // Validate and setup substitution model
        substitutionModel = (SubstitutionModel.Base) siteModel.substModelInput.get();
        if (!substitutionModel.canReturnComplexDiagonalization()) {
            throw new IllegalArgumentException("Substitution model must be able to return transition probabilities.");
        }
        // Validate and setup branch rate model
        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) {
            throw new IllegalArgumentException("Branch rate model must be specified in the current implementation.");
        }

        // Get state count and site count
        stateCount = substitutionModel.getStateCount();
        if (dataInput.get().getSiteCount() > 1) {
            throw new RuntimeException("Only one tissue per tip is allowed.");
        }

        // Setup origin
        originHeight = originInput.get().getValue();

        // Setup arrays
        // Initialize frequency array
        currentFreqs = new double[stateCount];
        // Setup matrices
        matrixDimensions = stateCount * stateCount;
        // Initialize probability arrays
        probabilities = new double[matrixDimensions];
        matrices = new double[matrixDimensions];
        // Setup tree node counts
        nodeCount = treeInput.get().getNodeCount();
        tipCount = treeInput.get().getLeafNodeCount();
        internalNodeCount = nodeCount - tipCount;
        // Initialize branch length arrays
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];
        // Setup origin arrays
        int partialsSize = stateCount;
        originPartials = new double[partialsSize];
        storedOriginPartials = new double[partialsSize];
        rootTransitionMatrix = new double[matrixDimensions];
        storedRootTransitionMatrix = new double[matrixDimensions];

        // Initialize other arrays
        patternLogLikelihoods = new double[1];
        matrixUpdateIndices = new int[1][nodeCount];
        branchLengths = new double[1][nodeCount];
        branchUpdateCount = new int[1];
        scaleBufferIndices = new int[internalNodeCount];
        storedScaleBufferIndices = new int[internalNodeCount];
        operations = new int[1][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
        operationCount = 0;

        // Setup BEAGLE
        setupBeagle();

        // Initialize tissue input states and output tag
        tag = tagInput.get();
        tipStates = new int[tipCount];
        for (Node node : treeInput.get().getExternalNodes()) {
            int code = data.getPattern(data.getTaxonIndex(node.getID()), 0);
            tipStates[node.getNr()] = dataType.getStatesForCode(code)[0];
        }

        // Initialize reconstructed tissue states arrays
        reconstructedStates = new int[nodeCount];
        storedReconstructedStates = new int[nodeCount];
        
        // Initialize the input observed tissue state for each tip in the tree
        treeTraits.addTrait(STATES_KEY, new TreeTrait.IA() {
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

        // Do not allow the tree root to be older than the origin 
        if (treeInput.get().getRoot().getHeight() >= originHeight) {
            return Double.NEGATIVE_INFINITY;
        }

        // Setup rescaling
        recomputeScaleFactors = false;
        if (everUnderflowed) {
            useScaleFactors = true;
            if (rescalingCountInner < RESCALE_TIMES) {
                recomputeScaleFactors = true;
                hasDirt = Tree.IS_FILTHY;
            }
            // Update rescaling counters
            rescalingCountInner++;
            rescalingCount++;
            if (rescalingCount > rescalingFrequency) {
                rescalingCount = 0;
                rescalingCountInner = 0;
            }
        }

        // Reset counters
        for (int i = 0; i < 1; i++) {
            branchUpdateCount[i] = 0;
        }
        operationCount = 0;

        // Calculate likelihood
        final Node root = treeInput.get().getRoot();
        traverse(root, true);

        boolean done;
        boolean firstRescaleAttempt = true;
        do {
            // Update partials in BEAGLE
            beagle.updatePartials(operations[0], operationCount, Beagle.NONE);
            int rootIndex = partialBufferHelper.getOffsetIndex(root.getNr());

            // Setup frequencies
            double[] frequencies = substitutionModel.getFrequencies();
            if (frequencies != currentFreqs) { beagle.setStateFrequencies(0, frequencies); }
            System.arraycopy(frequencies, 0, currentFreqs, 0, frequencies.length);

            // Calculate likelihood based on origin
            logP = calculateLikelihoodWithOrigin(root, rootIndex);

            if (Double.isNaN(logP) || Double.isInfinite(logP)) {
                everUnderflowed = true;
                logP = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt) {
                    useScaleFactors = true;
                    recomputeScaleFactors = true;

                    for (int i = 0; i < 1; i++) {
                        branchUpdateCount[i] = 0;
                    }

                    operationCount = 0;
                    traverse(root, false);
                    done = false;
                }
                done = true;
            }
            done = true;
            firstRescaleAttempt = false;
        } while (!done);

        areStatesRedrawn = false;   // Reset flag for ancestral state sampling
        return logP;
    }


    private int getScaleBufferIndex() {
        int cumulateScaleBufferIndex = Beagle.NONE;
        if (useScaleFactors) {
            if (recomputeScaleFactors) {
                scaleBufferHelper.flipOffset(internalNodeCount);
                cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                beagle.resetScaleFactors(cumulateScaleBufferIndex);
                beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, cumulateScaleBufferIndex);
            } else {
                cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
            }
        }
        return cumulateScaleBufferIndex;
    }

    private double calculateLikelihoodWithOrigin(Node root, int rootIndex) {
        // Get root partials
        double[] rootPartials = new double[stateCount];
        beagle.getPartials(rootIndex, Beagle.NONE, rootPartials);

        // Setup root transition matrix
        double br = branchRateModel.getRateForBranch(root);
        substitutionModel.getTransitionProbabilities(root, originHeight, root.getHeight(), br, probabilities);
        System.arraycopy(probabilities, 0, rootTransitionMatrix, 0, matrixDimensions);

        // Calculate origin partials
        calculateOriginPartials(rootPartials, rootTransitionMatrix, originPartials);

        // Scale origin partials if needed
        double originScaleFactorsSum = scaleOriginPartials();

        // Calculate final likelihood
        // Replace root partials with origin partials
        beagle.setPartials(partialBufferHelper.getOffsetIndex(root.getNr()), originPartials);

        // Calculate likelihood
        double[] sumLogLikelihoods = new double[1];
        int cumulateScaleBufferIndex = getScaleBufferIndex();
        beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0}, new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

        // Restore original root partials
        beagle.setPartials(partialBufferHelper.getOffsetIndex(root.getNr()), rootPartials);

        // Calculate final log likelihood
        double finalLogL = sumLogLikelihoods[0] + originScaleFactorsSum;
        return finalLogL;
    }

    private double scaleOriginPartials() {
        double originScaleFactorsSum = 0.0;
        if (useScaleFactors) {
            double scaleFactor = calculateScaleFactor(0, 0);
            if (scaleFactor < scalingThreshold) {
                applyScaleFactor(0, 0, scaleFactor);
                originScaleFactorsSum += Math.log(scaleFactor);
            }
        }
        return originScaleFactorsSum;
    }

    private double calculateScaleFactor(int pattern, int startIndex) {
        double scaleFactor = 0.0;
        for (int j = 0; j < stateCount; j++) {
            scaleFactor = Math.max(scaleFactor, originPartials[startIndex + j]);
        }
        return scaleFactor;
    }

    private void applyScaleFactor(int pattern, int startIndex, double scaleFactor) {
        for (int j = 0; j < stateCount; j++) {
            originPartials[startIndex + j] /= scaleFactor;
        }
    }


    /**
     * Traverse the tree to update transition probability matrices and subsequently calculate partial likelihoods for each node.
     *
     * @param node           node
     * @param flip           flip
     * @return boolean
     */
    private int traverse(Node node, boolean flip) {

        int nodeNum = node.getNr();

        // Decide if this node needs to be updated
        int update = (node.isDirty() | hasDirt);

        // Get the clock rate for the branch
        final double branchRate = branchRateModel.getRateForBranch(node);

        /* Calculate the branch length in number of substitutions to store it, where it is dependent on 
        the clock rate to convert realTimeLength * clockRate = numSubstitutionsLength */ 
        final double branchTime = node.getLength() * branchRate;

        // Update if its not the root and the node is dirty or the branch length has changed
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeNum])) {

            // Store the current node branch length in case it was changed, causing the update
            m_branchLengths[nodeNum] = branchTime;
            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            // Flip places a flag to calculate these values later in beagle updatePartials()
            if (flip) {
                matrixBufferHelper.flipOffset(nodeNum);
            }

            // Set which matrix to update
            final int eigenIndex = 0;
            final int updateCount = branchUpdateCount[eigenIndex];
            matrixUpdateIndices[eigenIndex][updateCount] = matrixBufferHelper.getOffsetIndex(nodeNum);

            // Get the new transition probability matrix and store it in beagle
            substitutionModel.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(), branchRate, probabilities);
            System.arraycopy(probabilities, 0, matrices,  0, matrixDimensions);
            int matrixIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
            beagle.setTransitionMatrix(matrixIndex, matrices, 1);

            branchLengths[eigenIndex][updateCount] = branchTime;
            branchUpdateCount[eigenIndex]++;

            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes to enforce post-order traversal
            Node child1 = node.getLeft();
            final int update1 = traverse(child1, flip);

            Node child2 = node.getRight();
            final int update2 = traverse(child2, flip);

            // If either child node was dirty, then update the parent node
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                int x = operationCount * Beagle.OPERATION_TUPLE_SIZE;

                // Flip places a flag to calculate these values later in beagle updatePartials()
                if (flip) {
                    partialBufferHelper.flipOffset(nodeNum);
                }

                final int[] operations = this.operations[0];
                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

                if (useScaleFactors) {
                    // get the index of this scaling buffer
                    int n = nodeNum - tipCount;

                    if (recomputeScaleFactors) {
                        // flip the indicator: can take either n or (internalNodeCount + 1) - n
                        scaleBufferHelper.flipOffset(n);

                        // store the index
                        scaleBufferIndices[n] = scaleBufferHelper.getOffsetIndex(n);

                        operations[x + 1] = scaleBufferIndices[n]; // Write new scaleFactor
                        operations[x + 2] = Beagle.NONE;

                    } else {
                        operations[x + 1] = Beagle.NONE;
                        operations[x + 2] = scaleBufferIndices[n]; // Read existing scaleFactor
                    }

                } else {
                    if (useAutoScaling) {
                        scaleBufferIndices[nodeNum - tipCount] = partialBufferHelper.getOffsetIndex(nodeNum);
                    }
                    operations[x + 1] = Beagle.NONE; // Not using scaleFactors
                    operations[x + 2] = Beagle.NONE;
                }

                // specify operations for beagle to perform later in updatePartials()
                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNr()); // source node 1
                operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNr()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNr()); // source node 2
                operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNr()); // source matrix 2

                operationCount++;

                update |= (update1 | update2);
            }
        }
        return update;
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

        hasDirt = Tree.IS_CLEAN;

        if (substitutionModel instanceof CalculationNode) {
            if (((CalculationNode) substitutionModel).isDirtyCalculation()) {
                hasDirt = Tree.IS_DIRTY;
                return true;
            }
        }
        
        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    /**
     * Stores the additional state other than model components
     */
    @Override
    public void store() {

        partialBufferHelper.storeState();
        eigenBufferHelper.storeState();
        matrixBufferHelper.storeState();

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
        }

        // Store origin partials
        System.arraycopy(originPartials, 0, storedOriginPartials, 0, originPartials.length);
        // Store root to origin branch transition matrix
        System.arraycopy(rootTransitionMatrix, 0, storedRootTransitionMatrix, 0, rootTransitionMatrix.length);

        // Store logP and reset isDirty to false
        super.store();

        // Store branch lengths
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);

        // Store ancestral states
        System.arraycopy(reconstructedStates, 0, storedReconstructedStates, 0, reconstructedStates.length);
        storedAreStatesRedrawn = areStatesRedrawn;
    }

    /**
     * Restores the state that was stored.
     */
    @Override
    public void restore() {
        
        partialBufferHelper.restoreState();
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.restoreState();
            System.arraycopy(storedScaleBufferIndices, 0, scaleBufferIndices, 0, storedScaleBufferIndices.length);
        }

        // restore origin partials
        System.arraycopy(storedOriginPartials, 0, originPartials, 0, originPartials.length);
        // Restore root to origin branch transition matrix
        System.arraycopy(storedRootTransitionMatrix, 0, rootTransitionMatrix, 0, rootTransitionMatrix.length);


        // Restore logP and reset isDirty to false
        logP = storedLogP;

        // Reset isDirty to false
        super.restore(); 

        // Restore branch lengths
        System.arraycopy(storedBranchLengths, 0, m_branchLengths, 0, m_branchLengths.length);

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

    // Sets the partials from a sequence in an alignment
    protected final void setStates(Beagle beagle, int nodeIndex, int taxon) {
        Alignment data = dataInput.get();
        int[] states = new int[1];
        states[0] = data.getDataType().getStatesForCode(data.getPattern(taxon, 0))[0];
        beagle.setTipStates(nodeIndex, states);
    }

    // Basic initialization of a beagle instance
    private void setupBeagle() {
        // Initialize buffer helpers
        partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);
        eigenBufferHelper = new BufferIndexHelper(1, 0);
        matrixBufferHelper = new BufferIndexHelper(nodeCount, 0);
        scaleBufferHelper = new BufferIndexHelper(internalNodeCount + 1, 0);

        // Load system properties
        preferredOrder = parseSystemPropertyIntegerArray("beagle.preferred.flags");

        // Initialize BEAGLE instance
        initializeBeagleInstance();

        // Setup tip states and pattern weights
        setupTipStatesAndWeights();
    }

    private void initializeBeagleInstance() {
        long preferenceFlags = 0;
        if (preferredOrder.size() > 0) {
            preferenceFlags = preferredOrder.get(0);
        }

        long requirementFlags = 0;
        requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();

        try {
            beagle = BeagleFactory.loadBeagleInstance(
                tipCount,
                partialBufferHelper.getBufferCount(),
                tipCount,
                stateCount,
                1,
                eigenBufferHelper.getBufferCount(),
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
    }

    private void setupTipStatesAndWeights() {
        Node[] nodes = treeInput.get().getNodesAsArray();
        for (int i = 0; i < tipCount; i++) {
            int taxon = dataInput.get().getTaxonIndex(nodes[i].getID());
            setStates(beagle, i, taxon);
        }

        double[] patternWeights = new double[1];
        patternWeights[0] = dataInput.get().getPatternWeight(0);
        beagle.setPatternWeights(patternWeights);
        beagle.setCategoryWeights(0, new double[]{1.0});
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

        public void storeState() {
            System.arraycopy(indexOffsets, 0, storedIndexOffsets, 0, indexOffsets.length);
        }

        public void restoreState() {
            int[] tmp = storedIndexOffsets;
            storedIndexOffsets = indexOffsets;
            indexOffsets = tmp;
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
            if (node.isRoot()) {
                // Handle a root node without origin
                beagle.getPartials(partialBufferHelper.getOffsetIndex(node.getNr()), Beagle.NONE, conditionalProbabilities);
                double[] rootFrequencies = substitutionModel.getFrequencies();
                for (int i = 0; i < stateCount; i++) {
                        conditionalProbabilities[i] *= rootFrequencies[i];
                    }
                reconstructedStates[node.getNr()] = Randomizer.randomChoicePDF(conditionalProbabilities);
            } else {
                // Handle an internal node
                double[] partialLikelihood = new double[stateCount];
                beagle.getPartials(partialBufferHelper.getOffsetIndex(node.getNr()), Beagle.NONE, partialLikelihood);

                if (node.isRoot()) {
                    // Handle root node (below origin) when origin is present since the origin is not physically in the tree
                    double[] oPs = new double[stateCount];
                    System.arraycopy(originPartials, 0, oPs, 0, originPartials.length);
                    double[] rootFrequencies = substitutionModel.getFrequencies();
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
            }

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