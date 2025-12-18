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
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.BeagleTreeLikelihood.PartialsRescalingScheme;
import beast.base.evolution.likelihood.GenericTreeLikelihood;

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.BeagleInfo;
import beagle.InstanceDetails;
import beagle.ResourceDetails;

/**
 * Uses the BEAGLE library to calculate the tree likelihood by modifying the BeagleTreeLikelihood
 * code to include the origin node branch and output a log likelihood identical to TideTree
 * for the case when editing occurs the entire duration of the experiment.
 *
 * @author Stephen Staklinski
 */
@Description("Uses the beagle library to calculate the tree likelihood by modifying the BeagleTreeLikelihood" +
            "code to include the origin node branch and output a log likelihood identical to TideTree" +
            "for the case when editing occurs the entire duration of the experiment.")
public class BeamBeagleTreeLikelihood extends GenericTreeLikelihood {

    /** Input for the origin parameter */
    public Input<RealParameter> originInput = new Input<>("origin", 
            "Start of the cell division process, usually start of the experiment.", 
            Input.Validate.OPTIONAL);

    /** Input for root frequencies */
    final public Input<Frequencies> rootFrequenciesInput = new Input<>("rootFrequencies", 
            "prior state frequencies at root, optional", 
            Input.Validate.OPTIONAL);

    /**
     * Enumeration defining different scaling options for numerical stability in BEAGLE calculations.
     * Each option represents a different strategy for handling scaling in likelihood computations.
     */
    public static enum Scaling {
        /** No scaling performed */
        none,
        
        /** Always perform scaling */
        always,
        
        /** Use default scaling behavior */
        _default
    }

    /** Input for scaling type */
    final public Input<Scaling> scaling = new Input<>("scaling", 
            "type of scaling to use, one of " + Arrays.toString(Scaling.values()) + 
            ". If not specified, the -beagle_scaling flag is used.", 
            Scaling._default, 
            Scaling.values());

    /** Origin parameter */
    protected RealParameter origin;

    /** Whether to use origin node */
    protected boolean useOrigin = false;

    /** BEAGLE instance counter */
    private static int instanceCount = 0;

    /** Resource order list */
    private static List<Integer> resourceOrder = null;

    /** Preferred flags list */
    private static List<Integer> preferredOrder = null;

    /** Required flags list */
    private static List<Integer> requiredOrder = null;

    /** Scaling order list */
    private static List<String> scalingOrder = null;

    /** Property name for resource order */
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";

    /** Property name for preferred flags */
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";

    /** Property name for required flags */
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";

    /** Threshold for scaling partials */
    private double scalingThreshold = 1.0E-100;

    /** Number of states in the model */
    int m_nStateCount;

    /** Number of nodes in the tree */
    int m_nNodeCount;

    /** Dimensions of the transition matrix */
    int matrixDimensions;

    /** Transition matrices */
    private double[] matrices;

    /** Current frequencies */
    private double[] currentFreqs;

    /** Number of eigen decompositions */
    private int eigenCount;

    /** Matrix update indices */
    private int[][] matrixUpdateIndices;

    /** Branch lengths */
    private double[][] branchLengths;

    /** Branch update counts */
    private int[] branchUpdateCount;

    /** Operations array */
    private int[][] operations;

    /** Operation list count */
    private int operationListCount;

    /** Operation counts */
    private int[] operationCount;

    /** Buffer index helper for partials */
    protected BufferIndexHelper partialBufferHelper;

    /** Buffer index helper for eigen values */
    private BufferIndexHelper eigenBufferHelper;

    /** Buffer index helper for matrices */
    protected BufferIndexHelper matrixBufferHelper;

    /** Buffer index helper for scaling */
    protected BufferIndexHelper scaleBufferHelper;

    /** Number of tip nodes */
    protected int tipCount;

    /** Number of internal nodes */
    protected int internalNodeCount;

    /** Number of patterns */
    protected int patternCount;

    /** Number of categories */
    protected int categoryCount;

    /** BEAGLE instance */
    protected Beagle beagle;

    /** Partial likelihoods for origin node */
    protected double[] originPartials;

    /** Stored partial likelihoods for origin node */
    protected double[] storedOriginPartials;

    /** Transition matrix for root branch */
    protected double[] rootTransitionMatrix;

    /** Stored transition matrix for root branch */
    protected double[] storedRootTransitionMatrix;

    /** Branch lengths */
    protected double[] m_branchLengths;

    /** Stored branch lengths */
    protected double[] storedBranchLengths;

    /** Dirt flag for tree updates */
    protected int hasDirt;

    /** Substitution model */
    protected SubstitutionModel substitutionModel;

    /** Site model */
    protected SiteModel.Base m_siteModel;

    /** Branch rate model */
    protected BranchRateModel.Base branchRateModel;

    /** Pattern log likelihoods */
    protected double[] patternLogLikelihoods;

    /** Probability array */
    protected double[] probabilities;

    /** Property name for scaling */
    private static final String SCALING_PROPERTY = "beagle.scaling";

    /** Property name for rescale frequency */
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";

    /** Default rescaling scheme */
    private static final PartialsRescalingScheme DEFAULT_RESCALING_SCHEME = 
            PartialsRescalingScheme.DYNAMIC;

    /** Current rescaling scheme */
    private PartialsRescalingScheme rescalingScheme = DEFAULT_RESCALING_SCHEME;

    /** Default rescale frequency */
    private static final int RESCALE_FREQUENCY = 10000;

    /** Current rescale frequency */
    private int rescalingFrequency = RESCALE_FREQUENCY;

    /** Number of times to rescale */
    private static final int RESCALE_TIMES = 1;

    /** Whether to use scale factors */
    protected boolean useScaleFactors = false;

    /** Whether to use auto scaling */
    private boolean useAutoScaling = false;

    /** Whether to recompute scale factors */
    private boolean recomputeScaleFactors = false;

    /** Whether underflow has occurred */
    private boolean everUnderflowed = false;

    /** Rescaling counter */
    private int rescalingCount = 0;

    /** Inner rescaling counter */
    private int rescalingCountInner = 0;

    /** Scale buffer indices */
    private int[] scaleBufferIndices;

    /** Stored scale buffer indices */
    private int[] storedScaleBufferIndices;

    /** Debug flags */
    private boolean storeRestoreDebug = false;
    private boolean transitionMatrixDebug = false;
    private boolean partialsDebug = false;
    private boolean debugInputData = false;

    @Override
    public void initAndValidate() {
        initialize();
    }

    /**
     * Initializes the BEAGLE instance and sets up necessary data structures.
     */
    private void initialize() {
        setupOrigin();
        setupModels();
        setupStateCounts();
        setupArrays();
        setupBeagle();
    }

    private void setupOrigin() {
        if (originInput.get() != null) {
            origin = originInput.get();
            useOrigin = true;
        }
    }

    private void setupModels() {
        // Validate and setup site model
        validateAndSetupSiteModel();
        
        // Validate and setup substitution model
        validateAndSetupSubstitutionModel();
        
        // Validate and setup branch rate model
        validateAndSetupBranchRateModel();
    }

    private void validateAndSetupSiteModel() {
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        
        // Validate site categories
        if (m_siteModel.getCategoryRates(null).length != 1) {
            throw new IllegalArgumentException("Site categories are not supported in the current implementation.");
        }
    }

    private void validateAndSetupSubstitutionModel() {
        substitutionModel = (SubstitutionModel.Base) m_siteModel.substModelInput.get();
        if (!substitutionModel.canReturnComplexDiagonalization()) {
            throw new IllegalArgumentException("Substitution model must be able to return transition probabilities.");
        }
    }

    private void validateAndSetupBranchRateModel() {
        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) {
            throw new IllegalArgumentException("Branch rate model must be specified in the current implementation.");
        }
    }

    private void setupStateCounts() {
        // Get state count and pattern count
        m_nStateCount = substitutionModel.getStateCount();
        patternCount = dataInput.get().getPatternCount();
        
        if (debugInputData) {
            System.out.println("There are " + m_nStateCount + " unique states.");
            System.out.println("There are " + patternCount + " unique site patterns.");
        }
        Log.warning.println("There are " + patternCount + " unique site patterns.");
    }

    private void setupArrays() {
        // Initialize frequency array
        currentFreqs = new double[m_nStateCount];
        
        // Setup category count and matrix dimensions
        categoryCount = m_siteModel.getCategoryCount();
        matrixDimensions = m_nStateCount * m_nStateCount;
        
        // Initialize probability arrays
        probabilities = new double[matrixDimensions];
        matrices = new double[matrixDimensions * categoryCount];
        
        // Setup tree node counts
        m_nNodeCount = treeInput.get().getNodeCount();
        tipCount = treeInput.get().getLeafNodeCount();
        internalNodeCount = m_nNodeCount - tipCount;
        
        // Initialize branch length arrays
        m_branchLengths = new double[m_nNodeCount];
        storedBranchLengths = new double[m_nNodeCount];
        
        // Set eigen count
        eigenCount = 1;
        
        // Setup origin-specific arrays if needed
        if (useOrigin) {
            setupOriginArrays();
        }
    }

    private void setupOriginArrays() {
        int partialsSize = patternCount * m_nStateCount;
        originPartials = new double[partialsSize];
        storedOriginPartials = new double[partialsSize];
        rootTransitionMatrix = new double[matrixDimensions];
        storedRootTransitionMatrix = new double[matrixDimensions];
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood
     */
    @Override
    public double calculateLogP() {
        if (debugEnabled()) {
            System.out.println("calculateLogP() called");
        }

        // Validate origin if used
        if (useOrigin && !validateOrigin()) {
            return Double.NEGATIVE_INFINITY;
        }

        // Initialize arrays if needed
        initializeArrays();

        // Setup rescaling
        setupRescaling();

        // Reset counters
        resetCounters();

        // Calculate likelihood
        return calculateLikelihood();
    }

    private boolean debugEnabled() {
        return transitionMatrixDebug || partialsDebug || storeRestoreDebug;
    }

    private boolean validateOrigin() {
        Double originHeight = origin.getValue();
        if (treeInput.get().getRoot().getHeight() >= originHeight) {
            if (debugEnabled()) {
                System.out.println("Tree height larger than origin.");
            }
            return false;
        }
        return true;
    }

    private void initializeArrays() {
        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        if (matrixUpdateIndices == null) {
            matrixUpdateIndices = new int[eigenCount][m_nNodeCount];
            branchLengths = new double[eigenCount][m_nNodeCount];
            branchUpdateCount = new int[eigenCount];
            scaleBufferIndices = new int[internalNodeCount];
            storedScaleBufferIndices = new int[internalNodeCount];
        }

        if (operations == null) {
            operations = new int[1][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
            operationCount = new int[1];
        }
    }

    private void setupRescaling() {
        recomputeScaleFactors = false;

        switch (rescalingScheme) {
            case ALWAYS:
                useScaleFactors = true;
                recomputeScaleFactors = true;
                break;
            case DYNAMIC:
                if (everUnderflowed) {
                    useScaleFactors = true;
                    if (rescalingCountInner < RESCALE_TIMES) {
                        recomputeScaleFactors = true;
                        hasDirt = Tree.IS_FILTHY;
                    }
                    updateRescalingCounters();
                }
                break;
            case DELAYED:
                if (everUnderflowed) {
                    useScaleFactors = true;
                    recomputeScaleFactors = true;
                    hasDirt = Tree.IS_FILTHY;
                    rescalingCount++;
                }
                break;
            default:
                break;
        }
    }

    private void updateRescalingCounters() {
        rescalingCountInner++;
        rescalingCount++;
        if (rescalingCount > RESCALE_FREQUENCY) {
            rescalingCount = 0;
            rescalingCountInner = 0;
        }
    }

    private void resetCounters() {
        for (int i = 0; i < eigenCount; i++) {
            branchUpdateCount[i] = 0;
        }
        operationListCount = 0;
        operationCount[0] = 0;
    }

    private double calculateLikelihood() {
        final Node root = treeInput.get().getRoot();
        traverse(root, true);

        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {
            logL = calculateSingleLikelihood(root);
            done = handleUnderflow(logL, root, firstRescaleAttempt);
            firstRescaleAttempt = false;
        } while (!done);

        logP = logL;
        logDebugInfo(logL);
        return logL;
    }

    private double calculateSingleLikelihood(Node root) {
        // Update partials in BEAGLE
        beagle.updatePartials(operations[0], operationCount[0], Beagle.NONE);
        int rootIndex = partialBufferHelper.getOffsetIndex(root.getNr());

        // Setup frequencies
        setupFrequencies();

        // Calculate likelihood based on origin
        return useOrigin && root.getHeight() != origin.getValue() ? 
                calculateLikelihoodWithOrigin(root, rootIndex) : 
                calculateStandardLikelihood(rootIndex);
    }

    private void setupFrequencies() {
        double[] frequencies = rootFrequenciesInput.get() == null ? 
                substitutionModel.getFrequencies() : 
                rootFrequenciesInput.get().getFreqs();

        int cumulateScaleBufferIndex = setupScaling();

        if (frequencies != currentFreqs) {
            beagle.setStateFrequencies(0, frequencies);
        }
        System.arraycopy(frequencies, 0, currentFreqs, 0, frequencies.length);
    }

    private int setupScaling() {
        int cumulateScaleBufferIndex = Beagle.NONE;
        if (useScaleFactors) {
            if (recomputeScaleFactors) {
                scaleBufferHelper.flipOffset(internalNodeCount);
                cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                beagle.resetScaleFactors(cumulateScaleBufferIndex);
                beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, 
                        cumulateScaleBufferIndex);
            } else {
                cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
            }
        }
        return cumulateScaleBufferIndex;
    }

    private double calculateLikelihoodWithOrigin(Node root, int rootIndex) {
        if (partialsDebug) {
            System.out.println("Using origin, so starting root to origin partials calculation.");
        }

        // Get root partials and transition matrix
        double[] rootPartials = getRootPartials(rootIndex);
        setupRootTransitionMatrix(root);

        // Calculate origin partials
        calculateOriginPartials(rootPartials, rootTransitionMatrix, originPartials);

        // Scale origin partials if needed
        double originScaleFactorsSum = scaleOriginPartials();

        // Calculate final likelihood
        return calculateFinalLikelihood(root, rootIndex, rootPartials, originScaleFactorsSum);
    }

    private double[] getRootPartials(int rootIndex) {
        double[] rootPartials = new double[patternCount * m_nStateCount];
        beagle.getPartials(rootIndex, Beagle.NONE, rootPartials);
        return rootPartials;
    }

    private void setupRootTransitionMatrix(Node root) {
        int rootNodeNum = root.getNr();
        double br = branchRateModel.getRateForBranch(root);

        if (transitionMatrixDebug) {
            logTransitionMatrixDebug(root, br);
        }

        substitutionModel.getTransitionProbabilities(root, origin.getValue(), root.getHeight(), br, probabilities);
        System.arraycopy(probabilities, 0, rootTransitionMatrix, 0, matrixDimensions);

        if (transitionMatrixDebug) {
            logTransitionMatrixDetails();
        }
    }

    private void logTransitionMatrixDebug(Node root, double br) {
        double len = origin.getValue() - root.getHeight();
        System.out.println("Root to origin length: " + len);
        System.out.println("Root to origin clockRate: " + br);
        double dist = len * br;
        System.out.println("Root to origin distance (length * clockRate): " + dist);
    }

    private void logTransitionMatrixDetails() {
        System.out.println();
        System.out.println("Root to origin transition matrix returned to likelihood:");
        for (int i = 0; i < m_nStateCount; i++) {
            System.out.println(Arrays.toString(Arrays.copyOfRange(probabilities, 
                    i * m_nStateCount, (i + 1) * m_nStateCount)));
        }
    }

    private double scaleOriginPartials() {
        if (!useScaleFactors) {
            return 0.0;
        }

        double[] originScaleFactors = new double[patternCount];
        int u = 0;
        for (int i = 0; i < patternCount; i++) {
            double scaleFactor = calculateScaleFactor(i, u);
            if (scaleFactor < scalingThreshold) {
                applyScaleFactor(i, u, scaleFactor);
                originScaleFactors[i] = Math.log(scaleFactor);
            } else {
                originScaleFactors[i] = 0.0;
            }
            u += m_nStateCount;
        }

        return Arrays.stream(originScaleFactors).sum();
    }

    private double calculateScaleFactor(int pattern, int startIndex) {
        double scaleFactor = 0.0;
        for (int j = 0; j < m_nStateCount; j++) {
            scaleFactor = Math.max(scaleFactor, originPartials[startIndex + j]);
        }
        return scaleFactor;
    }

    private void applyScaleFactor(int pattern, int startIndex, double scaleFactor) {
        for (int j = 0; j < m_nStateCount; j++) {
            originPartials[startIndex + j] /= scaleFactor;
        }
    }

    private double calculateFinalLikelihood(Node root, int rootIndex, double[] rootPartials, 
            double originScaleFactorsSum) {
        int rootNodeNum = root.getNr();
        int cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);

        // Replace root partials with origin partials
        beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), originPartials);

        // Calculate likelihood
        double[] sumLogLikelihoods = new double[1];
        beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, 
                new int[]{0}, new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

        // Restore original root partials
        beagle.setPartials(partialBufferHelper.getOffsetIndex(rootNodeNum), rootPartials);

        return sumLogLikelihoods[0] + originScaleFactorsSum;
    }

    private double calculateStandardLikelihood(int rootIndex) {
        double[] sumLogLikelihoods = new double[1];
        int cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
        beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, 
                new int[]{0}, new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);
        return sumLogLikelihoods[0];
    }

    private boolean handleUnderflow(double logL, Node root, boolean firstRescaleAttempt) {
        if (Double.isNaN(logL) || Double.isInfinite(logL)) {
            everUnderflowed = true;
            logL = Double.NEGATIVE_INFINITY;

            if (firstRescaleAttempt && 
                (rescalingScheme == PartialsRescalingScheme.DYNAMIC || 
                 rescalingScheme == PartialsRescalingScheme.DELAYED)) {
                return handleRescalingAttempt(root);
            }
            return true;
        }
        return true;
    }

    private boolean handleRescalingAttempt(Node root) {
        useScaleFactors = true;
        recomputeScaleFactors = true;

        for (int i = 0; i < eigenCount; i++) {
            branchUpdateCount[i] = 0;
        }

        operationCount[0] = 0;
        traverse(root, false);
        return false;
    }

    private void logDebugInfo(double logL) {
        if (partialsDebug) {
            System.out.println("Returning logL: " + logL);
            if (Double.isNaN(logL) || Double.isInfinite(logL)) {
                System.out.println("Likelihood calculation has underflow (is NaN or -Infinity) " +
                        "and partials scaling is not used in the current implementation, " +
                        "so returning -Infinity.");
            }
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

        if (transitionMatrixDebug || partialsDebug) {
            System.out.println();
            System.out.println();
            System.out.println("Traversing node: " + node.getNr());
        }

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

            if (transitionMatrixDebug) {
                System.out.println("Updating transition matrix for node: " + nodeNum);
                double len = node.getParent().getHeight() - node.getHeight();
                System.out.println("Node: " + node.getNr() + " real time length: " + len);
                System.out.println("Node: " + node.getNr() + " Clock Rate: " + branchRate);
                double dist = len * branchRate;
                System.out.println("Node: " + node.getNr() + " Branch dist (clock rate * length): " + branchTime);
            }

            // Get the new transition probability matrix and store it in beagle
            substitutionModel.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(), branchRate, probabilities);
            System.arraycopy(probabilities, 0, matrices,  0, matrixDimensions);
            int matrixIndex = matrixBufferHelper.getOffsetIndex(nodeNum);
            beagle.setTransitionMatrix(matrixIndex, matrices, 1);

            if (transitionMatrixDebug) {
                System.out.println();
                System.out.println("Node: " + node.getNr() + " transition matrix returned to likelihood:");
                for (int i = 0; i < m_nStateCount; i++) {
                    System.out.println(Arrays.toString(Arrays.copyOfRange(probabilities, i * m_nStateCount, (i + 1) * m_nStateCount)));
                }
            }

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

                int x = operationCount[operationListCount] * Beagle.OPERATION_TUPLE_SIZE;

                // Flip places a flag to calculate these values later in beagle updatePartials()
                if (flip) {
                    partialBufferHelper.flipOffset(nodeNum);
                }

                final int[] operations = this.operations[operationListCount];
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

                if (partialsDebug) {
                    System.out.println("Updating partials for node: " + nodeNum);
                    System.out.println("Operations: " + Arrays.toString(Arrays.copyOfRange(operations, x, x + Beagle.OPERATION_TUPLE_SIZE)));
                }

                operationCount[operationListCount]++;

                update |= (update1 | update2);
            }
        }
        return update;
    }

    /*
     * Calculate partials for the origin node of degree 1 that goes from the start of the experiment to the root
     */
    protected double[] calculateOriginPartials(double[] partials1, double[] matrices1, double[] partials3) {
        double sum1;
        int u = 0;
        int v = 0;

        for (int k = 0; k < patternCount; k++) {
            int w = 0;
            for (int i = 0; i < m_nStateCount; i++) {
                sum1 = 0.0;
                for (int j = 0; j < m_nStateCount; j++) {
                    sum1 += matrices1[w] * partials1[v + j];
                    w++;
                }
                partials3[u] = sum1;
                u++;
            }
            v += m_nStateCount;
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
        if (m_siteModel.isDirtyCalculation()) {
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

        if (useScaleFactors || useAutoScaling) { // Only store when actually used
            scaleBufferHelper.storeState();
            System.arraycopy(scaleBufferIndices, 0, storedScaleBufferIndices, 0, scaleBufferIndices.length);
//            storedRescalingCount = rescalingCount;
        }

        if (useOrigin) {
            // store origin partials
            System.arraycopy(originPartials, 0, storedOriginPartials, 0, originPartials.length);

            // Store root to origin branch transition matrix
            System.arraycopy(rootTransitionMatrix, 0, storedRootTransitionMatrix, 0, rootTransitionMatrix.length);

        }
        // Store logP and reset isDirty to false
        super.store();

        // store branch lengths
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);

        if (storeRestoreDebug) {
            System.out.println("Original logP: " + logP);
            System.out.println("Original branch lengths: " + Arrays.toString(m_branchLengths));
            double[] rootPartials = new double[patternCount * m_nStateCount];
            beagle.getPartials(partialBufferHelper.getOffsetIndex(treeInput.get().getRoot().getNr()), Beagle.NONE, rootPartials);
            System.out.println("Original root partials: " + Arrays.toString(rootPartials));
            System.out.println("Original origin partials: " + Arrays.toString(originPartials));
            System.out.println("Original root matrix: " + Arrays.toString(rootTransitionMatrix));
        }
    }

    /**
     * Restores the state that was stored.
     */
    @Override
    public void restore() {

        if (storeRestoreDebug) {
            System.out.println("New logP: " + logP);
            System.out.println("New branch lengths: " + Arrays.toString(m_branchLengths));
            double[] rootPartials = new double[patternCount * m_nStateCount];
            beagle.getPartials(partialBufferHelper.getOffsetIndex(treeInput.get().getRoot().getNr()), Beagle.NONE, rootPartials);
            System.out.println("New root partials: " + Arrays.toString(rootPartials));
            System.out.println("New origin partials: " + Arrays.toString(originPartials));
            System.out.println("New root matrix: " + Arrays.toString(rootTransitionMatrix));
        }
        
        partialBufferHelper.restoreState();
        eigenBufferHelper.restoreState();
        matrixBufferHelper.restoreState();

        if (useScaleFactors || useAutoScaling) {
            scaleBufferHelper.restoreState();
            int[] tmp2 = storedScaleBufferIndices;
            storedScaleBufferIndices = scaleBufferIndices;
            scaleBufferIndices = tmp2;
//            rescalingCount = storedRescalingCount;
        }

        if (useOrigin) {
            // restore origin partials
            double[] tmp3 = storedOriginPartials;
            storedOriginPartials = originPartials;
            originPartials = tmp3;

            // Restore root to origin branch transition matrix
            System.arraycopy(storedRootTransitionMatrix, 0, rootTransitionMatrix, 0, rootTransitionMatrix.length);
        }


        // Restore logP and reset isDirty to false
        logP = storedLogP;

        // Reset isDirty to false
        super.restore(); 

        // restore branch lengths
        double[] tmp = storedBranchLengths;
        storedBranchLengths = m_branchLengths;
        m_branchLengths = tmp;

        if (storeRestoreDebug) {
            System.out.println("Restored logP: " + logP);
            System.out.println("Restored branch lengths: " + Arrays.toString(m_branchLengths));
            double[] rootPartials2 = new double[patternCount * m_nStateCount];
            beagle.getPartials(partialBufferHelper.getOffsetIndex(treeInput.get().getRoot().getNr()), Beagle.NONE, rootPartials2);
            System.out.println("Restored root partials: " + Arrays.toString(rootPartials2));
            System.out.println("Restored origin partials: " + Arrays.toString(originPartials));
            System.out.println("Restored root matrix: " + Arrays.toString(rootTransitionMatrix));
        }

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


    private static List<String> parseSystemPropertyStringArray(String propertyName) {

        List<String> order = new ArrayList<>();

        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    String s = part.trim();
                    order.add(s);
                } catch (NumberFormatException nfe) {
                	Log.warning.println("Invalid getEigenDecompositionentry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }
    
    
    protected int getScaleBufferCount() {
        return internalNodeCount + 1;
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param nodeIndex     nodeIndex
     * @param taxon the taxon
     */
    protected final void setPartials(Beagle beagle,
                                     int nodeIndex, int taxon) {
        Alignment data = dataInput.get();

        double[] partials = new double[patternCount * m_nStateCount * categoryCount];

        int v = 0;
        for (int i = 0; i < patternCount; i++) {

        	double[] tipProbabilities = data.getTipLikelihoods(taxon,i);
            if (tipProbabilities != null) {
            	for (int state = 0; state < m_nStateCount; state++) {
            		partials[v++] = tipProbabilities[state];
            	}
            }
            else {
            	int stateCount = data.getPattern(taxon, i);
                boolean[] stateSet = data.getStateSet(stateCount);
                for (int state = 0; state < m_nStateCount; state++) {
                	 partials[v++] = (stateSet[state] ? 1.0 : 0.0);                
                }
            }
        }

        // if there is more than one category then replicate the partials for each
        int n = patternCount * m_nStateCount;
        int k = n;
        for (int i = 1; i < categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }

        beagle.setPartials(nodeIndex, partials);
    }

    public int getPatternCount() {
        return patternCount;
    }


    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param nodeIndex     nodeIndex
     * @param taxon         the taxon
     */
    protected final void setStates(Beagle beagle,
                                   int nodeIndex, int taxon) {
        Alignment data = dataInput.get();
        int i;

        int[] states = new int[patternCount];

        if (debugInputData) {
            System.out.println("Number of site patterns to set: " + patternCount);
        }

        for (i = 0; i < patternCount; i++) {
            int code = data.getPattern(taxon, i);
            int[] statesForCode = data.getDataType().getStatesForCode(code);
            if (statesForCode.length != 1) {
                System.err.println("Error: Invalid state code for taxon " + taxon + " at pattern " + i + ".");
                System.exit(1);
            }
            
            states[i] = statesForCode[0];
        
            if (debugInputData) {
                System.out.println("For pattern " + i + " the code is " + code + " and the stateForCode is " + Arrays.toString(statesForCode) + " so the state is " + states[i]);
            }

            }

        if (debugInputData) {
            System.out.println("Setting states for node: " + nodeIndex + " with states: " + Arrays.toString(states));
        }

        beagle.setTipStates(nodeIndex, states);
    }


    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {    	
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
        	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
            	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
	}


    /*
     * Basic initialization of a beagle instance taken from BeagleTreeLikelihood
     */
    private void setupBeagle() {
        // Initialize buffer helpers
        partialBufferHelper = new BufferIndexHelper(m_nNodeCount, tipCount);
        eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);
        matrixBufferHelper = new BufferIndexHelper(m_nNodeCount, 0);
        scaleBufferHelper = new BufferIndexHelper(internalNodeCount + 1, 0);

        // Load configuration from system properties
        loadSystemProperties();

        // Initialize BEAGLE instance
        initializeBeagleInstance();

        // Setup tip states and pattern weights
        setupTipStatesAndWeights();
    }

    private void loadSystemProperties() {
        if (resourceOrder == null) {
            resourceOrder = parseSystemPropertyIntegerArray(RESOURCE_ORDER_PROPERTY);
        }
        if (preferredOrder == null) {
            preferredOrder = parseSystemPropertyIntegerArray(PREFERRED_FLAGS_PROPERTY);
        }
        if (requiredOrder == null) {
            requiredOrder = parseSystemPropertyIntegerArray(REQUIRED_FLAGS_PROPERTY);
        }
        if (scalingOrder == null) {
            scalingOrder = parseSystemPropertyStringArray(SCALING_PROPERTY);
        }

        // Configure rescaling scheme
        configureRescalingScheme();
    }

    private void configureRescalingScheme() {
        rescalingScheme = DEFAULT_RESCALING_SCHEME;
        
        if (scalingOrder.size() > 0) {
            rescalingScheme = PartialsRescalingScheme.parseFromString(
                    scalingOrder.get(instanceCount % scalingOrder.size()));
        }

        if (scaling.get().equals(Scaling.always)) {
            rescalingScheme = PartialsRescalingScheme.ALWAYS;
        } else if (scaling.get().equals(Scaling.none)) {
            rescalingScheme = PartialsRescalingScheme.NONE;
        }

        // Set default behavior based on resource type
        if (rescalingScheme == PartialsRescalingScheme.DEFAULT) {
            int[] resourceList = getResourceList();
            rescalingScheme = (resourceList != null && resourceList[0] > 1) ? 
                    PartialsRescalingScheme.NONE : PartialsRescalingScheme.DYNAMIC;
        }
    }

    private int[] getResourceList() {
        if (resourceOrder.size() > 0) {
            return new int[]{resourceOrder.get(instanceCount % resourceOrder.size()), 0};
        }
        return null;
    }

    private void initializeBeagleInstance() {
        long preferenceFlags = getPreferenceFlags();
        long requirementFlags = getRequirementFlags();
        int[] resourceList = getResourceList();

        try {
            beagle = BeagleFactory.loadBeagleInstance(
                tipCount,
                partialBufferHelper.getBufferCount(),
                tipCount,
                m_nStateCount,
                patternCount,
                eigenBufferHelper.getBufferCount(),
                matrixBufferHelper.getBufferCount(),
                categoryCount,
                scaleBufferHelper.getBufferCount(),
                resourceList,
                preferenceFlags,
                requirementFlags
            );
        } catch (Exception e) {
            Log.warning.println("Error setting up BEAGLE. Check install.");
            System.exit(1);
        }

        logBeagleDetails();
    }

    private long getPreferenceFlags() {
        long flags = 0;
        int[] resourceList = getResourceList();
        
        if (resourceList != null && resourceList[0] > 0) {
            flags |= BeagleFlag.PROCESSOR_GPU.getMask();
        }
        
        if (preferredOrder.size() > 0) {
            flags = preferredOrder.get(instanceCount % preferredOrder.size());
        }
        
        if (rescalingScheme == PartialsRescalingScheme.AUTO) {
            flags |= BeagleFlag.SCALING_AUTO.getMask();
            useAutoScaling = true;
        }
        
        if (flags == 0 && resourceList == null && m_nStateCount == 4 && patternCount < 10000) {
            flags |= BeagleFlag.PROCESSOR_CPU.getMask();
        }
        
        return flags;
    }

    private long getRequirementFlags() {
        long flags = 0;
        if (requiredOrder.size() > 0) {
            flags = requiredOrder.get(instanceCount % requiredOrder.size());
        }
        flags |= BeagleFlag.EIGEN_COMPLEX.getMask();
        return flags;
    }

    private void logBeagleDetails() {
        InstanceDetails instanceDetails = beagle.getDetails();
        if (instanceDetails == null) {
            Log.warning.println("No external BEAGLE resources available");
            System.exit(1);
        }

        ResourceDetails resourceDetails = BeagleFactory.getResourceDetails(instanceDetails.getResourceNumber());
        if (resourceDetails == null) {
            Log.warning.println("Error retrieving BEAGLE resource for instance: " + instanceDetails);
            System.exit(1);
        }

        StringBuilder sb = new StringBuilder("Using BEAGLE version: " + BeagleInfo.getVersion() + " resource ");
        sb.append(resourceDetails.getNumber()).append(": ").append(resourceDetails.getName()).append("\n");
        
        if (resourceDetails.getDescription() != null) {
            Arrays.stream(resourceDetails.getDescription().split("\\|"))
                  .map(String::trim)
                  .filter(s -> !s.isEmpty())
                  .forEach(s -> sb.append("    ").append(s).append("\n"));
        }
        
        sb.append("    with instance flags: ").append(instanceDetails);
        Log.info.println(sb.toString());
        Log.warning.println("With " + patternCount + " unique site patterns.");
    }

    private void setupTipStatesAndWeights() {
        Node[] nodes = treeInput.get().getNodesAsArray();
        for (int i = 0; i < tipCount; i++) {
            int taxon = getTaxonIndex(nodes[i].getID(), dataInput.get());
            if (debugInputData) {
                System.out.println("Setting states for taxon: " + nodes[i].getID());
            }
            setStates(beagle, i, taxon);
        }

        double[] patternWeights = new double[patternCount];
        for (int i = 0; i < patternCount; i++) {
            patternWeights[i] = dataInput.get().getPatternWeight(i);
        }
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

    /**
     * Enumeration defining different schemes for rescaling partial likelihoods in BEAGLE.
     */
    public enum PartialsRescalingScheme {
        DEFAULT("default"),
        NONE("none"),
        DYNAMIC("dynamic"),
        ALWAYS("always"),
        DELAYED("delayed"),
        AUTO("auto");

        private final String text;

        PartialsRescalingScheme(String text) {
            this.text = text;
        }

        public String getText() {
            return text;
        }

        public static PartialsRescalingScheme parseFromString(String text) {
            for (PartialsRescalingScheme scheme : values()) {
                if (scheme.getText().compareToIgnoreCase(text) == 0) {
                    return scheme;
                }
            }
            return DEFAULT;
        }
    }
}