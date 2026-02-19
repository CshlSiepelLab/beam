package beam.likelihood;

import beast.base.core.Description;
import beast.base.evolution.likelihood.LikelihoodCore;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.List;

/**
 * Contains methods to calculate the partial likelihoods by using a simplified pruning
 * algorithm to save on computations given the irreversible assumptions of the substitution model.
 *
 * @author Stephen Staklinski
 */
@Description("Contains methods to calculate the partial likelihoods by using a simplified pruning " +
            "algorithm to save on computations given the irreversible assumptions of the substitution model.")
public class IrreversibleLikelihoodCore extends LikelihoodCore {

    /** Number of states in the model */
    protected int[] nrOfStatesPerSite;
    /** Number of nodes in the tree */
    protected int nrOfNodes;
    /** Number of nodes excluding the origin node */
    protected int numNodesNoOrigin;
    /** Number of sites in the alignment */
    protected int nrOfSites;
    /** Size of the partials array per site */
    protected int[] partialsSizes;
    /** Size of the transition matrix per site */
    protected int[] matrixSizes;
    /** Partial likelihoods for each node */
    protected double[][][][] partials;
    /** Transition matrices for each node */
    protected double[][][][] matrices;
    /** Current matrix indices for each node */
    protected int[] currentMatrixIndex;
    /** Stored matrix indices for each node */
    protected int[] storedMatrixIndex;
    /** Current partials indices for each node */
    protected int[] currentPartialsIndex;
    /** Stored partials indices for each node */
    protected int[] storedPartialsIndex;
    /** Ancestral states for each node and site */
    protected int[][] ancestralStates;
    /** Stored ancestral states for each node and site */
    protected int[][] storedAncestralStates;
    /** All possible states per site*/
    private int[][] allStatesPerSite;
    /** State representing unedited sequence */
    private static final int[] UNEDITED_STATE = new int[]{0};
    /** Whether to use scaling for numerical stability */
    protected boolean useScaling = false;
    /** Scaling factors for each node and site */
    protected double[][][] scalingFactors;
    /** Threshold for scaling partials */
    private static final double SCALING_THRESHOLD = 1.0E-100;
    /** State representing missing data */
    private int[] missingDataStatePerSite;
    /** Whether the substitution model is global */
    private boolean isGlobalModel;

    /**
     * Constructs a new IrreversibleLikelihoodCore with the specified parameters.
     */
    public IrreversibleLikelihoodCore(int nodeCount, int[] nrOfStates, int siteCount, int[] missingData, boolean isGlobal) {
        isGlobalModel = isGlobal;
        nrOfStatesPerSite = nrOfStates;
        nrOfNodes = nodeCount;
        numNodesNoOrigin = nrOfNodes - 1;
        nrOfSites = siteCount;
        missingDataStatePerSite = missingData;

        scalingFactors = new double[2][nodeCount][nrOfSites];

        matrixSizes = new int[nrOfSites];
        partialsSizes = new int[nrOfSites];
        for (int i = 0; i < siteCount; i++) {
            matrixSizes[i] = nrOfStatesPerSite[i] * nrOfStatesPerSite[i];
            partialsSizes[i] = nrOfStatesPerSite[i];
        }
        partials = new double[2][nodeCount][nrOfSites][];
        matrices = new double[2][nodeCount][nrOfSites][];
        // Initialize all arrays to 0
        for (int i = 0; i < nodeCount; i++) {
            for (int j = 0; j < nrOfSites; j++) {
                partials[0][i][j] = new double[partialsSizes[j]];
                partials[1][i][j] = new double[partialsSizes[j]];
                matrices[0][i][j] = new double[matrixSizes[j]];
                matrices[1][i][j] = new double[matrixSizes[j]];
                Arrays.fill(partials[0][i][j], 0.0);
                Arrays.fill(partials[1][i][j], 0.0);
                Arrays.fill(matrices[0][i][j], 0.0);
                Arrays.fill(matrices[1][i][j], 0.0);
                scalingFactors[0][i][j] = 0.0;
                scalingFactors[1][i][j] = 0.0;
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

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];
        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

        ancestralStates = new int[numNodesNoOrigin][nrOfSites];
        storedAncestralStates = new int[numNodesNoOrigin][nrOfSites];
        for (int i = 0; i < numNodesNoOrigin; i++) {
            for (int j = 0; j < nrOfSites; j++) {
                ancestralStates[i][j] = -1;
                storedAncestralStates[i][j] = -1;
            }
        }
    }

    /**
     * Initializes the partials at a node with the known states.
     *
     * @param leafIndex Index of the leaf node
     * @param states List of states for each site
     */
    public void setNodePartials(int leafIndex, List<Integer> states) {
        for (int i = 0; i < nrOfSites; i++) {
            // Set the known state to 1.0 for tips
            partials[currentPartialsIndex[leafIndex]][leafIndex][i][states.get(i)] = 1.0;
            // Set the ancestral state for tips
            ancestralStates[leafIndex][i] = (states.get(i) == missingDataStatePerSite[i]) ? -1 : states.get(i);
        }
    }

    /**
     * Sets possible ancestral states for a parent node based on its children's states.
     *
     * @param childIndex1 Index of the first child node
     * @param childIndex2 Index of the second child node
     * @param parentIndex Index of the parent node
     */
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
            // -1 initialization is the fallback here
        }
    }

    /**
     * Returns possible states for a given state.
     *
     * @param state The current state
     * @param siteNum The site number
     * @return Array of possible states
     */
    private int[] getPossibleStates(int state, int siteNum) {
        return state == 0 ? UNEDITED_STATE :
               state > 0 ? new int[]{0, state} :
               allStatesPerSite[siteNum];
    }

    /**
     * Calculates partial likelihoods at a node while first checking which partials
     * need to be calculated to save on computations when 0 values should just be propagated.
     *
     * @param childIndex1 Index of the first child node
     * @param childIndex2 Index of the second child node
     * @param parentIndex Index of the parent node
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

    /**
     * Calculates partial likelihoods and site log likelihoods at the cell division origin node.
     * Since this is the start of the experiment, the origin is known to be in the unedited state,
     * so we can only calculate the partials for that state and set the other to 0 since they
     * will not be used by the frequencies anyways.
     *
     * @param rootIndex Index of the root node
     * @param originIndex Index of the origin node
     * @return The log likelihood
     */
    public double calculateLogLikelihoods(int rootIndex, int originIndex) {
        double logP = 0.0;
        currentPartialsIndex[originIndex] = 1 - currentPartialsIndex[originIndex];

        final double[][] partials1 = partials[currentPartialsIndex[rootIndex]][rootIndex];
        final double[][] matrices1 = matrices[currentMatrixIndex[rootIndex]][rootIndex];
        final double[][] partials3 = partials[currentPartialsIndex[originIndex]][originIndex];

        for (int k = 0; k < nrOfSites; k++) {
            final int[] possibleStates = getPossibleStates(ancestralStates[rootIndex][k], k);
            double sum1 = 0.0;

            for (int j : possibleStates) {
                sum1 += matrices1[k][j] * partials1[k][j];
            }
            partials3[k][0] = sum1;

            if (useScaling) {
                scalePartials(originIndex, k, new int[]{0});
                logP += Math.log(partials3[k][0]) + getLogScalingFactor(k);
            } else {
                logP += Math.log(partials3[k][0]);
            }
        }

        return logP;
    }

    public void setNodeMatrices(int nodeIndex, double[][] siteMatrices) {
        currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
        for (int siteNum = 0; siteNum < nrOfSites; siteNum++) {
            System.arraycopy(siteMatrices[siteNum], 0, 
                            matrices[currentMatrixIndex[nodeIndex]][nodeIndex][siteNum], 0, 
                            matrixSizes[siteNum]);
        }
    }

    /**
     * Scales partials for numerical stability.
     *
     * @param nodeIndex Index of the node
     * @param siteNum Site number
     * @param possibleStates Array of possible states
     */
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
        scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][siteNum] = scaleFactor < SCALING_THRESHOLD ? Math.log(scaleFactor) : 0.0;
    }

    @Override
    public double getLogScalingFactor(int siteNum) {    
        return java.util.stream.IntStream.range(0, nrOfNodes)
            .mapToDouble(i -> scalingFactors[currentPartialsIndex[i]][i][siteNum])
            .sum();
    }

    /**
     * Sets whether to use scaling for numerical stability.
     *
     * @param status Whether to use scaling
     */
    public void setUseScaling(boolean status) {
        useScaling = status;
    }

    @Override
    public void store() {
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

    // Unused methods required by the interface
    @Override
    public void initialize(int nodeCount, int siteCount, int matrixCount, 
            boolean integrateCategories, boolean useAmbiguities) {}
    
    @Override
    public void finalize() throws java.lang.Throwable {}

    @Override
    public void createNodePartials(int nodeIndex) {}

    @Override
    public void setNodePartials(int nodeIndex, double[] partials) {}
    
    @Override
    public void getNodePartials(int nodeIndex, double[] partials) {}

    @Override
    public void setNodePartialsForUpdate(int nodeIndex) {}

    @Override
    public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {}

    @Override
    public void setNodeStates(int nodeIndex, int[] states) {}
    
    @Override
    public void getNodeStates(int nodeIndex, int[] states) {}
    
    @Override
    public void getNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {}

    @Override
    public void setNodeMatrixForUpdate(int nodeIndex) {}
    
    @Override
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {}
    
    @Override
    public void calculateLogLikelihoods(double[] partials, double[] frequencies, 
            double[] outLogLikelihoods) {}
    
    @Override
    protected void calculateIntegratePartials(double[] inPartials, double[] proportions, 
            double[] outPartials) {}
    
    @Override
    public void setUseScaling(double scale) {}

    @Override
    public void unstore() {}
}
