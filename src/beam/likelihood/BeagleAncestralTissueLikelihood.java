package beam.likelihood;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beagle.Beagle;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.UserDataType;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;
import beastclassic.evolution.tree.TreeTrait;
import beastclassic.evolution.tree.TreeTraitProvider;
import beastclassic.evolution.likelihood.LeafTrait;

/**
 * Ancestral State Tree Likelihood that extends to origin input for a branch above the root.
 *
 * @author Stephen Staklinski
 */
@Description("Ancestral State Tree Likelihood that extends the custom Beagle likelihood with origin input for the branch above the root.")
public class BeagleAncestralTissueLikelihood extends BeagleTissueLikelihood implements TreeTraitProvider {

    // Key for states trait
    public static final String STATES_KEY = "states";
    // Input for the label used for tissues written to posterior samples
    public Input<String> tagInput = new Input<>("tag", "label used to report trait in the output tree posterior samples", Validate.REQUIRED);
    // Input for leaf traits
    public Input<List<LeafTrait>> leafTraitsInput = new Input<>("leaftrait", "list of leaf traits", new ArrayList<>());

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
    // Number of possible states
    private int stateCount;
    // States at tip nodes
    private int[] tipStates;
    // Helper for tree traits
    protected TreeTraitProvider.Helper treeTraits = new Helper();


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (dataInput.get().getSiteCount() != 1 || m_siteModel.getCategoryCount() > 1) {
            throw new RuntimeException("Only single-site data and no site rate categories are supported, meaning only one tissue per tip is allowed.");
        }

        // Initialize data structures
        tag = tagInput.get();
        Alignment data = dataInput.get();
        dataType = data.getDataType();
        stateCount = dataType.getStateCount();
        TreeInterface treeModel = treeInput.get();
        int nodeCount = treeModel.getNodeCount();

        // Initialize tip states
        tipStates = new int[treeModel.getLeafNodeCount()];
        for (Node node : treeModel.getExternalNodes()) {
            String taxon = node.getID();
            int taxonIndex = data.getTaxonIndex(taxon);;
            int code = data.getPattern(taxonIndex, 0);
            tipStates[node.getNr()] = data.getDataType().getStatesForCode(code)[0];
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

    // Gets the state for a given node in the tree.
    public int getStateForNode(TreeInterface tree, Node node) {
        if (!areStatesRedrawn) {
            traverseSample(tree, tree.getRoot(), -1);
            areStatesRedrawn = true;
        }
        return reconstructedStates[node.getNr()];
    }


    @Override
    public double calculateLogP() {
        super.calculateLogP();
        areStatesRedrawn = false;
        return logP;
    }

    // Traverses the tree in pre-order to sample internal node states
    public void traverseSample(TreeInterface tree, Node node, int parentState) {
        int nodeNum = node.getNr();
        Node parent = node.getParent();
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
                    double[] oPs = new double[m_nStateCount];
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
            traverseSample(treeInput.get(), child1, reconstructedStates[node.getNr()]);
            Node child2 = node.getChild(1);
            traverseSample(treeInput.get(), child2, reconstructedStates[node.getNr()]);
            
        } else {
            // Handle leaf node where tissue is known
            reconstructedStates[node.getNr()] = tipStates[node.getNr()];
        }
    }

    @Override
    public void store() {
        super.store();
        System.arraycopy(reconstructedStates, 0, storedReconstructedStates, 0, reconstructedStates.length);
        storedAreStatesRedrawn = areStatesRedrawn;
    }

    @Override
    public void restore() {
        super.restore();
        System.arraycopy(storedReconstructedStates, 0, reconstructedStates, 0, storedReconstructedStates.length);
        areStatesRedrawn = storedAreStatesRedrawn;
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
