package beam.likelihood;

import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
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

import beam.likelihood.IrreversibleLikelihoodCore;
import beam.substitutionmodel.CrisprSubstitutionModel;


/**
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

    /** Input for the origin time of the cell division process */
    public Input<RealParameter> originInput = new Input<>("origin",
            "Start of the cell division process, usually start of the experiment.",
            Input.Validate.OPTIONAL);

    /** Alignment data */
    protected Alignment data;

    /** Height of the origin node */
    protected Double originHeight;

    /** Substitution model for mutations */
    protected CrisprSubstitutionModel substitutionModel;

    /** State representing missing data */
    protected int missingDataState;

    /** Number of sites in the alignment */
    protected int nrOfSites;

    /** Transition probability matrix */
    protected double[] probabilities;

    /** Core likelihood calculation engine */
    protected IrreversibleLikelihoodCore likelihoodCore;

    /**
     * Flag indicating the state of the tree:
     * <ul>
     *   <li>CLEAN=0: nothing needs to be recalculated for the node</li>
     *   <li>DIRTY=1: node partial needs to be recalculated</li>
     *   <li>FILTHY=2: indices for the node need to be recalculated</li>
     * </ul>
     */
    protected int hasDirt;

    @Override
    public void initAndValidate() {
        // Get and validate inputs
        substitutionModel = (CrisprSubstitutionModel) ((SiteModel.Base) siteModelInput.get()).substModelInput.get();
        data = substitutionModel.getData(); // We use the data from the substitution model since the missing data state (-1) is replaced there during the rate matrix setup

        if (data.getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("Number of taxa in alignment does not match number of leaves in tree");
        }
        if (originInput.get() == null || branchRateModelInput.get() == null) {
            throw new IllegalArgumentException("Origin time and branch rate model must be specified");
        }
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("Site model must be of type SiteModel.Base");
        }

        // Initialize model parameters
        originHeight = originInput.get().getValue();
        missingDataState = substitutionModel.getMissingState();
        nrOfSites = data.getPatternCount();

        // Initialize likelihood core
        int nrOfStates = substitutionModel.getStateCount();
        probabilities = new double[nrOfStates * nrOfStates];
        likelihoodCore = new IrreversibleLikelihoodCore(treeInput.get().getNodeCount() + 1, nrOfStates, nrOfSites, missingDataState);

        // Set initial states
        setStates(treeInput.get().getRoot());
    }

    /**
     * Sets the data at the tips in the likelihood core.
     *
     * @param node The node to process
     */
    protected void setStates(Node node) {
        if (node.isLeaf()) {
            int taxonIndex = data.getTaxonIndex(node.getID());
            if (taxonIndex == -1) {
                throw new RuntimeException("Could not find sequence " + node.getID() + " in the alignment");
            }

            // Get states directly without intermediate array
            List<Integer> seq = data.getCounts().get(taxonIndex);
            int[] states = new int[nrOfSites];
            for (int i = 0; i < nrOfSites; i++) {
                int state = seq.get(i);
                states[i] = state;
            }
            likelihoodCore.setNodePartials(node.getNr(), states);
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
            likelihoodCore.setUseScaling(true);
            traverse(root);
            if (logP == Double.NEGATIVE_INFINITY) {
                throw new RuntimeException("Likelihood is still negative infinity after scaling");
            }
        }

        return logP;
    }

    /**
     * Computes the transition pre-order, then does a post-order traversal
     * calculating the partial likelihoods by a modified pruning algorithm
     * taking advantage of the irreversibility of the substitution model
     * to reduce the number of ancestral states to propagate partial
     * likelihoods for.
     *
     * @param node The node to traverse
     * @return The update status of the node
     */
    protected int traverse(Node node) {
        final int nodeIndex = node.getNr();
        int update = (node.isDirty() | hasDirt);

        // Update transition probabilities if needed
        if (update != Tree.IS_CLEAN) {
            double parentHeight = node.isRoot() ? originHeight : node.getParent().getHeight();
            substitutionModel.getTransitionProbabilities(node, parentHeight, node.getHeight(), 
                    branchRateModelInput.get().getRateForBranch(node), probabilities);
            likelihoodCore.setNodeMatrix(nodeIndex, 0, probabilities);
            update |= Tree.IS_DIRTY;
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
                
                likelihoodCore.setPossibleAncestralStates(childIndex1, childIndex2, nodeIndex);
                likelihoodCore.calculatePartials(childIndex1, childIndex2, nodeIndex);

                // Calculate origin partials at root
                if (node.isRoot()) {
                    logP = likelihoodCore.calculateLogLikelihoods(nodeIndex, nodeIndex + 1);
                }

                update |= (update1 | update2);
            }
        }

        return update;
    }

    @Override
    public void store() {
        super.store();
        likelihoodCore.store();
    }

    @Override
    public void restore() {
        super.restore();
        likelihoodCore.restore();
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