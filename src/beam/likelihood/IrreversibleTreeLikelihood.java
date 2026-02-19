package beam.likelihood;

import java.util.Arrays;
import java.util.List;

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
            Input.Validate.REQUIRED);

    /** Alignment data */
    protected Alignment data;
    /** Height of the origin node */
    protected Double originHeight;
    /** Substitution model for mutations */
    protected CrisprSubstitutionModel substitutionModel;
    /** Number of sites in the alignment */
    protected int nrOfSites;
    /** Transition probability matrix */
    protected double[][] probabilities;
    /** Core likelihood calculation engine */
    protected IrreversibleLikelihoodCore likelihoodCore;
    /** Whether the substitution model is global */
    protected boolean isGlobalModel;
    /**
     * Flag indicating the state of the tree:
     * CLEAN=0: nothing needs to be recalculated for the node</li>
     * DIRTY=1: node partial needs to be recalculated</li>
     * FILTHY=2: indices for the node need to be recalculated</li>
     */
    protected int hasDirt;

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
        int[] nrOfStates = substitutionModel.getNrOfStatesPerSite();    // Reference, not copy
        isGlobalModel = substitutionModel.isGlobalModel();

        // Initialize transition probability matrices per site
        probabilities = new double[nrOfSites][];
        for (int i = 0; i < nrOfSites; i++) {
            probabilities[i] = new double[nrOfStates[i] * nrOfStates[i]];
        }

        // Initialize likelihood core
        likelihoodCore = new IrreversibleLikelihoodCore(treeInput.get().getNodeCount() + 1, 
                    nrOfStates, nrOfSites, substitutionModel.getMissingStates(), isGlobalModel);

        // Set initial states
        setStates(treeInput.get().getRoot());
    }

    // Sets the observed data in likelihood core for leaf nodes across all alignment sites
    protected void setStates(Node node) {
        if (node.isLeaf()) {
            int taxonIndex = data.getTaxonIndex(node.getID());
            if (taxonIndex == -1) throw new RuntimeException("Could not find sequence " + node.getID() + " in the alignment");
            likelihoodCore.setNodePartials(node.getNr(), data.getCounts().get(taxonIndex)); // Pass data as reference
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
            likelihoodCore.setNodeMatrices(nodeIndex, probabilities);
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