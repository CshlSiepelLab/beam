package beam.tree;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.BEASTObject;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import beam.likelihood.TissueLikelihood;

@Description("Logs tree annotated with tissues and including the origin branch length above the root.")
public class TreeWithTissuesAndOriginLogger extends BEASTObject implements Loggable {

    public Input<Tree> m_tree = new Input<Tree>("tree", "Tree to be logged", Validate.REQUIRED);

    public Input<TissueLikelihood> metadataInput = new Input<TissueLikelihood>("tissues", "Tissues to be logged with the tree nodes.", Validate.REQUIRED);

    TissueLikelihood tissues;

    @Override
    public void initAndValidate() {
        tissues = metadataInput.get();
    }

    @Override
    public void init(PrintStream out) {
        m_tree.get().init(out);
    }

    @Override
    public void log(long nSample, PrintStream out) {
        Tree tree = (Tree) m_tree.get().getCurrent();
        out.print("tree STATE_" + nSample + " = ");
        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), tree));
        out.print(";");
    }

    String toNewick(Node node, Tree tree) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), tree));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), tree));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }

        buf.append("[&");
        buf.append("location=");
        buf.append("\"").append(tissues.getStateStringForNode(tree, node)).append("\"");
        buf.append(']');
        
        // Use origin branch length for root, otherwise use normal branch length
        double branchLength = node.isRoot() 
            ? tissues.originHeight - node.getHeight() 
            : node.getLength();
        buf.append(":").append(branchLength);

        return buf.toString();
    }

    @Override
    public void close(PrintStream out) {
        m_tree.get().close(out);
    }
}