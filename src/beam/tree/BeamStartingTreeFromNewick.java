package beam.tree;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Stephen Staklinski
 */
@Description("This class uses the normal TreeParser abilities of reading in a starting tree" +
            " while also taking into account the origin tree height.")
public class BeamStartingTreeFromNewick extends TreeParser {

    // Optional file input to read a Newick string from.
    public Input<String> newickFileInput = new Input<>("newickFile",
            "Optional path to a file containing a Newick tree string. If provided, this value is read into the newick input string.",
            Input.Validate.OPTIONAL);

    // Input for the origin time of the cell division process
    public Input<RealParameter> originInput = new Input<>("origin",
            "Start of the cell division process, usually start of the experiment." +
            " If this is provided, the input tree will be rescaled so that the root height is close to this origin time." +
            " This is important for when a tree is input without branch lengths, so that here it starts with a tree" +
            " with a more reasonable height than the default of 0.01 total height. (Default: No input value, in which case" +
            " the tree is not rescaled.)", Input.Validate.OPTIONAL);

    @Override
    public void initAndValidate() {
        if (newickFileInput.get() != null) {
            // Make sure we do not have conflicting inputs of both file and string
            if (newickInput.get() != null) {
                throw new IllegalArgumentException("Both 'newick' and 'newickFile' inputs were provided. Please provide only one of these inputs.");
            }

            // Read in the provied Newick string from the file
            try {
                String fileNewick = Files.readString(Path.of(newickFileInput.get())).trim();
                if (fileNewick.endsWith(";")) {
                    fileNewick = fileNewick.substring(0, fileNewick.length() - 1);
                }
                newickInput.setValue(fileNewick, this);
            } catch (IOException e) {
                throw new IllegalArgumentException("Could not read Newick file at '" + newickFileInput.get() + "': " + e.getMessage(), e);
            }
        } else if (newickInput.get() == null) {
            throw new IllegalArgumentException("BeamStartingTreeFromNewick requires either 'newick' or 'newickFile' to be provided.");
        }

        super.initAndValidate();

        if (originInput.get() != null) {
            /*
             * Rescale the tree so that the root height matches close to the origin time
             * This is important for when a tree is input without branch lengths, so that here
             * it starts the MCMC with a tree with a more reasonable height than the default of
             * 0.01 total height. We offset the provided origin slightly to prevent errors
             * in which the tree height is exactly the same as the origin time.
             */
            System.out.println("Rescaling the input starting tree since the origin time was provided with the tree input.");
            System.out.println("Original input starting tree: ");
            System.out.println(getRoot().toNewick());

            double desired = originInput.get().getValue() - 0.1;
            Node root = getRoot();
            double current = root.getHeight();
            double factor = desired / current;

            System.out.println("Current root height: " + current);
            System.out.println("Target root height: " + desired);
            System.out.println("Scaling factor: " + factor);

            scaleHeights(root, factor);

            System.out.println("Rescaled starting tree: ");
            System.out.println(getRoot().toNewick());
        }
    }

    private void scaleHeights(Node node, double factor) {
        node.setHeight(node.getHeight() * factor);
        for (int i = 0; i < node.getChildCount(); i++) {
            scaleHeights(node.getChild(i), factor);
        }
    }
}
