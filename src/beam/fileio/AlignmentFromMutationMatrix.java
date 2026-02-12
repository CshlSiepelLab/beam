package beam.fileio;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * @author Stephen Staklinski
 */

@Description("Populate a alignment from a mutation matrix.")
public class AlignmentFromMutationMatrix extends Alignment {

    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file containing the mutation matrix.", Input.Validate.REQUIRED);

    public Input<String> sepInput = new Input<>("sep",
            "Separator used in the mutation matrix file (\",\" or \"\\t\"; default is tab).", "\t");

    public Input<Boolean> skipHeaderRowInput = new Input<>("skipHeaderRow",
            "If true, skip first row. (Default true.)", true);


    public AlignmentFromMutationMatrix() { }

    @Override
    public void initAndValidate() {

        // Guard against double-initialization
        if (!sequenceInput.get().isEmpty()) {
            super.initAndValidate();
            return;
        }

        final String sep = sepInput.get();

        try (BufferedReader reader = new BufferedReader(new FileReader(fileNameInput.get()))) {
            if (skipHeaderRowInput.get()) {
                reader.readLine(); // Skip header row
            }

            StringBuilder seqBuilder;
            String taxonName;
            String line;

            while ((line = reader.readLine()) != null) {
                if (line.isEmpty()) continue;
                line = line.trim();
                String[] fields = line.split(sep);

                // Row names are taxon names, and are in the first column
                taxonName = fields[0];  

                // Mutations are the remaining columns, and must be comma-separated for the Sequence datatype string
                String mutations = String.join(",", java.util.Arrays.copyOfRange(fields, 1, fields.length));    

                // Assign mutatins to taxon
                sequenceInput.setValue(new Sequence(taxonName, mutations), this);
            }
        } catch (IOException e) {
            throw new RuntimeException("Error reading from file '" + fileNameInput.get() + "': " + e.getMessage());
        }

        super.initAndValidate();

    }
}