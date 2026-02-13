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
 * Code based on classes from feast by Tim Vaughan
 */

@Description("Populate a alignment from a mutation matrix.")
public class AlignmentFromXSV extends Alignment {

    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file containing the mutation matrix.", Input.Validate.REQUIRED);

    public Input<String> sepInput = new Input<>("sep",
            "Separator used in the mutation matrix file (\",\" or \"\\t\"; default is tab).", "\t");

    public Input<Boolean> skipHeaderRowInput = new Input<>("skipHeaderRow",
            "If true, skip first row. (Default true.)", true);


    public AlignmentFromXSV() { }

    @Override
    public void initAndValidate() {
        // Separate two steps of initialization to allow for easier subclass extension with additional processing after reading in the sequences but before finalizing the alignment
        readSequencesFromFile();
        finalizeAlignment();
    }

    protected void finalizeAlignment() {
        super.initAndValidate();
    }

    protected void readSequencesFromFile() {
        // Guard against double-initialization
        if (!sequenceInput.get().isEmpty()) return;

        final String sep = sepInput.get();

        try (BufferedReader reader = new BufferedReader(new FileReader(fileNameInput.get()))) {
            if (skipHeaderRowInput.get()) {
                reader.readLine();
            }

            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                String[] fields = line.split(sep);
                // Row names are taxon names, and are in the first column
                String taxonName = fields[0].trim();  
                // Mutations are the remaining columns, and must be comma-separated for the Sequence datatype string
                String taxonData = String.join(",", java.util.Arrays.copyOfRange(fields, 1, fields.length));    
                // Assign mutations to taxon
                sequenceInput.setValue(new Sequence(taxonName, taxonData), this);
            }
        } catch (IOException e) {
            throw new RuntimeException("Error reading from file '" + fileNameInput.get() + "': " + e.getMessage());
        }
    }
}