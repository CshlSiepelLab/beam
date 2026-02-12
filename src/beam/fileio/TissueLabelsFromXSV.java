package beam.fileio;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.TraitSet;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * @author Stephen Staklinski
 * Code based on classes from feast by Tim Vaughan
 */
@Description("Initialize tissue labels from a CSV/TSV file.")
public class TissueLabelsFromXSV extends TraitSet {

    public Input<String> fileNameInput = new Input<>("fileName", "Name of CSV/TSV file to extract values from.",
            Input.Validate.REQUIRED);

    public Input<String> sepInput = new Input<>("sep",
            "Separator used in the mutation matrix file (\",\" or \"\\t\"; default is \",\").", ",");

    public Input<Boolean> skipHeaderRowInput = new Input<>("skipHeaderRow",
            "If true, skip first row. (Default false.)", false);

    public TissueLabelsFromXSV() {
        traitsInput.setRule(Input.Validate.OPTIONAL);
        traitNameInput.setValue("tissue", this);
    }

    @Override
    public void initAndValidate() {

        StringBuilder traitSB = new StringBuilder();

        String sep = sepInput.get();

        try (BufferedReader reader = new BufferedReader(new FileReader(fileNameInput.get()))) {
            if (skipHeaderRowInput.get()) {
                reader.readLine(); // Skip header row
            }

            String line;
            while ((line = reader.readLine()) != null) {
                String[] elements = line.split(sep);
                String taxonName = elements[0].trim();
                String traitValue = elements[1].trim();

                // Check that taxonName exists
                if (taxaInput.get().getTaxaNames().contains(taxonName)) {
                    if (traitSB.length()>0) {
                        traitSB.append(",");    // Separate entries with commas for the TraitSet string format
                    }
                    traitSB.append(taxonName).append("=").append(traitValue);
                } else {
                    throw new IllegalArgumentException("Taxon '" + taxonName + "' in tissue file is not found in the other inputs.");
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("Error reading from file '" + fileNameInput.get() + "': " + e.getMessage());
        }

        traitsInput.setValue(traitSB.toString(), this);

        super.initAndValidate();
    }
}