package beam.fileio;

import java.util.Set;
import java.util.LinkedHashSet;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Sequence;

import beam.fileio.AlignmentFromXSV;
import beam.datatype.TissueData;

/**
 * @author Stephen Staklinski
 */
@Description("Initialize tissue labels from a CSV/TSV file.")
public class TissueLabelsFromXSV extends AlignmentFromXSV {

    public TissueLabelsFromXSV() {}

    @Override
    public void initAndValidate() {

        // Validate that the datatype input is TissueData
        if (!(userDataTypeInput.get() instanceof TissueData)) {
            throw new IllegalArgumentException("TissueLabelsFromXSV requires a TissueData datatype.");
        }

        // Read in the data by standard XSV read in to an alignment
        readSequencesFromFile();

        // Check that the TissueData codeMap and stateCount are set, and if not, set them from the input data
        TissueData tissueData = (TissueData) userDataTypeInput.get();
        if (tissueData.getCodeMap() == null) {
            // Make sure the primary tissue is always first
            String primaryTissue = tissueData.primaryTissueInput.get();
            Set<String> uniqueStates = new LinkedHashSet<>();
            uniqueStates.add(primaryTissue);
            // Get other tissues in the order in which they appear in the stored data
            for (Sequence seq : sequenceInput.get()) {
                String state = seq.getData().trim();
                if (!state.equals(primaryTissue)) {
                    uniqueStates.add(state);
                }
            }

            // Build comma-separated codeMap
            tissueData.setCodeMap(String.join(",", uniqueStates));
            // Set state count to number of unique tissues
            tissueData.setStateCount(uniqueStates.size());

            System.out.println("Found " + tissueData.getStateCount() + " unique tissue labels in the input data. Internal ordering set based on input data to: " + tissueData.getCodeMap());
        }

        finalizeAlignment();
    }
}