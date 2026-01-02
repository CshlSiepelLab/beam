package beam.substitutionmodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.IntegerData;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;

/**
 * TideTree substitution model that can be used with the modified BEAGLE tree likelihood
 * under the assumption that editing happens during the entire experiment.
 *
 * @author Stephen Staklinski
 */
@Description("TideTree substitution model that can be used with the modified BEAGLE tree likelihood" +
            "under the assumption that editing happens during the entire experiment.")
public class BeamMutationSubstitutionModel extends SubstitutionModel.Base {

    // Input for the silencing rate throughout the experiment
    public final Input<RealParameter> silencingRateInput = new Input<>("silencingRate", "Rate at which barcodes are silenced throughout the entire experiment", Validate.REQUIRED);
    // Input for edit rates during the editing window
    public final Input<List<RealParameter>> editRatesInput = new Input<>("editRates", "Input rates at which edits are introduced into the genomic barcode during the editing window", new ArrayList<>(), Validate.OPTIONAL);
    // Input for edit rates during the editing window
    final public Input<Alignment> dataInput = new Input<>("data", "EditData from which to calculate edit rates if they are not provided", Validate.OPTIONAL);

    // Frequencies of states
    private double[] frequencies;
    // Edit rate parameter
    private RealParameter editRate_;
    // Silencing rate parameter
    private RealParameter silencingRate_;
    // Array of edit rates
    private Double[] editRates;
    // State representing missing data
    private int missingDataState;
    // Stored silencing rate for dirty state checking
    private double storedSilencingRate;

    @Override
    public void initAndValidate() {
        // Make sure that data is provided if edit rates are not
        if (editRatesInput.get().isEmpty() && dataInput.get() == null) {
            throw new RuntimeException("Either edit rates or data to calculate edit rates must be provided!");
        }

        // Edit rates can be provided
        Double editRateSum = 0.0;
        if (!editRatesInput.get().isEmpty()) {
            System.out.println("Using provided edit rates.");
            editRate_ = editRatesInput.get().get(0);
            editRates = new Double[editRate_.getDimension()];
            for (int i = 0; i < editRate_.getDimension(); i++) {
                editRates[i] = editRate_.getValue(i);
            }

            // Add edit rates to rate matrix
            for (double editRate : editRates) {
                if (editRate < 0) {
                    throw new RuntimeException("All edit rates must be positive!");
                }
                editRateSum += editRate;
            }
        }
        // Or data can be provided from which we will compute the edit rates
        else {
            System.out.println("Calculating empirical edit rates from data.");
            Alignment data = dataInput.get();
            DataType dataType = data.getDataType();
            int taxonCount = data.getTaxonCount();

            // Determine maximum state in the data, which is the missing data state
            int maxState = 0;
            for (int taxon = 0; taxon < taxonCount; taxon++) {
                List<Integer> seq = data.getCounts().get(taxon);
                for (int value : seq) {
                    maxState = Math.max(maxState, value);
                }
            }

            // Count occurrences of each state
            double[] stateCounts = new double[maxState - 1]; // Exclude missing data (last state)
            double total = 0.0;

            int maxValSeen = 0; // Make sure data is input in sequential order
            for (int taxon = 0; taxon < taxonCount; taxon++) {
                List<Integer> seq = data.getCounts().get(taxon);
                for (int value : seq) {
                    if (value == 0 || value == maxState) {
                        continue; // Skip unedited and missing data
                    }
                    stateCounts[value - 1] += 1;    // Use -1 index to skip unedited (0)
                    total += 1;

                    if (value > maxValSeen) {
                        if (value == maxValSeen + 1) {
                            maxValSeen = value;
                        }
                        else {
                            throw new RuntimeException("States in data must be sequentially ordered starting from 0!");
                        }
                    }
                }
            }

            // Convert counts → proportions for edit rates
            editRates = new Double[stateCounts.length];
            for (int i = 0; i < stateCounts.length; i++) {
                double editRate = stateCounts[i] / total;
                editRates[i] = editRate;
                editRateSum += editRate;
            }
        }

        System.out.println("Edit rates: " + Arrays.toString(editRates));

        if (Math.abs(editRateSum - 1.0) > 1e-5) {
            throw new RuntimeException("Sum of edit rates must be 1.0, but it is " + editRateSum + "!");
        }

        nrOfStates = editRates.length + 2;

        silencingRate_ = silencingRateInput.get();
        double silencingRate = silencingRate_.getValue();
        storedSilencingRate = silencingRate;

        if (silencingRate < 0) {
            throw new RuntimeException("Loss rate must be positive!");
        }

        // Missing data state is the last edit
        missingDataState = nrOfStates - 1;

        // Center root frequency on the unedited first state, regardless of input frequencies
        // as this is a property of the barcodes
        frequencies = new double[nrOfStates];
        frequencies[0] = 1;
    }

    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        // Calculate key parameters
        double silencingRate = silencingRate_.getValue();
        double delta = (startTime - endTime) * rate;
        double expOfDeltaLoss = Math.exp(-delta * silencingRate);
        double c = expOfDeltaLoss * (1 - Math.exp(-delta));

        // Initialize matrix to zeros
        Arrays.fill(matrix, 0.0);

        // Set absorbing state (bottom right corner)
        matrix[nrOfStates * nrOfStates - 1] = 1.0;

        // Set diagonal elements and loss probabilities in one pass
        for (int i = 0; i < nrOfStates - 1; i++) {
            int rowOffset = i * nrOfStates;
            matrix[rowOffset + i] = (i == 0) ? Math.exp(-delta * (1 + silencingRate)) : expOfDeltaLoss;
            matrix[rowOffset + (nrOfStates - 1)] = 1 - expOfDeltaLoss;
        }

        // Set edit probabilities (first row)
        for (int j = 1; j < nrOfStates - 1; j++) {
            matrix[j] = editRates[j - 1] * c;
        }
    }

    @Override
    public void store() {
        storedSilencingRate = silencingRate_.getValue();
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        // Since edit rates are fixed hyperparameters, the only thing that can change is the silencing rate
        return silencingRate_.getValue() != storedSilencingRate;
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        // Return true for BEAGLE compatibility
        return true;
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return null;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof IntegerData;
    }

    @Override
    public double[] getFrequencies() {
        return frequencies;
    }

    /**
     * Returns the state representing missing data.
     *
     * @return The missing data state
     */
    public int getMissingState() {
        return missingDataState;
    }
}

