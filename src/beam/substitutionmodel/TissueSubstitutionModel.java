package beam.substitutionmodel;

import java.lang.reflect.InvocationTargetException;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.inference.parameter.RealParameter;

/**
 * General tissue substitution model, implementing a GTR (default) or GTI tissue 
 * substitution model and tissue-specific frequencies parameterized by a primary 
 * tissue frequency (pi) and uniform distribution for remaining tissues.
 *
 * @author Stephen Staklinski
 */
@Description("General tissue substitution model implementation.")
public class TissueSubstitutionModel extends GeneralSubstitutionModel {

    // Input for the model structure
    public Input<String> modelInput = new Input<>("model",
            "String specifying the model structure (e.g., 'GTR' or 'GTI'). Default is GTR",
            "GTR", Validate.OPTIONAL);

    /** Input for the stationary frequency of the first state */
    public Input<RealParameter> piInput = new Input<>("pi",
            "Stationary frequency of the first state",
            Validate.REQUIRED);

    @Override
    public void initAndValidate() {

        if (!modelInput.get().equals("GTR") && !modelInput.get().equals("GTI")) {
            throw new IllegalArgumentException("Invalid model type: " + modelInput.get() + 
                ". Supported models are 'GTR' and 'GTI'.");
        }
        
        updateMatrix = true;
        frequencies = frequenciesInput.get();
        nrOfStates = frequencies.getFreqs().length;
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[ratesInput.get().getDimension()];

        validateInputRates();
        
        // Initialize the eigen system
        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException |
                 IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
    }

    /**
     * Validates that the number of input rates is correct for the GTR of GTI models.
     * The number of rates should be equal to ((nrOfStates * nrOfStates) - nrOfStates) / 2
     * for all off-diagonal rates on one side of the matrix for GTR or 
     * (nrOfStates * nrOfStates) - nrOfStates) for GTI.
     *
     * @throws IllegalArgumentException if the number of input rates is incorrect
     */
    private void validateInputRates() {
        // Calculate the expected number of rates based on the model type
        int expectedRates;
        if (modelInput.get().equals("GTR")) {
            // GTR has one rate for each off-diagonal tissue transition, but only on one side of the matrix
            expectedRates = ((nrOfStates * nrOfStates) - nrOfStates) / 2;
        } else {
            // GTI has one rate for each off-diagonal tissue transition
            expectedRates = nrOfStates * nrOfStates - nrOfStates;
        }

        // Validate the number of input rates against the expected number
        if (ratesInput.get().getDimension() != expectedRates) {
            throw new IllegalArgumentException(
                String.format("The number of input rates must be %d for all off-diagonal rates " +
                            "on one side of the matrix, but it is %d. " +
                            "Check the dimension of the input rate parameters.",
                    expectedRates, ratesInput.get().getDimension()));
        }
    }


    /**
     * Sets up the rate matrix for the GTR or GTI model.
     * The matrix is normalized to expect one substitution per unit time with respect
     * to the stationary frequencies. The frequencies are parameterized by pi specifying
     * the primary tissue frequency and 1-pi for the remaining tissues in a uniform distribution.
     */
    @Override
    public void setupRateMatrix() {
        // Calculate the stationary frequencies based on pi
         double pi = piInput.get().getValue();
         double[] piFreqs = new double[nrOfStates];
         piFreqs[0] = pi;
         double piUniformOthers = (1 - pi) / (nrOfStates - 1);
         for (int i = 1; i < nrOfStates; i++) {
             piFreqs[i] = piUniformOthers;
         }
        
        // Initialize rate matrix off-diagonal rates and calculate row sums in one pass
        double[] rowSums = new double[nrOfStates];
        int count = 0;
        if (modelInput.get().equals("GTR")) {
            for (int i = 0; i < nrOfStates; i++) {
                rateMatrix[i][i] = 0.0; // Set diagonal to zero
                for (int j = i + 1; j < nrOfStates; j++) {
                    double rate = relativeRates[count];
                    double rateI = rate * piFreqs[i];
                    double rateJ = rate * piFreqs[j];
                    rateMatrix[i][j] = rateJ;
                    rateMatrix[j][i] = rateI;
                    rowSums[i] += rateJ;
                    rowSums[j] += rateI;
                    count++;
                }
            }
        } else if (modelInput.get().equals("GTI")) {
            for (int i = 0; i < nrOfStates; i++) {
                rateMatrix[i][i] = 0.0; // Set diagonal to zero
                for (int j = 0; j < nrOfStates; j++) {
                    if (i != j) {
                        double rate = relativeRates[count];
                        double rateJ = rate * piFreqs[j];
                        rateMatrix[i][j] = rateJ;
                        rowSums[i] += rateJ;
                        count++;
                    }
                }
            }
        }
        
        // Set diagonal elements and calculate normalization factor
        double subst = 0.0;
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = -rowSums[i];
            subst += -rateMatrix[i][i] * piFreqs[i];
        }
        
        // Normalize the rate matrix to one expected substitution per unit time
        if (subst != 0.0) {
            for (int i = 0; i < nrOfStates; i++) {
                for (int j = 0; j < nrOfStates; j++) {
                    rateMatrix[i][j] = rateMatrix[i][j] / subst;
                }
            }
        }
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }
}
