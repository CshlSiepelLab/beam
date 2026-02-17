package beam.substitutionmodel;

import java.io.PipedInputStream;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.HashMap;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

/**
 * General GTR substitution model implementation for tissue-specific evolution.
 * This model implements a GTR substitution model with tissue-specific frequencies
 * parameterized by a primary tissue frequency (pi) and uniform distribution for
 * remaining tissues.
 *
 * @author Stephen Staklinski
 */
@Description("General GTR substitution model implementation.")
public class TissueSubstitutionModel extends GeneralSubstitutionModel {

    /** Input for the stationary frequency of the first state */
    public Input<RealParameter> piInput = new Input<>("pi",
            "Stationary frequency of the first state",
            Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        updateMatrix = true;
        frequencies = frequenciesInput.get();
        nrOfStates = frequencies.getFreqs().length;
        rateMatrix = new double[nrOfStates][nrOfStates];

        validateInputRates();
        initializeEigenSystem();
        relativeRates = new double[ratesInput.get().getDimension()];
    }

    /**
     * Validates that the number of input rates is correct for the GTR model.
     * The number of rates should be equal to ((nrOfStates * nrOfStates) - nrOfStates) / 2
     * for all off-diagonal rates on one side of the matrix.
     *
     * @throws IllegalArgumentException if the number of input rates is incorrect
     */
    private void validateInputRates() {
        int nrInputRates = ratesInput.get().getDimension();
        int expectedRates = ((nrOfStates * nrOfStates) - nrOfStates) / 2;
        
        if (nrInputRates != expectedRates) {
            throw new IllegalArgumentException(
                String.format("The number of input rates must be %d for all off-diagonal rates " +
                            "on one side of the matrix, but it is %d. " +
                            "Check the dimension of the input rate parameters.",
                    expectedRates, nrInputRates));
        }
    }

    /**
     * Initializes the eigen system for the substitution model.
     *
     * @throws IllegalArgumentException if there is an error creating the eigen system
     */
    private void initializeEigenSystem() {
        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException |
                 IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
    }

    /**
     * Sets up the rate matrix for the GTR model.
     * The matrix is normalized to expect one substitution per unit time with respect
     * to the stationary frequencies. The frequencies are parameterized by pi specifying
     * the primary tissue frequency and 1-pi for the remaining tissues in a uniform distribution.
     */
    @Override
    public void setupRateMatrix() {
        // Calculate tissue frequencies once
        double[] piFreqs = calculatePiFrequencies();
        double piUniformOthers = piFreqs[1]; // All non-primary tissues have same frequency
        
        // Initialize rate matrix and calculate row sums in one pass
        double[] rowSums = new double[nrOfStates];
        int count = 0;
        
        // First pass: set up off-diagonal rates and accumulate row sums
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
        
        // Second pass: set diagonal elements and calculate normalization factor
        double subst = 0.0;
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = -rowSums[i];
            subst += -rateMatrix[i][i] * piFreqs[i];
        }
        
        // Normalize the matrix if needed
        if (subst != 0.0) {
            double scale = 1.0 / subst;
            for (int i = 0; i < nrOfStates; i++) {
                double[] row = rateMatrix[i];
                for (int j = 0; j < nrOfStates; j++) {
                    row[j] *= scale;
                }
            }
        }
    }

    /**
     * Calculates the tissue frequencies based on the primary tissue frequency (pi).
     *
     * @return Array of tissue frequencies
     */
    private double[] calculatePiFrequencies() {
        double pi = piInput.get().getValue();
        double[] piFreqs = new double[nrOfStates];
        piFreqs[0] = pi;
        
        double piUniformOthers = (1 - pi) / (nrOfStates - 1);
        for (int i = 1; i < nrOfStates; i++) {
            piFreqs[i] = piUniformOthers;
        }
        
        return piFreqs;
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }
}
