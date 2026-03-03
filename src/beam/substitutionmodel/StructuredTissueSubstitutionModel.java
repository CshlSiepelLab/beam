package beam.substitutionmodel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.lang.reflect.InvocationTargetException;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.inference.parameter.RealParameter;

/**
 * Structured tissue substitution model implementing a user defined paramaterization of
 * the rate matrix (0 = fixed rate of 0; 1,2,...,N other values index the input rates as 
 * free parameters). Tissue-specific frequencies parameterized by a primary tissue frequency 
 * (pi) and uniform distribution for remaining tissues.
 * 
 * WARNING: Please make sure the order of states in input model structure file matches the order
 * of states in the input data, which can be input by the user directly to the TissueData alignment 
 * rather than relying on internal order determination which may not match expectations.
 *
 * @author Stephen Staklinski
 */
@Description("Structured tissue substitution model implementation.")
public class StructuredTissueSubstitutionModel extends GeneralSubstitutionModel {

    // Input for the model structure
    public Input<String> modelStructureFileInput = new Input<>("modelStructureFile",
            "Filepath for the model structure CSV file", Validate.REQUIRED);

    /** Input for the stationary frequency of the first state */
    public Input<RealParameter> piInput = new Input<>("pi",
            "Stationary frequency of the first state",
            Validate.REQUIRED);

    // Store model structure
    private int[][] rateIndices;

    @Override
    public void initAndValidate() {
        updateMatrix = true;
        frequencies = frequenciesInput.get();
        nrOfStates = frequencies.getFreqs().length;
        if (nrOfStates < 2) throw new IllegalArgumentException("Number of tissues must be at least 2.");
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[ratesInput.get().getDimension()];

        // Read model structure file
        rateIndices = readStructureFile(modelStructureFileInput.get());
        // Validate model structure dimensions match nrOfStates
        if (rateIndices.length != nrOfStates) throw new IllegalArgumentException("Model structure file must have " + nrOfStates + " rows.");
        for (int i = 0; i < nrOfStates; i++) {
            if (rateIndices[i].length != nrOfStates) throw new IllegalArgumentException("Model structure file must be square with dimension " + nrOfStates + ".");
        }
        
        // Initialize the eigen system
        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException |
                 IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
    }


    private int[][] readStructureFile(String filepath) {
        List<int[]> rows = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filepath))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                String[] tokens = line.split(",");
                int[] row = new int[tokens.length];
                for (int i = 0; i < tokens.length; i++) {
                    row[i] = Integer.parseInt(tokens[i].trim());
                }
                rows.add(row);
            }
        } catch (IOException e) {
            throw new IllegalArgumentException("Error reading model structure file: " + filepath, e);
        }

        int[][] matrix = new int[rows.size()][];
        for (int i = 0; i < rows.size(); i++) {
            matrix[i] = rows.get(i);
        }

        return matrix;
    }


    /**
     * Sets up the rate matrix based on the input structure.
     * The matrix is normalized to expect one substitution per unit time with respect
     * to the stationary frequencies. The frequencies are parameterized by pi specifying
     * the primary tissue frequency and 1-pi for the remaining tissues in a uniform distribution.
     */
    @Override
    public void setupRateMatrix() {
        // Calculate the stationary frequencies based on pi
         double pi = piInput.get().getValue();
         if (pi < 0.0 || pi > 1.0) throw new IllegalArgumentException("pi must be in [0,1], but it is " + pi);
         double[] piFreqs = new double[nrOfStates];
         piFreqs[0] = pi;
         double piUniformOthers = (1 - pi) / (nrOfStates - 1);
         for (int i = 1; i < nrOfStates; i++) {
             piFreqs[i] = piUniformOthers;
         }
        
        // Initialize rate matrix off-diagonal rates and calculate row sums in one pass
        double[] rowSums = new double[nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0.0; // Set diagonal to zero
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j) {
                    double rate = (rateIndices[i][j] == 0) ? 0.0 : relativeRates[rateIndices[i][j] - 1];
                    double rateJ = rate * piFreqs[j];
                    rateMatrix[i][j] = rateJ;
                    rowSums[i] += rateJ;
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
    protected boolean requiresRecalculation() {
        if (piInput.get() != null && piInput.get().somethingIsDirty()) {
            updateMatrix = true;
            return true;
        }
        return super.requiresRecalculation();   // Will always return true as currently setup
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }
}
