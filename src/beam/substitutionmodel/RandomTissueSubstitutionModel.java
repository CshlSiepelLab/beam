package beam.substitutionmodel;


import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;

/**
 * @author Stephen Staklinski
 **/

@Description("Random substitution model implementation that simply returns equillibrium distribution transition probabilities.")

public class RandomTissueSubstitutionModel extends GeneralSubstitutionModel {

    public Input<RealParameter> piInput = new Input<>("pi", "Stationary frequency of the first state", Validate.REQUIRED);

    /*
     * This class always returns transition probabilities equal to the stationary
     * distribution of the states. This is useful for testing an ancestral state
     * random sampling control model for model selection purposes between
     * informative vs. uninformative data.
     */

     @Override
     public void initAndValidate(){
         updateMatrix = true;
         frequencies = frequenciesInput.get();
         nrOfStates = frequencies.getFreqs().length;
         rateMatrix = new double[nrOfStates][nrOfStates];
 
         // Verify the number of input rates is correct
         int nrInputRates = ratesInput.get().getDimension();
         relativeRates = new double[ratesInput.get().getDimension()];
     }


    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        // always returns stationary distribution probabilities for random sampling
        // pi is the primary tissue frequency and the rest are uniform over the remaining tissues given 1-pi
        double pi = piInput.get().getValue();
        double[] piFreqs = new double[nrOfStates];
        piFreqs[0] = pi;
        double piUniformOthers = (1 - pi) / (nrOfStates - 1);
        for (int i = 1; i < nrOfStates; i++) {
            piFreqs[i] = piUniformOthers;
        }

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                matrix[i * nrOfStates + j] = piFreqs[j];
            }
        }
    }


    /** sets up rate matrix**/
    @Override
    public void setupRateMatrix() {
        /*
         * This does nothing because the rate matrix is not
         * used in the likelihood calculations.
         */
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }
}

