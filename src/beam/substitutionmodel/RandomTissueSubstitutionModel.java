package beam.substitutionmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.SubstitutionModel;

/**
* This class always returns transition probabilities equal to the stationary
* distribution of the states. This is useful for testing an ancestral state
* random sampling control model for model selection purposes between
* informative vs. uninformative data.
* 
* @author Stephen Staklinski
**/

@Description("Random substitution model implementation that simply returns equillibrium distribution transition probabilities.")

public class RandomTissueSubstitutionModel extends SubstitutionModel.Base {

    public Input<RealParameter> piInput = new Input<>("pi", "Stationary frequency of the first state", Validate.REQUIRED);

     @Override
     public void initAndValidate(){
        super.initAndValidate();
        nrOfStates = frequencies.getFreqs().length;
        if (nrOfStates < 2) throw new IllegalArgumentException("Number of tissues must be at least 2.");
     }


    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        // Always returns stationary distribution probabilities for random sampling
        // pi is the primary tissue frequency and the rest are uniform over the remaining tissues given 1-pi
        double pi = piInput.get().getValue();
        if (pi < 0.0 || pi > 1.0) throw new IllegalArgumentException("pi must be in [0,1], but it is " + pi);
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

    @Override
    protected boolean requiresRecalculation() {
        if (piInput.get() != null && piInput.get().somethingIsDirty()) {
            return true;
        }
        return false;
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return true;
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {return null;}

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType.getStateCount() != Integer.MAX_VALUE;
    }
}

