package beam.substitutionmodel;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.datatype.UserDataType;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;

import beam.datatype.TissueData;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * Logs the normalized rates from the rate matrix of a BEAM tissue substitution model.
 * 
 * @author Stephen Staklinski
 **/

@Description("Logger for BEAM tissue subsitution models.")
public class TissueSubstitutionModelLogger extends BEASTObject implements Loggable{

    public Input<GeneralSubstitutionModel> modelInput = new Input<>("model", "Beam general tissue substitution model.", Input.Validate.REQUIRED);
    public Input<TissueData> dataTypeInput = new Input<>("dataType", "User data type for the location data to generate more readable logs.", Input.Validate.REQUIRED);

    private int nrOfStates;
    protected GeneralSubstitutionModel model;

    public TissueSubstitutionModelLogger() {}

    @Override
    public void initAndValidate() {
        model = modelInput.get();
        nrOfStates = modelInput.get().getStateCount();
    }

    @Override
    public void init(PrintStream out) {
        String mainID = (getID() == null || getID().matches("\\s*")) ? "tissueSubstModelLogger" : getID();
        String relRatePrefix = mainID + ".relMigRate_";
        TissueData dataType = dataTypeInput.get();

        for (int i=0; i<nrOfStates; i++) {
            String iStr = dataTypeInput.get().getCode(i);
            for (int j=0 ; j<nrOfStates; j++) {
                if (j==i) { continue; }
                String jStr = dataTypeInput.get().getCode(j);
                out.print(relRatePrefix + iStr + "_" + jStr + "\t");
            }
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        model.setupRelativeRates();
        model.setupRateMatrix();
        double[][] rateMatrix = model.getRateMatrix();

        for (int i=0; i<nrOfStates; i++) {
            for (int j=0; j<nrOfStates; j++) {
                if (j==i) { continue; } // Skip logging the diagonal values
                out.print(rateMatrix[i][j] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) {}
}
