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
 * @author Stephen Staklinski
 **/

@Description("Logger for BEAM tissue subsitution models.")
public class BeamTissueSubstitutionModelLogger extends BEASTObject implements Loggable{

    public Input<GeneralSubstitutionModel> modelInput = new Input<>("model", "Beam general substitution model.", Input.Validate.REQUIRED);
    public Input<TissueData> dataTypeInput = new Input<>("dataType", "User data type for the location data to generate more readable logs.", Input.Validate.REQUIRED);

    private int nrOfStates;
    protected GeneralSubstitutionModel model;

    public BeamTissueSubstitutionModelLogger() { }

    @Override
    public void initAndValidate() {
        model = modelInput.get();
        nrOfStates = modelInput.get().getStateCount();
    }

    @Override
    public void init(PrintStream out) {
        String mainID = (getID() == null || getID().matches("\\s*")) ? "geoSubstModel" : getID();
        String relRatePrefix = mainID + ".relGeoRate_";
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
        // Logging normalized rates directly from the rate matrix, not the input rate parameters used to setup the rate matrix.
        model.setupRateMatrix();
        double[][] rateMatrix = model.getRateMatrix();

        for (int i=0; i<nrOfStates; i++) {
            for (int j=0; j<nrOfStates; j++) {
                if (j==i) { continue; }
                out.print(rateMatrix[i][j] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) { }
}
