package beam.datatype;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.DataType.Base;

/**
 * @author Stephen Staklinski
 */
@Description("This is a data class for tissue labels. It allows specifying a primary tissue and other tissues, " +
        "and automatically creates a code map for the tissues.")
public class TissueData extends Base {

    final public Input<String> primaryTissueInput = new Input<>("primaryTissue", "Primary tissue name", Validate.REQUIRED);
    final public Input<String> otherTissuesInput = new Input<>("otherTissues", "Optional comma-separated list of other tissue names. If not provided, " +
        "then other tissues will be read from the input data.", Validate.OPTIONAL);

    public void initAndValidate(){
        if (otherTissuesInput.get() != null) {
            codeMap = primaryTissueInput.get() + "," + otherTissuesInput.get();
            stateCount = codeMap.split(",").length;
            System.out.println("Internal ordering set based on input otherTissues: " + codeMap);
        } else {
            System.out.println("Only the primary tissue name was provided, so reading other tissue names in from the input data.");
            codeMap = null;
        }
    }

    @Override
    public String getTypeDescription() {
        return "TissueData";
    }

    @Override
    public int[] getStatesForCode(int code) {
        return new int[]{code};
    }

    public void setCodeMap(String newCodeMap) {
        codeMap = newCodeMap;
    }

    public void setStateCount(int newStateCount) {
        stateCount = newStateCount;
    }
}