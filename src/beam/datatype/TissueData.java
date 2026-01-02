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

    final public Input<String> primaryTissueInput = new Input<>("primaryTissue", "the name of the character", Validate.REQUIRED);
    final public Input<String> otherTissuesInput = new Input<>("otherTissues", "comma-separated list of other tissues", Validate.REQUIRED);

    public void initAndValidate(){
        codeMap = primaryTissueInput.get() + "," + otherTissuesInput.get();
        stateCount = codeMap.split(",").length;
    }

    @Override
    public String getTypeDescription() {
        return "tissueData";
    }

    @Override
    public int[] getStatesForCode(int code) {
        return new int[]{code};
    }
}