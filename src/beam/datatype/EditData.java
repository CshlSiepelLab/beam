package beam.datatype;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType.Base;

/**
 * @author Stephen Staklinski
 */
@Description("This is based on the EditData class from Sophie Seidel's tidetree. " +
        "It is a datatype for integer sequence representing the barcode states. " +
        "The state 0 is unedited, states 1, ..., N-1 represent particular indels, " +
        "state N represents the silenced state.")
public class EditData extends Base {

    public void initAndValidate(){
    }

    @Override
    public String getTypeDescription() {
        return "editData";
    }

    @Override
    public int[] getStatesForCode(int code) {
        return new int[]{code};
    }
}