package beam.substitutionmodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.IntegerData;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;

/**
 * CRISPR substitution models based on TideTree (global model across sites)
 * and LAML (sitewise model for heterogeneous sites).
 *
 * @author Stephen Staklinski
 */
@Description("CRISPR substitution models.")
public class CrisprSubstitutionModel extends SubstitutionModel.Base {
    private static final String MODEL_GLOBAL = "global";
    private static final String MODEL_SITEWISE = "sitewise";
    private static final String EDIT_RATE_MODE_EMPIRICAL = "empirical";
    private static final String EDIT_RATE_MODE_UNIFORM = "uniform";

    // Input for the silencing rate throughout the experiment
    public final Input<RealParameter> silencingRateInput = new Input<>("silencingRate", "Rate at which barcodes " +
            "are silenced throughout the entire experiment", Validate.REQUIRED);
    // Input for edit rates during the editing window
    public final Input<List<RealParameter>> editRatesInput = new Input<>("editRates", "Input rates at which edits " +
            "are introduced into the genomic barcode during the editing window", new ArrayList<>(), Validate.OPTIONAL);
    // Input for edit rates during the editing window
    public final Input<Alignment> dataInput = new Input<>("data", "EditData from which to calculate edit rates if " +
            "they are not provided", Validate.REQUIRED);
    // Input mode for calculating edit rates when they are not provided
    public final Input<String> editRateCalculationModeInput = new Input<>("editRateCalculationMode", "Method for " +
            "calculating edit rates from the data if they are not provided. Options are 'empirical' to calculate " +
            "edit rates proportional to the empirical frequency of each edit in the data, or 'uniform' to give " +
            "equal rates to all edits. Default is 'empirical'.", "empirical", Validate.OPTIONAL);
    // Input model structure for rate matrices
    public final Input<String> modelStructureInput = new Input<>("modelStructure", "Whether to use one rate matrix " +
            "for all sites (global) or a separate rate matrix for each site (sitewise). Default is global.", "global", 
            Validate.OPTIONAL);
    // Whether to print edit rate details out during initialiation
    public final Input<Boolean> printRates = new Input<>("printRates", "Whether to print " +
            "details about the intialization of edit rates. Default is false.", false, Validate.OPTIONAL);

    // Store a copy of data, since missing state will be replaced
    private Alignment data;
    // Number of taxa
    private int taxonCount;
    // Stored silencing rate
    private double storedSilencingRate;
    // Number of sites
    private int nrOfSites;
    // Array of max state found per site
    private int[] maxStates;
    // Array of edit rates
    private Double[][] editRates;
    // Sitewise number of states
    public int[] nrOfStatesPerSite;
    // Sitewise index for missing data state in the rate matrix
    private int[] missingDataStates;
    // Save model structure input for easy access
    private boolean is_global;
    // Save edit rate calculation mode input for easy access
    private boolean is_uniform;
    // Whether to print edit rate details out during initialiation
    private boolean printRates_;

    public CrisprSubstitutionModel() {
        // This model sets root frequencies internally (unedited state fixed to 1.0
        // in likelihood core), so make optional to avoid forcing unused XML config.
        frequenciesInput.setRule(Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        // Ensure valid inputs
        if (!MODEL_GLOBAL.equals(modelStructureInput.get()) && !MODEL_SITEWISE.equals(modelStructureInput.get())) {
            throw new RuntimeException("Invalid model structure: " + modelStructureInput.get() + ". Must be '" + MODEL_GLOBAL + "' or '" + MODEL_SITEWISE + "'.");
        }
        if (!EDIT_RATE_MODE_EMPIRICAL.equals(editRateCalculationModeInput.get()) && !EDIT_RATE_MODE_UNIFORM.equals(editRateCalculationModeInput.get())) {
            throw new RuntimeException("Invalid edit rate calculation mode: " + editRateCalculationModeInput.get() + ". Must be '" + EDIT_RATE_MODE_EMPIRICAL + "' or '" + EDIT_RATE_MODE_UNIFORM + "'.");
        }

        // Store model setup based on inputs
        is_global = MODEL_GLOBAL.equals(modelStructureInput.get());
        is_uniform = EDIT_RATE_MODE_UNIFORM.equals(editRateCalculationModeInput.get());
        printRates_ = printRates.get();

        // Read in mutation data and remap mutations to be sequential
        data = dataInput.get();
        taxonCount = data.getTaxonCount();
        nrOfSites = data.getSiteCount(); // Using all sites, not just unique patterns here
        reorderMutationsSequentially();

        // Setup the edit rates, either from input or calculated from the data
        initializeEditRates();
        nrOfStatesPerSite = new int[nrOfSites];
        for (int i = 0; i < nrOfSites; i++) {
            nrOfStatesPerSite[i] = maxStates[i] + 2;  // Number of states is number of edits + unedited state (0) + missing data state (-1)
        }

        // Missing data state is the last edit in the rate matrix, so remap input -1 characters to that state index
        missingDataStates = new int[nrOfSites];
        for (int i = 0; i < nrOfSites; i++) {
            missingDataStates[i] = nrOfStatesPerSite[i] - 1; // Missing data state index
        }
        replaceMissingDataStates();

        storedSilencingRate = silencingRateInput.get().getValue();
    }




    /* Mutations must be sequential to map correctly to rate matrix indices.
    * If not, remap them in input order so editRates[i] corresponds to the i-th mutation 
    * encountered (e.g., first mutation 34 → editRates[0], second 12 → editRates[1]).
    * This remapping can be done globally or sitewise.
    */
    public void reorderMutationsSequentially() {
        if (is_global) {
            reorderMutationsSequentiallyGlobal();
            return;
        }
        reorderMutationsSequentiallySitewise();
    }

    private void initializeEditRates() {
        maxStates = new int[nrOfSites];
        editRates = new Double[nrOfSites][];
        if (!editRatesInput.get().isEmpty()) {
            throw new RuntimeException("TODO: Not sure that the provided edit rates is setup currently in terms of mcmc store/restore. Will work in this soon.");
            // if (!is_global) {
            //     throw new RuntimeException("Sitewise model structure is not compatible with provided " +
            //             "edit rates in the current implementation. Please use the global model structure.");
            // }
            // getProvidedEditRates();
            // return;
        } else {
            calculateDynamicEditRates();
        }
        validateAndSumEditRates();
    }

    private void getProvidedEditRates() {
        if (printRates_) System.out.println("Using provided edit rates.");
        RealParameter editRate_ = editRatesInput.get().get(0);
        editRates[0] = new Double[editRate_.getDimension()];
        for (int i = 0; i < editRate_.getDimension(); i++) {
            editRates[0][i] = editRate_.getValue(i);
        }
        // Reuse the same array reference for all sites since we are using the global model
        for (int i = 1; i < nrOfSites; i++) {
            editRates[i] = editRates[0];
        }
    }

    private void calculateDynamicEditRates() {
        if (is_global) {
            calculateDynamicEditRatesGlobal();
        } else {
            calculateDynamicEditRatesSitewise();
        }
    }

    private void calculateDynamicEditRatesGlobal() {
        getMaxObservedStatesGlobal();
        if (is_uniform) {
            if (printRates_) System.out.println("Using global uniform edit rates.");
            createUniformEditRates();
        } else {
            if (printRates_) System.out.println("Calculating empirical frequency edit rates from data (global model).");
            collectEmpiricalStateCountsGlobal();
            normalizeCountsGlobal();
        }
    }

    private void calculateDynamicEditRatesSitewise() {
        getMaxObservedStatesSitewise();
        if (is_uniform) {
            if (printRates_) System.out.println("Using sitewise uniform edit rates.");
            createUniformEditRates();
        } else {
            if (printRates_) System.out.println("Calculating empirical frequency edit rates from data (sitewise model).");
            collectEmpiricalStateCountsSitewise();
            normalizeCountsSitewise();
        }
    }

    private void getMaxObservedStatesGlobal() {
        maxStates[0] = 0;
        for (int taxon = 0; taxon < taxonCount; taxon++) {
            List<Integer> seq = data.getCounts().get(taxon);
            for (int value : seq) {
                maxStates[0] = Math.max(maxStates[0], value);
            }
        }
        if (maxStates[0] <= 0) throw new RuntimeException("No edits observed in the data! Cannot calculate edit rates. Please check your input data.");
        // Other sites can just be copies of the same max state since we are using the global model
        for (int i = 1; i < nrOfSites; i++) {
            maxStates[i] = maxStates[0];
        }
    }

    private void getMaxObservedStatesSitewise() {
        for (int i = 0; i < nrOfSites; i++) {
            maxStates[i] = 0;
            for (int taxon = 0; taxon < taxonCount; taxon++) {
                maxStates[i] = Math.max(maxStates[i], data.getCounts().get(taxon).get(i));
            }
        }
    }

    private void collectEmpiricalStateCountsGlobal() {
        editRates[0] = new Double[maxStates[0]];
        Arrays.fill(editRates[0], 0.0);
        for (int taxon = 0; taxon < taxonCount; taxon++) {
            List<Integer> seq = data.getCounts().get(taxon);
            for (int value : seq) {
                if (value <= 0) continue;   // Skip unedited state (0) and missing data state (-1)
                editRates[0][value - 1] += 1;
            }
        }
        // Other sites can just be references
        for (int i = 1; i < nrOfSites; i++) {
            editRates[i] = editRates[0];
        }
    }

    private void collectEmpiricalStateCountsSitewise() {
        for (int i = 0; i < nrOfSites; i++) {
            if (maxStates[i] <= 0) continue;  // No edits observed at this site, so leave edit rates empty since they will never be used
            editRates[i] = new Double[maxStates[i]];
            Arrays.fill(editRates[i], 0.0);
            for (int taxon = 0; taxon < taxonCount; taxon++) {
                int value = data.getCounts().get(taxon).get(i);
                if (value <= 0) continue;   // Skip unedited state (0) and missing data state (-1)
                editRates[i][value - 1] += 1;
            }
        }
    }

    private void createUniformEditRates() {
        for (int i = 0; i < nrOfSites; i++) {
            if (maxStates[i] <= 0) continue;  // No edits observed at this site, so leave edit rates empty since they will never be used
            editRates[i] = new Double[maxStates[i]];
            Arrays.fill(editRates[i], 1.0 / maxStates[i]);
        }
    }

    private void normalizeCountsGlobal() {
        double total = 0.0;
        // Only count and normalize the first site since we are using the global model and other sites reference the same array
        for (int j = 0; j < maxStates[0]; j++) {
            total += editRates[0][j];
        }
        for (int j = 0; j < maxStates[0]; j++) {
            editRates[0][j] /= total;
        }
    }

    private void normalizeCountsSitewise() {
        // Normalize all sites independently since we are using the sitewise model
        for (int i = 0; i < nrOfSites; i++) {
            if (maxStates[i] <= 0) continue;  // No edits observed at this site, so skip normalization since edit rates will never be used
            double total = 0.0;
            for (int j = 0; j < maxStates[i]; j++) {
                total += editRates[i][j];
            }
            for (int j = 0; j < maxStates[i]; j++) {
                editRates[i][j] /= total;
            }
        }
    }

    private void validateAndSumEditRates() {
        for (int i = 0; i < nrOfSites; i++) {
            if (maxStates[i] <= 0) continue;  // No edits observed at this site, so skip validation since edit rates will never be used
            double sum = 0.0;
            for (double editRate : editRates[i]) {
                if (editRate < 0) {
                    throw new RuntimeException("All edit rates must be positive!");
                }
                sum += editRate;
            }
            if (Math.abs(sum - 1.0) > 1e-6) {
                throw new RuntimeException("Edit rates must sum to 1! Check input edit rates or internal calculations.");
            }
        }
        
        if (printRates_) {
            // Print edit rates for the user
            if (is_global) {
                System.out.println("Final edit rates (global model):");
                System.out.println(Arrays.toString(editRates[0]));
            } else {
                System.out.println("Final edit rates (sitewise model):");
                for (int i = 0; i < nrOfSites; i++) {
                    if (maxStates[i] <= 0) {
                        System.out.println("  Site " + (i + 1) + ": No edits observed, so no edit rates.");
                        continue; // No edits observed at this site, so do not print the array
                    }
                    System.out.println("  Site " + (i + 1) + ": ");
                    System.out.println(Arrays.toString(editRates[i]));
                }
            }
        }
    }

    private void replaceMissingDataStates() {
        // Overwrite the data to replace -1 states with the missing data state for proper use in the mutation model rate matrix
        for (int taxon = 0; taxon < taxonCount; taxon++) {
            for (int i = 0; i < nrOfSites; i++) {
                if (data.getCounts().get(taxon).get(i) == -1) {
                    data.getCounts().get(taxon).set(i, missingDataStates[i]); // Replace missing data state (-1) with the index for the missing data state in the rate matrix
                }
            }
        }
    }

    private void reorderMutationsSequentiallyGlobal() {
        boolean isSequential = true;
        int maxValSeen = 0;
        for (int taxon = 0; taxon < taxonCount; taxon++) {
            List<Integer> seq = data.getCounts().get(taxon);
            for (int value : seq) {
                if (value <= 0) continue;
                if (value == maxValSeen + 1) maxValSeen = value;
                else if (value > maxValSeen + 1) {
                    isSequential = false;
                    break;
                }
            }
            if (!isSequential) break;
        }

        if (isSequential) {
            return;
        }

        // Remap globally in input order so mutations become 1..K.
        Map<Integer, Integer> map = new HashMap<>();
        int nextSequentialMut = 1;
        for (int taxon = 0; taxon < taxonCount; taxon++) {
            List<Integer> seq = data.getCounts().get(taxon);
            for (int i = 0; i < seq.size(); i++) {
                int inputMut = seq.get(i);
                if (inputMut <= 0) continue;
                Integer sequentialMut = map.get(inputMut);
                if (sequentialMut == null) {
                    sequentialMut = nextSequentialMut++;
                    map.put(inputMut, sequentialMut);
                }
                seq.set(i, sequentialMut);
            }
        }
    }

    private void reorderMutationsSequentiallySitewise() {
        boolean isSequential = true;
        for (int site = 0; site < nrOfSites; site++) {
            int maxValSeen = 0;
            for (int taxon = 0; taxon < taxonCount; taxon++) {
                int value = data.getCounts().get(taxon).get(site);
                if (value <= 0) continue;
                if (value == maxValSeen + 1) maxValSeen = value;
                else if (value > maxValSeen + 1) {
                    isSequential = false;
                    break;
                }
            }
            if (!isSequential) break;
        }

        if (isSequential) {
            return;
        }

        // Remap each site independently in input order so mutations become 1..K per site.
        for (int site = 0; site < nrOfSites; site++) {
            Map<Integer, Integer> map = new HashMap<>();
            int nextSequentialMut = 1;
            for (int taxon = 0; taxon < taxonCount; taxon++) {
                int inputMut = data.getCounts().get(taxon).get(site);
                if (inputMut <= 0) continue;
                Integer sequentialMut = map.get(inputMut);
                if (sequentialMut == null) {
                    sequentialMut = nextSequentialMut++;
                    map.put(inputMut, sequentialMut);
                }
                data.getCounts().get(taxon).set(site, sequentialMut);
            }
        }
    }
    

    public double[] getCoreTransitionProbabilityValues(Node node, double startTime, double endTime, double rate) {
        // Calculate key parameters that are assumed to be consistent across sites
        double silencingRate = silencingRateInput.get().getValue();
        double delta = (startTime - endTime) * rate;
        double expOfDeltaLoss = Math.exp(-delta * silencingRate);

        double[] coreValues = new double[4];
        coreValues[0] = Math.exp(-delta * (1.0 + silencingRate));  // Top left/first state, in the matrix; Prob of unedited -> unedited (0 to 0)
        coreValues[1] = expOfDeltaLoss;             // Diagonal elements, except for first element in top row and last element in the final row; Prob of staying in the same edited state (i to i) from an edited starting state
        coreValues[2] = 1.0 - expOfDeltaLoss;           // Last column except for the bottom right/final entry; Prob of transitioning to the missing data state (i to -1) from any i>=0 (unedited or edited states)
        coreValues[3] = expOfDeltaLoss * (1.0 - Math.exp(-delta));           // First row edit probability scalar; Edit rate * scalar give the prob of editing from unedited to the specific edit state (0 to j) for j>0
        return coreValues;
    }

    @Override
    public void store() {
        storedSilencingRate = silencingRateInput.get().getValue();
        super.store();
    }

    @Override
    public boolean requiresRecalculation() {
        if (silencingRateInput.get().getValue() != storedSilencingRate) return true;
        return super.requiresRecalculation();
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        // Return true for BEAGLE compatibility
        return true;
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return null;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof IntegerData;
    }

    public int[] getMissingStates() {
        return missingDataStates;
    }

    public Alignment getData() {
        return data;
    }

    public int getSiteCount() {
        return nrOfSites;
    }

    public int[] getNrOfStatesPerSite() {
        return nrOfStatesPerSite;
    }

    public double getEditRate(int site, int editIndex) {
        return editRates[site][editIndex - 1];  // Edit index is 1-based for input but 0-base in the array
    }

    public double getSilencingRate() {
        return silencingRateInput.get().getValue();
    }

    // Required, but not used
    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {}
}
