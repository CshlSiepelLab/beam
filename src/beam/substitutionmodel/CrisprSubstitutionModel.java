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

    // Store a copy of data, since missing state will be replaced
    private Alignment data;
    // Number of taxa
    private int taxonCount;
    // Number of sites
    private int nrOfSites;
    // Silencing rate parameter
    private RealParameter silencingRate_;
    // Array of max state found per site
    private int[] maxStates;
    // Array of edit rates
    private Double[][] editRates;
    // Stored silencing rate for dirty state checking
    private double storedSilencingRate;
    // Sitewise number of states
    public int[] nrOfStatesPerSite;
    // Sitewise index for missing data state in the rate matrix
    private int[] missingDataStates;

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

        // Read in mutation data and remap mutations to be sequential
        data = dataInput.get();
        taxonCount = data.getTaxonCount();
        nrOfSites = data.getPatternCount();
        reorderMutationsSequentially();

        // Setup the edit rates, either from input or calculated from the data
        initializeEditRates();
        nrOfStatesPerSite = new int[nrOfSites];
        for (int i = 0; i < nrOfSites; i++) {
            nrOfStatesPerSite[i] = editRates[i].length + 2;  // Number of states is number of edits + unedited state (0) + missing data state (-1)
        }

        // Silencing rate is independent of edit rates, so set it up seperately
        silencingRate_ = silencingRateInput.get();  // Rate for the missing data state (-1)
        double silencingRate = silencingRate_.getValue();
        storedSilencingRate = silencingRate;
        if (silencingRate < 0) {
            throw new RuntimeException("Loss rate must be positive!");
        }

        // Missing data state is the last edit in the rate matrix, so remap input -1 characters to that state index
        missingDataStates = new int[nrOfSites];
        for (int i = 0; i < nrOfSites; i++) {
            missingDataStates[i] = nrOfStatesPerSite[i] - 1; // Missing data state index
        }
        replaceMissingDataStates();
    }




    /* Mutations must be sequential to map correctly to rate matrix indices.
    * If not, remap them in input order so editRates[i] corresponds to the i-th mutation 
    * encountered (e.g., first mutation 34 → editRates[0], second 12 → editRates[1]).
    * This remapping can be done globally or sitewise.
    */
    public void reorderMutationsSequentially() {
        if (MODEL_GLOBAL.equals(modelStructureInput.get())) {
            reorderMutationsSequentiallyGlobal();
            return;
        }
        reorderMutationsSequentiallySitewise();
    }

    private void initializeEditRates() {
        maxStates = new int[nrOfSites];
        editRates = new Double[nrOfSites][];
        if (!editRatesInput.get().isEmpty()) {
            if (MODEL_SITEWISE.equals(modelStructureInput.get())) {
                throw new RuntimeException("Sitewise model structure is not compatible with provided " +
                        "edit rates in the current implementation. Please use the global model structure.");
            }
            getProvidedEditRates();
            return;
        } else {
            calculateDynamicEditRates();
        }
        validateAndSumEditRates();
    }

    private void getProvidedEditRates() {
        System.out.println("Using provided edit rates.");
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
        if (MODEL_GLOBAL.equals(modelStructureInput.get())) {
            calculateDynamicEditRatesGlobal();
        } else {
            calculateDynamicEditRatesSitewise();
        }
    }

    private void calculateDynamicEditRatesGlobal() {
        getMaxObservedStatesGlobal();
        if (EDIT_RATE_MODE_UNIFORM.equals(editRateCalculationModeInput.get())) {
            System.out.println("Using global uniform edit rates.");
            createUniformEditRates();
        } else {
            System.out.println("Calculating empirical frequency edit rates from data (global model).");
            collectEmpiricalStateCountsGlobal();
            normalizeCountsGlobal();
        }
    }

    private void calculateDynamicEditRatesSitewise() {
        getMaxObservedStatesSitewise();
        if (EDIT_RATE_MODE_UNIFORM.equals(editRateCalculationModeInput.get())) {
            System.out.println("Using sitewise uniform edit rates.");
            createUniformEditRates();
        } else {
            System.out.println("Calculating empirical frequency edit rates from data (sitewise model).");
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
            if (maxStates[i] <= 0) {
                throw new RuntimeException("No edited states were found in the data for site " + (i + 1) + ".");
            }
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
        // Print edit rates for the user
        if (MODEL_GLOBAL.equals(modelStructureInput.get())) {
            System.out.println("Final edit rates (global model):");
            System.out.println(Arrays.toString(editRates[0]));
        } else {
            System.out.println("Final edit rates (sitewise model):");
            for (int i = 0; i < nrOfSites; i++) {
                System.out.println("  Site " + (i + 1) + ": ");
                System.out.println(Arrays.toString(editRates[i]));
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
        int siteCount = data.getPatternCount();
        boolean isSequential = true;
        for (int site = 0; site < siteCount; site++) {
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
        for (int site = 0; site < siteCount; site++) {
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



    // Required from extension of SubstitutionMode.Base, but not used since the sitewise version is used below
    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {}
    

    public void getTransitionProbabilitiesAllSites(Node node, double startTime, double endTime, double rate, double[][] matrix) {
        // Calculate key parameters that are assumed to be consistent across sites
        double silencingRate = silencingRate_.getValue();
        double delta = (startTime - endTime) * rate;
        double expOfDeltaLoss = Math.exp(-delta * silencingRate);
        double c = expOfDeltaLoss * (1 - Math.exp(-delta));

        int nrOfSitesToCalculate = nrOfSites;
        if (MODEL_GLOBAL.equals(modelStructureInput.get())) {
            nrOfSitesToCalculate = 1;   // Only calculate the first site in the global model since other sites reference the same matrix
        }

        for (int siteNum = 0; siteNum < nrOfSitesToCalculate; siteNum++) {
            // Initialize matrix to zeros
            Arrays.fill(matrix[siteNum], 0.0);
            int nrOfStates = nrOfStatesPerSite[siteNum];
            // Set absorbing state (bottom right corner)
            matrix[siteNum][matrix[siteNum].length - 1] = 1.0;
            // Set diagonal elements and loss probabilities in one pass
            for (int i = 0; i < nrOfStates - 1; i++) {
                int rowOffset = i * nrOfStates;
                matrix[siteNum][rowOffset + i] = (i == 0) ? Math.exp(-delta * (1 + silencingRate)) : expOfDeltaLoss;
                matrix[siteNum][rowOffset + (nrOfStates - 1)] = 1 - expOfDeltaLoss;
            }
            // Set edit probabilities (first row)
            for (int j = 1; j < nrOfStates - 1; j++) {
                double editRate = editRates[siteNum][j - 1];
                matrix[siteNum][j] = editRate * c;
            }
        }

        if (MODEL_GLOBAL.equals(modelStructureInput.get())) {
            // Reference the first site matrix for all other sites in the global model
            for (int siteNum = 1; siteNum < nrOfSites; siteNum++) {
                matrix[siteNum] = matrix[0];
            }
        }
    }

    @Override
    public void store() {
        storedSilencingRate = silencingRate_.getValue();
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        // Since edit rates are fixed hyperparameters, the only thing that can change is the silencing rate
        return silencingRate_.getValue() != storedSilencingRate;
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
}
