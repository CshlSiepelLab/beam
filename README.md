# Bayesian Evolutionary Analysis of Metastasis (BEAM)

<div style="text-align: left;">
  <img src="logo.jpg" alt="BEAM logo" width="250"/>
</div>

BEAM is a BEAST2 package for Bayesian cancer migration graph inference from CRISPR cell lineage tracing data. BEAM provides a joint model of lineage reconstruction and migration history inference from the raw data, avoiding the need to condition on a single phylogeny while instead inferring a distribution of migration graphs.


## Installation

1. Install BEAST2 from [beast2.org](https://www.beast2.org/)
2. Install BEAGLE from [beagle-dev.github.io](https://beagle-dev.github.io/)
3. Install required BEAST2 packages based on xml declaration. The example uses the following:
   ```bash
   packagemanager -add NS
   packagemanager -add BDSKY
   packagemanager -add BEAST_CLASSIC
   packagemanager -add feast
   ```
4. Install BEAM:
   - Download the latest release from [BEAM GitHub releases](https://github.com/StephenStaklinski/beam/releases)
   - Find your BEAST2 packages directory (check with `packagemanager -list`)
   - Create and set up BEAM in that directory
   ```bash
     mkdir beam && cd beam
     wget https://github.com/StephenStaklinski/beam/releases/download/v0.2.0/beam.v0.2.0.zip
     unzip beam.v0.2.0.zip
   ```

Verify installation by running `packagemanager -list`. Beam should appear in the list.

## Usage

### Basic Usage

Run BEAM using the BEAST2 command line interface:

```bash
beast -working examples/input.xml
```

