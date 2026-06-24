# SQOUnified R Package

The `SQOUnified` R package provides a suite of tools for assessing sediment quality based on multiple lines of evidence, including benthic indices, chemistry, and toxicity. The package integrates these assessments into a comprehensive site evaluation according to Sediment Quality Objectives (SQO).

This is based on the [California Sediment Quality Objectives (CASQO) Technical Manual](http://ftp.sccwrp.org/pub/download/DOCUMENTS/TechnicalReports/777_CASQO_TechnicalManual.pdf)

## Installation

To install the `SQOUnified` package, you can use the `devtools` package. Ensure all dependencies are installed beforehand.

```r
# Install necessary dependencies
install.packages(c("dbplyr","reshape2","vegan","dplyr","plyr","purrr","stringr","tidyr"))

# Install the package using devtools
devtools::install_github("SCCWRP/SQOUnified")
```

The package also comes with example data sets which are available when the package is loaded
```r
# load the package
library(SQOUnified)

# The package ships with sample data for each line of evidence:
#   benthic_sampledata           - bay/estuary benthic infauna
#   chem_sampledata              - sediment chemistry
#   tox_sampledata               - sediment toxicity
#   offshore_benthic_sampledata  - offshore (shelf) benthic infauna  } used together for the
#   offshore_station_sampledata  - offshore station info (depth, lat/long) } offshore BRI line of evidence
SQOUnified(benthic = benthic_sampledata, chem = chem_sampledata, tox = tox_sampledata)

# Include the offshore benthic line of evidence (offshore BRI) by also passing the
# offshore benthic infauna and the matching offshore station info
SQOUnified(
  benthic = benthic_sampledata,
  chem = chem_sampledata,
  tox = tox_sampledata,
  offshore_benthic = offshore_benthic_sampledata,
  offshore_stations = offshore_station_sampledata
)

### Other misc examples:
# Benthic integrated SQO (returns all individual indices as well)
benthic.sqo(benthic_sampledata)

# Offshore Benthic Response Index (Edition 14)
BRI.Offshore(offshore_benthic_sampledata, offshore_station_sampledata)

# Tox summary
tox.summary(tox_sampledata)

# Tox integrated SQO Score
tox.sqo(tox_sampledata)

# Chemical Score Index
CSI(chem_sampledata)

# Chemistry Integrated SQO Score
chem.sqo(chem_sampledata)
```

In addition to the small, self-contained `*_sampledata` sets above, the package also ships
larger compiled Southern California Bight regional monitoring datasets (`bight.unified.chem`,
`bight.unified.tox`, `bight.unified.benthic`, and `bight.unified.benthic.offshore`). These are
real, multi-survey-year data and can be used as more realistic inputs to the line-of-evidence
functions. See the [Example Datasets](#example-datasets) section below for details.

## Usage

### Input Data Requirements for Benthic, Chemistry, and Toxicity Functions

#### 1. Input Data for benthic functions (`benthic.sqo`, `BRI`,`RBI`,`IBI`,`RIVPACS`,`MAMBI`)

The benthic functions require a dataframe with the following columns:

- **StationID**: An alphanumeric identifier representing the sampling location.
- **Replicate**: A numeric value indicating the replicate sample number taken at the location.
- **SampleDate**: The date when the sample was collected.
- **Latitude**: Latitude in decimal degrees; use a negative sign for Western Hemisphere coordinates.
- **Longitude**: Longitude in decimal degrees; ensure a negative sign for Western coordinates.
- **Taxon**: The name of the organism, formatted according to SCAMIT standards. If no organisms are present, use `NoOrganismsPresent` with an abundance of `0`.
- **Abundance**: The number of each species observed in the sample.
- **Salinity**: The salinity level (in PSU) at the sampling location.
- **Stratum**: The habitat stratum (e.g., Bays, Estuaries).
- **Exclude**: Information on whether the sample should be excluded.


#### 2. Input Data for chemistry functions (`chem.sqo`, `chemdata_prep`, `CSI`, `LRM`)

The chemistry functions require a dataframe containing:

- **StationID**: Identifier for the station.
- **AnalyteName**: Name of the chemical analyte.
- **Result**: The measurement value for the analyte.
- **RL**: Reporting limit for the analyte.
- **MDL**: Method detection limit.
- **units**: (*optional*) - Metals should be in mg/dry kg (mg/kg dw) and all organic constituents should be in ug/dry kg (ug/kg dw).
- **fieldrep**: (*optional*) - Data is filtered to where fieldrep = 1 if this column is included
- **labrep**: (*optional*) - Data is filtered to where labrep = 1 if this column is included
- **sampletypecode**: (*optional*) Data is filtered to where sampletypecode = Result in order to prevent inclusion of QA/QC samples

**Notes**:
- Non-detect values should be marked as `-88`.
- All measurements should be expressed on a dry-weight basis: metals in `mg/dry kg` and organic compounds in `ug/dry kg`.


#### 3. Input Data for the offshore benthic function (`BRI.Offshore`)

The offshore Benthic Response Index requires depth in addition to the standard benthic columns, because tolerance values are applied per depth zone. Column names are case-insensitive. The data can be supplied either as a single dataframe that already contains all of these columns, or as two dataframes (benthic infauna + station info) that share `StationID` (this is how the bundled `offshore_benthic_sampledata` / `offshore_station_sampledata` are organized):

- **StationID**: An alphanumeric identifier representing the sampling location.
- **SampleDate**: The date when the sample was collected.
- **Replicate**: A numeric value indicating the replicate sample number.
- **Taxon**: The name of the organism, in SCAMIT Edition 14 naming conventions. If no organisms are present, use `NoOrganismsPresent` with an abundance of `0`.
- **Abundance**: The number of individuals counted for each taxon.
- **Depth**: Station depth in meters. Determines which depth-zone tolerance values are applied (shallow `<25m`, mid `35-110m`, deep `130-324m`, with overlap zones at `25-35m` and `110-130m`). Samples deeper than `324m` are flagged "BRI not Applicable".
- **Latitude**: Latitude in decimal degrees.
- **Longitude**: Longitude in decimal degrees; use a negative sign for Western Hemisphere coordinates.

**Note:** The BRI is calibrated for the Southern California Bight. Samples north of Point Conception (`lat > 34.45`) or south of the US-Mexico border (`lat < 32.52`) are retained but flagged with a geographic-range caution.


#### 4. Input Data for toxicity functions (`tox.sqo`, `tox.summary`)

The toxicity functions require a dataframe containing:

- **stationid**: Alphanumeric identifier of the location.
- **toxbatch**: Identifier linking results to the control sample.
- **species**: The genus and species of the tested organism.
- **sampletypecode**: Type of sample (e.g., Grab, CNEG). Control samples must be included.
- **matrix**: (*optional*) - Type of sample matrix (e.g., Whole Sediment, Sediment Water Interface). Be sure to not include Reference Toxicants
- **labrep**: Numeric identifier for the lab replicate, typically with five replicates per station and species pair.
- **result**: The percentage of survival or normal development observed in the test.



## All Package Functions

Below we will list all the functions of the SQOUnified package, grouped into 4 sections
- SQOUnified
- Toxicity
- Chemistry
- Benthic

### `SQOUnified` Function

The primary function `SQOUnified` computes integrated sediment quality scores based on the available lines of evidence: benthic, chemistry, and toxicity data. Every argument is optional and defaults to `NULL` - the function calculates scores for whatever data you provide, and produces an integrated site assessment for any station that has all three lines of evidence (Benthic, Chemistry, and Toxicity).

The benthic line of evidence can come from **either** bay/estuary benthic data (via `benthic`) **or** offshore (shelf) benthic data (via `offshore_benthic` + `offshore_stations`). The offshore path runs the offshore Benthic Response Index (`BRI.Offshore`) and maps its condition categories onto the standard benthic categories, so both feed the same integrated assessment.

**Example Usage:**

```r
# Load the package
library(SQOUnified)

# Example data files
benthic_data <- read.csv("path/to/your/benthic_data.csv")
chem_data <- read.csv("path/to/your/chem_data.csv")
tox_data <- read.csv("path/to/your/tox_data.csv")

# Run the function (bay/estuary benthic + chemistry + toxicity)
result <- SQOUnified(benthic = benthic_data, chem = chem_data, tox = tox_data)
print(result)

# To use offshore benthic data instead of (or in addition to) bay/estuary benthic data,
# supply the offshore infauna and its matching station info (which carries depth/lat/long)
offshore_benthic_data <- read.csv("path/to/your/offshore_benthic_data.csv")
offshore_station_data <- read.csv("path/to/your/offshore_station_data.csv")

result <- SQOUnified(
  chem = chem_data,
  tox = tox_data,
  offshore_benthic = offshore_benthic_data,
  offshore_stations = offshore_station_data
)
```

#### Parameters
- `benthic`: A dataframe containing bay/estuary benthic data for assessment.
- `chem`: A dataframe containing chemical concentration data.
- `tox`: A dataframe containing toxicity test results.
- `offshore_benthic`: (*optional*) A dataframe of offshore (shelf) benthic infauna for the offshore BRI line of evidence. Must be supplied together with `offshore_stations`.
- `offshore_stations`: (*optional*) A dataframe of offshore station info (`StationID`, `Depth`, `Latitude`, `Longitude`, `SampleDate`), joined to `offshore_benthic` on `StationID`.
- `logfile`: (*optional*) Path to the log file. Defaults to a timestamped `.Rmd` file under a `logs/` directory in the working directory.
- `verbose`: (*optional*) Logical; if `TRUE`, detailed intermediate tables are written to the log. Default `FALSE`.
- `logtitle`: (*optional*) Title used for the generated log. Default `'Unified SQO Logs'`.
- `knitlog`: (*optional*) Logical; if `TRUE`, the consolidated log is knit to a single self-contained HTML report. Default `FALSE`.

The function computes and logs the results for each input, generating an integrated score based on the criteria defined in the package. The output dataframe also includes all individual scores and indices:
- Integrated Chemistry Score
- CSI
- LRM
- Integrated Toxicity Score
- (Results from all the tox endpoint tests in the data, typically 'Mytilus NormDev' and 'Eohaustorius Survival')
- Integrated Benthic Score
- BRI
- RBI
- IBI
- RIVPACS
- MAMBI (Not used for calculation of SQO Score)
- Offshore BRI (when offshore data is provided; contributes to the Benthic line of evidence in the integrated assessment)



### Toxicity Functions

#### `tox.sqo`

The `tox.sqo` function assesses sediment toxicity by analyzing toxicity test results.

**Usage:**

```r
tox_results <- tox.sqo(tox_data)
```

- **Input:** `tox_data` - A dataframe containing toxicity test results.
- **Output:** A dataframe with SQO scores and categories for each station based on toxicity criteria.

#### `tox.summary`

The `tox.summary` function provides a summary of the toxicity results for easier interpretation and reporting.

**Usage:**

```r
tox_summary <- tox.summary(tox_data)
```

- **Input:** `tox_data` - A dataframe containing toxicity test results.
- **Output:** A summary dataframe providing an overview of toxicity conditions for all stations.


### Chemistry Functions

#### `chem.sqo`

The `chem.sqo` function evaluates the chemical condition of sediments using specific chemical metrics.

**Usage:**

```r
chem_results <- chem.sqo(chem_data)
```

- **Input:** `chem_data` - A dataframe containing sediment chemistry data.
- **Output:** A dataframe with SQO scores and categories for each station based on chemical criteria.

#### `chemdata_prep`

The `chemdata_prep` function reformats/prepares the data for getting the SQO score for the chemistry LOE

**Usage:**

```r
chem <- chemdata_prep(chem_data)
```

- **Input:** `chem_data` - A dataframe containing sediment chemistry data.
- **Output:** A dataframe prepared for SQO calculation for the chemistry LOE

**NOTE:** When you put your raw chemistry data in the chemistry functions, this chemdata_prep function is called from within, so no need to preprocess the data before so long as it has those required columns



#### `CSI`

The `CSI` function calculates the Chemical Score Index for the input chemistry data

**Usage:**

```r
csi.results <- CSI(chem_data)
```

- **Input:** `CSI` - A dataframe containing sediment chemistry data.
- **Output:** A dataframe with the Chemical Score Index results


#### `LRM`

The `LRM` function calculates the Chemical Score Index for the input chemistry data

**Usage:**

```r
lrm.results <- LRM(chem_data)
```

- **Input:** `LRM` - A dataframe containing sediment chemistry data.
- **Output:** A dataframe with the Logistic Regression Model results





### Benthic Functions
#### `benthic.sqo`

The `benthic.sqo` function computes benthic indices to evaluate the biological condition of the sediments.

**Usage:**

```r
benthic_results <- benthic.sqo(benthic_data)
```

- **Input:** `benthic_data` - A dataframe containing benthic community data.
- **Output:** A dataframe with SQO scores and categories for each station based on benthic criteria. It also contains the individual benthic index scores (BRI, RBI, IBI, RIVPACS and MAMBI)

**NOTE:** It will be the same for each of the benthic functions - BRI, RBI, IBI, RIVPACS and MAMBI


#### `BRI.Offshore`

The `BRI.Offshore` function computes the offshore Benthic Response Index using depth-zone-specific pollution tolerance values (p-codes) from SCAMIT Edition 14. Unlike the bay/estuary `BRI`, it accounts for depth-dependent community composition by applying separate tolerance scores in shallow (`<25m`), mid (`35-110m`), and deep (`130-324m`) depth zones; samples in the overlap ranges (`25-35m` and `110-130m`) receive the average of the two adjacent zone scores.

**Usage:**

```r
# The offshore infauna and the station info (depth/lat/long) can be passed as two dataframes
offshore_results <- BRI.Offshore(offshore_benthic_sampledata, offshore_station_sampledata)

# ...or as a single dataframe that already contains all required columns
offshore_results <- BRI.Offshore(combined_offshore_data)

# output_format controls the shape of the scores table: 'wide' (default) returns one row per
# sample; 'long' returns a tidy index/score/category layout (the form SQOUnified consumes)
offshore_results <- BRI.Offshore(combined_offshore_data, output_format = 'long')
```

- **Input:** Offshore benthic infauna with station info (see *Input Data for the offshore benthic function* above). Required columns: `StationID`, `SampleDate`, `Replicate`, `Taxon`, `Abundance`, `Depth`, `Latitude`, `Longitude`.
- **Output:** A named `list`. The main element, `bri_scores`, holds the BRI score, condition category (Reference, Marginal Deviation, Biodiversity Loss, Function Loss, Defaunation), class, and usage notes for each sample. The list also returns the intermediate tables so you can review how your data were processed:
  - `taxa_with_pcode` - submitted taxa that matched a p-code
  - `taxa_without_pcode` - submitted taxa with no p-code (with fuzzy-matched "did you mean" suggestions when the optional `fuzzyjoin` package is installed)
  - `all_taxa_by_sample` - every taxon joined to its p-code and depth-zone tolerance values, by sample

**NOTE:** When `SQOUnified` is called with `offshore_benthic` and `offshore_stations`, it calls `BRI.Offshore` internally (in `output_format = 'long'`), so there is no need to run this beforehand to feed the integrated assessment.



## Example Workflow

Here’s how you can use the `SQOUnified` package to perform a complete sediment quality assessment:

```r
# Load sample data
benthic_data <- read.csv("path/to/benthic_data.csv")
chem_data <- read.csv("path/to/chem_data.csv")
tox_data <- read.csv("path/to/tox_data.csv")

# Calculate individual SQO components, if you would like
benthic_results <- benthic.sqo(benthic_data)
chem_results <- chem.sqo(chem_data)
tox_results <- tox.sqo(tox_data)

# Get the tox summary table, on its own, if you would like
tox_summary <- tox.summary(tox_data)

# Compute the integrated SQO scores - a dataframe with all scores, including the overall site assessments
integrated_result <- SQOUnified(benthic = benthic_data, chem = chem_data, tox = tox_data)

# View the results
print(integrated_result)
```



## Example Datasets

The package bundles two tiers of example data, all lazy-loaded and available by name once the package is loaded (`library(SQOUnified)`).

### Small, self-contained sample data

These are compact, synthetic datasets designed to exercise each function quickly. They span a range of contamination/condition levels (a clean reference site through a heavily disturbed site).

| Dataset | Description | Used with |
|---|---|---|
| `benthic_sampledata` | Bay/estuary benthic infauna (one row per taxon per sample) | `benthic.sqo`, `BRI`, `RBI`, `IBI`, `RIVPACS`, `MAMBI`, `SQOUnified` |
| `chem_sampledata` | Sediment chemistry (long format, one row per analyte per station) | `chem.sqo`, `CSI`, `LRM`, `SQOUnified` |
| `tox_sampledata` | Sediment toxicity test results, including control samples | `tox.sqo`, `tox.summary`, `SQOUnified` |
| `offshore_benthic_sampledata` | Offshore (shelf) benthic infauna across a range of depths | `BRI.Offshore`, `SQOUnified` |
| `offshore_station_sampledata` | Offshore station info (`stationid`, `sampledate`, `latitude`, `longitude`, `depth`) | `BRI.Offshore`, `SQOUnified` |

#### About the offshore sample data

The offshore benthic line of evidence is split across **two** dataframes that are joined on `stationid`:

- `offshore_benthic_sampledata` holds the infauna - `stationid`, `replicate`, `sampledate`, `taxon`, `abundance`.
- `offshore_station_sampledata` holds the per-station environmental info needed for depth-zone scoring - `stationid`, `sampledate`, `latitude`, `longitude`, `depth`.

They are kept separate (rather than as one wide table) because the offshore BRI needs station `depth` to choose the correct depth-zone tolerance values, and because this mirrors how the data are typically stored. Both `BRI.Offshore` and `SQOUnified` accept them as two arguments and join them internally:

```r
# Standalone offshore BRI
BRI.Offshore(offshore_benthic_sampledata, offshore_station_sampledata)

# As the benthic line of evidence in the full integrated assessment
SQOUnified(
  chem = chem_sampledata,
  tox = tox_sampledata,
  offshore_benthic = offshore_benthic_sampledata,
  offshore_stations = offshore_station_sampledata
)
```

The example stations cover shallow, mid, and deep depth zones (as well as a heavily disturbed shallow site) so the depth-dependent behavior of `BRI.Offshore` can be observed.

### Compiled Bight regional monitoring data

The package also ships larger, real datasets compiled across multiple Southern California Bight regional monitoring survey years (1998-2023). These share a consistent schema and can be used as realistic inputs to the line-of-evidence functions.

| Dataset | Description | Used with |
|---|---|---|
| `bight.unified.chem` | Compiled sediment chemistry (long format) | `chem.sqo`, `SQOUnified` |
| `bight.unified.tox` | Compiled sediment toxicity results (with controls) | `tox.sqo`, `SQOUnified` |
| `bight.unified.benthic` | Compiled bay/estuary benthic infauna (with full WoRMS taxonomy) | `benthic.sqo`, `SQOUnified` |
| `bight.unified.benthic.offshore` | Compiled offshore (shelf) benthic infauna | `BRI.Offshore`, `SQOUnified` |

Each dataset is documented in the package help; use `?bight.unified.chem` (etc.) for the full column descriptions.



## Contributing

Contributions to `SQOUnified` are welcome, however, if you come across issues using the package, it would be preferable for one to submit an [issue](https://github.com/SCCWRP/SQOUnified/issues) rather than a pull request in most cases. If you have a request for a feature(s) then please do not hesitate at all to post an [issue](https://github.com/SCCWRP/SQOUnified/issues).


