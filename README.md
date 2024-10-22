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
- **fieldrep**: (*optional*) - Data is filtered to where fieldrep = 1 if this column is included
- **labrep**: (*optional*) - Data is filtered to where labrep = 1 if this column is included
- **sampletypecode**: (*optional*) Data is filtered to where sampletypecode = Result in order to prevent inclusion of QA/QC samples

**Notes**:
- Non-detect values should be marked as `-88`.
- All measurements should be expressed on a dry-weight basis: metals in `mg/dry kg` and organic compounds in `ug/dry kg`.


#### 3. Input Data for toxicity functions (`tox.sqo`, `tox.summary`)

The toxicity functions require a dataframe containing:

- **stationid**: Alphanumeric identifier of the location.
- **toxbatch**: Identifier linking results to the control sample.
- **species**: The genus and species of the tested organism.
- **sampletypecode**: Type of sample (e.g., Grab, CNEG). Control samples must be included.
- **matrix**: (*optional*) - Type of sample matrix (e.g., Whole Sediment, Sediment Water Interface). Be sure to not include Reference Toxicants
- **labrep**: Numeric identifier for the lab replicate, typically with five replicates per station and species pair.
- **result**: The percentage of survival or normal development observed in the test.


### `SQOUnified` Function

The primary function `SQOUnified` computes integrated sediment quality scores based on various lines of evidence: benthic, chemistry, and toxicity data. It requires input data for each category and outputs an integrated assessment.

**Example Usage:**

```r
# Load the package
library(SQOUnified)

# Example data files
benthic_data <- read.csv("path/to/your/benthic_data.csv")
chem_data <- read.csv("path/to/your/chem_data.csv")
tox_data <- read.csv("path/to/your/tox_data.csv")

# Run the function
result <- SQOUnified(benthic = benthic_data, chem = chem_data, tox = tox_data)
print(result)
```

#### Parameters
- `benthic`: A dataframe containing benthic data for assessment.
- `chem`: A dataframe containing chemical concentration data.
- `tox`: A dataframe containing toxicity test results.


The function will compute and log the results for each of these inputs, generating an integrated score based on the criteria defined in the package. The output dataframe also includes all individual scores and indices:
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
- **Output:** A dataframe with SQO scores and categories for each station based on benthic criteria.

#### `benthic.sqo`

The `benthic.sqo` function computes benthic indices to evaluate the biological condition of the sediments.

**Usage:**

```r
benthic_results <- benthic.sqo(benthic_data)
```

- **Input:** `benthic_data` - A dataframe containing benthic community data.
- **Output:** A dataframe with SQO scores and categories for each station based on benthic criteria. It also contains the individual benthic index scores (BRI, RBI, IBI, RIVPACS and MAMBI)

**NOTE:** It will be the same for each of the benthic functions - BRI, RBI, IBI, RIVPACS and MAMBI



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



## Contributing

Contributions to `SQOUnified` are welcome, however, if you come across issues using the package, it would be preferable for one to submit an issue rather than a pull request in most cases. If you have a request for a feature(s) then please do not hesitate at all to post an issue on the github page.


