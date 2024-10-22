# SQOUnified R Package
  <!-- badges: start -->
  [![Travis build status](https://travis-ci.com/SCCWRP/SQOUnified.svg?branch=master)](https://travis-ci.com/SCCWRP/SQOUnified)
  <!-- badges: end -->


The `SQOUnified` R package provides a suite of tools for assessing sediment quality based on multiple lines of evidence, including benthic indices, chemistry, and toxicity. The package integrates these assessments into a comprehensive site evaluation according to Sediment Quality Objectives (SQO).

This is based on the [California Sediment Quality Objectives (CASQO) Technical Manual](http://ftp.sccwrp.org/pub/download/DOCUMENTS/TechnicalReports/777_CASQO_TechnicalManual.pdf)

## Installation

To install the `SQOUnified` package, you can use the `devtools` package. Ensure all dependencies are installed beforehand.

```r
# Install necessary dependencies
install.packages(c("dplyr", "tidyr", "purrr", "data.table", "ggplot2", "readr"))

# Install the package using devtools
devtools::install_github("SCCWRP/SQOUnified")
```

## Usage

### 1. `SQOUnified` Function

The primary function `SQOUnified` computes integrated sediment quality scores based on various lines of evidence: benthic, chemistry, and toxicity data. It requires input data for each category and outputs an integrated assessment.

**Example Usage:**

```r
# Load the package
library(SQOUnified)

# Example data files
benthic_data <- read.csv("path/to/benthic_data.csv")
chem_data <- read.csv("path/to/chem_data.csv")
tox_data <- read.csv("path/to/tox_data.csv")

# Run the function
result <- SQOUnified(benthic = benthic_data, chem = chem_data, tox = tox_data)
print(result)
```

### Parameters
- `benthic`: A dataframe containing benthic data for assessment.
- `chem`: A dataframe containing chemical concentration data.
- `tox`: A dataframe containing toxicity test results.

The function will compute and log the results for each of these inputs, generating an integrated score based on the criteria defined in the package.

### 2. `chem.sqo`

The `chem.sqo` function evaluates the chemical condition of sediments using specific chemical metrics.

**Usage:**

```r
chem_results <- chem.sqo(chem_data)
```

- **Input:** `chem_data` - A dataframe containing sediment chemistry data.
- **Output:** A dataframe with SQO scores and categories for each station based on chemical criteria.

### 3. `benthic.sqo`

The `benthic.sqo` function computes benthic indices to evaluate the biological condition of the sediments.

**Usage:**

```r
benthic_results <- benthic.sqo(benthic_data)
```

- **Input:** `benthic_data` - A dataframe containing benthic community data.
- **Output:** A dataframe with SQO scores and categories for each station based on benthic criteria.

### 4. `tox.sqo`

The `tox.sqo` function assesses sediment toxicity by analyzing toxicity test results.

**Usage:**

```r
tox_results <- tox.sqo(tox_data)
```

- **Input:** `tox_data` - A dataframe containing toxicity test results.
- **Output:** A dataframe with SQO scores and categories for each station based on toxicity criteria.

### 5. `tox.summary`

The `tox.summary` function provides a summary of the toxicity results for easier interpretation and reporting.

**Usage:**

```r
tox_summary <- tox.summary(tox_data)
```

- **Input:** `tox_data` - A dataframe containing toxicity test results.
- **Output:** A summary dataframe providing an overview of toxicity conditions for all stations.

## Example Workflow

Here’s how you can use the `SQOUnified` package to perform a complete sediment quality assessment:

```r
# Load sample data
benthic_data <- read.csv("path/to/benthic_data.csv")
chem_data <- read.csv("path/to/chem_data.csv")
tox_data <- read.csv("path/to/tox_data.csv")

# Calculate individual SQO components
benthic_results <- benthic.sqo(benthic_data)
chem_results <- chem.sqo(chem_data)
tox_results <- tox.sqo(tox_data)

# Get a summary of the toxicity data
tox_summary <- tox.summary(tox_data)

# Compute the integrated SQO score
integrated_result <- SQOUnified(benthic = benthic_data, chem = chem_data, tox = tox_data)

# View the results
print(integrated_result)
```

## Data Requirements

Each function requires specific data formats:
- **benthic data**: Dataframe with information on benthic community structure.
- **chem data**: Dataframe with chemical concentrations.
- **tox data**: Dataframe with toxicity test results.

Ensure that the column names match the expected format specified in the package documentation for accurate computation.

## Contributing

Contributions to `SQOUnified` are welcome. If you encounter any issues or have feature requests, please open an issue or submit a pull request on the GitHub repository.

## License

This package is licensed under the [MIT License](LICENSE).
