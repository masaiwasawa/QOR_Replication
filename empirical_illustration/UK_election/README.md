# Children in Need Assessment Duration and Labour Majorities

This directory contains the replication files for the social-care
application in the paper

> **Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure**

The application combines publicly released council-year assessment
counts and duration quartiles with annual council-composition data. 
It examines differences in reconstructed mean assessment duration around
a 50 percent Labour seat-share threshold using prespecified fixed windows.

## Directory contents

```text
R/                  Reproduction scripts
data/raw/           Archived, unmodified public source files and manifest
data/derived/       Analysis data constructed from the raw files
output/             Estimation results and figures
output/paper_tables Manuscript-ready CSV and LaTeX tables
```

The archived source files in `data/raw/` are not manually edited.
Derived variables and cross-source joins are written only to
`data/derived/`.

## Software requirements

The analysis requires R (version 4.5.2 (2025-10-31)) and the `qonlyreg` package. 
Install the version of `qonlyreg` specified in the repository-level `README.md` before running this application.


The analysis also requires:

- `lpSolve` 5.6.23
- `ggplot2` 4.0.1

The optional download script verifies that the downloaded source files exactly match 
the versions used for the manuscript by comparing their SHA-256 checksums.
It uses the R package `openssl` when available or the system's checksum utility otherwise.

## Data sources

### Children in Need assessment duration

- Publisher: UK Department for Education (DfE)
- Dataset: [Characteristics of children in need](https://explore-education-statistics.service.gov.uk/data-catalogue/data-set/3fb5a93d-1eb3-4f3a-bdd0-5caea4556ce6)
- Archived file: `data/raw/children_in_need_assessment_duration_2014_2025.csv`
- Variables used: local-authority code and name, completed assessment
  count, and the 25th, 50th, and 75th percentiles of assessment
  duration


### Annual council composition

- Distributor: [electionresults.uk council data](https://electionresults.uk/councils/data)
- Original source: Open Council Data council-composition history
- Licence: CC BY-SA 4.0
- Source version (Git commit):
  [`29f587a2e761357c34cb828b0204861d8ec1c015`](https://github.com/fsargent/electionresults.uk/commit/29f587a2e761357c34cb828b0204861d8ec1c015)
- Downloaded file: `source-data/council/history2016-2025.csv`
- Archived file: `data/raw/council_composition_2016_2025.csv`

The file `data/raw/download_manifest.csv` records the source URL, source version, SHA-256 digest, and file size
for each archived input.

### Data retrieval

From the `R/` working directory, the script

```r
source("00_download_raw_data.R")
```

verifies the archived source files against their recorded SHA-256 digests. 
If a source file is missing, the script downloads the corresponding file from 
the recorded source URL and verifies its contents before saving it to `data/raw/`.

Because the archived source files are included in this replication package,
running this script is not required to reproduce the manuscript results.
Downloading missing files requires an internet connection and may take some
time depending on the connection and source availability.


## Reproduction steps

Open the `R/` directory and use it as the R working directory.

### Reproducing the manuscript results from the supplied analysis data

The supplied file

```text
data/derived/children_in_need_labour_rd_all_years.csv
```

is the analysis dataset used for the manuscript results. Run:

```r
source("02_run_qonlyreg_rd.R")
source("03_make_rd_quantile_plot.R")
source("04_make_paper_tables.R")
```

### Reconstructing the analysis data from the archived public files

To reconstruct the analysis dataset and reproduce all results, run:

```r
source("00_download_raw_data.R")
source("01_build_analysis_data.R")
source("02_run_qonlyreg_rd.R")
source("03_make_rd_quantile_plot.R")
source("04_make_paper_tables.R")
```

Existing output files with the same names are overwritten.


## Outputs

The principal outputs are:

- `output/qonlyreg_labour_fixed_window_point_estimates.csv`
- `output/qonlyreg_labour_fixed_window_simultaneous_region_summary.csv`
- `output/rd_quantiles_and_reconstructed_mean_main_window_4pp.png`
- `output/rd_quantiles_and_reconstructed_mean_main_window_4pp.pdf`
- `output/rd_quantiles_main_window_4pp.png`
- `output/rd_quantiles_main_window_4pp.pdf`
- `output/paper_tables/table_main_results.csv`
- `output/paper_tables/table_main_results.tex`

A complete list of manuscript and supporting outputs, together with
their intended use, is provided in:

```text
output/paper_tables/paper_output_manifest.csv
```

Additional audit files include:

- `data/derived/labour_rd_build_summary.csv`;
- `data/derived/unmatched_dfe_to_composition.csv`;
- `output/qonlyreg_labour_fixed_window_simultaneous_region_constraints.csv`.


