# Empirical Illustration: NAEP Educational Achievement

This folder reproduces the empirical illustration in the paper 

> **Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure**

using publicly released subgroup summaries from the National Assessment of Educational Progress (NAEP).

The analysis studies the relationship between economically disadvantaged status and mathematics achievement 
using only disclosed quantile information and known support restrictions.

## Directory contents

```text
NAEP_data/
├── README.md
├── NAEP_2024_math_grade8_eds.csv
├── NAEP_2024_math_grade8_eds_share.csv
├── payload_full.json
├── 0 data_retrieve_and make.R
├── 1 run_naep_analysis.R
├── naep_analysis_data.csv
├── naep_qonlyreg_results.csv
├── naep_qonlyreg_results.rds
└── table_naep_qonlyreg.tex
```

## Software requirements

The analysis requires R (version 4.5.2 (2025-10-31)) and the `qonlyreg` package. 
Install the version of `qonlyreg` specified in the repository-level `README.md` before running this application.

The main analysis script also requires:

- `dplyr`        1.1.4
- `tidyr`        1.3.1
- `ggplot2`      4.0.1
- `patchwork`    1.3.2
- `knitr`         1.50
- `kableExtra`   1.4.0


Rebuilding the analysis dataset with `0 data_retrieve_and make.R`
additionally requires 

- `httr`         1.4.7
- `jsonlite`     2.0.0
- `purrr`        1.2.0
- `tibble`       3.3.0

## Data source

The empirical illustration uses publicly released subgroup summaries
from the 2024 Grade 8 Mathematics Assessment of the National Assessment
of Educational Progress (NAEP).

Data were obtained from the NAEP Data Explorer:

https://www.nationsreportcard.gov/ndecore/xplore/NDE

The analysis uses subgroup summaries defined by the interaction between:

- U.S. jurisdiction;
- economically disadvantaged status.

No individual-level microdata are used.

All NAEP data used in this application were retrieved on 19 May 2026.

### NAEP Data Explorer exports

The following two files were manually exported from the NAEP Data Explorer.

1. Subgroup quantiles and mean scores, using:

   - Subject: Mathematics
   - Grade: 8
   - Scale: Composite
   - Jurisdictions: All
   - Group: Economically disadvantaged status
   - Statistics: Average scale scores and percentiles

2. Subgroup population shares, using:

   - Subject: Mathematics
   - Grade: 8
   - Scale: Composite
   - Jurisdictions: All
   - Group: Economically disadvantaged status
   - Statistics: Percentages

The exported files are included in this directory as:

- `NAEP_2024_math_grade8_eds.csv`
- `NAEP_2024_math_grade8_eds_share.csv`

### Effective total sample sizes

The file

- `payload_full.json`

contains the API request template obtained from the NAEP Data Explorer web interface.

The script

- `0 data_retrieve_and make.R`

uses this template together with
`NAEP_2024_math_grade8_eds.csv` and
`NAEP_2024_math_grade8_eds_share.csv` to construct:

- `naep_analysis_data.csv`

The script requires an internet connection and may take some time 
because it sends separate requests to the NAEP regression API for each jurisdiction.

The script is provided to document the construction of the analysis dataset. 
For ordinary replication of the manuscript results, use the supplied `naep_analysis_data.csv`; 
rerunning the data-construction script is not required.


## Reproduction steps

Set the working directory to this application folder.

The supplied file `naep_analysis_data.csv` is the analysis dataset used to produce the manuscript results. 
Run:

```r
source("1 run_naep_analysis.R")
```


### Output files

- `naep_qonlyreg_results.csv`

  Main empirical results in CSV format.

- `naep_qonlyreg_results.rds`

  R object containing the main empirical results.

- `table_naep_qonlyreg.tex`

  LaTeX table used in the manuscript.

Existing output files with the same names are overwritten.


## Notes

The published subgroup mean scores are used only for benchmark validation. 
The partial-identification analysis uses the known NAEP score bounds `[0, 500]`.

