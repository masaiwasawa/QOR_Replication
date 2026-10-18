# Healthcare-Associated Infections and Medical Expenditure Application

This directory contains the replication files for the healthcare-associated infections (HAI) and medical expenditure application in the paper

> **Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure**

The outcome is total medical expenditure per inpatient in Chinese yuan (CNY), expressed as the HAI-minus-non-HAI difference.

## Directory contents

```text
health care meta analysis/
├── README.md
├── health_care_meta_analysis.R
└── study_data.csv
```

Running the R script creates:

```text
table_liu_main.csv
table_liu_support_sensitivity.csv
```

## Software requirements

The analysis requires R (version 4.5.2 (2025-10-31)) and the `qonlyreg` package. 
Install the version of `qonlyreg` specified in the repository-level `README.md` before running this application.

No additional non-base R packages are required by the analysis script.

## Input data

The script reads:

```text
study_data.csv
```

The file contains the seven study or study-year comparisons used in the application. 
The arm-level sample sizes and reported summary statistics were manually transcribed from Table IV of:

- Liu, X., Spencer, A., Long, Y., Greenhalgh, C., Steeg, S., and Verma, A. (2022). 
“A systematic review and meta-analysis of disease burden of healthcare-associated infections in China: 
an economic burden perspective from general hospitals.” 
Journal of Hospital Infection, 123, 1–11. doi:10.1016/j.jhin.2022.02.005.
 
For the 2013–2015 comparisons reported by Li et al., the
paired-difference standard errors were obtained from the primary study:

- Li, H., Liu, X., Cui, D., Wang, Q., Mao, Z., Fang, L., Zhang, F., Yang, P., Wu, H., Ren, N., He, J., and Sun, J. (2017).
“Estimating the direct medical economic burden of health care-associated infections in public tertiary hospitals in Hubei Province, China.”
Asia Pacific Journal of Public Health, 29(5), 440–450. doi:10.1177/1010539517717366.
  
The input file includes:

- study identifiers;
- the reporting format for each comparison;
- HAI and control sample sizes;
- reported quartiles or reported means and standard deviations;
- the paired-difference standard errors and exchange-rate information used for the Li et al. comparisons.


## Running the analysis

Set the working directory to this application folder and run:

```r
source("health_care_meta_analysis.R")
```
Existing output files with the same names are overwritten.

Outputs:

-  `table_liu_main.csv`
- `table_liu_support_sensitivity.csv`
