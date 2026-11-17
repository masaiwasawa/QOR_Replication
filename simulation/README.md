# Simulation Replication Files

This directory contains the Monte Carlo simulation code for the paper

> **Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure**

The five scripts reproduce the simulation results reported in Section 5 and Supplemental Material D. They are mutually independent and may be run in any order.

## Directory contents

```text
simulation/
├── README.md
├── 01_main_point_estimation.R
├── 02_main_partial_identification.R
├── 03_supplementD4_compact_support.R
├── 04_supplementD5_laplace.R
└── 05_supplementD4_laplace_partial_identification.R
```

Each script creates its own output directory when run.

## Software requirements

The simulations require R and the following packages:

```r
           version
qonlyreg     0.1.0
dplyr        1.1.2
tidyr        1.3.0
purrr        1.0.1
ggplot2      4.0.3
doParallel  1.0.17
foreach      1.5.2
doRNG        1.8.6
patchwork    1.3.2
knitr         1.52
kableExtra   1.4.1
```

```text
R version:              4.5.3
Operating system:       Windows 11 Pro for Workstation
qonlyreg commit:        025f75a
```

## Scripts and manuscript outputs
Simulation designs, parameter settings, and disclosure grids are described in Section 5 and Supplemental Material D of the paper and are also specified directly in the R scripts.

### 01_main_point_estimation.R

Reproduces the point-estimation simulation in Section 5.1 and the
corresponding results in Supplemental Material D.

Main outputs:
- Figure 1
- Table D.1
- Figures D.2--D.5

Results are written to `results_smooth_qlm/`.

### 02_main_partial_identification.R

Reproduces the partial-identification simulation in Section 5.2 and
the corresponding results in Supplemental Material D.
The uniform-error design yields consistent finite-difference covariance
estimates even with fixed `K`.

Main outputs:
- Table 1
- Figure D.1
- Tables D.5--D.6
- Figures D.9--D.10

Results are written to `results_partial_qlm/`.

### 03_supplementD4_compact_support.R

Reproduces the point-estimation simulation under asymmetric truncated-normal
errors reported in Supplemental Material D-4.1.

Main outputs:
- Table D.2
- Figure D.6

Results are written to `results_compact_truncnorm_qlm/`.

### 04_supplementD5_laplace.R

Reproduces the Laplace-error simulation reported in Supplemental Material D-5.

Main outputs:
- Table D.4
- Figure D.8

Results are written to `results_nonsmooth/`.

### 05_supplementD4_laplace_partial_identification.R

Reproduces the partial-identification simulation under standardized,
asymmetric truncated-normal errors reported in Supplemental Material D-4.2.
The data-generating process is the same as in
`03_supplementD4_compact_support.R`.
This experiment examines the same finite-difference procedure when its
covariance estimator is not consistent for fixed `K`.

Main outputs:

- Table D.3: `table_region_coverage_all_grids.tex`
- Figure D.7: `plot_length_beta1.png`

Results are written to `results_compact_truncnorm_qlm_partial/`.

## Running the simulations

Set the working directory to this `simulation/` directory and run the scripts individually. For example:

```r
source("01_main_point_estimation.R")
source("02_main_partial_identification.R")
source("03_supplementD4_compact_support.R")
source("04_supplementD5_laplace.R")
source("05_supplementD4_laplace_partial_identification.R")
```

Because the simulations are computationally intensive and use parallel workers, running them one at a time is recommended.

The scripts are independent; it is not necessary to run them in numerical order. Each script creates its own output directory and may overwrite existing files with the same names.

## Reproducing figures and tables without rerunning all replications

Each script saves both the full replication-level results and a compact summary object. To regenerate figures and tables without rerunning the Monte Carlo loop, the script can be modified to load the saved object at the summary stage, for example:

```r
results <- readRDS("results_smooth_qlm/results.rds")
```

or, when only the summarized results are needed:

```r
summary_tbl <- readRDS("results_smooth_qlm/summary_tbl.rds")
```

Equivalent paths apply to the other four output directories. The archived summary files are the quickest route for checking the manuscript tables and figures.

## Reproducibility notes

- All stochastic simulations use seed `123` and reproducible parallel streams from `doRNG`.
- The exact numerical results depend on the archived version of `qonlyreg` and the R package environment. 
- The number of parallel workers affects runtime but should not affect results when `doRNG` is used as in these scripts.
- `meta_info.txt` records the number of replications, cores, start and end times, and elapsed runtime for each completed run.
- No external data files are required for these simulations.
