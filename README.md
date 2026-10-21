# Replication Files for *Quantile-Only Regression*

This repository contains the code and data used to reproduce the simulation and empirical results in

> **Masamune Iwasawa and Hayato Nakanishi, “Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure.”**

The repository is organized into separate directories for the Monte Carlo simulations and empirical applications. Detailed instructions are provided in the README file within each directory.

## Repository structure

```text
QOR_Replication/
├── README.md
├── simulation/
│   ├── README.md
│   ├── 01_main_point_estimation.R
│   ├── 02_main_partial_identification.R
│   ├── 03_supplementD4_compact_support.R
│   └── 04_supplementD5_laplace.R
└── applications/
    ├── README.md
    ├── NAEP_data/
    ├── UK_election/
    ├── Health care meta analysis/
    └── vitamin_d_meta_analysis/
```

- `simulation/` contains the Monte Carlo code for Section 5 and Supplemental Material D.
- `applications/` contains the data and code for the NAEP, social care, healthcare-associated infection expenditure, and vitamin D meta-analysis applications.

## Required version of `qonlyreg`

The manuscript results were produced using `qonlyreg` version `0.1.0` at Git commit `025f75a`.

Install this exact version from GitHub before running the replication files:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "masaiwasawa/qonlyreg@025f75a",
  upgrade = "never"
)
```

The installed version can be checked with:

```r
packageVersion("qonlyreg")
utils::packageDescription("qonlyreg")$RemoteSha
```

The expected package version is `0.1.0`, and the recorded SHA should begin with `025f75a`.

## Reproducing the results

### Simulations

Follow the instructions in simulation/00_README.md. 
The four simulation scripts are mutually independent and may be run in any order. 
Because they are computationally intensive and use parallel workers, running them one at a time is recommended.

### Empirical applications

Follow the instructions in the README file in each application directory. 
Each application directory documents its data sources, required packages, script order, and correspondence between generated outputs and the manuscript tables and figures.

Scripts should generally be run with the relevant subdirectory as the working directory because the replication files use relative paths.

## Data sources

The empirical applications use publicly available statistical releases or summary information reported in published studies.

- Where automated retrieval is feasible, the corresponding download and data-construction scripts are supplied.
- Where inputs were transcribed from published tables, the application directory documents the original source and the construction of the analysis data.
- No confidential microdata are required to reproduce the reported results.

## Software environment

The exact R and package versions used for each component are reported in the corresponding subdirectory README or `sessionInfo.txt` file. 
Random-number seeds are set in the simulation and analysis scripts where stochastic computation is used.

The simulation results were generated with reproducible parallel random-number streams. 
The number of available processor cores affects runtime but should not affect the numerical results when the scripts are run with the recorded seeds and package versions.

## Citation

Until a final journal citation is available, please cite the arXiv preprint as:

> Iwasawa, Masamune, and Hayato Nakanishi. 2026. “Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure.” *arXiv preprint* [arXiv:2610.27224](https://arxiv.org/abs/2610.27224).

BibTeX:

```bibtex
@article{iwasawa2026quantile,
  title   = {Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure},
  author  = {Iwasawa, Masamune and Nakanishi, Hayato},
  journal = {arXiv preprint arXiv:2610.27224},
  year    = {2026},
  url     = {https://arxiv.org/abs/2610.27224}
}
```

## Authors

- **Masamune Iwasawa**, Doshisha University  
  Email: `miwasawa@mail.doshisha.ac.jp`
- **Hayato Nakanishi**, Kanagawa University  
  Email: `hnakanishi@kanagawa-u.ac.jp`
