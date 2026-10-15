# Self-thresholding angular-radial models of extremal dependence: replication package

Manuele Leonelli and Dani Gamerman. Paper submitted to *Bayesian Analysis*. Repository: https://github.com/manueleleonelli/star.

```
star/            source of the R package star (version 0.2.0)
replication/     scripts that reproduce every table, figure and number of the paper and the supplement,
                 with the data and the saved results of the runs reported in the paper
```

## Setup

R (4.2 or later) with a C++ compiler. From `replication/`, run `source("00_setup.R")`, which installs Rcpp and
`star` from `star/`; the package alone can also be installed with
`remotes::install_github("manueleleonelli/star", subdir = "star")`. The comparison with the geometric model of
Wadsworth and Campbell (column WC of Tables 2, 4 and S7) also needs their package `geometricMVE` (version 0.1.6), which
is not redistributed here. The scripts are written for Linux or macOS, where they fork with `parallel::mclapply`; on
Windows they run on one core.

Every script runs from its own folder (`setwd()` to it, then `source()` the script, or `Rscript script.R`). Paths are
relative to that folder.

## Two ways to reproduce

1. **From the saved results (minutes).** Each folder contains the results of the runs reported in the paper (`results/`,
   chunk folders such as `htres/` or `res_all/`), and every fitting script skips work whose output exists. Running the
   table and figure scripts listed below therefore rebuilds the tables, figures and quoted numbers without any MCMC.
2. **From scratch (days of computing).** Delete the saved results of a folder and rerun its scripts in the order given.
   All runs are seeded. Costs are indicated below; the main simulation studies dominate (30 to 40 core-hours).

Simulated datasets are reproduced exactly on any platform (for example, the empirical, conditional extremes and
geometric estimates rebuilt from the seeds coincide with the saved ones), and a seeded chain is reproduced exactly on
the same machine. Across machines, floating-point differences in the likelihood can change individual accept or reject
decisions, so rerun chains agree with the saved ones up to Monte Carlo error rather than digit for digit.

## Where each result comes from

Folders are relative to `replication/`. Arrows give the order in which scripts are run.

### Main paper

| Output | Script(s) | Notes |
|---|---|---|
| Figures 1, 2, 3 | `figures/fig_illustrations.R` | Writes the TikZ files `fig-thresholds.tex`, `fig-scalefam.tex`, `fig-geometry.tex` (seconds). |
| Table 1 | none | Closed-form values from Section 3.2. |
| Figure 4 | `theory/fig_geometry3.R` | Uses the trivariate code shipped with `star` (`inst/trivariate/star3.R`). |
| Section 4.2 (prior sensitivity) | `simulations/sensitivity/run_sensitivity.R` → `simulations/tables/make_tables_supp.R` | Full rerun 1.5 to 2 hours. |
| Section 5, Tables 2 and 3 | `simulations/main/run_all.R`; `simulations/competitors/run_ht.R`, `run_wc.R`, `emp_known.R`, `emp_s3.R`, `ht_boot_known.R` → `assemble.R` → `simulations/tables/make_tables_main.R` | Main studies 30 to 40 core-hours; `STAR_SCALE = "test"` checks the pipeline in about 15 minutes. |
| Section 6.1, Figure 5 (a) to (e), Table 4 (STAR columns) | `applications/beijing/beijing_stations.R` → `app1_showcase.R` | Fits all twelve stations (about 3 minutes each); the Wanliu fit is the one reported. |
| Table 4 (HT and WC columns) and the comparison with the conditional extremes and geometric models in Section 6.1 | `applications/beijing/competitors_wanliu.R` | Run after `app1_showcase.R`; needs `geometricMVE`. Block bootstrap (B = 200) about 15 minutes on two cores, skipped when `results/competitors_wanliu_boot.rds` exists. Writes the rows of Table 4 to `tab_beijing.tex`. |
| Figure 5 (f) and the return level sets of Section 6.1 | `applications/beijing/wanliu_cut60.R` → `sec6_cut.R` → `ba_figure.R` | Cut posterior with sixty marginal draws; `ba_figure.R` writes the figure in the journal layout. |
| Convergence, extremal index and block bootstrap in Section 6.1 | `applications/beijing/app1_checks.R` | Bootstrap about 40 minutes on two cores. |
| Section 6.1, class probabilities under the cut and joint posteriors | `applications/beijing/wanliu_cut60.R`, `wanliu_joint.R` → `sec6_cut.R` | |
| Section 6.2, Table 5, Figure 6 | `applications/beijing3/data3.R` → `app3.R` → `tables.R`; `proj_app3.R` → `fig_app3.R` | Trivariate fit about 6 minutes; the figure also reads the bivariate Wanliu fit from `applications/beijing/results`. |

### Supplementary material

| Output | Script(s) | Notes |
|---|---|---|
| Constants of Proposition 3 (Section S2.6), numerical check | `theory/theory_check.R` | Under a minute; not reported in the paper. |
| Table S1 | `theory/geomtab.R` | Seconds. |
| Section S3.6, acceptance rates | `simulations/sampler_check/acceptance_rates.R` | About 10 minutes on two cores. |
| Table S2 | `simulations/main/s10_diagnostics.R` → `simulations/tables/make_tables_supp.R` | |
| Section S4.1, Table S3 | `simulations/identification/run_all.R` (01 to 04), `05_undercoverage.R`, `sensitivity/run_sensitivity.R` → `simulations/tables/make_tables_supp.R` | `05_undercoverage.R` takes 3 to 4 hours on 8 cores. |
| Tables S4 to S9, Sections S4.2 to S4.6 | as Tables 2 and 3 (`make_tables_main.R`); Table S8 from `simulations/competitors/ht_sens_run.R` | |
| Table S10, Section S4.7 | `simulations/sensitivity/run_sensitivity.R`, `simulations/main/s7_slab.R` → `simulations/tables/make_tables_supp.R` | |
| Table S11, Section S4.8 | `simulations/cut_vs_joint/run_all.R` → `tab_all.R` → `tab_s11.R`; `run_own.R` | About 2 hours per setting on two cores. |
| Table S12, Section S4.9 | `applications/beijing3/calib.R` → `tables.R` | Twenty fits of about 4 minutes. |
| Table S13, Section S5.1 | `applications/beijing/beijing_stations.R`; `clean_counts.R` for the cleaning counts | |
| Table S14, Section S5.2 | `applications/duck/app2_wavesurge.R` → `tables_duck.R` | Downloads the NOAA series on a first run (internet needed); fit 20 to 30 minutes. |

The trivariate implementation described in Section S3.7 is `star/inst/trivariate/star3.R`, loaded in R with
`source(system.file("trivariate", "star3.R", package = "star"))`.

## Data

* **Beijing** (`replication/data/beijing/`): the Beijing Multi-Site Air Quality data of Zhang et al. (2017), hourly
  records at twelve stations from March 2013 to February 2017, from the UCI Machine Learning Repository
  (https://archive.ics.uci.edu/dataset/501/beijing+multi+site+air+quality+data, licence CC BY 4.0). The files
  `PRSA_Data_<station>_20130301-20170228.csv` are included unchanged, renamed `<station>.csv` and compressed with xz
  (`<station>.csv.xz`, read directly by `read.csv`).
* **Duck, North Carolina**: hourly significant wave height at NDBC buoy 44100 and water level with its astronomical
  prediction at NOS tide gauge 8651370, 2010 to 2023, downloaded by `applications/duck/app2_wavesurge.R` from the
  National Data Buoy Center and NOAA CO-OPS. The daily maxima used in the paper are saved in
  `applications/duck/results/wavesurge.rds`.

## Software

The scripts were checked with R 4.3.3 on Linux, with `star` 0.2.0 and `geometricMVE` 0.1.6. During development the
package was called `bstar`, and the saved results were produced under that name; `star` 0.2.0 is version 0.1.4 of
`bstar` with the new names.
