## 00_setup.R -- installs what the replication scripts need. Run once from this folder: source("00_setup.R").
## R >= 4.2 and a C++ compiler (Rtools on Windows, Xcode command line tools on macOS) are required.
if (!requireNamespace("Rcpp", quietly = TRUE)) install.packages("Rcpp")
install.packages("../star", repos = NULL, type = "source")      # the package source in the folder star/
library(star); cat("star", as.character(packageVersion("star")), "installed\n")
## Optional: coda (effective sample sizes in beijing3/calib.R).
if (!requireNamespace("coda", quietly = TRUE)) message("optional package coda not installed")
## Only for the geometric competitor (WC) in simulations/competitors/run_wc.R: the geometricMVE package (version 0.1.6)
## of Wadsworth and Campbell, not redistributed here; install it from the authors' source.
if (!requireNamespace("geometricMVE", quietly = TRUE)) message("geometricMVE not installed: needed only by simulations/competitors/run_wc.R")
## The scripts use the base package parallel (mclapply); on Windows they run on one core.
