# run_all.R -- driver for the identification study (Study 0)
#
# FIRST set the working directory to this folder, for example
#   setwd(".../replication/simulations/identification")
# or in RStudio: Session > Set Working Directory > To Source File Location.
# Then source this file. Adjust REPS / NS / ITER inside each script for
# paper-quality runs; defaults are sized for a workstation afternoon.

if (!file.exists("00_utils.R"))
  stop("Set the working directory to the identification folder first, e.g.\n",
       '  setwd(".../replication/simulations/identification")')
source("01_information.R")   # ~3 min per configuration, no MCMC
source("02_contraction.R")   # heaviest component
source("03_ridge.R")
