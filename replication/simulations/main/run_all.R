# run_all.R -- one simulation run, 50 replicates, every fit under
# the star defaults (v_p = 0.75, 30000 iterations). Single provenance for the
# simulation section of the paper.
#
#   setwd(".../replication/simulations/main")
#   Sys.setenv(STAR_SCALE = "beach", STAR_CORES = "8")   # your core count
#   source("run_all.R")
#
# FULLY CHECKPOINTED at 10-replicate granularity: every study writes
# results/<study>_<cell>_c##.rds and re-sourcing resumes where it stopped, so
# the run can be spread across several sessions. Set STAR_SCALE = "test" for a
# ~15 min end-to-end validation before committing the machine.
#
# COST: at beach scale the full run is heavy (order 30-40 core-hours, Study 2
# the largest single piece). On 8 cores plan for a long weekend, or run the
# studies one at a time. To trade a little threshold-block ESS at the largest n
# for roughly a third less time, lower ITER in 00_common.R from 30000 to 20000.

if (!file.exists("00_common.R")) stop('setwd to the simulations/main folder first')
t0 <- Sys.time()
# Study 3 must precede Studies 6 and 8 (they reuse results/s3_truths.rds).
source("s1_calibration.R")        # calibration under the model
source("s3_misspecification.R")   # writes s3_truths.rds
source("s2_classification.R")     # classification operating characteristics
source("s4_ablation.R")           # self-thresholding ablation (well specified)
source("s5_gauss_sweep.R")        # P(AD) across Gaussian rho (money figure)
source("s6_hill.R")               # frequentist Hill-type eta benchmark
source("s7_slab.R")               # slab sensitivity at the boundary
source("s8_misspec_ablation.R")   # ablation under misspecification
source("s9_asymmetric.R")         # asymmetric structure (within + Khoudraji)
source("s10_diagnostics.R")       # R-hat / ESS per block
source("s11_theory_checks.R")     # Proposition 3 in the fitted object
cat(sprintf("\nSimulation run complete in %.1f h. Results in results/*.rds\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))
