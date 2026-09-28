MAIN SIMULATION STUDIES (Section 5 of the paper and Sections S3.6, S4.2 to S4.6 of the supplement)

Every fit uses the star defaults (v_p = 0.75, 30,000 iterations, 10,000 burn-in, thinning 5), 50 replicates per cell.

  setwd(".../replication/simulations/main")
  Sys.setenv(STAR_SCALE = "beach", STAR_CORES = "8")     # your core count
  source("run_all.R")

STAR_SCALE = "test" runs the whole pipeline at toy size in about 15 minutes, as a check of the installation.
Results are written to results/<study>_<cell>_c##.rds (checkpoints) and results/s#_*.rds (collated); re-sourcing
resumes where a run stopped. The collated results of the paper's run are included in results/, so the tables can be
rebuilt without rerunning: see ../tables/make_tables_main.R and ../tables/make_tables_supp.R.

s1  calibration under the model                                   Table S4, Table 3 (left)
s2  classification under the three-component prior                Table S5, Table 3 (right)
s3  misspecification (Gaussian, inverted logistic, Frank, t)      Table S6, Table 2 (right); writes s3_truths.rds
s4  self-thresholding ablation, well specified                    Table 3 (centre)
s5  Gaussian correlation sweep, P(AD | data)                      starred rows of Table S5
s6  Hill-type estimator of eta                                    Hill column of Table S6
s7  slab sensitivity at the boundary                              Table S10 (c)
s8  ablation under misspecification                               Table S9, Table 3 (centre)
s9  asymmetric structure                                          Section S4.2
s10 convergence diagnostics                                       Table S2
s11 Proposition 3 in the fitted object                            Section S4.2

Cost at full scale: about 30 to 40 core-hours, Study 2 the largest. Study 3 must precede Studies 6 and 8.
