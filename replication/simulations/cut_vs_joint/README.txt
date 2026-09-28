Code for Section S4.8 of the Supplementary Material (R, package star). Run each script from this folder.

joint.R          margins, plug-in cut, cut posterior (segmented sampler cutx_mcmc) and joint posterior samplers
run_joint.R      simulation design (model at chibar = 0.6 pushed through known margins); read by run_all.R and tab_all.R
run_all.R        plug-in cut, cut posterior and joint posterior on the same datasets: GPD margins (Table S11) and a
                 lognormal body below the 0.9 quantile, n = 2000 and 400, twenty replicates; thinned draws saved in res_all/
tab_all.R        functionals of the saved draws, with bias, coverage and width over replicates
tab_s11.R        body of Table S11
outside_check.R  chi(q) read through the copula of the fit against the reading of the paper, for the Gaussian,
                 inverted logistic, Frank and Student t copulas (sources ../main/00_common.R; override with SIMS_DIR)
run_own.R        the same copula samples fitted on the model's own margins (P(AD | data) quoted in S4.8)
The Wanliu scripts of S4.8 (wanliu_joint.R, wanliu_cut60.R, sec6_cut.R) are in ../../applications/beijing, next to
results/fit_Wanliu.rds.
