# s10_diagnostics.R -- MCMC convergence diagnostics: 4 chains per dataset,
# Gelman-Rubin R-hat and batch-means ESS per parameter, at representative
# configurations. Reports the worst R-hat and worst ESS across parameters.
source("00_common.R"); source("00_extra.R")
configs <- list(mid_n1000 = list(gen=function() rstar(1000, theta_model(0.6))$E),
                mid_n4000 = list(gen=function() rstar(4000, theta_model(0.6))$E),
                gauss09_n4000 = list(gen=function() sim_gauss(4000, 0.9)))
R <- if (SCALE=="test") 3 else 12
NCH <- if (SCALE=="test") 2 else 4
one_rep <- function(seed, cname) {
  set.seed(seed); E <- configs[[cname]]$gen()
  chains <- lapply(seq_len(NCH), function(ch) {
    f <- star_mcmc(E, prior=star_prior("beta"), n_iter=SET_NEW$ITER,
                    burnin=SET_NEW$BURN, thin=SET_NEW$THIN, verbose=FALSE)
    f$draws[, c("chibar","alpha","u1","u2","p","z0")] })
  Rh <- gelman_rubin(chains)
  ess <- apply(do.call(rbind, chains), 2, ess_bm)[c("chibar","alpha","u1","u2","p","z0")]
  data.frame(cfg=cname, seed=seed, max_rhat=max(Rh), rhat_p=Rh["p"],
             min_ess=min(ess, na.rm=TRUE), ess_p=ess["p"], ess_chibar=ess["chibar"])
}
res <- list()
for (cname in names(configs))
  res[[cname]] <- run_cell("s10", cname, R, function(r) one_rep(1e5+41*r+nchar(cname), cname), chunk = 4)
res <- do.call(rbind, res); saveRDS(res, "results/s10_diagnostics.rds")
cat("\n===== Study 10: convergence diagnostics =====\n")
print(aggregate(cbind(max_rhat, rhat_p, min_ess, ess_p, ess_chibar) ~ cfg, res, mean), digits=3, row.names=FALSE)
cat("\ntarget: max R-hat < 1.01, ESS > 400; p and the threshold block are the\n",
    "binding cases, so ess_p is the number to watch.\n")
