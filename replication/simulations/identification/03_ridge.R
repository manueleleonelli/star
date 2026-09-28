# 03_ridge.R -- Component C: the anchoring ridge experiment
# One dataset, two priors: the anchoring prior of the paper versus a weak
# one. Reports posterior correlations on the sampling scale and pairs plots;
# discharges the todo of Section 4.4 (ridge without anchoring, removal with).

source("00_utils.R")
set.seed(303)

n <- 2000
th0 <- theta_true(chibar = 0.6, p = 2)
dat <- rstar(n, th0)

ITER <- 12000; BURN <- 4000; THIN <- 4

pri_anchor <- star_prior("beta", kappa0 = 0.8, s_u = 1)
pri_weak   <- star_prior("beta", kappa0 = 0.5, s_u = 5)

for (nm in c("anchor", "weak")) {
  pri <- if (nm == "anchor") pri_anchor else pri_weak
  fit <- star_mcmc(dat$E, prior = pri, n_iter = ITER, burnin = BURN,
                    thin = THIN, verbose = FALSE)
  D <- with(as.data.frame(fit$draws),
            cbind(log_u1 = log(u1), log_u2 = log(u2), log_p = log(p),
                  logit_chibar = qlogis(pmin(pmax(chibar, 1e-6), 1 - 1e-6)),
                  log_alpha = log(alpha), b0 = b0, z0 = z0))
  Cor <- cor(D)
  cat("\n=====", nm, "prior =====\n")
  cat("posterior correlations with |corr| > 0.5:\n")
  big <- which(abs(Cor) > 0.5 & upper.tri(Cor), arr.ind = TRUE)
  if (nrow(big) == 0) cat("  none\n")
  for (k in seq_len(nrow(big)))
    cat(sprintf("  %-14s %-14s %6.2f\n", colnames(D)[big[k,1]],
                colnames(D)[big[k,2]], Cor[big[k,1], big[k,2]]))
  cat(sprintf("posterior sd of log p: %.3f | of log u1: %.3f\n",
              sd(D[, "log_p"]), sd(D[, "log_u1"])))
  pdf(sprintf("ridge_pairs_%s.pdf", nm), 7, 7); pairs(D, pch = ".", col = "grey40"); dev.off()
  saveRDS(list(Cor = Cor, sd = apply(D, 2, sd)), sprintf("ridge_%s.rds", nm))
}
