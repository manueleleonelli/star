# 01_information.R -- Component A: local information analysis
# Numerical observed information per observation (a refinement would average the
# Hessian over 20 replicate datasets, or use expected information via score
# outer products, to stabilize boundary configurations) at the truth, on a large
# sample, on the unconstrained sampling scale. Identification = positive
# definite information; weak identification = small eigenvalues, whose
# eigenvectors NAME the flat directions. No MCMC involved.
#
# Output: eigenvalues, leading flat directions, implied asymptotic
# correlations and per-parameter asymptotic sd at n = 2000.

source("00_utils.R")
set.seed(101)

N_INFO <- 25000                    # sample size per replicate dataset
N_H <- 8                           # Hessians averaged over N_H datasets:
                                   # stabilizes the flat directions that a
                                   # single observed information misses
configs <- list(AI_mid = theta_true(chibar = 0.6, p = 2),
                AI_rect = theta_true(chibar = 0.6, p = 8),
                near_AD = theta_true(chibar = 0.95, p = 2))

q0 <- -log(1 - 0.8)
for (nm in names(configs)) {
  cat(sprintf("configuration %s: simulating n = %g and computing the Hessian (~2-3 min)...\n", nm, N_INFO))
  th0 <- configs[[nm]]
  x0 <- pack_theta(th0, q0)
  hvec <- rep(2e-3, length(x0))
  hvec[match(c("log_o1", "log_o2", "log_p"), par_names)] <- 0.05   # wide steps
  Hs <- lapply(seq_len(N_H), function(k) {
    dat <- rstar(N_INFO, th0)
    num_hessian(make_ll(dat$r, dat$w, q0), x0, h = hvec) / N_INFO })
  H <- Reduce(`+`, Hs) / N_H
  Ihat <- -(H + t(H)) / 2
  if (min(eigen(Ihat, symmetric = TRUE, only.values = TRUE)$values) < 0)
    cat("note: information not positive definite even after averaging;\n",
        "     increase N_H, or interpret the negative directions as flat.\n")          # information per observation
  ev <- eigen(Ihat, symmetric = TRUE)
  cat("\n=====", nm, "=====\n")
  cat("eigenvalues of per-observation information:\n")
  print(signif(ev$values, 3))
  cat("\nflattest direction (smallest eigenvalue) loadings:\n")
  v <- ev$vectors[, length(par_names)]
  print(sort(round(setNames(v, par_names), 2)[abs(v) > 0.15], decreasing = TRUE))
  cat("\nsecond flattest direction loadings:\n")
  v2 <- ev$vectors[, length(par_names) - 1]
  print(sort(round(setNames(v2, par_names), 2)[abs(v2) > 0.15], decreasing = TRUE))
  V <- solve(Ihat)                          # asymptotic covariance (per obs)
  Cor <- cov2cor(V)
  dimnames(Cor) <- list(par_names, par_names)
  big <- which(abs(Cor) > 0.6 & upper.tri(Cor), arr.ind = TRUE)
  cat("\nasymptotic correlations with |corr| > 0.6:\n")
  for (k in seq_len(nrow(big)))
    cat(sprintf("  %-12s %-12s %6.2f\n", par_names[big[k, 1]],
                par_names[big[k, 2]], Cor[big[k, 1], big[k, 2]]))
  cat("\nasymptotic sd at n = 2000 (natural-ish scale):\n")
  print(signif(setNames(sqrt(diag(V) / 2000), par_names), 2))
  saveRDS(list(I = Ihat, eig = ev, Cor = Cor), sprintf("info_%s.rds", nm))
}
