# s1_calibration.R -- Study 1: frequentist calibration under the model
# Truths chibar in {0.15, 0.45, 0.75, 0.95} x n in {N1, N2} (+ one asymmetric
# cell); bias, RMSE and 90% coverage for chibar, eta, sigma(1/2), rho at three
# angles, p and chi(0.99). Atoms are the subject of Study 2.

source("00_common.R")
set.seed(11)

cells <- expand.grid(chibar = c(0.15, 0.45, 0.75, 0.95), n = c(SET$N1, SET$N2))
cells$u1 <- 2.2; cells$u2 <- 2.2
cells <- rbind(cells, data.frame(chibar = 0.45, n = SET$N2, u1 = 2.8, u2 = 1.7))

one_rep <- function(seed, chibar, n, u1, u2) {
  set.seed(seed)
  th0 <- theta_model(chibar = chibar, u1 = u1, u2 = u2)
  dat <- rstar(n, th0)
  fit <- star_mcmc(dat$E, prior = star_prior("beta"), n_iter = SET$ITER,
                    burnin = SET$BURN, thin = SET$THIN, verbose = FALSE)
  m <- fit_metrics(fit)
  tq <- -log(1 - 0.99); gl <- seq(0.5/201, 1-0.5/201, length.out = 201)
  tr <- c(chibar = chibar, eta = (1 + chibar)/2, sigma_half = 1 + chibar,
          setNames(rho_star(c(0.15, 0.5, 0.85), u1, u2, 2), paste0("rho", c(0.15,0.5,0.85))),
          p = 2,
          chi_q = star:::joint_surv_cpp(tq, tq, gl, rep(1/201, 201), th0$phi, th0$alpha,
                                 th0$b, u1, u2, 2, chibar, th0$z) / 0.01,
          lam03 = { wg <- seq(0.002,0.998,length.out=499)
                    sg <- sigma_star(wg, chibar, th0$z)
                    min(pmax(0.3/wg, 0.7/(1-wg))/sg) })
  m$truth <- tr[m$quantity]; m$err <- m$med - m$truth
  m$cover <- as.numeric(m$truth >= m$lo & m$truth <= m$hi)
  cbind(m, chibar0 = chibar, n = n, u1 = u1, u2 = u2, seed = seed)
}

res <- list()
for (i in seq_len(nrow(cells))) {
  cl <- cells[i, ]
  cell_id <- sprintf("chi%02.0f_n%d_u%.1f", 100*cl$chibar, cl$n, cl$u1)
  res[[i]] <- run_cell("s1", cell_id, SET$R, function(r)
    one_rep(1e4*i + r, cl$chibar, cl$n, cl$u1, cl$u2))
}
res <- do.call(rbind, res); saveRDS(res, "results/s1_all.rds")

cat("\n===== Study 1: bias / RMSE / coverage =====\n")
agg <- aggregate(cbind(err, cover) ~ quantity + chibar0 + n, res,
                 function(x) c(m = mean(x)))
agg2 <- aggregate((err)^2 ~ quantity + chibar0 + n, res, function(x) sqrt(mean(x)))
names(agg2)[4] <- "rmse"
out <- merge(agg, agg2)
print(out[order(out$quantity, out$chibar0, out$n), ], digits = 2, row.names = FALSE)

