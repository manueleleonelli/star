# s4_ablation.R -- Study 4: the value of self-thresholding, by ablation.
# Same model, same priors; the only difference is whether the threshold curve
# is estimated or fixed at the level-kappa calibrated curve. Reports the
# sensitivity of eta and chi(0.999) across fixed kappa and interval coverage.

source("00_common.R")
set.seed(44)

th0 <- theta_model(chibar = 0.6)
tq3 <- -log(1 - 0.999); gl <- seq(0.5/201, 1-0.5/201, length.out = 201)
truth <- c(eta = 0.8,
           chi999 = star:::joint_surv_cpp(tq3, tq3, gl, rep(1/201, 201), th0$phi,
                    th0$alpha, th0$b, th0$u1, th0$u2, th0$p, th0$chibar,
                    th0$z) / 0.001)

calib_u <- function(r, w, kappa) {   # u with p = 2 giving exceedance 1 - kappa
  f <- function(u) mean(r > rho_star(w, u, u, 2)) - (1 - kappa)
  uniroot(f, c(0.3, 12))$root
}

one_rep <- function(seed, variant) {
  set.seed(seed)
  dat <- rstar(SET$N2, th0); r <- rowSums(dat$E); w <- dat$E[,1]/r
  fx <- if (variant == "full") NULL else {
    u <- calib_u(r, w, as.numeric(variant)); list(u1 = u, u2 = u, p = 2) }
  fit <- star_mcmc(dat$E, prior = star_prior("beta"), n_iter = SET$ITER,
                    burnin = SET$BURN, thin = SET$THIN, fix_thr = fx,
                    verbose = FALSE)
  et <- eta_draws(fit)
  cq <- sapply(round(seq(1, nrow(fit$draws), length.out = 120)), function(s) {
    th <- star:::theta_row(fit, s)
    star:::joint_surv_cpp(tq3, tq3, gl, rep(1/201, 201), th$phi, th$alpha, th$b,
                   th$u1, th$u2, th$p, th$chibar, th$z) / 0.001 })
  data.frame(variant = variant, seed = seed,
             eta_med = median(et), eta_w = diff(quantile(et, c(.05, .95))),
             eta_cover = as.numeric(truth["eta"] >= quantile(et, .05) &
                                    truth["eta"] <= quantile(et, .95)),
             chi999_med = median(cq),
             chi999_cover = as.numeric(truth["chi999"] >= quantile(cq, .05) &
                                       truth["chi999"] <= quantile(cq, .95)))
}

res <- list()
for (v in c("full", "0.8", "0.9", "0.95")) {
  res[[v]] <- run_cell("s4", paste0("v", gsub("\\.", "", v)), SET$R,
                       function(r) one_rep(4e4 + 97*r + nchar(v), v))
}
res <- do.call(rbind, res); saveRDS(res, "results/s4_all.rds")

cat("\n===== Study 4: self-thresholding ablation (truth eta = 0.8) =====\n")
print(aggregate(cbind(eta_med, eta_w, eta_cover, chi999_med, chi999_cover)
                ~ variant, res, mean), digits = 3, row.names = FALSE)
cat("\nreading: the fixed-kappa rows should show (i) point estimates moving\n",
    "with kappa, the sensitivity the full model removes, and (ii) narrower\n",
    "intervals with lower coverage, threshold uncertainty conditioned away.\n")
