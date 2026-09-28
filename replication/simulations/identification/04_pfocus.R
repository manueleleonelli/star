# 04_pfocus.R -- Component D: coordinate versus curve identification for the
# threshold. Hypothesis: the CURVE rho(.) contracts regularly even where the
# coordinates (u1, u2, p) trade off, so the function is the identified
# estimand. Reports posterior quantiles of p, pointwise widths of rho(w),
# and the sup-norm error of the fitted curve. Package v0.1.1 places a
# log-normal prior on p over (0, Inf), so p = 1 is interior; the p = 1
# undercoverage seen under v0.1.0 was a boundary artifact of the old
# Uniform(0, 1] prior on 1/p.

source("00_utils.R")
library(parallel)
set.seed(404)

REPS <- 20
NS <- c(1000, 4000)
PS <- c(1, 2, 8)
ITER <- 8000; BURN <- 3000; THIN <- 5
CORES <- getOption("mc.cores", 8)
wgrid <- seq(0.02, 0.98, length.out = 97)

one_rep <- function(seed, n, ptrue) {
  set.seed(seed)
  th0 <- theta_true(chibar = 0.6, p = ptrue)
  dat <- rstar(n, th0)
  fit <- star_mcmc(dat$E, prior = star_prior("beta"), n_iter = ITER,
                    burnin = BURN, thin = THIN, verbose = FALSE)
  S <- nrow(fit$draws)
  R <- sapply(seq_len(S), function(s) {
    th <- star:::theta_row(fit, s); rho_star(wgrid, th$u1, th$u2, th$p) })
  rho0 <- rho_star(wgrid, th0$u1, th0$u2, th0$p)
  ci <- apply(R, 1, quantile, c(0.05, 0.95))
  data.frame(n = n, ptrue = ptrue, seed = seed,
             p_med = median(fit$draws[, "p"]),
             p_lo = quantile(fit$draws[, "p"], 0.05),
             p_hi = quantile(fit$draws[, "p"], 0.95),
             rho_width_mean = mean(ci[2, ] - ci[1, ]),
             rho_cover = mean(rho0 >= ci[1, ] & rho0 <= ci[2, ]),
             rho_sup_err = max(abs(rowMeans(R) - rho0)))
}

grid <- expand.grid(seed = seq_len(REPS), n = NS, ptrue = PS)
cat(sprintf("Component D: %d fits; progress by chunk...\n", nrow(grid)))
res <- run_par(seq_len(nrow(grid)), function(i)
  one_rep(grid$seed[i] + 7 * grid$n[i] + 991 * grid$ptrue[i],
          grid$n[i], grid$ptrue[i]), cores = CORES)
res <- do.call(rbind, res)
saveRDS(res, "pfocus.rds")

cat("\nmedian over replicates, by (p, n):\n")
agg <- aggregate(cbind(p_med, p_lo, p_hi, rho_width_mean, rho_cover, rho_sup_err)
                 ~ ptrue + n, res, median)
print(within(agg, { p_CI_width <- p_hi - p_lo }), digits = 3)
cat("\nreading: if p_CI_width stalls while rho_width_mean and rho_sup_err\n",
    "shrink like 1/sqrt(n), the curve is identified and the coordinates trade off.\n")
