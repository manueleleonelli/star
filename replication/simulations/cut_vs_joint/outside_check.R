## outside_check.R -- outside the model class, on exact exponential margins: does reading chi(q) through the copula of
## the fitted model (its own marginal quantiles) do better or worse than the paper's reading at exponential quantiles?
## Generators and truths of Section 5 (s3_misspecification.R); continuous prior on chibar as there.
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
wd0 <- getwd(); setwd(Sys.getenv("SIMS_DIR", "../main"))   # folder with 00_common.R of the simulation scripts
Sys.setenv(STAR_SCALE = "beach"); source("00_common.R")
setwd(wd0)
gens <- list(gauss05 = list(f = function(n) sim_gauss(n, 0.5), gen = "gauss", par = 0.5),
             invlog05 = list(f = function(n) sim_invlogistic(n, 0.5), gen = "invlogistic", par = 0.5),
             frank5 = list(f = function(n) sim_frank(n, 5), gen = "frank", par = 5),
             t05_4 = list(f = function(n) sim_t(n, 0.5, 4), gen = "t", par = c(0.5, 4)))
QS <- c(0.95, 0.99, 0.999); N <- 2000; R <- 10
set.seed(33); truths <- lapply(gens, function(g) sapply(QS, function(q) true_chi_q(g$gen, g$par, q, nmc = 2e7)))
gl <- seq(0.5 / 201, 1 - 0.5 / 201, length.out = 201); glw <- rep(1 / 201, 201)
js6 <- function(t, a, b) star:::joint_surv_cpp(a, b, gl, glw, t$phi, t$alpha, t$b, t$u1, t$u2, t$p, t$chibar, t$z)
mqp <- function(t, pe, j) uniroot(function(v) log(max(if (j == 1) js6(t, v, 0) else js6(t, 0, v), 1e-300)) - log(pe), c(0.001, 80), tol = 1e-8)$root
dir.create("res_out", showWarnings = FALSE)
jobs <- do.call(c, lapply(names(gens), function(g) lapply(1:R, function(r) list(g = g, r = r))))
invisible(mclapply(jobs, function(j) { f <- sprintf("res_out/%s_%02d.rds", j$g, j$r); if (file.exists(f)) return(NULL)
  set.seed(3e4 * match(j$g, names(gens)) + j$r); E <- gens[[j$g]]$f(N)
  fit <- star_mcmc(E, prior = star_prior("beta"), n_iter = 12000, burnin = 4000, thin = 5, verbose = FALSE)
  idx <- round(seq(1, nrow(fit$draws), length.out = 120))
  V <- t(sapply(idx, function(s) { t <- star:::theta_row(fit, s)
    c(std = sapply(QS, function(q) js6(t, -log(1 - q), -log(1 - q)) / (1 - q)),
      cop = sapply(QS, function(q) js6(t, mqp(t, 1 - q, 1), mqp(t, 1 - q, 2)) / (1 - q))) }))
  out <- list(g = j$g, r = j$r, med = apply(V, 2, median), lo = apply(V, 2, quantile, .05), hi = apply(V, 2, quantile, .95), truth = truths[[j$g]])
  saveRDS(out, f); NULL }, mc.cores = 2, mc.preschedule = FALSE))
res <- lapply(list.files("res_out", full.names = TRUE), readRDS)
for (g in names(gens)) { Rg <- Filter(function(z) z$g == g, res); if (!length(Rg)) next
  M <- t(sapply(Rg, `[[`, "med")); L <- t(sapply(Rg, `[[`, "lo")); H <- t(sapply(Rg, `[[`, "hi")); tr <- rep(truths[[g]], 2)
  cat(sprintf("\n%s (%d replicates), truth chi(0.95, 0.99, 0.999) = %s\n", g, length(Rg), paste(sprintf("%.3f", truths[[g]]), collapse = ", ")))
  cat("  reading at exponential quantiles (paper): bias", sprintf("%+.3f", colMeans(M[, 1:3, drop = FALSE]) - tr[1:3]), " rmse", sprintf("%.3f", sqrt(colMeans(sweep(M[, 1:3, drop = FALSE], 2, tr[1:3])^2))),
      " cover", sprintf("%.2f", colMeans(sweep(L[, 1:3, drop = FALSE], 2, tr[1:3], "<=") & sweep(H[, 1:3, drop = FALSE], 2, tr[1:3], ">="))), "\n")
  cat("  copula of the fit:                        bias", sprintf("%+.3f", colMeans(M[, 4:6, drop = FALSE]) - tr[4:6]), " rmse", sprintf("%.3f", sqrt(colMeans(sweep(M[, 4:6, drop = FALSE], 2, tr[4:6])^2))),
      " cover", sprintf("%.2f", colMeans(sweep(L[, 4:6, drop = FALSE], 2, tr[4:6], "<=") & sweep(H[, 4:6, drop = FALSE], 2, tr[4:6], ">="))), "\n") }
