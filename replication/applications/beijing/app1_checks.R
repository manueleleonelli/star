## app1_checks.R -- checks behind Section 6 of the paper for the Wanliu fit: convergence over four
## chains, serial dependence (extremal index, alternate-day fits, moving-block bootstrap) and the
## tails of the margins implied by the fitted model. Run from this folder after beijing_stations.R:
## source("app1_checks.R"). The fits are resumable (results/wanliu_*.rds, results/boot/*.rds); the
## bootstrap takes about 40 minutes on two cores and is skipped with STAR_BOOT=0.
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
suppressMessages(source("app_common.R")); library(parallel)
CORES <- as.integer(Sys.getenv("STAR_CORES", "2"))
gl <- seq(0.5/201, 1-0.5/201, length.out = 201); glw <- rep(1/201, 201)
x <- readRDS("results/fit_Wanliu.rds"); a <- x$data; E <- x$E; n <- nrow(E); r <- rowSums(E); w <- E[, 1] / r
season <- format(a$date, "%Y")
fit_one <- function(EE, seed, iter = 30000) { set.seed(seed)
  star_mcmc(EE, prior = star_prior("threecomp", p_min = 1), n_iter = iter, burnin = iter / 3, thin = 5, verbose = FALSE) }

## ---- 1. three further chains and the alternate-day fits ------------------------------------------
jobs <- list(list(tag = "chain2", E = E, seed = 2027), list(tag = "chain3", E = E, seed = 2028), list(tag = "chain4", E = E, seed = 2029),
             list(tag = "odd", E = E[seq(1, n, 2), ], seed = 2026), list(tag = "even", E = E[seq(2, n, 2), ], seed = 2026))
todo <- Filter(function(j) !file.exists(sprintf("results/wanliu_%s.rds", j$tag)), jobs)
invisible(mclapply(todo, function(j) saveRDS(list(tag = j$tag, E = j$E, fit = fit_one(j$E, j$seed)), sprintf("results/wanliu_%s.rds", j$tag)),
                   mc.cores = CORES, mc.preschedule = FALSE))

## ---- 2. convergence: rank-normalized split R-hat (bulk and tail) and effective sample sizes -----
derived <- function(fit, EE) { D <- fit$draws; rr <- rowSums(EE); ww <- EE[, 1] / rr
  t(sapply(seq_len(nrow(D)), function(k) { t <- star:::theta_row(fit, k)
    js <- function(l1, l2) star:::joint_surv_cpp(-log(1-l1), -log(1-l2), gl, glw, t$phi, t$alpha, t$b, t$u1, t$u2, t$p, t$chibar, t$z)
    c(exc = mean(rr > rho_star(ww, t$u1, t$u2, t$p)), P9090 = js(.9, .9), P9595 = js(.95, .95), P9990 = js(.99, .9), P9099 = js(.9, .99), P9999 = js(.99, .99),
      AD = as.numeric(fit$component[k] == "one"), AI0 = as.numeric(fit$component[k] == "zero")) })) }
fits <- c(list(x$fit), lapply(2:4, function(k) readRDS(sprintf("results/wanliu_chain%d.rds", k))$fit))
V <- lapply(fits, function(fit) cbind(fit$draws, loglik = fit$loglik[seq_len(nrow(fit$draws))], derived(fit, E)))
rnorm_rank <- function(v) qnorm((rank(v, ties.method = "average") - 3/8) / (length(v) + 1/4))
rhat_basic <- function(M) { N <- nrow(M); B <- N * var(colMeans(M)); W <- mean(apply(M, 2, var)); sqrt(((N - 1)/N * W + B/N) / W) }
rhat_rank <- function(M) { h <- floor(nrow(M) / 2); S <- cbind(M[1:h, , drop = FALSE], M[(nrow(M) - h + 1):nrow(M), , drop = FALSE])
  if (var(as.vector(S)) == 0) return(NA)
  max(rhat_basic(matrix(rnorm_rank(as.vector(S)), nrow(S))), rhat_basic(matrix(rnorm_rank(abs(as.vector(S) - median(S))), nrow(S)))) }
ess1 <- function(v) { if (var(v) == 0) return(NA); ac <- acf(v, lag.max = min(length(v) - 1, 1000), plot = FALSE)$acf[-1]; s <- 0
  for (k in seq(1, length(ac) - 1, by = 2)) { pr <- ac[k] + ac[k + 1]; if (pr < 0) break; s <- s + pr }; length(v) / (1 + 2 * s) }
conv <- t(sapply(colnames(V[[1]]), function(v) { M <- sapply(V, function(X) X[, v]); c(rhat = rhat_rank(M), ess = sum(apply(M, 2, ess1)), mean = mean(M)) }))
print(signif(conv, 4)); cat(sprintf("max R-hat %.4f, min effective sample size %.0f over four chains\n", max(conv[, "rhat"], na.rm = TRUE), min(conv[, "ess"], na.rm = TRUE)))

## ---- 3. serial dependence ------------------------------------------------------------------------
l1 <- function(v) { m <- mean(v); sum(sapply(split(v, season), function(u) sum((u[-1] - m) * (u[-length(u)] - m)))) / sum((v - m)^2) }
cat(sprintf("lag-one autocorrelation within seasons: PM2.5 %.2f, O3 %.2f\n", l1(a$PM2.5), l1(a$O3)))
## intervals estimator of the extremal index (Ferro and Segers, 2003), interexceedance times within seasons
theta_fs <- function(ex) { Ti <- unlist(lapply(split(which(ex), season[ex]), diff)); N1 <- length(Ti); if (N1 < 2) return(NA)
  th <- if (max(Ti) <= 2) 2 * sum(Ti)^2 / (N1 * sum(Ti^2)) else 2 * sum(Ti - 1)^2 / (N1 * sum((Ti - 1) * (Ti - 2))); min(1, th) }
D <- x$fit$draws; th <- lapply(round(seq(1, nrow(D), length.out = 400)), function(k) star:::theta_row(x$fit, k))
thc <- sapply(th, function(t) theta_fs(r > rho_star(w, t$u1, t$u2, t$p)))
cat(sprintf("extremal index, exceedances of the threshold curve: posterior median %.2f, 90%% (%.2f, %.2f)\n", median(thc), quantile(thc, .05), quantile(thc, .95)))
for (q in c(0.8, 0.9)) { t <- -log(1 - q); je <- E[, 1] > t & E[, 2] > t; cat(sprintf("extremal index, joint exceedances of the %.1f quantiles: %.2f (%d exceedances)\n", q, theta_fs(je), sum(je))) }
## alternate days
H <- lapply(c("odd", "even"), function(tg) readRDS(sprintf("results/wanliu_%s.rds", tg)))
chi90 <- function(fit) { D <- fit$draws; mean(sapply(round(seq(1, nrow(D), length.out = 400)), function(k) { t <- star:::theta_row(fit, k)
  star:::joint_surv_cpp(-log(.1), -log(.1), gl, glw, t$phi, t$alpha, t$b, t$u1, t$u2, t$p, t$chibar, t$z) / 0.1 })) }
cat(sprintf("chi(0.9), posterior mean: all days %.3f, odd days %.3f, even days %.3f\n", chi90(x$fit), chi90(H[[1]]$fit), chi90(H[[2]]$fit)))

## ---- 4. moving-block (ten-day blocks within seasons) against iid bootstrap of the posterior means ---
if (Sys.getenv("STAR_BOOT", "1") != "0") {
  dir.create("results/boot", showWarnings = FALSE); B <- as.integer(Sys.getenv("STAR_NBOOT", "60")); L <- 10
  starts <- unlist(lapply(split(seq_len(n), season), function(ix) ix[1:(length(ix) - L + 1)]))
  draw_idx <- function(type) if (type == "iid") sample.int(n, n, replace = TRUE) else head(as.vector(outer(0:(L - 1), sample(starts, ceiling(n / L), replace = TRUE), "+")), n)
  summ <- function(fit, EE) { D <- fit$draws; rr <- rowSums(EE); ww <- EE[, 1] / rr
    v <- t(sapply(round(seq(1, nrow(D), length.out = 300)), function(k) { t <- star:::theta_row(fit, k)
      js <- function(l1, l2) star:::joint_surv_cpp(-log(1-l1), -log(1-l2), gl, glw, t$phi, t$alpha, t$b, t$u1, t$u2, t$p, t$chibar, t$z)
      c(chi90 = js(.9, .9) / .1, chi95 = js(.95, .95) / .05, chi99 = js(.99, .99) / .01, p9990 = js(.99, .9), exc = mean(rr > rho_star(ww, t$u1, t$u2, t$p)),
        AI = as.numeric(fit$component[k] != "one"), chibar = t$chibar) }))
    rbind(mean = colMeans(v), sd = apply(v, 2, sd)) }
  bj <- c(lapply(1:B, function(b) list(type = "blk", b = b)), lapply(1:B, function(b) list(type = "iid", b = b)))
  bj <- Filter(function(j) !file.exists(sprintf("results/boot/%s_%03d.rds", j$type, j$b)), bj)
  invisible(mclapply(bj, function(j) { set.seed(5000 + j$b + 1000 * (j$type == "iid")); aa <- a[draw_idx(j$type), ]; EE <- standardize(aa$PM2.5, aa$O3)$E
    fit <- tryCatch(star_mcmc(EE, prior = star_prior("threecomp", p_min = 1), n_iter = 10000, burnin = 4000, thin = 3, verbose = FALSE), error = function(e) NULL)
    if (!is.null(fit)) saveRDS(list(type = j$type, b = j$b, S = summ(fit, EE)), sprintf("results/boot/%s_%03d.rds", j$type, j$b)) }, mc.cores = CORES, mc.preschedule = FALSE))
  Lb <- lapply(list.files("results/boot", full.names = TRUE), readRDS); ty <- sapply(Lb, `[[`, "type")
  M <- t(sapply(Lb, function(z) z$S["mean", ])); SD <- t(sapply(Lb, function(z) z$S["sd", ]))
  sb <- apply(M[ty == "blk", ], 2, sd); si <- apply(M[ty == "iid", ], 2, sd); ps <- colMeans(SD)
  cat(sprintf("bootstrap replicates: %d block, %d iid\n", sum(ty == "blk"), sum(ty == "iid")))
  print(round(rbind(sd_block = sb, sd_iid = si, block_over_iid = sb / si, posterior_sd = ps, block_over_posterior = sb / ps), 3))
}

## ---- 5. tails of the margins implied by the fitted model, against the exponential -----------------
wf <- seq(0.0005, 0.9995, length.out = 1999); dw <- wf[2] - wf[1]; xs <- 2:7
surv <- function(t, rr, lam_, rho_, sig_) { pi_ <- pgamma(rho_, shape = t$alpha, rate = lam_)
  ifelse(rr > rho_, (1 - pi_) * exp(-(rr - rho_) / sig_), 1 - pgamma(rr, shape = t$alpha, rate = lam_)) }
Mg <- t(sapply(th[seq(1, 400, 2)], function(t) { K <- length(t$phi); fW <- rowSums(sapply(1:K, function(k) t$phi[k] * dbeta(wf, k, K - k + 1))); fW <- fW / sum(fW * dw)
  lam_ <- lambda_star(wf, t$b); rho_ <- rho_star(wf, t$u1, t$u2, t$p); sig_ <- sigma_star(wf, t$chibar, t$z)
  c(sapply(xs, function(x0) sum(fW * surv(t, x0 / wf, lam_, rho_, sig_)) * dw), sapply(xs, function(x0) sum(fW * surv(t, x0 / (1 - wf), lam_, rho_, sig_)) * dw)) }))
cat("P(E > x) / exp(-x) implied by the model at x =", xs, "\n  PM2.5:", round(colMeans(Mg[, seq_along(xs)]) / exp(-xs), 2),
    "\n  O3:   ", round(colMeans(Mg[, length(xs) + seq_along(xs)]) / exp(-xs), 2), "\n")
