## calib.R -- small calibration check of the trivariate implementation (Section S4.9 and Table S12 of the supplement)
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
source(system.file("trivariate", "star3.R", package = "star")); library(parallel)
K <- 4; ix <- basis_index(K); phi <- rep(1, nrow(ix)); phi[ix[, 3] == 0] <- 1.6; phi <- phi / sum(phi)   # tilted towards the (1,2) edge
th0 <- list(phi = phi, alpha = 3, b = c(0.1, -0.1, 0), u = c(2.2, 2.2, 2.2), p = 2, chib = c(0.6, 0.2, 0.2), nu = c(0.5, 0.5, 0.5))
dir.create("results/calib", recursive = TRUE, showWarnings = FALSE); saveRDS(th0, "results/th0.rds"); n <- 2000; NREP <- as.integer(Sys.getenv("NREP", "20"))
truth <- c(eta12 = eta3(th0, 1:2), eta13 = eta3(th0, c(1, 3)), eta23 = eta3(th0, 2:3), eta123 = eta3(th0),
           chi12_99 = chi3(th0, K, .99, 1:2), chi13_99 = chi3(th0, K, .99, c(1, 3)), chi123_99 = chi3(th0, K, .99),
           chi123_9 = chi3(th0, K, .9), exc = exc_frac3(th0, K), rho_c = rho3(matrix(1/3, 1, 3), th0$u, th0$p),
           u1 = th0$u[1], u2 = th0$u[2], u3 = th0$u[3], p = th0$p, chib12 = th0$chib[1], chib13 = th0$chib[2], chib23 = th0$chib[3])
one <- function(r) {
  f <- sprintf("results/calib/rep%02d.rds", r); if (file.exists(f)) return(readRDS(f))
  set.seed(1000 + r); E <- rstar3(n, th0, K)$E; t0 <- Sys.time()
  fit <- mcmc3(E, n_iter = 20000, burnin = 8000, thin = 5, K = K, seed = 2000 + r, verbose = FALSE)
  mins <- as.numeric(Sys.time() - t0, units = "mins")
  D <- fit$draws; idx <- round(seq(1, nrow(D), length.out = 200))
  S <- t(sapply(idx, function(s) { th <- unpack(D[s, ], K)
    c(eta12 = eta3(th, 1:2), eta13 = eta3(th, c(1, 3)), eta23 = eta3(th, 2:3), eta123 = eta3(th),
      chi12_99 = chi3(th, K, .99, 1:2), chi13_99 = chi3(th, K, .99, c(1, 3)), chi123_99 = chi3(th, K, .99),
      chi123_9 = chi3(th, K, .9), exc = exc_frac3(th, K), rho_c = rho3(matrix(1/3, 1, 3), th$u, th$p),
      u1 = th$u[1], u2 = th$u[2], u3 = th$u[3], p = th$p, chib12 = th$chib[1], chib13 = th$chib[2], chib23 = th$chib[3]) }))
  out <- list(rep = r, median = apply(S, 2, median), lo = apply(S, 2, quantile, .05), hi = apply(S, 2, quantile, .95), acc = fit$acc, mins = mins,
              ess_lp = tryCatch(coda::effectiveSize(fit$logpost), error = function(e) NA))
  saveRDS(out, f); out
}
dir.create("calib", showWarnings = FALSE)
res <- mclapply(seq_len(NREP), one, mc.cores = 2, mc.preschedule = FALSE)
med <- t(sapply(res, `[[`, "median")); lo <- t(sapply(res, `[[`, "lo")); hi <- t(sapply(res, `[[`, "hi"))
tab <- data.frame(truth = truth[colnames(med)], mean_median = colMeans(med), bias = colMeans(med) - truth[colnames(med)],
                  rmse = sqrt(colMeans(sweep(med, 2, truth[colnames(med)])^2)),
                  coverage = colMeans(sweep(lo, 2, truth[colnames(med)], "<=") & sweep(hi, 2, truth[colnames(med)], ">=")),
                  width = colMeans(hi - lo))
print(round(tab, 4)); cat("minutes per fit:", round(mean(sapply(res, `[[`, "mins")), 1), "\n"); print(round(colMeans(t(sapply(res, `[[`, "acc"))), 2))
saveRDS(list(tab = tab, res = res, truth = truth), "results/calib3.rds")
