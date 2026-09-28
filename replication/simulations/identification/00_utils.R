# 00_utils.R -- shared utilities for the identification study (Study 0)
# Requires: library(star)

## portable parallel apply: forks on unix, socket cluster on Windows
run_par <- function(X, FUN, cores = getOption("mc.cores", 4)) {
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { library(star); source("00_utils.R") })
    parallel::clusterExport(cl, ls(globalenv()), envir = globalenv())
    parallel::parLapply(cl, X, FUN)
  } else {
    parallel::mclapply(X, FUN, mc.cores = cores)
  }
}

library(star)

## ---- default true model and basis sizes -----------------------------------
K_W <- 5; K_lam <- 2; K_nu <- 2
theta_true <- function(chibar = 0.6, p = 2, u = 2.2, nubar = 0.5)
  list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = rep(0, K_lam + 1),
       u1 = u, u2 = u, p = p, chibar = chibar, z = rep(log(nubar), K_nu + 1))

## ---- packing on the sampling (unconstrained) scale ------------------------
## order: psi (K_W-1), la, b (K_lam+1), lo1, lo2, lz, lchi, z (K_nu+1)
par_names <- c(paste0("psi", 1:(K_W - 1)), "log_alpha", paste0("b", 0:K_lam),
               "log_o1", "log_o2", "log_p", "logit_chibar",
               paste0("z", 0:K_nu))

pack_theta <- function(theta, q0) {
  psi <- log(theta$phi[-K_W] / theta$phi[K_W])
  c(psi, log(theta$alpha), theta$b,
    log(theta$u1 - q0), log(theta$u2 - q0),
    log(theta$p), qlogis(theta$chibar), theta$z)
}

unpack_theta <- function(x, q0) {
  i <- 0; take <- function(k) { out <- x[(i + 1):(i + k)]; i <<- i + k; out }
  psi <- take(K_W - 1); la <- take(1); b <- take(K_lam + 1)
  lo1 <- take(1); lo2 <- take(1); lp <- take(1); lchi <- take(1)
  z <- take(K_nu + 1)
  e <- exp(c(psi, 0)); phi <- e / sum(e)
  list(phi = phi, alpha = exp(la), b = b,
       u1 = q0 + exp(lo1), u2 = q0 + exp(lo2),
       p = exp(lp), chibar = plogis(lchi), z = z)
}

## log-likelihood as a function of the packed vector, for fixed data
make_ll <- function(r, w, q0) function(x) loglik_star(r, w, unpack_theta(x, q0))

## ---- numerical Hessian (central differences) ------------------------------
num_hessian <- function(f, x, h = 1e-3) {
  ## h may be a vector of per-coordinate steps. The log-likelihood is
  ## DISCONTINUOUS in the threshold coordinates (observations cross the
  ## curve), so those coordinates need wide steps that average across many
  ## crossings; the result is the curvature of the smoothed envelope, the
  ## object the posterior actually responds to.
  d <- length(x); if (length(h) == 1) h <- rep(h, d)
  H <- matrix(NA_real_, d, d); f0 <- f(x)
  for (i in 1:d) {
    ei <- rep(0, d); ei[i] <- h[i]
    H[i, i] <- (f(x + ei) - 2 * f0 + f(x - ei)) / h[i]^2
    if (i < d) for (j in (i + 1):d) {
      ej <- rep(0, d); ej[j] <- h[j]
      H[i, j] <- H[j, i] <-
        (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) /
        (4 * h[i] * h[j])
    }
  }
  H
}

## ---- functionals of interest and their true values ------------------------
functionals <- function(theta, wpts = c(0.15, 0.5, 0.85), M = 401) {
  g <- seq(0.5 / M, 1 - 0.5 / M, length.out = M)
  exc <- 1 - sum(fW_star(g, theta$phi) *
                 pgamma(rho_star(g, theta$u1, theta$u2, theta$p),
                        shape = theta$alpha,
                        rate = lambda_star(g, theta$b))) / M
  c(setNames(rho_star(wpts, theta$u1, theta$u2, theta$p), paste0("rho", wpts)),
    sigma_half = sigma_star(0.5, theta$chibar, theta$z),
    nu_half = nu_star(0.5, theta$z),
    chibar = theta$chibar, alpha = theta$alpha, p = theta$p,
    u1 = theta$u1, u2 = theta$u2, exc_frac = exc)
}

fit_functionals <- function(fit, wpts = c(0.15, 0.5, 0.85)) {
  S <- nrow(fit$draws)
  t(sapply(seq_len(S), function(s) {
    th <- star:::theta_row(fit, s); functionals(th, wpts)
  }))
}

## width / coverage summary for one replicate
summarize_rep <- function(fit, truth, level = 0.9) {
  Fm <- fit_functionals(fit)
  lo <- apply(Fm, 2, quantile, (1 - level) / 2)
  hi <- apply(Fm, 2, quantile, 1 - (1 - level) / 2)
  data.frame(quantity = colnames(Fm), post_sd = apply(Fm, 2, sd),
             width = hi - lo, cover = as.numeric(truth >= lo & truth <= hi),
             post_med = apply(Fm, 2, median), truth = truth)
}
