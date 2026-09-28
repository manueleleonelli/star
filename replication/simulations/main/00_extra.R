# 00_extra.R -- helpers for the clean 50-replicate run. Source AFTER 00_common.R.
# All studies share the settings in SET (from 00_common.R); SET_NEW is kept as
# an alias so the study scripts read uniformly.
SET_NEW <- SET

## Hill-type estimator of eta from the structure variable min(E1, E2).
## On exponential margins P(min > t) ~ L e^{-t/eta}, so the mean excess of
## min(E1,E2) above a high threshold estimates eta directly; the normal
## approximation gives a simple interval.
hill_eta <- function(E, kappa = 0.9, level = 0.9) {
  Tmin <- pmin(E[, 1], E[, 2])
  u <- quantile(Tmin, kappa); exc <- Tmin[Tmin > u] - u; k <- length(exc)
  eta <- mean(exc); se <- eta / sqrt(k); z <- qnorm(1 - (1 - level) / 2)
  c(eta = eta, lo = eta - z * se, hi = eta + z * se, k = k)
}

## Non-exchangeable, asymmetric, asymptotically dependent generator via the
## Khoudraji device with a Gumbel base: C(u,v) = (uv)^{...} with unequal
## exponents a != b makes the extremal dependence angularly asymmetric.
sim_khoudraji <- function(n, theta_gumbel = 2, a = 0.9, b = 0.5) {
  UI <- cbind(runif(n), runif(n))          # independence part
  VG <- sim_gumbel_uv(n, theta_gumbel)     # Gumbel part (AD)
  U <- pmax(UI[, 1]^(1 / (1 - a)), VG[, 1]^(1 / a))
  V <- pmax(UI[, 2]^(1 / (1 - b)), VG[, 2]^(1 / b))
  cbind(-log(1 - U), -log(1 - V))
}

## Within-model asymmetric truth: angle-varying shape and unequal anchors.
theta_asym <- function(chibar = 0.6)
  list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0),
       u1 = 2.8, u2 = 1.7, p = 2, chibar = chibar,
       z = c(log(0.15), log(0.8), log(4)))   # increasing nu(.) across the angle

## Convergence diagnostics from a list of chains (matrices, rows = draws).
gelman_rubin <- function(chains) {          # classic potential scale reduction
  m <- length(chains); nn <- nrow(chains[[1]]); p <- ncol(chains[[1]])
  out <- numeric(p)
  for (j in seq_len(p)) {
    xj <- sapply(chains, function(C) C[, j])
    means <- colMeans(xj); B <- nn * var(means)
    W <- mean(apply(xj, 2, var))
    out[j] <- sqrt(((nn - 1) / nn * W + B / nn) / W)
  }
  setNames(out, colnames(chains[[1]]))
}
ess_bm <- function(x) {                     # batch-means effective sample size
  S <- length(x); B <- max(1, floor(sqrt(S))); nb <- floor(S / B)
  if (nb < 2) return(NA_real_)
  bm <- colMeans(matrix(x[seq_len(B * nb)], nrow = B))
  tau <- max(1, B * var(bm) / var(x)); S / tau
}
