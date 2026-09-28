## Heffernan-Tawn (2004) conditional extremes benchmark for chi(q) = P(E1 > t_q, E2 > t_q) / (1 - q).
## Applied pipeline: semiparametric margins (empirical below the kappa quantile, GPD above),
## Laplace scale, conditional model Y_{-j} | Y_j > u  =  a Y_j + Y_j^b Z with the empirical
## residual distribution, both conditioning directions averaged, chi(q) by simulation.
lap_q <- function(p) ifelse(p < 0.5, log(2 * p), -log(2 * (1 - p)))
lap_surv <- function(y) ifelse(y < 0, 1 - 0.5 * exp(y), 0.5 * exp(-y))

gpd_fit <- function(exc) {
  nll <- function(par) {
    s <- exp(par[1]); xi <- par[2]; z <- 1 + xi * exc / s
    if (any(z <= 0)) return(1e10)
    length(exc) * log(s) + (1 + 1 / xi) * sum(log(z))
  }
  nll0 <- function(par) { s <- exp(par[1]); length(exc) * log(s) + sum(exc) / s }
  o <- optim(c(log(mean(exc)), 0.1), nll, method = "L-BFGS-B", lower = c(-10, -0.49), upper = c(10, 1.5))
  c(sigma = exp(o$par[1]), xi = o$par[2])
}
## semiparametric marginal cdf: returns list(F = function, u, sigma, xi)
semipar <- function(x, kappa = 0.9) {
  n <- length(x); u <- quantile(x, kappa, names = FALSE); g <- gpd_fit(x[x > u] - u)
  Fx <- function(t) {
    p <- ecdf(x)(t) * n / (n + 1)
    hi <- t > u
    p[hi] <- 1 - (1 - kappa) * pmax(1 + g["xi"] * (t[hi] - u) / g["sigma"], 0)^(-1 / g["xi"])
    pmin(p, 1 - 0.5 / (n + 1))               # keep the Laplace transform finite beyond a finite GPD endpoint
  }
  list(F = Fx, u = u, sigma = g[["sigma"]], xi = g[["xi"]])
}
ht_fit <- function(y1, y2, u, start = NULL) {
  idx <- y1 > u; x <- y1[idx]; y <- y2[idx]
  nll <- function(par) {
    a <- par[1]; b <- par[2]; m <- par[3]; s <- exp(par[4]); xb <- x^b
    -sum(dnorm(y, a * x + m * xb, s * xb, log = TRUE))
  }
  best <- NULL
  starts <- if (is.null(start)) list(c(0.5, 0.2, 0, 0), c(0.9, 0.0, 0, 0), c(0.1, 0.5, 0, 0)) else list(start)
  for (st in starts) {
    o <- try(optim(st, nll, method = "L-BFGS-B", lower = c(-1, -5, -20, -10), upper = c(1, 1, 20, 10)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o
  }
  a <- best$par[1]; b <- best$par[2]
  list(a = a, b = b, z = (y - a * x) / x^b, par = best$par)
}
## P(Y1 > v1, Y2 > v2) from the conditional model fitted on Y1 > u, by simulation (v1 >= u)
ht_joint <- function(fit, v1, v2, nsim) {
  x <- v1 + rexp(nsim)                       # Laplace tail above v1 > 0 is unit exponential
  y <- fit$a * x + x^fit$b * sample(fit$z, nsim, replace = TRUE)
  lap_surv(v1) * mean(y > v2)
}
## chi(q) for a bivariate sample E on the data scale, at t_q = -log(1 - q)
ht_chi <- function(E, qs, kappa = 0.9, nsim = 20000, starts = NULL, fits = FALSE) {
  m1 <- semipar(E[, 1], kappa); m2 <- semipar(E[, 2], kappa)
  y1 <- lap_q(m1$F(E[, 1])); y2 <- lap_q(m2$F(E[, 2]))
  u <- lap_q(kappa)
  f12 <- ht_fit(y1, y2, u, starts[[1]]); f21 <- ht_fit(y2, y1, u, starts[[2]])
  if (fits) return(list(f12$par, f21$par))
  sapply(qs, function(q) {
    tq <- -log(1 - q); v1 <- lap_q(m1$F(tq)); v2 <- lap_q(m2$F(tq))
    v1 <- max(v1, u); v2 <- max(v2, u)
    0.5 * (ht_joint(f12, v1, v2, nsim) + ht_joint(f21, v2, v1, nsim)) / (1 - q)
  })
}
ht_boot <- function(E, qs, B, seed, kappa = 0.9, nsim = 10000, level = 0.9) {
  set.seed(seed); n <- nrow(E); st <- ht_chi(E, qs, kappa, fits = TRUE)   # bootstrap fits start from the full-sample fit
  b <- replicate(B, ht_chi(E[sample.int(n, n, replace = TRUE), , drop = FALSE], qs, kappa, nsim, starts = st))
  b <- matrix(b, nrow = length(qs)); a <- (1 - level) / 2
  cbind(lo = apply(b, 1, quantile, a, na.rm = TRUE), hi = apply(b, 1, quantile, 1 - a, na.rm = TRUE))
}
