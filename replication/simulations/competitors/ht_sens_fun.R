## HT sensitivity: margins (semiparametric at 0.9 vs true) x conditioning threshold kappa_d in {0.8,0.85,0.9,0.95},
## on the same datasets as Sections S3.2, S3.4 and S3.6 (regenerated from stored seeds). Point estimates only.
suppressMessages(library(star)); library(parallel)
source("ht.R"); source("gens.R")
theta_model <- function(chibar = 0.6, p = 2, u1 = 2.2, u2 = 2.2, nubar = 0.5)
  list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0), u1 = u1, u2 = u2, p = p, chibar = chibar, z = rep(log(nubar), 3))
KD <- c(0.8, 0.85, 0.9, 0.95); gl <- seq(0.5/201, 1 - 0.5/201, length.out = 201); gw <- rep(1/201, 201)
## true marginal cdfs of the model, by quadrature on a grid, log-survivor interpolated
true_margins <- function(th) {
  tg <- seq(0, 25, by = 0.01)
  S1 <- sapply(tg, function(t) star:::joint_surv_cpp(t, 0, gl, gw, th$phi, th$alpha, th$b, th$u1, th$u2, th$p, th$chibar, th$z))
  S2 <- sapply(tg, function(t) star:::joint_surv_cpp(0, t, gl, gw, th$phi, th$alpha, th$b, th$u1, th$u2, th$p, th$chibar, th$z))
  f <- function(S) { l <- approxfun(tg, log(pmax(S, 1e-300)), rule = 2); function(x) 1 - exp(l(x)) }
  list(F1 = f(S1), F2 = f(S2))
}
exp_margins <- list(F1 = function(x) 1 - exp(-x), F2 = function(x) 1 - exp(-x))
## HT chi(q) with given marginal cdfs (F1, F2) and conditioning threshold kd
ht_chi_F <- function(E, qs, F1, F2, kd, nsim = 20000, starts = NULL, fits = FALSE) {
  cl <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)
  y1 <- lap_q(cl(F1(E[, 1]))); y2 <- lap_q(cl(F2(E[, 2]))); u <- lap_q(kd)
  f12 <- ht_fit(y1, y2, u, starts[[1]]); f21 <- ht_fit(y2, y1, u, starts[[2]])
  if (fits) return(list(f12$par, f21$par))
  sapply(qs, function(q) { tq <- -log(1 - q); v1 <- max(lap_q(F1(tq)), u); v2 <- max(lap_q(F2(tq)), u)
    0.5 * (ht_joint(f12, v1, v2, nsim) + ht_joint(f21, v2, v1, nsim)) / (1 - q) })
}
one <- function(E, qs, Ftrue, sd) {
  m1 <- semipar(E[, 1], 0.9); m2 <- semipar(E[, 2], 0.9)
  do.call(rbind, lapply(KD, function(kd) { set.seed(sd + 2e6)
    s <- ht_chi_F(E, qs, m1$F, m2$F, kd); set.seed(sd + 3e6); k <- ht_chi_F(E, qs, Ftrue$F1, Ftrue$F2, kd)
    data.frame(kd = kd, q = qs, semi = s, known = k) }))
}
