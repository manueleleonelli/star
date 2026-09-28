## theory_check.R -- numerical check of Proposition 3: joint and marginal tail decay and the
## copula tail dependence coefficient at the boundary chibar = 1 (formula in Section S2).
## Integrals are computed on the scale exp(s) to avoid underflow. Needs star for rho_star.
suppressMessages(library(star))
check <- function(chibar, nu, alpha, u, p, fw, lamf) {
  tau <- function(w) pmin(w, 1 - w) / pmax(w, 1 - w); sig <- function(w) 1 + chibar * tau(w)^(1 + nu)
  rho <- function(w) rho_star(pmin(pmax(w, 1e-12), 1 - 1e-12), u[1], u[2], p)
  A <- function(w) 1 - pgamma(rho(w), alpha, rate = lamf(w)); B <- function(w) fw(w) * A(w) * exp(rho(w) / sig(w))
  Isc <- function(level, s, brk = 0.5) { f <- function(w) B(w) * exp(s - level(w) / sig(w))
    sum(sapply(list(c(1e-12, brk), c(brk, 1 - 1e-12)), function(ab) integrate(f, ab[1], ab[2], subdivisions = 20000, rel.tol = 1e-10, abs.tol = 0)$value)) }
  eta <- (1 + chibar) / 2; c0 <- 2 * eta^2 * B(0.5) / (1 + 3 * chibar + 2 * chibar * nu)
  kap <- if (chibar == 1) B(0.5) * (1 + nu) / (nu * (2 + nu)) else 0; c1 <- B(1 - 1e-12) + kap; c2 <- B(1e-12) + kap
  dl <- log(c1 / c2) / 4; I <- integrate(function(y) exp(-2 * abs(y) - 2 * (1 + nu) * abs(y + dl)), -Inf, Inf)$value
  cat(sprintf("chibar = %.2f: c = %.4f, c1 = %.4f, c2 = %.4f, copula chi = %.4f\n", chibar, c0, c1, c2, if (chibar == 1) B(0.5) * I / sqrt(c1 * c2) else 0))
  for (t in c(25, 100, 400)) cat(sprintf("  t = %3.0f: t e^{t/eta} P(E1>t,E2>t) = %.4f, t e^t P(E1>t) = %.4f, t e^t P(E2>t) = %.4f\n", t,
    t * Isc(function(w) t / pmin(w, 1 - w), t / eta), t * Isc(function(w) t / w, t), t * Isc(function(w) t / (1 - w), t)))
  for (s in c(25, 100, 400)) {
    t1 <- uniroot(function(t) log(Isc(function(w) t / w, s)), c(s - 2 * log(s) - 5, s + 5), tol = 1e-12)$root
    t2 <- uniroot(function(t) log(Isc(function(w) t / (1 - w), s)), c(s - 2 * log(s) - 5, s + 5), tol = 1e-12)$root
    cat(sprintf("  1 - u = exp(-%3.0f): copula P(U1>u, U2>u)/(1-u) = %.4f\n", s, Isc(function(w) pmax(t1 / w, t2 / (1 - w)), s, brk = t1 / (t1 + t2)))) }
}
lin <- function(w) 0.5 + w; lam <- function(w) exp(0.3 - 0.6 * w)
check(1, 0.7, 1.6, c(2.5, 3.5), 1.8, lin, lam)     # boundary: copula chi converges to the formula
check(0.5, 0.7, 1.6, c(2.5, 3.5), 1.8, lin, lam)   # interior: copula chi converges to zero
