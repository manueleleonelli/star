## s1_numerical_checks.R -- the numerical checks quoted after Proposition S3 in Section S1.2 of the supplement.
## (1) Order d - 1 of Laplace's method at d = 3 with sigma = 1 and B = 1. By symmetry the joint integral
##     int_S exp(-t max_j 1/w_j) dw equals 3 int_0^{1/3} exp(-t/w) (1 - 3w) dw, and the marginal integral
##     int_S exp(-t/w_1) dw equals int_0^1 exp(-t/w) (1 - w) dw; both are computed on the scale exp(3t) or exp(t).
##     t^2 e^{3t} times the joint integral tends to 1/9, t e^{3t} times it to zero, and t^2 e^t times the marginal
##     integral to 1.
## (2) The 0.99 quantiles of the three margins of the trivariate model at the parameter values of the calibration
##     check of Section S4.9, by quadrature on a fine grid of the simplex, against 4.61 for the standard exponential.
## Run from this folder: Rscript s1_numerical_checks.R (under a minute).

cat("joint integral, sigma = 1, B = 1 (t^2 e^{3t} I tends to 1/9 = 0.111):\n")
for (t in c(20, 40, 80, 160, 320)) {
  v <- integrate(function(w) 3 * exp(-t * (1 / w - 3)) * (1 - 3 * w), 0, 1 / 3, rel.tol = 1e-12, subdivisions = 10000)$value
  cat(sprintf("  t = %3d: t^2 e^{3t} I = %.3f, t e^{3t} I = %.5f\n", t, t^2 * v, t * v))
}
cat("marginal integral (t^2 e^t I tends to 1):\n")
for (t in c(20, 40, 80, 160, 320)) {
  v <- integrate(function(w) exp(-t * (1 / w - 1)) * (1 - w), 0, 1, rel.tol = 1e-12, subdivisions = 10000)$value
  cat(sprintf("  t = %3d: t^2 e^t I = %.2f\n", t, t^2 * v))
}

## (2) margins of the trivariate model at the calibration truth
source(system.file("trivariate", "star3.R", package = "star"))
K <- 4; ix <- basis_index(K); phi <- rep(1, nrow(ix)); phi[ix[, 3] == 0] <- 1.6; phi <- phi / sum(phi)
th0 <- list(phi = phi, alpha = 3, b = c(0.1, -0.1, 0), u = c(2.2, 2.2, 2.2), p = 2,
            chib = c(0.6, 0.2, 0.2), nu = c(0.5, 0.5, 0.5))
G <- simplex_grid(1500); W <- G$W
fw <- fW_grid(th0, K, G); lam <- lambda3(W, th0$b); rho <- rho3(W, th0$u, th0$p); sg <- sigma3(W, th0$chib, th0$nu)
surv_R <- function(r) ifelse(r <= rho, pgamma(r, th0$alpha, rate = lam, lower.tail = FALSE),
                             pgamma(rho, th0$alpha, rate = lam, lower.tail = FALSE) * exp(-(r - rho) / sg))
q99 <- sapply(1:3, function(j) uniroot(function(t) sum(fw * surv_R(t / W[, j])) / sum(fw) - 0.01, c(2, 6), tol = 1e-8)$root)
cat(sprintf("0.99 quantiles of the three margins: %.2f %.2f %.2f (standard exponential %.2f)\n", q99[1], q99[2], q99[3], -log(0.01)))
