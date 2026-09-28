## gens.R -- copula generators on standard exponential margins (Gaussian, inverted logistic, Frank, Student t), as in ../main/00_common.R
sim_gauss <- function(n, rho) {
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  cbind(-log(pnorm(-z1)), -log(pnorm(-z2)))
}
rstable_pos <- function(n, a) {
  U <- runif(n, 0, pi); E <- rexp(n)
  (sin(a * U) / sin(U)^(1 / a)) * (sin((1 - a) * U) / E)^((1 - a) / a)
}
sim_gumbel_uv <- function(n, theta) {
  S <- rstable_pos(n, 1 / theta)
  cbind(exp(-(rexp(n) / S)^(1 / theta)), exp(-(rexp(n) / S)^(1 / theta)))
}
sim_invlogistic <- function(n, alphaIL) { uv <- sim_gumbel_uv(n, 1 / alphaIL); -log(uv) }
sim_frank <- function(n, th) {
  u <- runif(n); t <- runif(n)
  v <- -1 / th * log(1 + t * (1 - exp(-th)) / (t * (exp(-th * u) - 1) - exp(-th * u)))
  cbind(-log(1 - u), -log(1 - v))
}
sim_t <- function(n, rho, df) {
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  s <- sqrt(rchisq(n, df) / df)
  cbind(-log(pt(-z1 / s, df)), -log(pt(-z2 / s, df)))
}
