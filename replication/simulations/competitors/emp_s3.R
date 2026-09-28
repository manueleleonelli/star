## Empirical and Hill benchmarks on the exact datasets of Study 3 (regenerated from the stored seeds)
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
gens <- list(gauss05 = function(n) sim_gauss(n, 0.5), invlog05 = function(n) sim_invlogistic(n, 0.5),
             frank5 = function(n) sim_frank(n, 5), t05_4 = function(n) sim_t(n, 0.5, 4))
hill_eta <- function(E, kappa = 0.9, level = 0.9) {
  Tmin <- pmin(E[, 1], E[, 2]); u <- quantile(Tmin, kappa); exc <- Tmin[Tmin > u] - u; k <- length(exc)
  eta <- mean(exc); se <- eta / sqrt(k); z <- qnorm(1 - (1 - level) / 2)
  c(eta = eta, lo = eta - z * se, hi = eta + z * se)
}
QS <- c(0.95, 0.99, 0.999)
s3 <- readRDS("../main/results/s3_all.rds"); tr <- readRDS("../main/results/s3_truths.rds")
out <- do.call(rbind, lapply(seq_len(nrow(s3)), function(i) {
  g <- as.character(s3$gen[i]); seed <- s3$seed[i]; n <- s3$n[i]
  set.seed(seed); E <- gens[[g]](n)
  U <- rank(E[, 1]) / (n + 1); V <- rank(E[, 2]) / (n + 1)
  emp_rank  <- sapply(QS, function(q) mean(U > q & V > q) / (1 - q))
  emp_known <- sapply(QS, function(q) { t <- -log(1 - q); mean(E[, 1] > t & E[, 2] > t) / (1 - q) })
  h <- hill_eta(E)
  data.frame(gen = g, seed = seed,
             er95 = emp_rank[1] - tr[[g]]$chi[1], er99 = emp_rank[2] - tr[[g]]$chi[2], er999 = emp_rank[3] - tr[[g]]$chi[3],
             ek95 = emp_known[1] - tr[[g]]$chi[1], ek99 = emp_known[2] - tr[[g]]$chi[2], ek999 = emp_known[3] - tr[[g]]$chi[3],
             hill_err = h["eta"] - tr[[g]]$eta, hill_cover = as.numeric(tr[[g]]$eta >= h["lo"] & tr[[g]]$eta <= h["hi"]),
             model_95 = s3$chi95_err[i], model_99 = s3$chi99_err[i], model_999 = s3$chi999_err[i],
             model_eta = s3$eta_err[i], model_eta_cover = s3$eta_cover[i])
}))
saveRDS(out, "results/emp_s3.rds")
rmse <- function(x) sqrt(mean(x^2))
agg <- function(f) aggregate(. ~ gen, out[, setdiff(names(out), "seed")], f)
cat("=== bias ===\n"); print(agg(mean), digits = 2)
cat("=== RMSE ===\n"); print(agg(rmse), digits = 2)
