# 00_common.R -- shared infrastructure for the main simulation studies
# SCALE: "test" (minutes, validates the pipeline), "beach" (overnight on ~8
# cores), "paper" (final runs). Set via Sys.setenv(STAR_SCALE = "beach").

library(star); library(parallel)
SCALE <- Sys.getenv("STAR_SCALE", "beach")
CORES <- as.integer(Sys.getenv("STAR_CORES", max(1, detectCores() - 1)))
SET <- switch(SCALE,
  test  = list(R = 2,  ITER = 1500,  BURN = 600,   THIN = 3, N1 = 400,  N2 = 800),
  beach = list(R = 50, ITER = 30000, BURN = 10000, THIN = 5, N1 = 1000, N2 = 4000),
  paper = list(R = 200, ITER = 30000, BURN = 10000, THIN = 5, N1 = 1000, N2 = 4000))
dir.create("results", showWarnings = FALSE)

run_par <- function(X, FUN, cores = CORES) {
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(cores); on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, { library(star); source("00_common.R") })
    clusterExport(cl, ls(globalenv()), envir = globalenv())
    parLapply(cl, X, FUN)
  } else mclapply(X, FUN, mc.cores = cores)
}

## checkpointed cell runner: results/<study>_<cell>_<chunk>.rds, resume-safe
run_cell <- function(study, cell, nrep, FUN, chunk = 10) {
  chunks <- split(seq_len(nrep), ceiling(seq_len(nrep) / chunk))
  out <- list()
  for (k in seq_along(chunks)) {
    f <- sprintf("results/%s_%s_c%02d.rds", study, cell, k)
    if (file.exists(f)) { out[[k]] <- readRDS(f); next }
    res <- run_par(chunks[[k]], FUN)
    bad <- !vapply(res, is.data.frame, logical(1))
    if (any(bad)) stop("worker error in ", study, " ", cell, ": ",
                       conditionMessage(attr(res[[which(bad)[1]]], "condition")))
    res <- do.call(rbind, res); saveRDS(res, f); out[[k]] <- res
    cat(sprintf("[%s] %s chunk %d/%d done (%s)\n", study, cell, k,
                length(chunks), format(Sys.time(), "%H:%M")))
  }
  do.call(rbind, out)
}

## ---- data generators on standard exponential margins ----------------------
theta_model <- function(chibar = 0.6, p = 2, u1 = 2.2, u2 = 2.2, nubar = 0.5)
  list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0),
       u1 = u1, u2 = u2, p = p, chibar = chibar, z = rep(log(nubar), 3))

sim_gauss <- function(n, rho) {
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  cbind(-log(pnorm(-z1)), -log(pnorm(-z2)))
}
rstable_pos <- function(n, a) {          # positive stable, LT exp(-s^a)
  U <- runif(n, 0, pi); E <- rexp(n)
  (sin(a * U) / sin(U)^(1 / a)) * (sin((1 - a) * U) / E)^((1 - a) / a)
}
sim_gumbel_uv <- function(n, theta) {    # Gumbel copula uniforms
  S <- rstable_pos(n, 1 / theta)
  cbind(exp(-(rexp(n) / S)^(1 / theta)), exp(-(rexp(n) / S)^(1 / theta)))
}
sim_logistic <- function(n, alphaL) {    # AD logistic on exponential margins
  uv <- sim_gumbel_uv(n, 1 / alphaL); -log(1 - uv)
}
sim_invlogistic <- function(n, alphaIL) { # AI inverted logistic
  uv <- sim_gumbel_uv(n, 1 / alphaIL); -log(uv)
}
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

## ---- true tail quantities per generator ------------------------------------
## chi(q) = P(U > q, V > q) / (1 - q); closed forms where available, else MC
true_chi_q <- function(gen, par, q, nmc = 2e7) {
  if (gen == "frank") {
    th <- par; Cuu <- -1/th * log(1 + (exp(-th*q) - 1)^2 / (exp(-th) - 1))
    return((1 - 2*q + Cuu) / (1 - q))
  }
  if (gen %in% c("logistic", "invlogistic")) {
    theta <- 1 / par
    if (gen == "logistic") {
      Cuu <- exp(-(2 * (-log(q))^theta)^(1/theta))
      return((1 - 2*q + Cuu) / (1 - q))
    } else {
      t <- -log(1 - q)                  # exponential quantile at level q
      return(exp(-(2 * t^theta)^(1/theta)) / (1 - q))
    }
  }
  E <- switch(gen, gauss = sim_gauss(nmc, par), t = sim_t(nmc, par[1], par[2]))
  tq <- -log(1 - q)
  mean(E[, 1] > tq & E[, 2] > tq) / (1 - q)
}
true_eta <- function(gen, par) switch(gen,
  gauss = (1 + par) / 2, invlogistic = 2^(-par), frank = 0.5,
  logistic = 1, t = 1, model = (1 + par) / 2)
true_lambda_adf <- function(gen, par, om) {   # ADF at omega
  wg <- seq(0.001, 0.999, length.out = 1999)
  sig <- switch(gen,
    gauss = (1 - par^2) / (1 - 2 * par * sqrt(wg * (1 - wg))),
    invlogistic = (wg^(1/par) + (1 - wg)^(1/par))^(-par),
    frank = rep(1, length(wg)),
    logistic = NULL, t = NULL)
  if (is.null(sig)) return(max(om, 1 - om))
  min(pmax(om / wg, (1 - om) / (1 - wg)) / sig)
}

## ---- per-fit metrics -------------------------------------------------------
fit_metrics <- function(fit, wpts = c(0.15, 0.5, 0.85), q_chi = 0.99, level = 0.9) {
  Fm <- star:::theta_row; S <- nrow(fit$draws)
  gl <- seq(0.5/201, 1 - 0.5/201, length.out = 201); glw <- rep(1/201, 201)
  idx <- round(seq(1, S, length.out = min(150, S)))
  tq <- -log(1 - q_chi)
  vals <- t(sapply(idx, function(s) {
    th <- Fm(fit, s)
    c(chibar = th$chibar, eta = (1 + th$chibar)/2,
      sigma_half = 1 + th$chibar,
      setNames(rho_star(wpts, th$u1, th$u2, th$p), paste0("rho", wpts)),
      p = th$p,
      chi_q = star:::joint_surv_cpp(tq, tq, gl, glw, th$phi, th$alpha, th$b,
                             th$u1, th$u2, th$p, th$chibar, th$z) / (1 - q_chi),
      lam03 = { wg <- seq(0.002, 0.998, length.out = 499)
                sg <- sigma_star(wg, th$chibar, th$z)
                min(pmax(0.3/wg, 0.7/(1 - wg)) / sg) })
  }))
  a <- (1 - level)/2
  data.frame(quantity = colnames(vals),
             med = apply(vals, 2, median),
             lo = apply(vals, 2, quantile, a), hi = apply(vals, 2, quantile, 1 - a))
}
