## wanliu_joint.R -- the Wanliu application under the cut (stored fit) and under the joint posterior
source("joint.R")
marg_prior <- function(lsig, xi) dnorm(lsig, 0, 10, log = TRUE) + dnorm(xi, 0, 0.5, log = TRUE)   # vaguer scale prior for data in microg/m3
x <- readRDS("results/fit_Wanliu.rds"); a <- x$data; Y <- cbind(a$PM2.5, a$O3); n <- nrow(Y)
fc <- x$fit; gx <- x$margins$gx; gy <- x$margins$gy
prior <- star_prior("threecomp", p_min = 1)
if (!file.exists("results/fit_wanliu_joint.rds")) { set.seed(2026); t0 <- Sys.time()
  fj <- joint_mcmc(Y, prior, n_iter = 30000, burnin = 10000, thin = 5, verbose = TRUE); fj$minutes <- as.numeric(Sys.time() - t0, units = "mins")
  saveRDS(fj, "results/fit_wanliu_joint.rds") } else fj <- readRDS("results/fit_wanliu_joint.rds")
cat("joint minutes:", round(fj$minutes, 1), " acceptance:", round(fj$acc_rates, 2), "\n")
gl <- seq(0.5 / 201, 1 - 0.5 / 201, length.out = 201); glw <- rep(1 / 201, 201)
q90 <- function(v) quantile(v, c(0.05, 0.5, 0.95), names = FALSE)
Qsemi <- function(p, xdat, u, p_u, sigma, xi) { out <- quantile(xdat, pmin(p, 1 - p_u), names = FALSE, type = 6); ab <- p > 1 - p_u
  out[ab] <- u + sigma / xi * (((1 - p[ab]) / p_u)^(-xi) - 1); out }
summ <- function(fit, joint) {
  D <- fit$draws; idx <- round(seq(1, nrow(D), length.out = 400)); th <- lapply(idx, function(s) star:::theta_row(fit, s))
  cp <- c(zero = mean(fit$component == "zero"), slab = mean(fit$component == "slab"), one = mean(fit$component == "one"))
  chi <- sapply(c(0.9, 0.95, 0.99), function(q) { tt <- -log(1 - q); sapply(th, function(t) star:::joint_surv_cpp(tt, tt, gl, glw, t$phi, t$alpha, t$b, t$u1, t$u2, t$p, t$chibar, t$z)) / (1 - q) })
  exc <- sapply(th, function(t) { wg <- gl; rho <- rho_star(wg, t$u1, t$u2, t$p); lam <- lambda_star(wg, t$b)
    sum(fW_star(wg, t$phi) * pgamma(rho, shape = t$alpha, rate = lam, lower.tail = FALSE)) / length(wg) })
  ## one-season return level and threshold on the diagonal, mapped to the original scale draw by draw
  orig <- t(sapply(seq_along(idx), function(k) { t <- th[[k]]; w0 <- 0.5
    rho <- rho_star(w0, t$u1, t$u2, t$p); sg <- sigma_star(w0, t$chibar, t$z); lam <- lambda_star(w0, t$b); A <- pgamma(rho, shape = t$alpha, rate = lam, lower.tail = FALSE)
    pr <- 1 / 183; r1 <- if (A > pr) rho + sg * log(A / pr) else qgamma(1 - pr, shape = t$alpha, rate = lam)
    if (joint) { m <- fit$mdraws[idx[k], ]; m1 <- c(fit$S[[1]]$u, fit$S[[1]]$p_u, exp(m[1]), m[2]); m2 <- c(fit$S[[2]]$u, fit$S[[2]]$p_u, exp(m[3]), m[4]) }
    else { m1 <- c(gx$u, gx$p_u, gx$sigma, gx$xi); m2 <- c(gy$u, gy$p_u, gy$sigma, gy$xi) }
    f <- function(e, j) { mm <- if (j == 1) m1 else m2; Qsemi(1 - exp(-e), Y[, j], mm[1], mm[2], mm[3], mm[4]) }
    c(rl1 = f(r1 / 2, 1), rl2 = f(r1 / 2, 2), th1 = f(rho / 2, 1), th2 = f(rho / 2, 2)) }))
  xi <- if (joint) cbind(xi1 = fit$mdraws[idx, 2], xi2 = fit$mdraws[idx, 4]) else NULL
  list(cp = cp, chi = apply(chi, 2, q90), eta = q90((1 + sapply(th, `[[`, "chibar")) / 2), exc = q90(exc), orig = apply(orig, 2, q90), xi = if (!is.null(xi)) apply(xi, 2, q90) else NULL)
}
sc <- summ(fc, FALSE); sj <- summ(fj, TRUE)
pr <- function(lab, s) { cat("\n====", lab, "====\n"); cat("class probs (0, slab, 1):", round(s$cp, 3), "\n")
  cat("chi(0.9, 0.95, 0.99) 5/50/95%:\n"); print(round(s$chi, 3)); cat("eta:", round(s$eta, 3), "  exceedance fraction:", round(s$exc, 3), "\n")
  cat("diagonal one-season return level (PM2.5, O3) and threshold, 5/50/95%:\n"); print(round(s$orig)); if (!is.null(s$xi)) { cat("xi draws:\n"); print(round(s$xi, 3)) } }
pr("cut (stored fit)", sc); pr("joint posterior", sj)
cat("\nGPD MLE shapes:", round(c(gx$xi, gy$xi), 3), "  scales:", round(c(gx$sigma, gy$sigma), 1), "\n")
## GPD scales under the joint posterior (posterior medians over all draws), against the MLE above
cat("GPD scales, joint posterior median:", round(apply(exp(fj$mdraws[, c("lsig1", "lsig2")]), 2, median), 1), "\n")
saveRDS(list(cut = sc, joint = sj), "results/wanliu_compare.rds")
