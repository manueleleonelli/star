## app3.R -- trivariate Beijing application: PM2.5, PM10 and O3 at Wanliu (Section 6.2, Table 5)
source(system.file("trivariate", "star3.R", package = "star")); K <- 4
## semiparametric margins as in the main paper: empirical below the 0.9 quantile, GPD above
gpd_fit <- function(x, u) { y <- x[x > u] - u
  nll <- function(p) { s <- exp(p[1]); xi <- p[2]; if (any(1 + xi * y / s <= 0)) return(1e10); length(y) * log(s) + (1 + 1/xi) * sum(log1p(xi * y / s)) }
  o <- optim(c(log(sd(y)), 0.1), nll); list(sigma = exp(o$par[1]), xi = o$par[2], u = u, p_u = mean(x > u)) }
to_exp <- function(x, g) { F <- ecdf(x)(x) * length(x) / (length(x) + 1); ab <- x > g$u
  F[ab] <- 1 - g$p_u * pmax(1 + g$xi * (x[ab] - g$u) / g$sigma, 1e-12)^(-1 / g$xi); -log(1 - F) }
from_exp <- function(e, x, g) { F <- 1 - exp(-e); out <- quantile(x, pmin(F, 1 - g$p_u), names = FALSE, type = 6)
  ab <- F > 1 - g$p_u; out[ab] <- g$u + g$sigma / g$xi * (((1 - F[ab]) / g$p_u)^(-g$xi) - 1); out }
d <- readRDS("results/daily6_wanliu.rds"); ok <- is.finite(d$PM2.5) & is.finite(d$PM10) & is.finite(d$O3); x <- d[ok, c("PM2.5", "PM10", "O3")]
n <- nrow(x); g <- lapply(x, function(v) gpd_fit(v, quantile(v, 0.9))); E <- sapply(1:3, function(j) to_exp(x[[j]], g[[j]])); colnames(E) <- names(x)
cat("n =", n, " GPD shapes:", round(sapply(g, `[[`, "xi"), 2), "\n")
ITER <- as.integer(Sys.getenv("ITER", "40000")); BURN <- as.integer(Sys.getenv("BURN", "15000"))
if (!file.exists("results/fit_app3.rds")) {
  t0 <- Sys.time(); fit <- mcmc3(E, n_iter = ITER, burnin = BURN, thin = 5, K = K, seed = 2026, verbose = TRUE)
  fit$minutes <- as.numeric(Sys.time() - t0, units = "mins"); saveRDS(fit, "results/fit_app3.rds") } else fit <- readRDS("results/fit_app3.rds")
cat("minutes:", round(fit$minutes, 1), " acceptance:", round(fit$acc, 2), "\n")
D <- fit$draws; R <- rowSums(E); W <- E / R
## ---- summaries -------------------------------------------------------------------------------
q3 <- function(v) c(lo = unname(quantile(v, .05)), med = unname(median(v)), hi = unname(quantile(v, .95)))
P <- t(apply(D, 1, function(z) { th <- unpack(z, K); c(alpha = th$alpha, u = th$u, p = th$p, chib = th$chib, nu = th$nu[1], rho_c = rho3(matrix(1/3, 1, 3), th$u, th$p)) }))
cat("\nparameters (5%, median, 95%):\n"); print(round(t(apply(P, 2, q3)), 3))
idx <- round(seq(1, nrow(D), length.out = 400)); TH <- lapply(idx, function(s) unpack(D[s, ], K))
subsets <- list("12" = 1:2, "13" = c(1, 3), "23" = 2:3, "123" = 1:3)
ETA <- sapply(subsets, function(C) sapply(TH, eta3, C = C)); cat("\neta_C:\n"); print(round(apply(ETA, 2, q3), 3))
## empirical (Hill-type) eta: mean excess of min_{j in C} E_j above its 0.9 quantile
eta_emp <- sapply(subsets, function(C) { T <- apply(E[, C, drop = FALSE], 1, min); u <- quantile(T, 0.9); mean(T[T > u] - u) })
cat("empirical eta:", round(eta_emp, 3), "\n")
QS <- c(0.9, 0.95, 0.99)
CHI <- lapply(subsets, function(C) sapply(QS, function(q) sapply(TH, chi3, K = K, q = q, C = C)))
chi_emp <- sapply(subsets, function(C) sapply(QS, function(q) { t <- -log(1 - q); mean(apply(E[, C, drop = FALSE] > t, 1, all)) / (1 - q) }))
cat("\nchi_C(q) model median (empirical):\n"); for (nm in names(subsets)) cat(nm, ":", paste(sprintf("%.3f (%.3f)", apply(CHI[[nm]], 2, median), chi_emp[, nm]), collapse = "  "), "\n")
EXC <- sapply(TH, exc_frac3, K = K); cat("\nexceedance fraction:", round(q3(EXC), 3), " observed tail fraction at posterior median threshold:",
  round(mean(R > rho3(W, apply(P[, 2:4], 2, median), median(P[, "p"]))), 3), "\n")
## joint events on the original scale (days per warm season of 183 days), with observed counts over four seasons
events <- list("all three > q0.9" = c(.9, .9, .9), "PM2.5, PM10 > q0.95, O3 > q0.9" = c(.95, .95, .9), "all three > q0.95" = c(.95, .95, .95),
               "PM2.5, PM10 > q0.99" = c(.99, .99, 0), "PM2.5, O3 > q0.99" = c(.99, 0, .99), "all three > q0.99" = c(.99, .99, .99))
EV <- t(sapply(events, function(lv) { t <- ifelse(lv > 0, -log(1 - lv), 0); p <- sapply(TH, joint_surv3, K = K, t = t)
  obs <- sum(apply(sweep(E, 2, t, ">") | matrix(rep(t == 0, each = n), n), 1, all)); c(q3(183 * p), obs = obs / 4, obs4 = obs) }))
cat("\njoint events, days per season (5%, median, 95%), observed per season:\n"); print(round(EV, 3))
## marginal quantile levels of the threshold on the axes
axes <- sapply(TH, function(th) 1 - exp(-th$u)); cat("threshold on the axes, marginal levels:", round(apply(axes, 1, median), 3), "\n")
saveRDS(list(E = E, x = x, g = g, P = P, ETA = ETA, eta_emp = eta_emp, CHI = CHI, chi_emp = chi_emp, EXC = EXC, EV = EV, axes = axes, TH = TH, n = n), "results/app3_summary.rds")
