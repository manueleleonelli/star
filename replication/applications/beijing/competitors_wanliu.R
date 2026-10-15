## competitors_wanliu.R -- the alternatives of the simulation study applied to the Wanliu pair of Section 6.1:
## the conditional extremes model (HT) and the geometric model of Wadsworth and Campbell (WC), on the same
## standardized data as STAR, with the empirical estimator alongside. Outputs chi(q), eta (WC) and the joint
## exceedance probabilities of Table 4, with 90% moving-block bootstrap intervals (ten-day blocks within seasons,
## margins fixed as for STAR). Run from this folder after beijing_stations.R and app1_showcase.R; needs geometricMVE.
## About 15 minutes on two cores (B = 200 bootstrap samples).
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
suppressMessages(library(star))
source("../../simulations/competitors/ht.R"); source("../../simulations/competitors/wc.R")
x <- readRDS("results/fit_Wanliu.rds"); E <- x$E; n <- nrow(E); season <- format(x$data$date, "%Y")
S <- readRDS("results/beijing_showcase.rds")
B <- as.integer(Sys.getenv("STAR_NBOOT", "200")); CORES <- as.integer(Sys.getenv("STAR_CORES", "2"))
QS <- c(0.9, 0.95, 0.99, 0.999)
REG <- list(c(.9, .9), c(.95, .95), c(.99, .9), c(.9, .99), c(.99, .99))       # rows of Table 4

## HT on the standardized scale: Laplace transform of the exponential margins, conditioning at the 0.9 quantile,
## both conditioning directions averaged
ht_all <- function(E, nsim = 50000) {
  y1 <- lap_q(1 - exp(-E[, 1])); y2 <- lap_q(1 - exp(-E[, 2])); u <- lap_q(0.9)
  f12 <- ht_fit(y1, y2, u); f21 <- ht_fit(y2, y1, u)
  pj <- function(l1, l2) { v1 <- max(lap_q(l1), u); v2 <- max(lap_q(l2), u); 0.5 * (ht_joint(f12, v1, v2, nsim) + ht_joint(f21, v2, v1, nsim)) }
  c(setNames(sapply(QS, function(q) pj(q, q) / (1 - q)), paste0("chi", QS)),
    setNames(sapply(REG, function(r) pj(r[1], r[2])), sapply(REG, paste, collapse = "_")), a12 = f12$a, a21 = f21$a)
}
## WC: truncated-Gamma model above the 0.95 rolling quantile of R | W; a region that straddles the threshold is
## estimated by the empirical frequency below the threshold plus the model above it
wc_all <- function(E, nsim = 2e5) {
  r <- rowSums(E); w <- E[, 1] / r; r0 <- QR.2d(r = r, w = w, method = "empirical", tau = 0.95)$r0w
  below <- !(r > r0 & !is.na(r0))
  wf <- wc_fit(E); f <- wf$fit
  xs <- if (f$add) sim.2d(w = wf$wexc, r0w = wf$r0w, nsim = nsim, par = f$mle, add.gauge = TRUE,
                          gauge1 = gauge_rvad, gauge2 = gauge_gaussian, par1ind = 2, par2ind = 3)
        else sim.2d(w = wf$wexc, r0w = wf$r0w, nsim = nsim, par = f$mle, gfun = f$gfun)
  pj <- function(l1, l2) { t1 <- -log(1 - l1); t2 <- -log(1 - l2)
    mean(below & E[, 1] > t1 & E[, 2] > t2) + wf$pexc * mean(xs[, 1] > t1 & xs[, 2] > t2) }
  c(setNames(sapply(QS, function(q) pj(q, q) / (1 - q)), paste0("chi", QS)),
    setNames(sapply(REG, function(r) pj(r[1], r[2])), sapply(REG, paste, collapse = "_")), eta = wc_eta(wf), gauge = match(wf$which, c("rvad", "gauss", "invlog", "square", "mix")))
}
emp_all <- function(E) c(setNames(sapply(QS, function(q) mean(E[, 1] > -log(1 - q) & E[, 2] > -log(1 - q)) / (1 - q)), paste0("chi", QS)),
                         setNames(sapply(REG, function(r) mean(E[, 1] > -log(1 - r[1]) & E[, 2] > -log(1 - r[2]))), sapply(REG, paste, collapse = "_")))

set.seed(2026); est <- list(HT = ht_all(E), WC = wc_all(E), emp = emp_all(E))
## moving-block bootstrap, ten-day blocks within seasons, margins fixed
L <- 10; starts <- unlist(lapply(split(seq_len(n), season), function(ix) ix[1:(length(ix) - L + 1)]))
f <- "results/competitors_wanliu_boot.rds"
if (!file.exists(f)) {
  bt <- mclapply(1:B, function(b) { set.seed(7000 + b); idx <- head(as.vector(outer(0:(L - 1), sample(starts, ceiling(n / L), replace = TRUE), "+")), n)
    EE <- E[idx, ]; list(HT = tryCatch(ht_all(EE), error = function(e) NULL), WC = tryCatch(suppressWarnings(wc_all(EE)), error = function(e) NULL)) },
    mc.cores = CORES, mc.preschedule = FALSE)
  saveRDS(bt, f)
} else bt <- readRDS(f)
ci <- function(m) { M <- do.call(rbind, Filter(Negate(is.null), lapply(bt, `[[`, m))); apply(M, 2, quantile, c(0.05, 0.95), na.rm = TRUE) }
CI <- list(HT = ci("HT"), WC = ci("WC"))
saveRDS(list(est = est, CI = CI, B = B), "results/competitors_wanliu.rds")

## ---- report: chi(q), eta and Table 4 ----
fmt <- function(v, d = 3) formatC(v, format = "f", digits = d)
cat(sprintf("selected WC gauge: %s; HT dependence a (12, 21): %.2f, %.2f; 90%% interval of a (O3 given PM2.5 large): (%.2f, %.2f)\n",
            c("rvad", "gauss", "invlog", "square", "mix")[est$WC["gauge"]], est$HT["a12"], est$HT["a21"], CI$HT[1, "a12"], CI$HT[2, "a12"]))
cat("\nchi(q)        STAR mean (90%)           HT (90% block bootstrap)   WC (90% block bootstrap)   empirical\n")
star_chi <- sapply(QS, function(q) { j <- which.min(abs(S$qs - q)); c(S$chi_m[j], S$chi_s[1, j], S$chi_s[3, j]) })
for (k in seq_along(QS)) { nm <- paste0("chi", QS[k])
  cat(sprintf("q = %-6s %s (%s, %s)       %s (%s, %s)        %s (%s, %s)        %s\n", QS[k], fmt(star_chi[1, k]), fmt(star_chi[2, k]), fmt(star_chi[3, k]),
      fmt(est$HT[nm]), fmt(CI$HT[1, nm]), fmt(CI$HT[2, nm]), fmt(est$WC[nm]), fmt(CI$WC[1, nm]), fmt(CI$WC[2, nm]), fmt(est$emp[nm]))) }
cat(sprintf("\neta: WC %s (%s, %s); STAR posterior median and empirical as in Section 6.1\n", fmt(est$WC["eta"], 2), fmt(CI$WC[1, "eta"], 2), fmt(CI$WC[2, "eta"], 2)))
cat("\nTable 4 rows: STAR mean (90%) | HT (90%) | WC (90%) | empirical\n")
for (k in seq_along(REG)) { nm <- paste(REG[[k]], collapse = "_"); tb <- S$tab4[k, ]
  cat(sprintf("%-10s %.4f (%.4f, %.4f) | %.4f (%.4f, %.4f) | %.4f (%.4f, %.4f) | %.4f\n", nm, tb["mean"], tb["lo"], tb["hi"],
      est$HT[nm], CI$HT[1, nm], CI$HT[2, nm], est$WC[nm], CI$WC[1, nm], CI$WC[2, nm], est$emp[nm])) }

## ---- rows of Table 4: expected days per 183-day warm season, 90% intervals, and the observed rate over four seasons ----
lab <- c("0.90", "0.95", "0.99", "0.90", "0.99"); lab2 <- c("0.90", "0.95", "0.90", "0.99", "0.99")
rows <- sapply(seq_along(REG), function(k) { nm <- paste(REG[[k]], collapse = "_"); tb <- S$tab4[k, ]
  dg <- if (183 * tb["mean"] >= 1) 1 else 2; d <- function(v) formatC(183 * v, format = "f", digits = dg)
  sprintf("$%s$ & $%s$ & $%s$ & $(%s,\\,%s)$ & $%s$ & $(%s,\\,%s)$ & $%s$ & $(%s,\\,%s)$ & $%s$ \\\\", lab[k], lab2[k],
          d(tb["mean"]), d(tb["lo"]), d(tb["hi"]), d(est$HT[nm]), d(CI$HT[1, nm]), d(CI$HT[2, nm]),
          d(est$WC[nm]), d(CI$WC[1, nm]), d(CI$WC[2, nm]), d(est$emp[nm])) })
writeLines(rows, "tab_beijing.tex"); cat("\n", rows, sep = "\n")
