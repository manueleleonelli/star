## beijing_stations.R -- twelve Beijing stations (Zhang et al. 2017): cross-station cleaning of hourly
## PM2.5 and O3, warm-season daily maxima, one fit per station (Section S5.1 and Table S13 of the supplement).
## Run from this folder after installing star: source("beijing_stations.R"). About 3 minutes per
## station on one core; resumable (results/fit_<station>.rds).
##
## Cleaning. An hourly value is set to missing when it exceeds 300 microg/m3 and
## three times the largest reading at the other eleven stations in the same hour (sensor faults).
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
dir <- "../../data/beijing"    # the twelve station files of Zhang et al. (2017), PRSA_Data_*.csv renamed by station and compressed with xz
fs <- list.files(dir, pattern = "\\.csv\\.xz$"); L <- lapply(file.path(dir, fs), read.csv); names(L) <- sub("\\.csv\\.xz$", "", fs)
key <- function(d) paste(d$year, d$month, d$day, d$hour); k0 <- key(L[[1]])
grab <- function(v) { M <- sapply(L, function(d) d[[v]][match(k0, key(d))]); colnames(M) <- names(L); M }
clean <- function(M, lev = 300, fac = 3) {
  for (j in seq_len(ncol(M))) { om <- suppressWarnings(apply(M[, -j, drop = FALSE], 1, max, na.rm = TRUE))
    bad <- !is.na(M[, j]) & M[, j] > lev & is.finite(om) & M[, j] > fac * om
    cat(sprintf("  %-14s flagged %d\n", colnames(M)[j], sum(bad))); M[bad, j] <- NA }
  M }
cat("O3\n"); O3 <- clean(grab("O3")); cat("PM2.5\n"); PM <- clean(grab("PM2.5"))
date <- as.Date(sprintf("%d-%02d-%02d", L[[1]]$year, L[[1]]$month, L[[1]]$day))
daily <- lapply(names(L), function(st) {
  a <- merge(aggregate(list(PM2.5 = PM[, st]), list(date = date), max, na.rm = TRUE), aggregate(list(O3 = O3[, st]), list(date = date), max, na.rm = TRUE), by = "date")
  a <- a[is.finite(a$PM2.5) & is.finite(a$O3) & as.integer(format(a$date, "%m")) %in% 4:9, ]; a })
names(daily) <- names(L)
for (st in names(daily)) { a <- daily[[st]]; cat(sprintf("%-14s n=%d  maxPM %4.0f maxO3 %4.0f  spearman %.2f  q99 PM %.0f O3 %.0f\n", st, nrow(a), max(a$PM2.5), max(a$O3), cor(a$PM2.5, a$O3, method = "spearman"), quantile(a$PM2.5, .99), quantile(a$O3, .99))) }
saveRDS(daily, "results/daily12.rds")

## ---- fits (same pipeline as app1_showcase.R) --------------------------------------------
suppressMessages(source("app_common.R")); library(parallel)
CORES <- as.integer(Sys.getenv("STAR_CORES", "2"))
ITER <- as.integer(Sys.getenv("STAR_ITER", "30000")); BURN <- round(ITER/3)
one <- function(st) {
  f <- sprintf("results/fit_%s.rds", st); if (file.exists(f)) return(invisible(NULL))
  a <- daily[[st]]; s <- standardize(a$PM2.5, a$O3); E <- s$E
  set.seed(2026); t0 <- Sys.time()
  fit <- star_mcmc(E, prior = star_prior("threecomp", p_min = 1), n_iter = ITER, burnin = BURN, thin = 5, verbose = FALSE)
  rep <- report(fit, E, st, days = 183)
  saveRDS(list(station = st, data = a, E = E, margins = s[c("gx", "gy")], fit = fit, report = rep, minutes = as.numeric(Sys.time() - t0, units = "mins")), f)
  cat("done", st, round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
}
order_st <- c("Dongsi", "Nongzhanguan", "Tiantan", "Guanyuan", "Wanshouxigong", "Aotizhongxin", "Wanliu", "Gucheng", "Changping", "Shunyi", "Huairou", "Dingling")
todo <- order_st[!file.exists(sprintf("results/fit_%s.rds", order_st))]
invisible(mclapply(todo, function(st) tryCatch(one(st), error = function(e) cat("ERR", st, conditionMessage(e), "\n")), mc.cores = CORES, mc.preschedule = FALSE))

## ---- summary across stations (Table S13 and the text of Section S5.1) -----------------------
## 400 equally spaced posterior draws per station, as in app1_showcase.R; predictive means are exact
## the station fits only (results/ also holds fit_wanliu_joint.rds from wanliu_joint.R)
fs <- file.path("results", sprintf("fit_%s.rds", names(daily))); gl <- seq(0.5/201, 1-0.5/201, length.out = 201); glw <- rep(1/201, 201)
wg <- seq(0.005, 0.995, length.out = 199)
rows <- lapply(fs, function(f) { x <- readRDS(f); r <- x$report; E <- x$E; n <- nrow(E); fit <- x$fit; D <- fit$draws
  R <- rowSums(E); W <- E[, 1] / R
  idx <- round(seq(1, nrow(D), length.out = 400)); th <- lapply(idx, function(s) star:::theta_row(fit, s))
  sl <- D[, "chibar"][fit$component == "slab"]
  ## predictive check of joint exceedance counts at the level-q quantiles
  pc <- function(l) { t <- -log(1 - l); p <- sapply(th, function(tt) star:::joint_surv_cpp(t, t, gl, glw, tt$phi, tt$alpha, tt$b, tt$u1, tt$u2, tt$p, tt$chibar, tt$z))
    set.seed(3); N <- rbinom(length(p), n, p); c(obs = sum(E[,1] > t & E[,2] > t), mean = n * mean(p), lo = quantile(N, .05, names = FALSE), hi = quantile(N, .95, names = FALSE)) }
  a <- pc(0.9); b <- pc(0.95)
  ## days outside the one-season return level set (posterior median radius), near the axes and in between
  rho <- t(sapply(th, function(t) rho_star(wg, t$u1, t$u2, t$p))); sig <- t(sapply(th, function(t) sigma_star(wg, t$chibar, t$z)))
  lam <- t(sapply(th, function(t) lambda_star(wg, t$b))); piw <- t(sapply(seq_along(th), function(i) pgamma(rho[i, ], shape = th[[i]]$alpha, rate = lam[i, ])))
  fW <- t(sapply(th, function(t) { K <- length(t$phi); rowSums(sapply(1:K, function(k) t$phi[k] * dbeta(wg, k, K - k + 1))) }))
  rl <- t(sapply(seq_along(th), function(i) { rb <- qgamma(1 - 1/183, shape = th[[i]]$alpha, rate = lam[i, ]); tl <- (1 - piw[i, ]) > 1/183
    rb[tl] <- rho[i, tl] + sig[i, tl] * log((1 - piw[i, tl]) * 183); rb }))
  rlm <- apply(rl, 2, median); out <- R > rlm[sapply(W, function(v) which.min(abs(wg - v)))]; ax <- wg < 0.25 | wg > 0.75
  pe <- rowMeans(sapply(seq_along(th), function(i) { sv <- ifelse(rlm > rho[i, ], (1 - piw[i, ]) * exp(-(rlm - rho[i, ]) / sig[i, ]), 1 - pgamma(rlm, shape = th[[i]]$alpha, rate = lam[i, ]))
    f <- fW[i, ] / sum(fW[i, ]); c(axis = sum((f * sv)[ax]), mid = sum((f * sv)[!ax])) }))
  ## empirical estimate of eta: mean excess of min(E1, E2) over its 0.9 quantile
  Tm <- pmin(E[, 1], E[, 2]); hill <- mean(Tm[Tm > quantile(Tm, 0.9)] - quantile(Tm, 0.9))
  data.frame(station = x$station, n = n, rho_s = round(cor(x$data$PM2.5, x$data$O3, method = "spearman"), 2), xi_pm = round(x$margins$gx$xi, 2), xi_o3 = round(x$margins$gy$xi, 2),
             P0 = round(r$cp[1], 2), Pslab = round(r$cp[2], 2), P1 = round(r$cp[3], 2), chibar_slab = round(median(sl), 2), eta = round(median((1 + D[, "chibar"]) / 2), 2), eta_emp = round(hill, 2),
             exc = round(r$exc[2], 2), p = round(r$p[2], 1), chi90 = round(r$chi[1, 1], 3), chi90_emp = round(r$chi_emp[1], 3), chi95 = round(r$chi[1, 2], 3), chi95_emp = round(r$chi_emp[2], 3), chi99 = round(r$chi[1, 3], 3),
             obs90 = a[["obs"]], mean90 = a[["mean"]], pred90 = sprintf("%.1f (%.0f,%.0f)", a["mean"], a["lo"], a["hi"]), obs95 = b[["obs"]], mean95 = b[["mean"]], pred95 = sprintf("%.1f (%.0f,%.0f)", b["mean"], b["lo"], b["hi"]),
             out_axis = sum(out & (W < 0.25 | W > 0.75)), exp_axis = n * pe[["axis"]], out_mid = sum(out & W >= 0.25 & W <= 0.75), exp_mid = n * pe[["mid"]], minutes = round(x$minutes, 1)) })
tab <- do.call(rbind, rows); rownames(tab) <- NULL
print(tab[, c("station","n","rho_s","P0","Pslab","P1","chibar_slab","exc","chi90","chi90_emp","chi99","obs95","pred95")], row.names = FALSE)
print(tab[, c("station","eta","eta_emp","obs90","pred90","out_axis","exp_axis","out_mid","exp_mid")], row.names = FALSE, digits = 3)
cat(sprintf("totals: joint exceedances at 0.95 %d against %.1f, at 0.9 %d against %.1f; outside the one-season sets near the axes %d against %.1f, in between %d against %.1f\n",
            sum(tab$obs95), sum(tab$mean95), sum(tab$obs90), sum(tab$mean90), sum(tab$out_axis), sum(tab$exp_axis), sum(tab$out_mid), sum(tab$exp_mid)))
saveRDS(tab, "results/summary12.rds")
## LaTeX rows of Table S13 (tab:stations), written to tab_stations.tex
f2 <- function(v) sprintf("$%.2f$", v); f3 <- function(v) sprintf("$%.3f$", v); fi <- function(v) sprintf("$%d$", as.integer(v))
pr95 <- sub("^(\\S+) \\((\\d+),(\\d+)\\)$", "$\\1$ $(\\2,\\\\,\\3)$", tab$pred95)
rows_tex <- paste(tab$station, fi(tab$n), f2(tab$P0), f2(tab$Pslab), f2(tab$P1), f2(tab$eta), f2(tab$eta_emp), f2(tab$exc), f3(tab$chi90), f3(tab$chi90_emp), f3(tab$chi99),
                  fi(tab$obs95), pr95, fi(tab$out_axis), sprintf("$%.1f$", tab$exp_axis), sep = " & ")
writeLines(paste(rows_tex, "\\\\"), "tab_stations.tex"); cat("wrote tab_stations.tex\n")
