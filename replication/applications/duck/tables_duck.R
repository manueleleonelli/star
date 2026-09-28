## tables_duck.R -- Table S14 (tab:duck) and the numbers of Section S5.2 of the supplement from the stored fit
## results/wavesurge.rds written by app2_wavesurge.R; no download and no refitting. Run from this folder:
## source("tables_duck.R"). Writes the LaTeX rows of Table S14 to tab_duck.tex.
suppressMessages(library(star))
x <- readRDS("results/wavesurge.rds"); r <- x$report; E <- x$E; n <- nrow(E)

## ---- Section S5.2 ----
qs <- c(0.9, 0.95, 0.99)
chi_emp <- sapply(qs, function(q) mean(E[, 1] > -log(1 - q) & E[, 2] > -log(1 - q)) / (1 - q))
cat(sprintf("paired days n = %d, years %s to %s\n", n, format(min(x$data$date), "%Y"), format(max(x$data$date), "%Y")))
cat("empirical chi(0.9, 0.95, 0.99):", sprintf("%.2f", chi_emp), "\n")
cat("GPD shapes: wave", sprintf("%.2f", x$margins$gx$xi), " surge", sprintf("%.2f", x$margins$gy$xi), "\n")
cat("class probabilities (0, slab, 1):", sprintf("%.2f", r$cp), "\n")
cat("chibar posterior 5/50/95%:", sprintf("%.3f", r$chibar), "\n")
cat("Lambda(omega) at omega =", c(0.2, 0.35, 0.5, 0.65, 0.8), ":", sprintf("%.3f", r$adf$summary[1, ]), " boundary max(omega, 1 - omega):", c(0.8, 0.65, 0.5, 0.65, 0.8), "\n")
cat(sprintf("exceedance fraction %.2f, 90%% (%.2f, %.2f)\n", r$exc[2], r$exc[1], r$exc[3]))
cat("fitted chi(0.9, 0.95, 0.99), posterior mean:", sprintf("%.2f", r$chi[1, ]), "\n")
tb <- r$regions
cat(sprintf("empirical over fitted joint exceedance probability: %.2f at (0.9, 0.9), %.2f at (0.99, 0.99)\n", tb[1, "empirical"] / tb[1, "mean"], tb[5, "empirical"] / tb[5, "mean"]))
cat(sprintf("days per year at (0.9, 0.9): %.1f against %.1f; at (0.99, 0.99): %.1f against %.1f\n", tb[1, "days"], 365 * tb[1, "empirical"], tb[5, "days"], 365 * tb[5, "empirical"]))

## ---- Table S14: 365 days per year, 90% intervals from the stored report (300 posterior draws) ----
f2 <- function(v) sprintf("$%.2f$", v); f3 <- function(v) sprintf("%.3f", v)
rows <- sprintf("%s & %s & $%s$ & $(%s,\\,%s)$ & $%s$ & $%.1f$ & $%.1f$ \\\\", f2(tb[, "lev1"]), f2(tb[, "lev2"]), f3(tb[, "mean"]), f3(tb[, 4]), f3(tb[, 5]),
                f3(tb[, "empirical"]), tb[, "days"], 365 * tb[, "empirical"])
writeLines(rows, "tab_duck.tex"); cat("Table S14 rows (tab_duck.tex):\n"); writeLines(rows)
