## tables.R -- Table 5 of the main paper (from results/app3_summary.rds) and Table S12 of the supplement (from results/calib3.rds)
f3 <- function(x, d = 3) formatC(x, format = "f", digits = d)
## ---- calibration ----
if (file.exists("results/calib3.rds")) {
  cl <- readRDS("results/calib3.rds"); tab <- cl$tab
  rows <- c(eta12 = "$\\eta_{12}$", eta13 = "$\\eta_{13}$", eta23 = "$\\eta_{23}$", eta123 = "$\\eta_{123}$",
            chi12_99 = "$\\chi_{12}(0.99)$", chi13_99 = "$\\chi_{13}(0.99)$", chi123_9 = "$\\chi_{123}(0.9)$", chi123_99 = "$\\chi_{123}(0.99)$",
            exc = "exceedance fraction", rho_c = "$\\rho$ at the barycentre", u1 = "$u_1$", p = "$p$", chib12 = "$\\bar{\\chi}_{12}$", chib13 = "$\\bar{\\chi}_{13}$")
  out <- c("\\begin{table}[t]", "\\centering",
    "\\caption{Calibration check of the trivariate implementation: twenty replicates at $n=2000$ from the model with $\\bar{\\chi}=(0.6,0.2,0.2)$, $\\nu=0.5$, $u_j=2.2$, $p=2$ and $\\alpha=3$; truth, mean of the posterior medians, bias, root mean squared error, coverage and mean width of $90\\%$ intervals. Chains of $20{,}000$ iterations with $8{,}000$ of burn-in.}",
    "\\label{tab:calib3}", "\\footnotesize", "\\begin{tabular}{l cccccc}", "\\toprule",
    "quantity & truth & mean median & bias & RMSE & coverage & width \\\\", "\\midrule")
  for (k in names(rows)) { r <- tab[k, ]; d <- if (grepl("chi", k) && grepl("_", k)) 4 else 3
    out <- c(out, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", rows[k], f3(r$truth, d), f3(r$mean_median, d), f3(r$bias, d), f3(r$rmse, d), f3(r$coverage, 2), f3(r$width, d))) }
  out <- c(out, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(out, "tab_calib3.tex"); cat("minutes per fit:", round(mean(sapply(cl$res, `[[`, "mins")), 1), "\n"); print(round(tab, 4))
}
## ---- application ----
if (file.exists("results/app3_summary.rds")) {
  S <- readRDS("results/app3_summary.rds"); q3 <- function(v) quantile(v, c(.05, .5, .95))
  sub <- c("12" = "PM$_{2.5}$, PM$_{10}$", "13" = "PM$_{2.5}$, O$_3$", "23" = "PM$_{10}$, O$_3$", "123" = "all three")
  out <- c("\\begin{table}[t]", "\\centering",
    "\\caption{Trivariate fit at Wanliu. Left: residual tail dependence coefficients $\\eta_C$ of the three pairs and the triple, posterior median and $90\\%$ interval, with the empirical mean-excess estimate. Centre: $\\chi_C(q)$, posterior median with the empirical value in parentheses. Right: joint exceedances of marginal quantiles in days per warm season, posterior median and $90\\%$ interval against the observed rate over four seasons.}",
    "\\label{tab:app3}", "\\footnotesize", "\\setlength{\\tabcolsep}{3.5pt}", "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{l c c c ccc}", "\\toprule",
    "subset $C$ & $\\eta_C$ & $90\\%$ interval & empirical & $\\chi_C(0.9)$ & $\\chi_C(0.95)$ & $\\chi_C(0.99)$ \\\\", "\\midrule")
  for (nm in names(sub)) { e <- q3(S$ETA[, nm]); ch <- apply(S$CHI[[nm]], 2, median); em <- S$chi_emp[, nm]
    out <- c(out, sprintf("%s & %s & (%s, %s) & %s & %s (%s) & %s (%s) & %s (%s) \\\\", sub[nm], f3(e[2]), f3(e[1], 2), f3(e[3], 2), f3(S$eta_emp[nm]),
      f3(ch[1]), f3(em[1]), f3(ch[2]), f3(em[2]), f3(ch[3]), f3(em[3]))) }
  out <- c(out, "\\midrule", "\\multicolumn{7}{l}{joint exceedances, days per warm season} \\\\",
    "event & median & $90\\%$ interval & observed & & & \\\\", "\\midrule")
  EV <- S$EV
  evl <- c("all three above their $0.9$ quantiles", "PM$_{2.5}$, PM$_{10}$ above $0.95$, O$_3$ above $0.9$", "all three above $0.95$",
           "PM$_{2.5}$, PM$_{10}$ above $0.99$", "PM$_{2.5}$, O$_3$ above $0.99$", "all three above $0.99$")
  for (i in seq_len(nrow(EV))) out <- c(out, sprintf("%s & %s & (%s, %s) & %s & & & \\\\", evl[i], f3(EV[i, "med"], 2), f3(EV[i, "lo"], 2), f3(EV[i, "hi"], 2), f3(EV[i, "obs"], 2)))
  out <- c(out, "\\bottomrule", "\\end{tabular}}", "\\end{table}")
  writeLines(out, "tab_app3.tex")
  P <- S$P; cat("\nparameter medians and intervals:\n"); print(round(t(apply(P, 2, q3)), 3))
  cat("threshold marginal levels on axes:", round(apply(S$axes, 1, median), 3), "\n"); cat("exceedance fraction:", round(q3(S$EXC), 3), " days:", round(q3(S$EXC) * S$n), "\n")
  cat("GPD shapes:", round(sapply(S$g, `[[`, "xi"), 2), "\n")
}
