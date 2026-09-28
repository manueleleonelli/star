## tab_s11.R -- LaTeX body of Table S11 (three treatments of the margins, GPD margins, n = 2000 and 400)
source("tab_all.R")
f <- function(x, d) { s <- formatC(x, format = "f", digits = d); sub("^-", "$-$", s) }
rows <- list(c("chi95_std", "$\\chi(0.95)$", 1, 3), c("chi99_std", "$\\chi(0.99)$", 1, 3), c("pstar_dir", "$10^4p^\\ast$", 1e4, 2),
             c("mstar_dir", "$10^3P(Y_1>y^\\ast)$", 1e3, 2), c("xi1", "$\\xi_1$", 1, 3))
out <- c()
for (n in c(2000, 400)) { A <- lapply(c(cut = "cut", cutx = "cutx", joint = "joint"), function(a) agg("gp", n, a))
  out <- c(out, sprintf("\\multicolumn{11}{l}{$n=%d$, %d exceedances per margin; $P(\\mathrm{AD}\\mid\\text{data})$ %s, %s and %s} \\\\", n, n / 10, f(A$cut$pad, 2), f(A$cutx$pad, 2), f(A$joint$pad, 2)))
  for (r in rows) { v <- r[1]; sc <- as.numeric(r[3]); d <- as.integer(r[4]); tr <- truth[sub("_.*", "", v)] * sc
    cells <- unlist(lapply(A, function(z) c(f(z$bias[v] * sc, d), f(z$cover[v], 2), f(z$width[v] * sc, d))))
    out <- c(out, paste0(r[2], " & ", f(tr, d), " & ", paste(cells, collapse = " & "), " \\\\")) }
  if (n == 2000) out <- c(out, "\\midrule") }
writeLines(out, "tab_s11_body.tex"); cat(out, sep = "\n")
## also the LN design at n = 400, and the numbers quoted in the text
for (d in c("gp", "ln")) for (n in c(400, 2000)) { A <- lapply(c(cut = "cut", cutx = "cutx", joint = "joint"), function(a) agg(d, n, a)); if (is.null(A$cut)) next
  cat(sprintf("\n%s %d: p* cover %.2f %.2f %.2f | p* rel bias %+.0f%% %+.0f%% %+.0f%% | mstar rel bias %+.0f%% %+.0f%% %+.0f%% | mstar skl %+.0f%% %+.0f%% | pstar skl %+.0f%% (cov %.2f) %+.0f%% (cov %.2f) | pmix cover %.2f %.2f %.2f | xi1 %+.3f %+.3f %+.3f xi2 %+.3f %+.3f %+.3f | P(AD) %.2f %.2f %.2f\n", d, n,
    A$cut$cover["pstar_dir"], A$cutx$cover["pstar_dir"], A$joint$cover["pstar_dir"], 100 * A$cut$bias["pstar_dir"] / truth["pstar"], 100 * A$cutx$bias["pstar_dir"] / truth["pstar"], 100 * A$joint$bias["pstar_dir"] / truth["pstar"],
    100 * A$cut$bias["mstar_dir"] / 0.005, 100 * A$cutx$bias["mstar_dir"] / 0.005, 100 * A$joint$bias["mstar_dir"] / 0.005, 100 * A$cut$bias["mstar_skl"] / 0.005, 100 * A$cutx$bias["mstar_skl"] / 0.005,
    100 * A$cut$bias["pstar_skl"] / truth["pstar"], A$cut$cover["pstar_skl"], 100 * A$cutx$bias["pstar_skl"] / truth["pstar"], A$cutx$cover["pstar_skl"],
    A$cut$cover["pmix_dir"], A$cutx$cover["pmix_dir"], A$joint$cover["pmix_dir"], A$cut$bias["xi1"], A$cutx$bias["xi1"], A$joint$bias["xi1"], A$cut$bias["xi2"], A$cutx$bias["xi2"], A$joint$bias["xi2"], A$cut$pad, A$cutx$pad, A$joint$pad)) }
