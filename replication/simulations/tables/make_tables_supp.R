## make_tables_supp.R -- Supplement Tables S2 (tab:diag), S3 (tab:ident) and
## S10 (tab:sens), the TikZ coordinates of Figure S1 (fig:ident), and the
## numbers quoted in main Section 4.2 and supplement Sections S3.6, S4.1 and
## S4.7. Works from saved results only; no MCMC is rerun.
##
## Run from simulations/tables. Inputs:
##   ../main/results/s10_diagnostics.rds, ../main/results/s7_slab.rds
##   ../identification/contraction.rds, pfocus.rds, ridge_anchor.rds,
##     ridge_weak.rds, results05/undercoverage.rds
##   ../sensitivity/results/kappa.rds, prior.rds, basis.rds, alphaw.rds,
##     ridge300.rds
## Outputs: tab_diag.tex, tab_ident.tex, tab_sens.tex, fig_ident_coords.tex
## (the contents of each table float, and the plot lines of the figure).
## Optional check against the paper: set STAR_PAPER_DIR to the folder holding
## supplement-ba.tex; every table row and figure line is then compared with it.

ID <- "../identification"; SE <- "../sensitivity/results"; MA <- "../main/results"

## ---- formatting helpers --------------------------------------------------
num <- function(x, d) {                        # fixed decimals, no negative zero
  s <- formatC(x, format = "f", digits = d)
  sub("^-(0\\.?0*)$", "\\1", s)
}
f <- function(x, d) paste0("$", num(x, d), "$")
row <- function(label, cells) paste0(label, " & ", paste(cells, collapse = " & "), " \\\\")
tabular <- function(spec, head, body) c(sprintf(r"(\begin{tabular}{%s})", spec), head,
                                        r"(\midrule)", body, r"(\bottomrule)", r"(\end{tabular})")
emit <- function(lines, file) { writeLines(lines, file); cat(lines, sep = "\n"); cat("\n") }

## =========================================================================
## Table S2 (tab:diag): convergence diagnostics, four chains per dataset,
## twelve datasets per configuration. Input: s10_diagnostics.rds
## =========================================================================
dg <- readRDS(file.path(MA, "s10_diagnostics.rds"))
cfg_lab <- c(mid_n1000     = r"(STAR, $\bar{\chi}=0.6$, $n=1000$)",
             mid_n4000     = r"(STAR, $\bar{\chi}=0.6$, $n=4000$)",
             gauss09_n4000 = r"(Gaussian, $\varrho=0.9$, $n=4000$)")
body_diag <- unname(sapply(names(cfg_lab), function(cf) {
  x <- dg[dg$cfg == cf, ]
  row(cfg_lab[cf], c(f(mean(x$max_rhat), 3), f(max(x$max_rhat), 3), f(mean(x$min_ess), 0),
                     f(min(x$min_ess), 0), f(mean(x$ess_p), 0), f(mean(x$ess_chibar), 0)))
}))
tab_diag <- tabular("l cc cc c c", c(r"(\toprule)",
  r"( & \multicolumn{2}{c}{largest $\hat{R}$} & \multicolumn{2}{c}{smallest ESS} & ESS of $p$ & ESS of $\bar{\chi}$ \\)",
  r"(\cmidrule(lr){2-3}\cmidrule(lr){4-5})",
  r"(Configuration & mean & max & mean & min & mean & mean \\)"), body_diag)
cat("\n===== Table S2 (tab:diag) =====\n"); emit(tab_diag, "tab_diag.tex")

## =========================================================================
## Table S3 (tab:ident). (a) contraction slopes, input contraction.rds;
## (b) coverage of rho(1/2) at n = 8000, input results05/undercoverage.rds;
## (c) anchoring prior, n = 2000 from ridge_{anchor,weak}.rds (03_ridge.R)
## and n = 300 from ridge300.rds (run_sensitivity.R)
## =========================================================================
ct <- readRDS(file.path(ID, "contraction.rds"))
wmed <- aggregate(width ~ quantity + n, ct, median)
slope <- sapply(split(wmed, wmed$quantity), function(a) unname(coef(lm(log(width) ~ log(n), a))[2]))
qa <- c("alpha", "u1", "u2", "exc_frac", "rho0.15", "rho0.85", "rho0.5", "chibar", "nu_half", "p")
body_a <- paste(paste(f(slope[qa], 2), collapse = " & "), r"(\\)")

uc <- readRDS(file.path(ID, "results05/undercoverage.rds"))
um <- aggregate(cbind(rho_cover, ess_logp) ~ iter + v_p, uc, mean)
body_b <- sapply(c(1.5, 0.75), function(v) {
  s <- um[um$v_p == v & um$iter == 8000, ]; l <- um[um$v_p == v & um$iter == 30000, ]
  row(f(v, if (v == 1.5) 1 else 2), c(f(s$rho_cover, 2), f(s$ess_logp, 0), f(l$rho_cover, 2), f(l$ess_logp, 0)))
})

## n = 2000: "corr." is the largest absolute posterior correlation, always the
## shape-rate pair (log_alpha, b0) of the Gamma body; "across" is the largest
## absolute correlation over all remaining pairs. Restricting instead to pairs
## in different sampler blocks (ang, body, thr, tail) gives 0.34 and 0.28:
## the weak-prior 0.30 is the log_u1-log_p pair inside the threshold block.
blk <- c(log_u1 = "thr", log_u2 = "thr", log_p = "thr", logit_chibar = "tail",
         log_alpha = "body", b0 = "body", z0 = "tail")
rid <- lapply(c(anchor = "anchor", weak = "weak"), function(nm) {
  x <- readRDS(file.path(ID, sprintf("ridge_%s.rds", nm))); C <- abs(x$Cor); diag(C) <- 0
  top <- which(C == max(C), arr.ind = TRUE)[1, ]
  stopifnot(setequal(rownames(C)[top], c("log_alpha", "b0")))
  C2 <- C; C2["log_alpha", "b0"] <- C2["b0", "log_alpha"] <- 0
  B <- outer(blk[rownames(C)], blk[colnames(C)], "!=")
  c(corr = max(C), across = max(C2), across_blocks = max(C * B), sd = unname(x$sd["log_u1"]))
})
r300 <- aggregate(cbind(max_corr, sd_log_u1, exc) ~ prior, readRDS(file.path(SE, "ridge300.rds")), mean)
body_c <- sapply(c(anchor = "anchoring", weak = "weak"), function(pl) {
  nm <- if (pl == "anchoring") "anchor" else "weak"; a <- rid[[nm]]; b <- r300[r300$prior == nm, ]
  row(pl, c(f(a["corr"], 2), f(a["across"], 2), f(a["sd"], 2), f(b$max_corr, 2), f(b$sd_log_u1, 2), f(b$exc, 2)))
})

tab_ident <- c(
  r"(\begin{tabular}{cccccccccc})", r"(\multicolumn{10}{l}{(a) contraction slopes}\\)", r"(\toprule)",
  r"($\alpha$ & $u_1$ & $u_2$ & exc. & $\rho(0.15)$ & $\rho(0.85)$ & $\rho(1/2)$ & $\bar{\chi}$ & $\nu(1/2)$ & $p$ \\)",
  r"(\midrule)", body_a, r"(\bottomrule)", r"(\end{tabular})", "", r"(\vspace{0.8em})",
  r"(\begin{tabular}{c cc cc})", r"(\multicolumn{5}{l}{(b) coverage of $\rho(1/2)$ at $n=8000$}\\)", r"(\toprule)",
  r"( & \multicolumn{2}{c}{$8000$ iterations} & \multicolumn{2}{c}{$30{,}000$ iterations} \\)",
  r"(\cmidrule(lr){2-3}\cmidrule(lr){4-5})", r"($v_p$ & coverage & ESS & coverage & ESS \\)",
  r"(\midrule)", body_b, r"(\bottomrule)", r"(\end{tabular}\hspace{1.5em})",
  r"(\begin{tabular}{l ccc ccc})", r"(\multicolumn{7}{l}{(c) anchoring prior}\\)", r"(\toprule)",
  r"( & \multicolumn{3}{c}{$n=2000$} & \multicolumn{3}{c}{$n=300$} \\)",
  r"(\cmidrule(lr){2-4}\cmidrule(lr){5-7})", r"(prior & corr. & across & sd & corr. & sd & exc. \\)",
  r"(\midrule)", body_c, r"(\bottomrule)", r"(\end{tabular})")
cat("\n===== Table S3 (tab:ident) =====\n"); emit(tab_ident, "tab_ident.tex")
cat(sprintf("n = 2000, largest correlation between different sampler blocks: anchoring %.3f, weak %.3f\n",
            rid$anchor["across_blocks"], rid$weak["across_blocks"]))

## =========================================================================
## Table S10 (tab:sens). (a) paired shifts from the default fit, inputs
## kappa.rds and basis.rds, and the alpha(w) error, input alphaw.rds;
## (b) atom weights and slab, input prior.rds; (c) concentrated slab, input
## s7_slab.rds. In (a) each cell is the mean within-dataset shift with the
## largest absolute shift over the five datasets in brackets; chi(0.99) is a
## relative shift in percent.
## =========================================================================
shifts <- function(D, key, ref, cells) do.call(rbind, lapply(cells, function(cc) {
  m <- merge(D[D[[key]] == cc, ], D[D[[key]] == ref, ], by = "seed")
  dd <- sapply(c("exc", "rho_half", "chibar"), function(v) m[[paste0(v, ".x")]] - m[[paste0(v, ".y")]])
  dd <- cbind(dd, chi99 = 100 * (m$chi99.x / m$chi99.y - 1))
  data.frame(setting = cc, t(c(mean = colMeans(dd), max = apply(abs(dd), 2, max))))
}))
kap <- readRDS(file.path(SE, "kappa.rds")); bas <- readRDS(file.path(SE, "basis.rds"))
SH <- rbind(shifts(kap, "cell", "k08", c("k07", "k09", "s05", "s2")),
            shifts(bas, "basis", "default", c("small", "large")))
sh_lab <- c(k07 = r"($\kappa_0=0.7$)", k09 = r"($\kappa_0=0.9$)", s05 = r"($s_u=0.5$)", s2 = r"($s_u=2$)",
            small = r"($(K_W,K_\lambda,K_\nu)=(3,1,1)$)", large = r"($(K_W,K_\lambda,K_\nu)=(8,4,4)$)")
dig <- c(exc = 3, rho_half = 2, chibar = 3, chi99 = 1)
cell_sh <- function(r, v) paste0(f(r[[paste0("mean.", v)]], dig[v]), " $(", num(r[[paste0("max.", v)]], dig[v]), ")$")
body_sa <- sapply(seq_len(nrow(SH)), function(i) row(sh_lab[SH$setting[i]], sapply(names(dig), function(v) cell_sh(SH[i, ], v))))
alw <- readRDS(file.path(SE, "alphaw.rds")); af <- alw$fits; tr <- alw$truth
err_aw <- c(exc = mean(af$exc) - tr[["exc"]], rho_half = mean(af$rho_half) - tr[["rho_half"]],
            chibar = mean(af$chibar) - 0.6,                      # true chibar of theta0(0.6)
            chi99 = mean(100 * (af$chi99 / tr[["chi99"]] - 1)))
body_sa <- c(body_sa, row(r"(body shape $\alpha(w)$, error)", mapply(f, err_aw, dig)))

pr <- aggregate(P1 ~ prior + chibar0, readRDS(file.path(SE, "prior.rds")), mean)
pr_lab <- c(spike10 = r"($\pi_0=\pi_1=0.1$)", default = r"($\pi_0=\pi_1=0.25$)",
            spike40 = r"($\pi_0=\pi_1=0.4$)", slab22 = r"($\mathrm{Beta}(2,2)$ slab)")
body_sb <- sapply(names(pr_lab), function(p) row(pr_lab[p], f(pr$P1[pr$prior == p][order(pr$chibar0[pr$prior == p])], 2)))

s7 <- aggregate(P1 ~ gen + n + slab, readRDS(file.path(MA, "s7_slab.rds")), mean)
g_lab <- c(gauss09 = r"(Gaussian, $\varrho=0.9$)", m_chi095 = r"(STAR, $\bar{\chi}=0.95$)", m_chi1 = r"(STAR, $\bar{\chi}=1$)")
P7 <- function(g, n, s) s7$P1[s7$gen == g & s7$n == n & s7$slab == s]
body_sc <- sapply(names(g_lab), function(g) row(g_lab[g], f(c(P7(g, 2000, "uniform"), P7(g, 2000, "concentrated"),
                                                              P7(g, 8000, "uniform"), P7(g, 8000, "concentrated")), 2)))

tab_sens <- c(
  r"(\begin{tabular}{l cccc})", r"(\multicolumn{5}{l}{(a) shifts in the summaries}\\)", r"(\toprule)",
  r"(Setting & exc. & $\rho(1/2)$ & $\bar{\chi}$ & $\chi(0.99)$, \% \\)",
  r"(\midrule)", body_sa, r"(\bottomrule)", r"(\end{tabular})", "", r"(\vspace{0.8em})",
  r"(\begin{tabular}{l cc})", r"(\multicolumn{3}{l}{(b) atom weights and slab}\\)", r"(\toprule)",
  r"(Prior & $\bar{\chi}=0.5$ & $\bar{\chi}=0.95$ \\)",
  r"(\midrule)", body_sb, r"(\bottomrule)", r"(\end{tabular}\hspace{1.5em})",
  r"(\begin{tabular}{l cc cc})", r"(\multicolumn{5}{l}{(c) concentrated slab}\\)", r"(\toprule)",
  r"( & \multicolumn{2}{c}{$n=2000$} & \multicolumn{2}{c}{$n=8000$} \\)",
  r"(\cmidrule(lr){2-3}\cmidrule(lr){4-5})",
  r"(Generator & uniform & $\mathrm{Beta}(2,0.8)$ & uniform & $\mathrm{Beta}(2,0.8)$ \\)",
  r"(\midrule)", body_sc, r"(\bottomrule)", r"(\end{tabular})")
cat("\n===== Table S10 (tab:sens) =====\n"); emit(tab_sens, "tab_sens.tex")

## =========================================================================
## Figure S1 (fig:ident), TikZ coordinates (log10 n, log10 width).
## Left: median 90% interval width over 30 replicates, input contraction.rds.
## Right: median pointwise width of rho(.) (solid) and median(p_hi) minus
## median(p_lo) (dashed), input pfocus.rds. The published coordinates were
## entered from the console summaries of 02_contraction.R (widths printed at
## three significant digits) and 04_pfocus.R (print(..., digits = 3), which
## gives each column a common number of decimals), so the widths are rounded
## as displayed there before taking log10.
## =========================================================================
coords <- function(x, y) paste(sprintf("(%s,%s)", num(x, 3), num(y, 3)), collapse = " ")
curve <- function(style, x, y) sprintf(r"(\draw[%s] plot coordinates {%s};)", style, coords(x, y))
label <- function(col, x, y, txt) sprintf(r"(\node[%s,anchor=west] at (%s,%s) {\tiny %s};)", col, num(x, 2), num(y, 3), txt)

ns_l <- sort(unique(wmed$n))
left <- list(p = c("blue", "$p$"), chibar = c("violet", r"($\bar{\chi}$)"),
             rho0.5 = c("red!80!black", r"($\rho(1/2)$)"), rho0.15 = c("orange", r"($\rho(0.15)$)"),
             exc_frac = c("black!60", r"($\text{exc.\ frac.}$)"))
yl <- lapply(names(left), function(q) { a <- wmed[wmed$quantity == q, ]; log10(signif(a$width[order(a$n)], 3)) })
fig_left <- c(mapply(function(q, y) curve(paste0("thick,", left[[q]][1]), log10(ns_l), y), names(left), yl),
              mapply(function(q, y) label(left[[q]][1], 3.95, tail(y, 1), left[[q]][2]), names(left), yl))

pf <- aggregate(cbind(p_lo, p_hi, rho_width_mean) ~ ptrue + n, readRDS(file.path(ID, "pfocus.rds")), median)
pf$p_CI_width <- pf$p_hi - pf$p_lo
shown3 <- function(x) as.numeric(format(x, digits = 3))     # as print(df, digits = 3) displays a column
pf$rho_disp <- shown3(pf$rho_width_mean); pf$pw_disp <- shown3(pf$p_CI_width)
ns_r <- sort(unique(pf$n)); right <- c(`1` = "orange", `2` = "red!80!black", `8` = "blue")
fig_right <- unlist(lapply(names(right), function(pt) {
  a <- pf[pf$ptrue == as.numeric(pt), ]; a <- a[order(a$n), ]
  c(curve(paste0("thick,", right[pt]), log10(ns_r), log10(a$rho_disp)),
    curve(paste0("thick,dashed,", right[pt]), log10(ns_r), log10(a$pw_disp)))
}))
fig_right <- c(fig_right, sapply(names(right), function(pt) {
  a <- pf[pf$ptrue == as.numeric(pt) & pf$n == max(ns_r), ]
  label(right[pt], 3.66, log10(a$pw_disp), sprintf("$p=%s$", pt)) }))
fig_ident <- c("% (a) contraction of 90% intervals", unname(fig_left),
               "% (b) curve (solid) against exponent (dashed)", unname(fig_right))
cat("\n===== Figure S1 (fig:ident), coordinates =====\n"); emit(fig_ident, "fig_ident_coords.tex")

## sensitivity to the display rounding: the same coordinates from unrounded
## widths differ by at most about 0.001; log-log slopes of the solid curves
raw_l <- unlist(lapply(names(left), function(q) { a <- wmed[wmed$quantity == q, ]; log10(a$width[order(a$n)]) }))
pr_o <- pf[order(pf$ptrue, pf$n), ]
cat(sprintf("max |unrounded - displayed| log10 width: left %.4f, right solid %.4f, right dashed %.4f\n",
            max(abs(raw_l - unlist(yl))), max(abs(log10(pr_o$rho_width_mean) - log10(pr_o$rho_disp))),
            max(abs(log10(pr_o$p_CI_width) - log10(pr_o$pw_disp)))))
sl_r <- sapply(split(pr_o, pr_o$ptrue), function(a) diff(log10(a$rho_width_mean)) / diff(log10(a$n)))
cat("log-log slopes of the solid curves (p = 1, 2, 8):", num(sl_r, 2), "\n")

## =========================================================================
## Numbers quoted in the text: main Section 4.2, supplement Sections S3.6,
## S4.1 and S4.7. Inputs as above. 'paper' is the value as printed.
## =========================================================================
claims <- list()
cl <- function(where, what, paper, value, ok) claims[[length(claims) + 1]] <<-
  data.frame(where = where, quantity = what, paper = paper, reproduced = value, ok = ok)
K <- SH[SH$setting %in% c("k07", "k09", "s05", "s2"), ]
wd <- mean(kap$chibar_w[kap$cell == "k08"])
for (w in c("main 4.2", "S4.7")) {
  S <- if (w == "main 4.2") K else SH                     # S4.7 adds the basis sizes
  cl(w, "max |mean shift|, exceedance fraction", "at most 0.025", num(max(abs(S$mean.exc)), 4), num(max(abs(S$mean.exc)), 3) == "0.025")
  if (w == "main 4.2")
    cl(w, "max |mean shift|, rho(1/2)", "at most 0.15", num(max(abs(S$mean.rho_half)), 3), num(max(abs(S$mean.rho_half)), 2) == "0.15")
  cl(w, "max |mean shift|, chibar", "at most 0.02", num(max(abs(S$mean.chibar)), 4), num(max(abs(S$mean.chibar)), 2) == "0.02")
  cl(w, "chibar 90% interval width (default fits; all fits)", "near 0.61",
     paste(num(wd, 3), num(mean(kap$chibar_w), 3), sep = "; "), num(wd, 2) == "0.61")
  cl(w, "max |mean shift|, chi(0.99), %", "under 4", num(max(abs(S$mean.chi99)), 2), max(abs(S$mean.chi99)) < 4)
}
dslab <- sapply(names(g_lab), function(g) sapply(c(2000, 8000), function(n) c(P7(g, n, "uniform") - P7(g, n, "concentrated"),
  round(P7(g, n, "uniform"), 2) - round(P7(g, n, "concentrated"), 2))))
cl("S4.7", "drop in P(AD) under Beta(2,0.8) slab (unrounded; from table entries)", "0.10 to 0.18",
   sprintf("%s to %s; %s to %s", num(min(dslab[c(1, 3), ]), 3), num(max(dslab[c(1, 3), ]), 3),
           num(min(dslab[c(2, 4), ]), 2), num(max(dslab[c(2, 4), ]), 2)),
   num(min(dslab[c(2, 4), ]), 2) == "0.10" && num(max(dslab[c(2, 4), ]), 2) == "0.18")
g8 <- c(P7("gauss09", 8000, "uniform"), P7("gauss09", 8000, "concentrated"))
cl("S4.7", "Gaussian 0.9, P(AD) at n = 8000, both slabs", "above one half", paste(num(g8, 2), collapse = ", "), all(g8 > 0.5))
n7 <- table(readRDS(file.path(MA, "s7_slab.rds"))[, c("gen", "n", "slab")])
nds <- c(length(unique(kap$seed)), length(unique(bas$seed)), length(unique(readRDS(file.path(SE, "prior.rds"))$seed)))
cl("S4.7", "replicates per cell in (c); datasets in (a), (b)", "fifty; five",
   sprintf("%s; %s", paste(unique(as.vector(n7)), collapse = ","), paste(unique(nds), collapse = ",")),
   all(n7 == 50) && all(nds == 5))
reg <- slope[c("alpha", "u1", "u2", "exc_frac", "rho0.15", "rho0.85")]
cl("S4.1", "slopes of alpha, anchors, exc. frac., flanks", "-0.50 to -0.67",
   sprintf("%s to %s", num(max(reg), 2), num(min(reg), 2)), num(max(reg), 2) == "-0.50" && num(min(reg), 2) == "-0.67")
cl("S4.1", "slower: rho(1/2), chibar, nu(1/2)", "slower than -1/2",
   paste(num(slope[c("rho0.5", "chibar", "nu_half")], 2), collapse = ", "), all(slope[c("rho0.5", "chibar", "nu_half")] > max(reg)))
cl("S4.1", "coverage of rho(1/2), v_p = 1.5 (8000; 30000 iterations)", "near 0.8",
   paste(num(um$rho_cover[um$v_p == 1.5][order(um$iter[um$v_p == 1.5])], 2), collapse = ", "),
   all(abs(um$rho_cover[um$v_p == 1.5] - 0.8) < 0.05))
cl("S4.1", "coverage of rho(1/2), v_p = 0.75, 30000 iterations", "0.90",
   num(um$rho_cover[um$v_p == 0.75 & um$iter == 30000], 2), num(um$rho_cover[um$v_p == 0.75 & um$iter == 30000], 2) == "0.90")
cl("S4.1", "largest correlation across blocks, anchoring; weak (n = 2000)", "below 0.35",
   sprintf("%s; %s", num(rid$anchor["across"], 3), num(rid$weak["across"], 3)), max(rid$anchor["across"], rid$weak["across"]) < 0.35)
cl("S4.1", "exc. frac. at n = 300, anchoring; weak", "threshold kept from drifting up",
   sprintf("%s; %s", num(r300$exc[r300$prior == "anchor"], 3), num(r300$exc[r300$prior == "weak"], 3)),
   r300$exc[r300$prior == "anchor"] > r300$exc[r300$prior == "weak"])
cl("S4.1 / Fig. S1", "replicates per n in contraction study", "30", paste(unique(table(ct$n) / length(unique(ct$quantity))), collapse = ","),
   all(table(ct$n) / length(unique(ct$quantity)) == 30))
cl("S3.6", "datasets per configuration in Table S2", "twelve", paste(table(dg$cfg)[names(cfg_lab)], collapse = ","), all(table(dg$cfg) == 12))
for (cf in c("mid_n4000", "gauss09_n4000")) {
  x <- dg[dg$cfg == cf, ]
  cl("S3.6", sprintf("%s: datasets where p has the smallest ESS; the largest R-hat", cf), "p is the binding parameter",
     sprintf("%d/%d; %d/%d", sum(abs(x$ess_p - x$min_ess) < 1e-8), nrow(x), sum(abs(x$rhat_p - x$max_rhat) < 1e-8), nrow(x)),
     all(abs(x$ess_p - x$min_ess) < 1e-8))
}
cl("S3.6", "Robbins-Monro acceptance rates (four blocks; jump)", "0.21-0.26; near 0.4",
   "not stored in these results", NA)
claims <- do.call(rbind, claims)
cat("\n===== Numbers quoted in the text =====\n")
cat(sprintf("%-15s %-68s | paper: %-31s | reproduced: %-28s | %s\n", claims$where, claims$quantity, claims$paper,
            claims$reproduced, ifelse(is.na(claims$ok), "n/a", ifelse(claims$ok, "ok", "MISMATCH"))), sep = "")

## =========================================================================
## Optional comparison with the paper source (STAR_PAPER_DIR/supplement-ba.tex)
## =========================================================================
paper <- file.path(Sys.getenv("STAR_PAPER_DIR", "."), "supplement-ba.tex")
if (file.exists(paper)) {
  L <- readLines(paper)
  norm <- function(s) gsub("\\s+", " ", trimws(gsub("\\\\new\\{([^}]*)\\}", "\\1", s)))
  block <- function(label, env) {
    i <- grep(sprintf("\\label{%s}", label), L, fixed = TRUE)
    a <- max(grep(sprintf("\\begin{%s}", env), L[1:i], fixed = TRUE))
    b <- i - 1 + min(grep(sprintf("\\end{%s}", env), L[i:length(L)], fixed = TRUE)); L[a:b]
  }
  rows_of <- function(x) {                   # body rows: between \midrule and \bottomrule
    inb <- cumsum(grepl("^\\s*\\\\midrule", x)) > cumsum(grepl("^\\s*\\\\bottomrule", x))
    norm(x[inb & !grepl("^\\s*\\\\midrule", x)])
  }
  nums <- function(s) as.numeric(unlist(regmatches(s, gregexpr("-?[0-9]*\\.?[0-9]+", s))))
  compare <- function(name, mine, theirs) {
    mine <- norm(mine); cat(sprintf("\n%s: %d lines generated, %d in the paper\n", name, length(mine), length(theirs)))
    for (k in seq_along(mine)) {
      same <- k <= length(theirs) && mine[k] == theirs[k]
      d <- if (k <= length(theirs) && length(nums(mine[k])) == length(nums(theirs[k])))
        max(abs(nums(mine[k]) - nums(theirs[k]))) else NA
      cat(sprintf("  line %2d: %s", k, if (same) "identical" else sprintf("DIFFERS (max abs diff %s)", format(d))), "\n")
      if (!same) cat("    generated:", mine[k], "\n    paper:    ", if (k <= length(theirs)) theirs[k] else "<none>", "\n")
    }
  }
  cat("\n===== Comparison with", paper, "=====\n")
  compare("Table S2", body_diag, rows_of(block("tab:diag", "table")))
  compare("Table S3", c(body_a, body_b, body_c), rows_of(block("tab:ident", "table")))
  compare("Table S10", c(body_sa, body_sb, body_sc), rows_of(block("tab:sens", "table")))
  fl <- norm(block("fig:ident", "figure"))
  fl <- fl[grepl("^\\\\draw\\[thick.*plot coordinates|^\\\\node\\[[^]]*anchor=west", fl)]
  compare("Figure S1", fig_ident[!grepl("^%", fig_ident)], fl)
} else cat("\n(paper source not found; set STAR_PAPER_DIR to compare with supplement-ba.tex)\n")
