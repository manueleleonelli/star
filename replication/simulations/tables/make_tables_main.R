# make_tables_main.R -- simulation tables and in-text numbers of the paper,
# rebuilt from saved results only (no MCMC or simulation is rerun).
# Run from this folder:  Rscript make_tables_main.R
#
# Inputs:  ../main/results/*.rds and ../competitors/results/*.rds
# Outputs: tab_<name>.tex in this folder, each holding the tabular rows between
#          \midrule and \bottomrule of one tabular; everything is also printed.
#
# Blocks (table, label, inputs)
#   1. main Table 2   tab:key        s1_all, s3_all, ht_all, wc_all, emp_known, emp_s3
#   2. main Table 3   tab:rest       s1_all, s4_all, s8_misspec_ablation, s2_all,
#                                    ht_all, emp_known
#   3. Table S4       tab:calib      s1_all
#   4. Table S5       tab:class      s2_all; rows marked * from s5_gauss_sweep
#   5. Table S6       tab:misspec    s3_truths, s3_all; Hill column from s6_hill
#   6. Table S7       tab:empirical  s1_all, s3_all, ht_all, wc_all, emp_known,
#                                    emp_s3, htkboot_s1 (HT coverage, true margins)
#   7. Table S8       tab:htthresh   htsens_s1, htsens_s3, htsens_s4, s1_all, s3_all, s4_all
#   8. Table S9       tab:misspecabl s8_misspec_ablation
#   9. numbers quoted in the text of main Section 5 and supplement Sections
#      S4.2 to S4.6 (s9_asymmetric and s11_theory_checks enter only here)
# Tables S2 (s10) and S10 (s7) and the text of S4.1 and S4.7 are built by
# make_tables_supp.R. Where the text quotes a difference or ratio, both the
# unrounded value and the one from the rounded table entries are printed.
#
# Conventions: RMSE is over the 50 replicates of a cell; coverages and class
# probabilities are means over replicates; errors are posterior median (or
# point estimate) minus truth. "emp." is the empirical estimator on the known
# exponential margins; HT in Table 2 (right) uses the true margins and in
# Table 2 (left) and Table 3 the semiparametric ones. Bold marks the smallest
# value after rounding, ties included.

MR <- "../main/results"; CR <- "../competitors/results"
rd <- function(dir, f) readRDS(file.path(dir, f))
rmse <- function(x) sqrt(mean(x^2))
QS <- c(0.95, 0.99, 0.999)

## ---- formatting helpers ----------------------------------------------------
num <- function(x, d) {                       # rounded string, no negative zero
  s <- sprintf("%.*f", d, x)
  sub("^-(0\\.?0*)$", "\\1", s)
}
cell <- function(x, d, bold = FALSE)
  ifelse(rep_len(bold, length(x)), sprintf("$\\mathbf{%s}$", num(x, d)), sprintf("$%s$", num(x, d)))
is_min <- function(x, d) { r <- as.numeric(num(x, d)); r == min(r) }
bcells <- function(x, d) cell(x, d, is_min(x, d))   # bold the rounded minimum
row <- function(...) {
  s <- paste(paste(c(...), collapse = " & "), "\\\\")
  gsub(" {2,}", " ", s)
}
emit <- function(name, lines) {
  f <- sprintf("tab_%s.tex", name)
  writeLines(lines, f)
  cat(sprintf("\n%% ---- %s ----\n", f)); cat(lines, sep = "\n")
}

## ---- inputs ----------------------------------------------------------------
s1  <- rd(MR, "s1_all.rds");   s2 <- rd(MR, "s2_all.rds");  s3 <- rd(MR, "s3_all.rds")
tr3 <- rd(MR, "s3_truths.rds"); s4 <- rd(MR, "s4_all.rds"); s5 <- rd(MR, "s5_gauss_sweep.rds")
s6  <- rd(MR, "s6_hill.rds");  s8 <- rd(MR, "s8_misspec_ablation.rds")
s9  <- rd(MR, "s9_asymmetric.rds"); s11 <- rd(MR, "s11_theory_checks.rds")
ht  <- rd(CR, "ht_all.rds");   wc <- rd(CR, "wc_all.rds");  ek <- rd(CR, "emp_known.rds")
es3 <- rd(CR, "emp_s3.rds");   hkb <- rd(CR, "htkboot_s1.rds")
hs1 <- rd(CR, "htsens_s1.rds"); hs3 <- rd(CR, "htsens_s3.rds"); hs4 <- rd(CR, "htsens_s4.rds")

## ---- Study 1 cells (within the model, q = 0.99) ------------------------------
key1 <- function(d) sprintf("%.2f_%d_%.1f", d$chibar0, as.integer(d$n), d$u1)
by1  <- function(d, v, f) tapply(d[[v]], key1(d), f)
k1   <- function(cb, n, u1 = 2.2) sprintf("%.2f_%d_%.1f", cb, as.integer(n), u1)
CB   <- c(0.15, 0.45, 0.75, 0.95)
s1q  <- s1[s1$quantity == "chi_q", ]
hs1$es <- hs1$semi - hs1$truth; hs1$ek <- hs1$known - hs1$truth
W <- list(                                   # per-cell summaries at q = 0.99
  star_rmse = by1(s1q, "err", rmse),       star_cov = by1(s1q, "cover", mean),
  ht_rmse   = by1(ht$s1, "ht_err", rmse),  ht_cov   = by1(ht$s1, "ht_cover", mean),
  htk_cov   = by1(hkb, "cover", mean),
  wc_rmse   = by1(wc$s1, "wc_err", rmse),
  emp_rmse  = by1(ek$s1, "emp_err", rmse), emp_cov  = by1(ek$s1, "emp_cover", mean),
  truth     = by1(s1q, "truth", mean))

## ---- Study 3 generators (copulas, n = 4000) in long format ------------------
long3 <- function(d, cols, gen = d$gen)
  do.call(rbind, lapply(1:3, function(j)
    data.frame(gen = gen, seed = d$seed, q = QS[j], err = d[[cols[j]]])))
key3 <- function(d) paste(d$gen, d$q)
by3  <- function(d, v = "err", f = rmse) tapply(d[[v]], key3(d), f)
k3   <- function(g, q) paste(g, q)
hs3$es <- hs3$semi - hs3$truth; hs3$ek <- hs3$known - hs3$truth
V <- list(                                   # RMSE of chi(q) per generator and q
  star = by3(long3(s3, c("chi95_err", "chi99_err", "chi999_err"))),
  htk  = by3(ht$s3, "htk_err"),              # HT, true (exponential) margins
  hte  = by3(ht$s3, "ht_err"),               # HT, semiparametric margins
  wc   = by3(wc$s3, "wc_err"),
  emp  = by3(long3(es3, c("ek95", "ek99", "ek999"))))
G3 <- c("gauss05", "invlog05", "frank5", "t05_4")
GL3 <- c(gauss05 = "Gaussian, $\\varrho=0.5$",
         invlog05 = "Inverted logistic, $\\alpha_{\\mathrm{IL}}=0.5$",
         frank5 = "Frank, $\\theta=5$", t05_4 = "Student $t$, $(0.5,4)$")

## ---- Study 4 (ablation within the model, chibar = 0.6, q = 0.999) -----------
tr4 <- unique(hs4$truth); f4 <- s4[s4$variant == "full", ]
VAR <- c("full", "0.8", "0.9", "0.95")
A4  <- sapply(c("chi999_med", "chi999_cover", "eta_w", "eta_med"), function(v)
  tapply(s4[[v]], s4$variant, mean))[VAR, ]
A8  <- aggregate(cbind(chi95_cover, chi99_cover, chi999_cover, chi999_err, chi999_w)
                 ~ gen + variant, s8, mean)
a8  <- function(g, v, col) A8[A8$gen == g & A8$variant == v, col]

## ---- Study 2 and Study 5 (classification) -----------------------------------
A2 <- aggregate(cbind(P1, correct_05, chibar_mean) ~ gen + n, s2, mean)
a2 <- function(g, n, col) A2[A2$gen == g & A2$n == n, col]
s5$correct_05 <- as.numeric(s5$P1 <= 0.5)       # all Gaussian truths are AI
A5 <- aggregate(cbind(P1, correct_05, chibar_mean) ~ rho + n, s5, mean)
a5 <- function(r, n, col) A5[A5$rho == r & A5$n == n, col]

######################################################################
## 1. Table 2 (tab:key): RMSE of chi(q) at n = 4000
######################################################################
key_left <- lapply(CB, function(cb) {
  k <- k1(cb, 4000)
  bcells(c(W$star_rmse[k], W$ht_rmse[k], W$wc_rmse[k], W$emp_rmse[k]), 3) })
key_right <- lapply(QS, function(q) unlist(lapply(c("gauss05", "invlog05"), function(g) {
  k <- k3(g, q); bcells(c(V$star[k], V$htk[k], V$wc[k], V$emp[k]), 3) })))
tab_key <- sapply(1:4, function(i)
  row(cell(CB[i], 2), key_left[[i]], "",
      if (i <= 3) c(cell(QS[i], if (QS[i] == 0.999) 3 else 2), key_right[[i]]) else rep("", 9)))
emit("key", tab_key)

######################################################################
## 2. Table 3 (tab:rest): coverage, threshold ablation, classification
######################################################################
rest_left <- c(lapply(CB, function(cb) { k <- k1(cb, 4000)
  c(sprintf("%s, $0.99$", cell(cb, 2)),
    cell(c(W$star_cov[k], W$ht_cov[k], W$emp_cov[k]), 2)) }),
  list(c("$0.6$, $0.999$", cell(c(mean(f4$chi999_cover), mean(ht$s4$ht_cover),
                                   mean(ek$s4$emp_cover)), 2))))
VLAB <- c(full = "estimated", "0.8" = "$\\kappa=0.8$", "0.9" = "$\\kappa=0.9$",
          "0.95" = "$\\kappa=0.95$")
rest_centre <- lapply(VAR, function(v)
  c(VLAB[[v]], cell(A4[v, "chi999_med"], 3),
    cell(c(A4[v, "chi999_cover"], a8("gauss05", v, "chi999_cover"),
           a8("frank5", v, "chi999_cover")), 2)))
G2 <- c("m_chi0", "m_chi05", "m_chi095", "m_chi1", "gauss05", "gauss09",
        "invlog05", "logist05", "t05_4")
AD2 <- c(m_chi1 = TRUE, logist05 = TRUE, t05_4 = TRUE)
GL2 <- c(m_chi0 = "\\new{STAR}, $\\bar{\\chi}=0$", m_chi05 = "\\new{STAR}, $\\bar{\\chi}=0.5$",
         m_chi095 = "\\new{STAR}, $\\bar{\\chi}=0.95$", m_chi1 = "\\new{STAR}, $\\bar{\\chi}=1$",
         gauss05 = "Gaussian, $\\varrho=0.5$", gauss09 = "Gaussian, $\\varrho=0.9$",
         invlog05 = "inverted logistic, $0.5$", logist05 = "logistic, $0.5$",
         t05_4 = "Student $t$, $(0.5,4)$")
cls <- function(g) if (isTRUE(AD2[g])) "AD" else "AI"
tab_rest <- sapply(seq_along(G2), function(i) {
  g <- G2[i]
  row(if (i <= 5) rest_left[[i]] else rep("", 4), "",
      if (i <= 4) rest_centre[[i]] else rep("", 5), "",
      GL2[[g]], cls(g), cell(a2(g, 8000, "P1"), 2)) })
emit("rest", tab_rest)

######################################################################
## 3. Table S4 (tab:calib): calibration under the model
######################################################################
CELLS <- data.frame(cb = c(CB, CB, 0.45), n = c(rep(1000, 4), rep(4000, 5)),
                    u1 = c(rep(2.2, 8), 2.8))
QTY <- list(chibar = c(2, 2), chi_q = c(3, 3), lam03 = c(3, 3), rho0.5 = c(2, 2))
tab_calib <- sapply(seq_len(nrow(CELLS)), function(i) {
  cl <- CELLS[i, ]; k <- k1(cl$cb, cl$n, cl$u1)
  blocks <- unlist(lapply(names(QTY), function(qn) {
    d <- s1[s1$quantity == qn, ]; dg <- QTY[[qn]]
    c(cell(by1(d, "err", mean)[k], dg[1]), cell(by1(d, "err", rmse)[k], dg[2]),
      cell(by1(d, "cover", mean)[k], 2)) }))
  lab <- if (cl$u1 == 2.2) cell(cl$cb, 2) else sprintf("%s, asym.", cell(cl$cb, 2))
  row(lab, cell(cl$n, 0), blocks) })
emit("calib", tab_calib)

######################################################################
## 4. Table S5 (tab:class): three-component prior as a classifier
######################################################################
CL <- list(                                   # label, class, source
  list("\\new{STAR}, $\\bar{\\chi}=0$", "AI", "m_chi0"),
  list("\\new{STAR}, $\\bar{\\chi}=0.5$", "AI", "m_chi05"),
  list("\\new{STAR}, $\\bar{\\chi}=0.95$", "AI", "m_chi095"),
  list("\\new{STAR}, $\\bar{\\chi}=1$", "AD", "m_chi1"),
  list("Gaussian, $\\varrho=0.3^{*}$", "AI", 0.3),
  list("Gaussian, $\\varrho=0.5$", "AI", "gauss05"),
  list("Gaussian, $\\varrho=0.7^{*}$", "AI", 0.7),
  list("Gaussian, $\\varrho=0.9$", "AI", "gauss09"),
  list("Inverted logistic, $\\alpha_{\\mathrm{IL}}=0.5$", "AI", "invlog05"),
  list("Logistic, $\\alpha_{\\mathrm{L}}=0.5$", "AD", "logist05"),
  list("Student $t$, $(0.5,4)$", "AD", "t05_4"))
tab_class <- sapply(CL, function(r) {
  src <- r[[3]]
  get <- if (is.numeric(src)) function(n, col) a5(src, n, col) else function(n, col) a2(src, n, col)
  vals <- unlist(lapply(c("P1", "correct_05", "chibar_mean"), function(col)
    c(get(2000, col), get(8000, col))))
  row(r[[1]], r[[2]], cell(vals, 2)) })
emit("class", tab_class)

######################################################################
## 5. Table S6 (tab:misspec): misspecification at n = 4000
######################################################################
GL6 <- GL3; GL6["invlog05"] <- "Inverted logistic, $0.5$"
M3 <- aggregate(cbind(eta_err, chi99_err, chi999_err, lam03_err) ~ gen, s3, mean)
H6 <- tapply(s6$err, s6$gen, mean)            # Hill-type estimator of eta (Study 6)
tab_misspec <- sapply(G3, function(g) {
  t <- tr3[[g]]; m <- M3[M3$gen == g, ]
  tv <- c(t$eta, t$chi[2], t$chi[3], t$lam03)
  dg <- ifelse(abs(tv) < 0.01, 3, 2)             # 3 decimals for truths below 0.01
  row(GL6[[g]], mapply(cell, tv, dg), "",
      cell(c(m$eta_err, H6[[g]], m$chi99_err, m$chi999_err, m$lam03_err), 2)) })
emit("misspec", tab_misspec)

######################################################################
## 6. Table S7 (tab:empirical): comparison with competing estimators
######################################################################
tab_emp_top <- unlist(lapply(CB, function(cb) sapply(c(1000, 4000), function(n) {
  k <- k1(cb, n)
  row(cell(cb, 2), cell(n, 0), cell(W$truth[k], 3),
      cell(c(W$star_rmse[k], W$ht_rmse[k], W$wc_rmse[k], W$emp_rmse[k]), 3),
      cell(c(W$star_cov[k], W$ht_cov[k], W$htk_cov[k], W$emp_cov[k]), 2)) })))
emit("empirical_top", tab_emp_top)
tab_emp_bot <- unlist(lapply(G3, function(g) c(
  sapply(QS, function(q) { k <- k3(g, q)
    row(if (q == QS[1]) GL3[[g]] else "", cell(q, if (q == 0.999) 3 else 2),
        cell(tr3[[g]]$chi[QS == q], 3),
        cell(c(V$star[k], V$htk[k], V$hte[k], V$wc[k], V$emp[k]), 3)) }),
  if (g != tail(G3, 1)) "\\addlinespace")))
emit("empirical_bottom", tab_emp_bot)

######################################################################
## 7. Table S8 (tab:htthresh): HT conditioning threshold x margins
######################################################################
KD <- c(0.8, 0.85, 0.9, 0.95)
ES1 <- lapply(c("es", "ek"), function(v) sapply(KD, function(kd) by1(hs1[hs1$kd == kd, ], v, rmse)))
hs4$es <- hs4$semi - hs4$truth; hs4$ek <- hs4$known - hs4$truth
ES4 <- lapply(c("es", "ek"), function(v) sapply(KD, function(kd) rmse(hs4[[v]][hs4$kd == kd])))
tab_ht_top <- c(
  unlist(lapply(seq_len(nrow(CELLS)), function(i) {
    cl <- CELLS[i, ]; k <- k1(cl$cb, cl$n, cl$u1)
    lab <- if (cl$u1 == 2.2) cell(cl$cb, 2) else sprintf("%s, asym.", cell(cl$cb, 2))
    row(lab, cell(cl$n, 0), "$0.99$",
        bcells(c(W$star_rmse[k], ES1[[1]][k, ], ES1[[2]][k, ]), 3)) })),
  row("$0.60$", "$4000$", "$0.999$",
      bcells(c(rmse(f4$chi999_med - tr4), ES4[[1]], ES4[[2]]), 3)))
emit("htthresh_top", tab_ht_top)
ES3 <- lapply(c("es", "ek"), function(v) sapply(KD, function(kd) by3(hs3[hs3$kd == kd, ], v)))
tab_ht_bot <- unlist(lapply(G3, function(g) c(
  sapply(QS, function(q) { k <- k3(g, q)
    row(if (q == QS[1]) GL3[[g]] else "", cell(q, if (q == 0.999) 3 else 2),
        bcells(c(V$star[k], ES3[[1]][k, ], ES3[[2]][k, ]), 3)) }),
  if (g != tail(G3, 1)) "\\addlinespace")))
emit("htthresh_bottom", tab_ht_bot)

######################################################################
## 8. Table S9 (tab:misspecabl): ablation under misspecification
######################################################################
GL9 <- c(gauss05 = "Gaussian, $\\varrho=0.5$", frank5 = "Frank, $\\theta=5$")
VL9 <- c(full = "full", VLAB[-1])
tab_mabl <- unlist(lapply(names(GL9), function(g) c(
  sapply(VAR, function(v)
    row(if (v == "full") GL9[[g]] else "", VL9[[v]],
        cell(c(a8(g, v, "chi95_cover"), a8(g, v, "chi99_cover"), a8(g, v, "chi999_cover")), 2),
        cell(c(a8(g, v, "chi999_err"), a8(g, v, "chi999_w")), 3))),
  if (g == "gauss05") "\\midrule")))
emit("misspecabl", tab_mabl)

######################################################################
## 9. Numbers quoted in the text (computed value, then value in the paper)
######################################################################
say <- function(what, x, paper) cat(sprintf("  %-66s %-22s [paper: %s]\n", what,
  paste(x, collapse = ", "), paper))
rng <- function(x, d = 2) paste(num(range(x), d), collapse = " to ")
r3 <- function(x) as.numeric(num(x, 3))     # value as printed in a 3-decimal table
cat("\n==== replicates per cell (text: fifty; Study 11: fifteen fits) ====\n")
nrep <- c(table(key1(s1q)), table(s2$gen, s2$n), table(s3$gen), table(s4$variant),
          table(s5$rho, s5$n), table(s6$gen), table(s8$gen, s8$variant),
          nrow(s9$within), nrow(s9$khoud))
say("Studies 1 to 6, 8 and 9, range over cells", rng(nrep, 0), "50")
say("Study 11 fits", length(unique(s11$seed)), "15")

cat("\n==== Main Section 5 ====\n")
cc <- sapply(CB, function(cb) k1(cb, 4000))
say("q=0.99, n=4000: chibar where STAR has smallest RMSE",
    CB[sapply(cc, function(k) which.min(c(W$star_rmse[k], W$ht_rmse[k], W$wc_rmse[k], W$emp_rmse[k])) == 1)],
    "all but the weakest (0.45, 0.75, 0.95)")
say("STAR coverage, chibar 0.15 to 0.75 ('near nominal')", rng(W$star_cov[cc[1:3]]), "near 0.90")
say("HT coverage, chibar 0.15 to 0.95 ('seven in ten')", rng(W$ht_cov[cc]), "about 0.7")
say("emp. coverage, smallest ('as few as three')", num(min(W$emp_cov[cc]), 2), "0.30")
say("q=0.999 coverage STAR, HT, emp.", num(c(mean(f4$chi999_cover), mean(ht$s4$ht_cover),
    mean(ek$s4$emp_cover)), 2), "0.90, 0.58, 0.04")
gi <- c("gauss05", "invlog05")
say("Gauss/IL q=0.99: STAR smallest RMSE",
    sapply(gi, function(g) { k <- k3(g, 0.99); V$star[k] < min(V$htk[k], V$wc[k], V$emp[k]) }), "TRUE, TRUE")
beh <- sapply(gi, function(g) sapply(c(0.95, 0.999), function(q) { k <- k3(g, q)
  c(raw = V$star[k] - V$htk[k], rounded = r3(V$star[k]) - r3(V$htk[k])) }))
say("STAR minus HT(true) at q=0.95/0.999, max (unrounded)", num(max(beh[c(1, 3), ]), 4), "0.011")
say("STAR minus HT(true), max of rounded table entries", num(max(beh[c(2, 4), ]), 3), "0.011")
emp_ratio <- sapply(gi, function(g) { k <- k3(g, 0.999)
  V$emp[k] / c(STAR = V$star[k], HT = V$htk[k], WC = V$wc[k]) })
say("q=0.999: emp. RMSE / (STAR, HT, WC RMSE), range", rng(emp_ratio, 1), "three to five times")
say("best HT level, Gauss/IL (est. and true margins, by q)",
    unlist(lapply(gi, function(g) sapply(QS, function(q) { k <- k3(g, q)
      paste(KD[which.min(ES3[[1]][k, ])], KD[which.min(ES3[[2]][k, ])], sep = "/") }))), "low")
say("best HT level, chibar 0.75/0.95, n=4000 (est./true)",
    sapply(cc[3:4], function(k) paste(KD[which.min(ES1[[1]][k, ])], KD[which.min(ES1[[2]][k, ])], sep = "/")),
    "high")
say("Frank: STAR mean error chi(0.99), chi(0.999) ('overestimated')",
    num(unlist(M3[M3$gen == "frank5", c("chi99_err", "chi999_err")]), 3), "positive")
say("t: STAR mean error chi(0.99), chi(0.999) ('underestimates')",
    num(unlist(M3[M3$gen == "t05_4", c("chi99_err", "chi999_err")]), 3), "negative")
WE <- tapply(wc$s3$eta_err[wc$s3$q == 0.99], wc$s3$gen[wc$s3$q == 0.99], mean)
say("mean eta error Gauss, IL: STAR", num(M3$eta_err[match(gi, M3$gen)], 2), "inflated")
say("mean eta error Gauss, IL: WC (geometric)", num(WE[gi], 2), "not inflated")
say("chi(0.999) median, kappa 0.8 -> 0.95: relative change (raw)",
    num(1 - A4["0.95", "chi999_med"] / A4["0.8", "chi999_med"], 2), "forty percent")
say("  same from the rounded medians", num(1 - r3(A4["0.95", "chi999_med"]) / r3(A4["0.8", "chi999_med"]), 2),
    "forty percent")
say("Gaussian coverage chi(0.999): kappa 0.8, 0.9, estimated, 0.95",
    num(c(a8("gauss05", "0.8", "chi999_cover"), a8("gauss05", "0.9", "chi999_cover"),
          a8("gauss05", "full", "chi999_cover"), a8("gauss05", "0.95", "chi999_cover")), 2),
    "0.10, 0.48, 0.78, high level as well")
say("Frank coverage chi(0.999), largest over variants", num(max(A8$chi999_cover[A8$gen == "frank5"]), 2),
    "no level covers")
p8 <- sapply(G2, function(g) a2(g, 8000, "P1"))
say("n=8000 P(AD): atoms chibar=0 and 1, interior slab 0.5",
    num(p8[c("m_chi0", "m_chi1", "m_chi05")], 2), "0.00, 0.81, 0.06")
say("n=8000 P(AD), AD truths", rng(p8[c("m_chi1", "logist05", "t05_4")]), "0.77 to 0.97")
say("n=8000 P(AD), strongly dependent AI truths",
    rng(p8[c("gauss05", "gauss09", "invlog05", "m_chi095")]), "0.60 to 0.82")

cat("\n==== Supplement S4.2 (calibration and further checks) ====\n")
s1c <- s1[s1$quantity %in% c("chibar", "chi_q") & s1$u1 == 2.8, ]
say("asym. cell coverage, chibar and chi(0.99)", num(tapply(s1c$cover, s1c$quantity, mean)[c("chibar", "chi_q")], 2),
    "0.82 to 0.84")
say("rho(1/2) bias, all cells", rng(by1(s1[s1$quantity == "rho0.5", ], "err", mean)), "biased upwards")
say("rho(1/2) coverage, smallest", num(min(by1(s1[s1$quantity == "rho0.5", ], "cover", mean)), 2), "0.78")
say("chibar=0.95 coverage of chi(0.99), n=1000 and 4000",
    num(W$star_cov[c(k1(0.95, 1000), k1(0.95, 4000))], 2), "0.54, 0.78")
say("within-model asymmetry rho(0.3)-rho(0.7): fitted, true",
    num(c(mean(s9$within$rho_asym_hat), mean(s9$within$rho_asym_true)), 2), "-0.91, -0.92")
say("within-model asymmetry sigma(0.3)-sigma(0.7): true, fitted",
    num(c(mean(s9$within$sig_asym_true), mean(s9$within$sig_asym_hat)), 2), "0.11, 0.02")
say("Khoudraji-Gumbel P(AD | data)", num(mean(s9$khoud$P1), 2), "0.96")
s11$t <- -log(1 - s11$q)
fit11 <- summary(lm(log(chi_num / chi_theory) ~ t, data = s11))$coefficients
say("log{chi(q)/(c/t)} ~ t, 15 fits x 4 levels: intercept (SE)",
    sprintf("%s (%s)", num(fit11[1, 1], 3), num(fit11[1, 2], 3)), "0.003 (0.037)")
say("  slope", num(fit11[2, 1], 3), "-0.057")
say("  -(1/eta - 1) at eta = 0.95 (eta not stored in s11)", num(-(1 / 0.95 - 1), 3), "-0.053")
say("  eta implied by the slope, 1/(1 - slope)", num(1 / (1 - fit11[2, 1]), 3), "0.95")
say("  q range", paste(range(s11$q)), "0.9 to 0.9999")

cat("\n==== Supplement S4.3 (classification) ====\n")
say("P(AD), chibar=0 and 0.5, largest", num(max(A2$P1[A2$gen %in% c("m_chi0", "m_chi05")]), 2), "at most 0.23")
say("correct decisions, AD truths", rng(A2$correct_05[A2$gen %in% c("m_chi1", "logist05", "t05_4")]),
    "86 to 100 percent")
say("P(AD), strongly dependent AI (both n)",
    rng(A2$P1[A2$gen %in% c("gauss05", "gauss09", "invlog05", "m_chi095")]), "0.60 to 0.82")
say("Gaussian family P(AD) at rho=0.3 and 0.9 (both n)",
    num(c(a5(0.3, 2000, "P1"), a5(0.3, 8000, "P1"), a2("gauss09", 2000, "P1"), a2("gauss09", 8000, "P1")), 2),
    "about 1/2 to about 4/5")

cat("\n==== Supplement S4.4 (misspecification) ====\n")
say("Gauss, IL: largest |error| in chi(0.99), chi(0.999)",
    num(max(abs(unlist(M3[M3$gen %in% gi, c("chi99_err", "chi999_err")]))), 3), "within 0.03")
say("Gauss, IL: eta error", num(M3$eta_err[match(gi, M3$gen)], 2), "0.16, 0.22")
say("Gauss, IL: Hill eta error", num(H6[gi], 2), "smaller")
say("t: chi(0.999) error; Hill eta error", num(c(M3$chi999_err[M3$gen == "t05_4"], H6[["t05_4"]]), 2),
    "-0.16, -0.19")

cat("\n==== Supplement S4.5 (competing estimators) ====\n")
c8 <- with(CELLS[1:8, ], k1(cb, n))
say("within model: cells where STAR has smallest RMSE (of 8)",
    sum(sapply(c8, function(k) which.min(c(W$star_rmse[k], W$ht_rmse[k], W$wc_rmse[k], W$emp_rmse[k])) == 1)),
    "all but one")
say("HT coverage over the 8 cells, either margins", rng(c(W$ht_cov[c8], W$htk_cov[c8])), "0.52 to 0.84")
dbl <- sapply(gi, function(g) { k <- k3(g, 0.999)
  c(est = ES3[[1]][k, 4] / ES3[[1]][k, 1], true = ES3[[2]][k, 4] / ES3[[2]][k, 1]) })
say("Gauss, IL at q=0.999: HT RMSE ratio kappa 0.95 / 0.8", num(dbl, 1), "roughly doubles")
strong <- c(k1(0.75, 1000), k1(0.95, 1000), k1(0.75, 4000), k1(0.95, 4000))
rat_raw <- sapply(strong, function(k) c(ES1[[1]][k, ], ES1[[2]][k, ]) / W$star_rmse[k])
rat_rnd <- sapply(strong, function(k) r3(c(ES1[[1]][k, ], ES1[[2]][k, ])) / r3(W$star_rmse[k]))
say("chibar >= 0.75: HT RMSE / STAR RMSE, unrounded", rng(rat_raw, 1), "1.2 to 2.6")
say("  same from the rounded table entries", rng(rat_rnd, 1), "1.2 to 2.6")

cat("\n==== Supplement S4.6 (value of self-thresholding) ====\n")
say("chi(0.999) median: kappa 0.8, 0.95, estimated", num(A4[c("0.8", "0.95", "full"), "chi999_med"], 3),
    "0.021, 0.013, 0.021")
say("eta 90% width: kappa 0.8, 0.95, estimated", num(A4[c("0.8", "0.95", "full"), "eta_w"], 2),
    "0.23, 0.31, 0.23")
say("chi(0.999) coverage over variants", rng(A4[, "chi999_cover"]), "0.90 to 0.96")
