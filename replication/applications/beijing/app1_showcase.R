## app1_showcase.R -- every posterior summary reported in Section 6.1 of the paper, Table 4 and the
## six-panel Figure 5 (fig-beijing.tex; panel (f) is replaced by sec6_cut.R), from the stored Wanliu fit
## results/fit_Wanliu.rds produced by beijing_stations.R. Run from this folder: source("app1_showcase.R").
suppressMessages(library(star))
x <- readRDS("results/fit_Wanliu.rds"); fit <- x$fit; E <- x$E; D <- fit$draws; n <- nrow(E)
a <- x$data; gx <- x$margins$gx; gy <- x$margins$gy
r <- rowSums(E); w <- E[, 1] / r
q90 <- function(v) quantile(v, c(0.05, 0.5, 0.95), names = FALSE)
gl <- seq(0.5/201, 1 - 0.5/201, length.out = 201); glw <- rep(1/201, 201)
idx <- round(seq(1, nrow(D), length.out = 400))
th <- lapply(idx, function(s) star:::theta_row(fit, s))
wg <- seq(0.005, 0.995, length.out = 199)

## ---- threshold, exceedances, tail scale ----
rho <- t(sapply(th, function(t) rho_star(wg, t$u1, t$u2, t$p)))
sig <- t(sapply(th, function(t) sigma_star(wg, t$chibar, t$z)))
lam <- t(sapply(th, function(t) lambda_star(wg, t$b)))
piw <- t(sapply(seq_along(th), function(i) pgamma(rho[i, ], shape = th[[i]]$alpha, rate = lam[i, ])))
exc <- sapply(th, function(t) mean(r > rho_star(w, t$u1, t$u2, t$p)))
cat("exceedance fraction", round(q90(exc), 3), " effective tail n", round(q90(exc) * n), "\n")
cat("p", round(q90(D[, "p"]), 2), " u1", round(q90(D[, "u1"]), 2), " u2", round(q90(D[, "u2"]), 2), " alpha", round(q90(D[, "alpha"]), 2), "\n")
cat("rho(1/2)", round(q90(rho[, which.min(abs(wg - 0.5))]), 2), " rho on axes (w=0.005)", round(q90(rho[, 1]), 2), " q95 of exp", round(-log(0.05), 2), "\n")
## marginal exponential quantile of the axis anchors
cat("anchor u1 as exp quantile level", round(q90(1 - exp(-D[, "u1"])), 3), "\n")
cat("sigma(1/2) = 1 + chibar", round(q90(1 + D[, "chibar"]), 2), "\n")

## ---- classes and coefficients ----
cp <- class_probs(fit); cat("class probs", round(cp, 3), "\n")
cb <- D[, "chibar"]; sl <- cb[fit$component == "slab"]
cat("chibar marginal", round(q90(cb), 3), " slab-conditional", round(q90(sl), 3), " eta", round(q90((1 + cb) / 2), 3), "\n")
cat("P(chibar > 0.5)", round(mean(cb > 0.5), 3), " P(AI) =", round(1 - cp[3], 3), "\n")

## ---- chi(q) curve, chibar(q), empirical ----
qs <- 1 - 10^(-seq(0.7, 3, length.out = 47)); tq <- -log(1 - qs)
js <- t(sapply(th, function(t) sapply(tq, function(tt) star:::joint_surv_cpp(tt, tt, gl, glw, t$phi, t$alpha, t$b, t$u1, t$u2, t$p, t$chibar, t$z))))
chi <- sweep(js, 2, 1 - qs, "/")
chibarq <- 2 * log(1 - matrix(qs, nrow(js), ncol(js), byrow = TRUE)) / log(js) - 1
chi_s <- apply(chi, 2, q90); chibar_s <- apply(chibarq, 2, q90)
qe <- c(0.8, 0.85, 0.9, 0.95, 0.975, 0.99); te <- -log(1 - qe)
chi_emp <- sapply(te, function(tt) mean(E[, 1] > tt & E[, 2] > tt)) / (1 - qe)
for (qq in c(0.9, 0.95, 0.99, 0.999)) { j <- which.min(abs(qs - qq)); cat(sprintf("chi(%.3f): mean %.3f  90%% (%.3f, %.3f)  chibar(q) %.2f (%.2f, %.2f)\n", qq, mean(chi[, j]), chi_s[1, j], chi_s[3, j], mean(chibarq[, j]), chibar_s[1, j], chibar_s[3, j])) }
cat("empirical chi at", qe, ":", round(chi_emp, 3), "\n")

## ---- angular dependence function, model and Hill-type empirical ----
om <- seq(0.05, 0.95, by = 0.025); grid <- seq(0.001, 0.999, length.out = 999)
adf <- t(sapply(th, function(t) { sg <- sigma_star(grid, t$chibar, t$z); sapply(om, function(o) min(pmax(o / grid, (1 - o) / (1 - grid)) / sg)) }))
adf_s <- apply(adf, 2, q90); adf_m <- colMeans(adf)
hill_adf <- function(o, kappa = 0.9) { Tm <- pmin(E[, 1] / o, E[, 2] / (1 - o)); u <- quantile(Tm, kappa); 1 / mean(Tm[Tm > u] - u) }
adf_emp <- sapply(om, hill_adf)
for (o in c(0.2, 0.35, 0.5)) { j <- which.min(abs(om - o)); cat(sprintf("Lambda(%.2f): mean %.3f  90%% (%.3f, %.3f)  Hill %.3f\n", o, adf_m[j], adf_s[1, j], adf_s[3, j], adf_emp[j])) }

## ---- angular density of exceedances ----
fW <- t(sapply(th, function(t) { K <- length(t$phi); rowSums(sapply(1:K, function(k) t$phi[k] * dbeta(wg, k, K - k + 1))) }))
fWe <- fW * (1 - piw); fWe <- fWe / rowSums(fWe) / (wg[2] - wg[1])
fWe_m <- colMeans(fWe); cat("exceedance-angle density: mean at w=0.5", round(fWe_m[which.min(abs(wg-0.5))], 2), " mass in (0.35,0.65)", round(mean(rowSums(fWe[, wg > 0.35 & wg < 0.65]) * (wg[2]-wg[1])), 3), " mass w<0.2", round(mean(rowSums(fWe[, wg < 0.2]) * (wg[2]-wg[1])), 3), " mass w>0.8", round(mean(rowSums(fWe[, wg > 0.8]) * (wg[2]-wg[1])), 3), "\n")
fW_m <- colMeans(fW); cat("body angular density mass in (0.35,0.65)", round(mean(rowSums(fW[, wg > 0.35 & wg < 0.65]) * (wg[2]-wg[1])), 3), "\n")

## ---- return level sets on the original scale ----
Qsemi <- function(p, xdat, g) {           # semiparametric marginal quantile at probability p
  out <- quantile(xdat, p, names = FALSE, type = 6)
  ab <- p > 1 - g$p_u
  out[ab] <- g$u + g$sigma / g$xi * (((1 - p[ab]) / g$p_u)^(-g$xi) - 1)
  out
}
rls <- function(i, pr) {                   # radius with P(R > r | W = w) = pr, for draw i
  t <- th[[i]]; rb <- qgamma(1 - pr, shape = t$alpha, rate = lam[i, ])
  tail <- (1 - piw[i, ]) > pr
  rb[tail] <- rho[i, tail] + sig[i, tail] * log((1 - piw[i, tail]) / pr)
  rb
}
to_orig <- function(rr) cbind(Qsemi(1 - exp(-rr * wg), a$PM2.5, gx), Qsemi(1 - exp(-rr * (1 - wg)), a$O3, gy))
Ts <- c(1, 5, 20); pr <- 1 / (183 * Ts)
rl <- lapply(pr, function(p) t(sapply(seq_along(th), function(i) rls(i, p))))
rl_med <- lapply(rl, function(m) apply(m, 2, median)); rl_lo <- lapply(rl, function(m) apply(m, 2, quantile, 0.05)); rl_hi <- lapply(rl, function(m) apply(m, 2, quantile, 0.95))
rho_med <- apply(rho, 2, median)
for (k in seq_along(Ts)) { j <- which.min(abs(wg - 0.5)); o <- to_orig(rl_med[[k]]); ol <- to_orig(rl_lo[[k]]); oh <- to_orig(rl_hi[[k]])
  cat(sprintf("T=%d seasons: diagonal point (PM2.5, O3) = (%.0f, %.0f), 90%% (%.0f, %.0f) to (%.0f, %.0f); axis intercepts PM2.5 %.0f O3 %.0f\n", Ts[k], o[j,1], o[j,2], ol[j,1], ol[j,2], oh[j,1], oh[j,2], o[1,1], o[199,2])) }
o <- to_orig(rho_med); cat("threshold curve on original scale: diagonal", round(o[which.min(abs(wg-0.5)), ]), " axes", round(c(o[1,1], o[199,2])), "\n")
cat("original-scale quantiles PM2.5 0.9/0.95/0.99:", round(quantile(a$PM2.5, c(.9,.95,.99))), " O3:", round(quantile(a$O3, c(.9,.95,.99))), " max PM2.5", max(a$PM2.5), " max O3", max(a$O3), "\n")
cat("empirical count outside the T=1 median set:", sum(sapply(1:n, function(i) { j <- which.min(abs(wg - w[i])); r[i] > rl_med[[1]][j] })), " expected", n / 183, "\n")
cat("empirical count above threshold median:", sum(sapply(1:n, function(i) { j <- which.min(abs(wg - w[i])); r[i] > rho_med[j] })), "\n")

## ---- posterior predictive check: joint exceedance counts ----
## N ~ Binomial(n, p(theta)) with p the joint survival of the level-q quantiles; the predictive mean is exact
pcheck <- function(q) { tt <- -log(1 - q)
  p <- sapply(th, function(t) star:::joint_surv_cpp(tt, tt, gl, glw, t$phi, t$alpha, t$b, t$u1, t$u2, t$p, t$chibar, t$z))
  set.seed(3); N <- rbinom(length(p), n, p); obs <- sum(E[, 1] > tt & E[, 2] > tt)
  c(q = q, obs = obs, mean = n * mean(p), lo = quantile(N, 0.05, names = FALSE), hi = quantile(N, 0.95, names = FALSE), p_le = mean(N <= obs)) }
pc <- t(sapply(c(0.8, 0.85, 0.9, 0.95, 0.975, 0.99), pcheck))
cat("posterior predictive check of joint exceedance counts:\n"); print(round(pc, 3))
## ---- days outside the one-season return level set: observed against posterior predictive ----
dw <- wg[2] - wg[1]
jw <- sapply(w, function(v) which.min(abs(wg - v))); out1 <- r > rl_med[[1]][jw]
p_out <- sapply(seq_along(th), function(i) { sv <- ifelse(rl_med[[1]] > rho[i, ], (1 - piw[i, ]) * exp(-(rl_med[[1]] - rho[i, ]) / sig[i, ]), 1 - pgamma(rl_med[[1]], shape = th[[i]]$alpha, rate = lam[i, ]))
  sum(fW[i, ] * sv) * dw / sum(fW[i, ] * dw) })
set.seed(3); N_out <- rbinom(length(p_out), n, p_out)
cat(sprintf("days outside the one-season set: observed %d (%d with w < 0.25, %d with w > 0.75), predictive mean %.1f, 90%% (%.0f, %.0f), P(N >= obs) %.2f\n", sum(out1), sum(out1 & w < 0.25), sum(out1 & w > 0.75),
            n * mean(p_out), quantile(N_out, 0.05), quantile(N_out, 0.95), mean(N_out >= sum(out1))))

## ---- Table 4: joint exceedance probabilities at mixed levels ----
regs <- list(c(.9, .9), c(.95, .95), c(.99, .9), c(.9, .99), c(.99, .99))
tab4 <- t(sapply(regs, function(rg) { t1 <- -log(1 - rg[1]); t2 <- -log(1 - rg[2])
  v <- sapply(th, function(t) star:::joint_surv_cpp(t1, t2, gl, glw, t$phi, t$alpha, t$b, t$u1, t$u2, t$p, t$chibar, t$z))
  c(lev1 = rg[1], lev2 = rg[2], mean = mean(v), lo = quantile(v, .05, names = FALSE), hi = quantile(v, .95, names = FALSE), empirical = mean(E[, 1] > t1 & E[, 2] > t2), days = 183 * mean(v)) }))
cat("Table 4:\n"); print(signif(tab4, 3))
s <- list(wg = wg, rho = rho, sig = sig, piw = piw, exc = exc, cp = cp, cb = cb, sl = sl, qs = qs, chi_s = chi_s, chi_m = colMeans(chi), chibar_s = chibar_s,
          qe = qe, chi_emp = chi_emp, om = om, adf_s = adf_s, adf_m = adf_m, adf_emp = adf_emp, fWe_m = fWe_m, fWe_s = apply(fWe, 2, q90),
          Ts = Ts, rl_med = rl_med, rl_lo = rl_lo, rl_hi = rl_hi, rho_med = rho_med, to_orig = to_orig, a = a, E = E, fig = x$report$fig, pc = pc, tab4 = tab4, N_out = N_out)
saveRDS(s, "results/beijing_showcase.rds")

## ---- Figure 5 as TikZ -------------------------------------------------------------
out <- c("% Drawing for Figure \\ref{fig:beijing}, included in main.tex with \\input{fig-beijing}")
add <- function(...) out <<- c(out, sprintf(...))
coords <- function(x, y, d = 3) paste(sprintf(paste0("(%.", d, "f,%.", d, "f)"), x, y), collapse = " ")
fig <- s$fig; ps <- fig$pts_scaled; pr <- fig$pts_raw
Qsemi <- function(p, xdat, g) {
  out <- quantile(xdat, p, names = FALSE, type = 6); ab <- p > 1 - g$p_u
  out[ab] <- g$u + g$sigma / g$xi * (((1 - p[ab]) / g$p_u)^(-g$xi) - 1); out
}
bj <- x; gx <- bj$margins$gx; gy <- bj$margins$gy; a <- bj$data; wg <- s$wg
to_orig <- function(rr) cbind(Qsemi(1 - exp(-rr * wg), a$PM2.5, gx), Qsemi(1 - exp(-rr * (1 - wg)), a$O3, gy))
s$to_orig <- to_orig

## Every panel is drawn in a frame of the same physical size, PW x PH cm, with the data area
## starting at (L, B) cm and measuring W x H cm, so that axes and captions align across panels.
PW <- 5.0; PH <- 4.9; L <- 0.95; B <- 1.1; W <- 3.85; H <- 3.5
panel_begin <- function(x0, y0, sx, sy) {
  add("\\begin{tikzpicture}")
  add("\\path[use as bounding box] (0,0) rectangle (%.2f,%.2f);", PW, PH)
  add("\\begin{scope}[shift={(%.4f,%.4f)}]\\begin{scope}[xscale=%.5f,yscale=%.5f]", L - sx * x0, B - sy * y0, sx, sy)
}
panel_end <- function(caption, last = FALSE, newrow = FALSE) {
  add("\\end{scope}\\end{scope}")
  add("\\node[anchor=base] at (%.2f,0.15) {\\scriptsize %s};", PW / 2, caption)
  add(if (newrow) "\\end{tikzpicture}\\\\[4pt]%%" else if (last) "\\end{tikzpicture}" else "\\end{tikzpicture}\\hspace{0.25cm}%%")
}
dot <- function(x, y, col = "gray!65", size = 0.35) for (i in seq_along(x)) add("\\node[circle,fill=%s,inner sep=%.2fpt] at (%.3f,%.3f) {};", col, size, x[i], y[i])

## ---- (a) scaled sample and posterior limit set --------------------------------
## draws at the atom chibar = 0 all have the triangle x1 + x2 = 1 as boundary: drawn once, in red
zero <- x$fit$component[round(seq(1, nrow(D), length.out = length(fig$sig)))] == "zero"
sa <- H / 1.55
panel_begin(0, 0, sa, sa)
add("\\begin{scope}\\clip (-0.03,-0.03) rectangle (1.58,1.58);")
dot(ps[, 1], ps[, 2])
add("\\end{scope}")
add("\\draw[dashed,gray] (1,0) -- (1,1) -- (0,1);")
for (cv in fig$sig[!zero]) add("\\draw[blue!70,opacity=0.3,semithick] plot coordinates {%s};", coords(cv[, 1], cv[, 2]))
if (any(zero)) add("\\draw[red!75!black,thick] (0,1) -- (1,0);")
add("\\draw[->,thick] (0,0) -- (1.55,0) node[below left=-1pt] {\\scriptsize $x_1$};")
add("\\draw[->,thick] (0,0) -- (0,1.55) node[left=-1pt] {\\scriptsize $x_2$};")
add("\\draw (1,%.4f) -- (1,%.4f) node[below] {\\scriptsize 1};", 0.05 / sa, -0.05 / sa)
add("\\draw (%.4f,1) -- (%.4f,1) node[left] {\\scriptsize 1};", 0.05 / sa, -0.05 / sa)
panel_end("(a) scaled sample, posterior limit set")

## ---- (b) data and posterior threshold curve ----------------------------------
sb <- H / 8.3
panel_begin(0, 0, sb, sb)
add("\\begin{scope}\\clip (-0.1,-0.1) rectangle (8.4,8.4);")
dot(pr[, 1], pr[, 2])
for (cv in fig$rho) add("\\draw[red!75!black,opacity=0.18,semithick] plot coordinates {%s};", coords(cv[, 1], cv[, 2]))
add("\\end{scope}")
add("\\draw[->,thick] (0,0) -- (8.3,0) node[below left=-1pt] {\\scriptsize $E_1$};")
add("\\draw[->,thick] (0,0) -- (0,8.3) node[left=-1pt] {\\scriptsize $E_2$};")
add("\\foreach \\t in {2,4,6} { \\draw (\\t,%.4f) -- (\\t,%.4f) node[below] {\\scriptsize \\t}; \\draw (%.4f,\\t) -- (%.4f,\\t) node[left] {\\scriptsize \\t}; }", 0.05 / sb, -0.05 / sb, 0.05 / sb, -0.05 / sb)
panel_end("(b) data, posterior threshold curve")

## ---- (c) posterior of chibar under the three-component prior ------------------
br <- seq(0, 1, by = 0.05); h <- hist(s$sl, breaks = br, plot = FALSE)$counts / length(s$cb) / 0.05
ytop <- max(1, max(h) + 0.05); sxc <- W / 1.2; syc <- H / (ytop + 0.1)
panel_begin(0, 0, sxc, syc)
for (k in seq_along(h)) add("\\fill[blue!30] (%.3f,0) rectangle (%.3f,%.3f);", br[k], br[k + 1], h[k])
xx <- rep(br, each = 2); xx <- xx[-c(1, length(xx))]; add("\\draw[blue!70!black,semithick] plot coordinates {%s};", coords(xx, rep(h, each = 2)))
add("\\draw[gray,dashed] (0,0.5) -- (1,0.5);")
add("\\draw[gray,dashed] (-0.04,0.25) -- (0.08,0.25); \\draw[gray,dashed] (0.92,0.25) -- (1.04,0.25);")
add("\\node[gray,anchor=south west,inner sep=1pt] at (0.62,0.51) {\\scriptsize prior};")
add("\\draw[red!75!black,line width=2.2pt] (0,0) -- (0,%.3f);", s$cp[1])
add("\\draw[red!75!black,line width=2.2pt] (1,0) -- (1,%.3f);", s$cp[3])
add("\\draw[->,thick] (0,0) -- (1.17,0) node[below left=-1pt] {\\scriptsize $\\bar{\\chi}$};")
add("\\draw[->,thick] (0,0) -- (0,%.4f);", ytop + 0.1)
add("\\foreach \\t/\\l in {0/0,0.5/0.5,1/1} \\draw (\\t,%.4f) -- (\\t,%.4f) node[below] {\\scriptsize \\l};", 0.05 / syc, -0.05 / syc)
add("\\foreach \\t in {0.25,0.5,0.75,1} \\draw (%.4f,\\t) -- (%.4f,\\t) node[left] {\\scriptsize \\t};", 0.05 / sxc, -0.05 / sxc)
add("\\node[red!75!black,anchor=north east,inner sep=1pt] at (1.2,%.4f) {\\scriptsize atoms $%.2f$, $%.2f$};", ytop + 0.1, s$cp[1], s$cp[3])
add("\\node[blue!70!black,anchor=north east,inner sep=1pt] at (1.2,%.4f) {\\scriptsize slab $%.2f$};", ytop + 0.1 - 0.36 / syc, s$cp[2])
panel_end("(c) posterior of $\\bar{\\chi}$", newrow = TRUE)

## ---- (d) chi(q): posterior mean, 90% band, empirical, independence ------------
xq <- -log10(1 - s$qs); sxd <- W / 2.55; syd <- H / 0.38
panel_begin(0.7, 0, sxd, syd)
add("\\fill[blue!15] plot coordinates {%s} -- plot coordinates {%s} -- cycle;", coords(xq, s$chi_s[1, ]), coords(rev(xq), rev(s$chi_s[3, ])))
add("\\draw[dashed,gray] plot coordinates {%s};", coords(xq, 1 - s$qs))
add("\\draw[blue!70!black,semithick] plot coordinates {%s};", coords(xq, s$chi_m))
for (k in seq_along(s$qe)) add("\\node[circle,fill,inner sep=0.9pt] at (%.3f,%.3f) {};", -log10(1 - s$qe[k]), s$chi_emp[k])
add("\\draw[->,thick] (0.7,0) -- (3.25,0) node[right=-1pt] {\\scriptsize $q$};")
add("\\draw[->,thick] (0.7,0) -- (0.7,0.38) node[left=-1pt] {\\scriptsize $\\chi(q)$};")
add("\\foreach \\t/\\l in {0.699/.8,1/.9,1.301/.95,2/.99,3/.999} \\draw (\\t,%.4f) -- (\\t,%.4f) node[below] {\\scriptsize \\l};", 0.05 / syd, -0.05 / syd)
add("\\foreach \\t in {0.1,0.2,0.3} \\draw (%.4f,\\t) -- (%.4f,\\t) node[left] {\\scriptsize \\t};", 0.7 + 0.05 / sxd, 0.7 - 0.05 / sxd)
panel_end("(d) $\\chi(q)$, posterior and empirical")

## ---- (e) posterior predictive check: joint exceedance counts ------------------
pc <- s$pc; sq <- function(v) sqrt(v); sxe <- W / 6.8; sye <- H / 8.5
panel_begin(0.3, 0, sxe, sye)
for (k in seq_len(nrow(pc))) {
  add("\\draw[blue!70!black,line width=3pt,opacity=0.35] (%d,%.3f) -- (%d,%.3f);", k, sq(pc[k, "lo"]), k, sq(pc[k, "hi"]))
  add("\\draw[blue!70!black,semithick] (%.2f,%.3f) -- (%.2f,%.3f);", k - 0.22, sq(pc[k, "mean"]), k + 0.22, sq(pc[k, "mean"]))
  add("\\node[circle,fill,inner sep=1pt] at (%d,%.3f) {};", k, sq(pc[k, "obs"]))
}
add("\\draw[->,thick] (0.3,0) -- (7.1,0) node[right=-1pt] {\\scriptsize $q$};")
add("\\draw[->,thick] (0.3,0) -- (0.3,8.5) node[left=-1pt] {\\scriptsize count};")
for (k in seq_len(nrow(pc))) add("\\draw (%d,%.4f) -- (%d,%.4f) node[below] {\\scriptsize %s};", k, 0.05 / sye, k, -0.05 / sye, sub("^0", "", format(pc[k, "q"])))
for (v in c(1, 5, 10, 20, 40, 60)) add("\\draw (%.4f,%.3f) -- (%.4f,%.3f) node[left] {\\scriptsize %d};", 0.3 + 0.05 / sxe, sq(v), 0.3 - 0.05 / sxe, sq(v), v)
panel_end("(e) joint exceedances of level $q$")

## ---- (f) return level sets on the original scale (units of 100 micrograms per cubic metre) ----
sxf <- W / 6.3; syf <- H / 5.3
panel_begin(0, 0, sxf, syf)
o5l <- s$to_orig(s$rl_lo[[2]]) / 100; o5h <- s$to_orig(s$rl_hi[[2]]) / 100
add("\\begin{scope}\\clip (-0.05,-0.05) rectangle (6.3,5.3);")
add("\\fill[blue!15] plot coordinates {%s} -- plot coordinates {%s} -- cycle;", coords(o5l[, 1], o5l[, 2]), coords(rev(o5h[, 1]), rev(o5h[, 2])))
dot(s$a$PM2.5 / 100, s$a$O3 / 100)
ot <- s$to_orig(s$rho_med) / 100
add("\\draw[red!75!black,densely dotted,semithick] plot coordinates {%s};", coords(ot[, 1], ot[, 2]))
o1 <- s$to_orig(s$rl_med[[1]]) / 100; o5 <- s$to_orig(s$rl_med[[2]]) / 100
add("\\draw[blue!70!black,semithick] plot coordinates {%s};", coords(o1[, 1], o1[, 2]))
add("\\draw[blue!70!black,semithick,dashed] plot coordinates {%s};", coords(o5[, 1], o5[, 2]))
add("\\end{scope}")
j <- which.min(abs(s$wg - 0.5))
add("\\node[blue!70!black,anchor=north west,inner sep=1pt] at (%.3f,%.3f) {\\scriptsize 1};", o1[j, 1] + 0.02, o1[j, 2] - 0.02)
add("\\node[blue!70!black,anchor=south west,inner sep=1pt] at (%.3f,%.3f) {\\scriptsize 5};", o5[j, 1] + 0.05, o5[j, 2] + 0.02)
add("\\draw[->,thick] (0,0) -- (6.3,0) node[below left=-1pt] {\\scriptsize PM$_{2.5}$};")
add("\\draw[->,thick] (0,0) -- (0,5.3) node[left=-1pt] {\\scriptsize O$_3$};")
add("\\foreach \\t/\\l in {2/200,4/400} \\draw (\\t,%.4f) -- (\\t,%.4f) node[below] {\\scriptsize \\l};", 0.05 / syf, -0.05 / syf)
add("\\foreach \\t/\\l in {1/100,2/200,3/300,4/400} \\draw (%.4f,\\t) -- (%.4f,\\t) node[left] {\\scriptsize \\l};", 0.05 / sxf, -0.05 / sxf)
panel_end("(f) return level sets, $\\mu$g/m$^3$", last = TRUE)
writeLines(out, "fig-beijing.tex")
cat("lines", length(out), "\n")
## numbers for the text
j <- which.min(abs(s$wg - 0.5)); for (k in 1:3) { o <- s$to_orig(s$rl_med[[k]]); cat(sprintf("T=%d: diag (%.0f,%.0f); PM axis end %.0f (O3 %.0f); O3 axis end %.0f (PM %.0f)\n", s$Ts[k], o[j,1], o[j,2], o[199,1], o[199,2], o[1,2], o[1,1])) }
cat("threshold on original scale: diag", round(100 * ot[j, ]), " ends", round(100 * c(ot[199, 1], ot[1, 2])), "\n")
