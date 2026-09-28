## sec6_cut.R -- return level sets of Section 6.1 under the cut posterior (margins drawn from their posterior, theta
## refitted on each standardized sample; 60 x 40 draws from results/wanliu_cut60_fits.rds, see wanliu_cut60.R), the
## numbers for the text, and panel (f) of Figure 5, which replaces the last tikzpicture of fig-beijing.tex (written by
## app1_showcase.R) in fig-beijing-new.tex. ba_figure.R then gives the Bayesian Analysis layout.
suppressMessages(library(star))
x <- readRDS("results/fit_Wanliu.rds"); a <- x$data; Y <- cbind(a$PM2.5, a$O3); n <- nrow(Y); gx <- x$margins$gx; gy <- x$margins$gy
source("joint.R")                                                    # std_setup, to_exp, gpd_ll, marg_prior (unused here)
S <- list(std_setup(Y[, 1]), std_setup(Y[, 2]))
P <- readRDS("results/wanliu_cut60_fits.rds")
wg <- seq(0.005, 0.995, length.out = 199); jd <- which.min(abs(wg - 0.5)); dw <- wg[2] - wg[1]
Qsemi <- function(p, xdat, u, p_u, sigma, xi) { out <- quantile(xdat, pmin(p, 1 - p_u), names = FALSE, type = 6); ab <- p > 1 - p_u
  out[ab] <- u + sigma / xi * (((1 - p[ab]) / p_u)^(-xi) - 1); out }
Tv <- function(y, j, ps) { Sj <- S[[j]]; ls <- ps[2 * j - 1]; xi <- ps[2 * j]; out <- -log(1 - ecdf(Sj$y)(y) * length(Sj$y) / (length(Sj$y) + 1))
  ab <- y > Sj$u; z <- 1 + xi * (y[ab] - Sj$u) / exp(ls); out[ab] <- ifelse(z > 0, -log(Sj$p_u) + log(pmax(z, 1e-300)) / xi, 60); out }
## draws: theta and psi
D <- do.call(c, lapply(P$fits, function(f) lapply(seq_len(nrow(f$draws)), function(s) list(th = star:::theta_row(f, s), psi = unname(f$m), comp = f$component[s]))))
comp <- sapply(D, `[[`, "comp"); cat(sprintf("cut posterior class probabilities: zero %.3f slab %.3f one %.3f (draws %d)\n", mean(comp == "zero"), mean(comp == "slab"), mean(comp == "one"), length(D)))
gl <- seq(0.5 / 201, 1 - 0.5 / 201, length.out = 201); glw <- rep(1 / 201, 201)
chi <- sapply(D, function(d) { t <- d$th; sapply(c(0.9, 0.95, 0.99), function(q) { tt <- -log(1 - q); star:::joint_surv_cpp(tt, tt, gl, glw, t$phi, t$alpha, t$b, t$u1, t$u2, t$p, t$chibar, t$z) / (1 - q) }) })
cat("cut posterior chi(0.9, 0.95, 0.99) means:", round(rowMeans(chi), 3), " 5%:", round(apply(chi, 1, quantile, .05), 3), " 95%:", round(apply(chi, 1, quantile, .95), 3), "\n")
cat("P(xi_O3 < 0) over the marginal draws:", mean(P$md2[, 2] < 0), "\n")
rT <- function(t, w, pr) { lam <- lambda_star(w, t$b); rho <- rho_star(w, t$u1, t$u2, t$p); sg <- sigma_star(w, t$chibar, t$z); piw <- pgamma(rho, shape = t$alpha, rate = lam)
  rb <- qgamma(1 - pr, shape = t$alpha, rate = lam); tl <- (1 - piw) > pr; rb[tl] <- rho[tl] + sg[tl] * log((1 - piw[tl]) / pr); rb }
to_o <- function(rr, ps) cbind(Qsemi(1 - exp(-rr * wg), Y[, 1], S[[1]]$u, S[[1]]$p_u, exp(ps[1]), ps[2]), Qsemi(1 - exp(-rr * (1 - wg)), Y[, 2], S[[2]]$u, S[[2]]$p_u, exp(ps[3]), ps[4]))
C1 <- lapply(D, function(d) to_o(rT(d$th, wg, 1 / 183), d$psi)); C5 <- lapply(D, function(d) to_o(rT(d$th, wg, 1 / (5 * 183)), d$psi))
Cth <- lapply(D, function(d) to_o(rho_star(wg, d$th$u1, d$th$u2, d$th$p), d$psi))
## the diagonal point of each draw's set (w = 1/2 on that draw's standardized scale): posterior quantiles for the text
dq <- function(L, p) c(quantile(sapply(L, function(m) m[jd, 1]), p), quantile(sapply(L, function(m) m[jd, 2]), p))
## curves and bands for the figure: pointwise quantiles of the distance along rays from the origin, concentrations scaled by their 0.9 quantiles
sc <- c(Qsemi(0.9, Y[, 1], S[[1]]$u, S[[1]]$p_u, gx$sigma, gx$xi), Qsemi(0.9, Y[, 2], S[[2]]$u, S[[2]]$p_u, gy$sigma, gy$xi)); ang <- seq(0.01, pi / 2 - 0.01, length.out = 400)
ray <- function(C) { z1 <- C[, 1] / sc[1]; z2 <- C[, 2] / sc[2]; a <- atan2(z2, z1); r <- sqrt(z1^2 + z2^2); o <- order(a); approx(a[o], r[o], ang, rule = 2, ties = mean)$y }
cq <- function(L, p) { R <- sapply(L, ray); rr <- apply(R, 1, quantile, p); cbind(sc[1] * rr * cos(ang), sc[2] * rr * sin(ang)) }
m1 <- cq(C1, .5); l1 <- cq(C1, .05); h1 <- cq(C1, .95); m5 <- cq(C5, .5); l5 <- cq(C5, .05); h5 <- cq(C5, .95); mth <- cq(Cth, .5)
d1 <- rbind(dq(C1, .5), dq(C1, .05), dq(C1, .95)); d5 <- rbind(dq(C5, .5), dq(C5, .05), dq(C5, .95)); dth <- dq(Cth, .5)
cat(sprintf("per-draw diagonal points: threshold (%.0f, %.0f); one-season (%.0f, %.0f), 90%% (%.0f, %.0f) to (%.0f, %.0f); five-season (%.0f, %.0f), 90%% (%.0f, %.0f) to (%.0f, %.0f)\n", dth[1], dth[2], d1[1,1], d1[1,2], d1[2,1], d1[2,2], d1[3,1], d1[3,2], d5[1,1], d5[1,2], d5[2,1], d5[2,2], d5[3,1], d5[3,2]))
jd <- which.min(abs(atan2(d1[1, 2] / sc[2], d1[1, 1] / sc[1]) - ang))   # the ray through the median diagonal point, for labels and checks

cat(sprintf("threshold curve on the diagonal: (%.0f, %.0f)\n", mth[jd, 1], mth[jd, 2]))
cat(sprintf("one-season set, diagonal: (%.0f, %.0f), 90%% from (%.0f, %.0f) to (%.0f, %.0f)\n", m1[jd, 1], m1[jd, 2], l1[jd, 1], l1[jd, 2], h1[jd, 1], h1[jd, 2]))
cat(sprintf("five-season set, diagonal: (%.0f, %.0f), 90%% from (%.0f, %.0f) to (%.0f, %.0f)\n", m5[jd, 1], m5[jd, 2], l5[jd, 1], l5[jd, 2], h5[jd, 1], h5[jd, 2]))
cat(sprintf("median five-season set: largest PM2.5 %.0f, largest O3 %.0f; upper band: largest O3 %.0f\n", max(m5[, 1]), max(m5[, 2]), max(h5[, 2])))
cat(sprintf("median one-season set: largest PM2.5 %.0f, largest O3 %.0f\n", max(m1[, 1]), max(m1[, 2])))
## days outside the median one-season set: curve and data mapped to the plug-in standardized scale, compared in polar form
psh <- c(log(gx$sigma), gx$xi, log(gy$sigma), gy$xi)
zd1 <- Y[, 1] / sc[1]; zd2 <- Y[, 2] / sc[2]; ad <- atan2(zd2, zd1); rdd <- sqrt(zd1^2 + zd2^2); rm1 <- sqrt((m1[, 1] / sc[1])^2 + (m1[, 2] / sc[2])^2)
out <- rdd > approx(ang, rm1, ad, rule = 2)$y
Ed <- cbind(Tv(Y[, 1], 1, psh), Tv(Y[, 2], 2, psh)); rd <- rowSums(Ed); wd <- Ed[, 1] / rd
cat(sprintf("days outside the median one-season set: %d (w < 0.25: %d, w > 0.75: %d)\n", sum(out), sum(out & wd < 0.25), sum(out & wd > 0.75)))
## predictive: for each draw, the probability of a new day outside the median set, on that draw's standardized scale
p_out <- sapply(D, function(d) { t <- d$th; E <- cbind(Tv(m1[, 1], 1, d$psi), Tv(m1[, 2], 2, d$psi)); rk <- rowSums(E); wk <- E[, 1] / rk; ok <- order(wk)
  rr <- approx(wk[ok], rk[ok], wg, rule = 2)$y; lam <- lambda_star(wg, t$b); rho <- rho_star(wg, t$u1, t$u2, t$p); sg <- sigma_star(wg, t$chibar, t$z); piw <- pgamma(rho, shape = t$alpha, rate = lam)
  sv <- ifelse(rr > rho, (1 - piw) * exp(-(rr - rho) / sg), 1 - pgamma(rr, shape = t$alpha, rate = lam))
  K <- length(t$phi); fW <- rowSums(sapply(1:K, function(k) t$phi[k] * dbeta(wg, k, K - k + 1))); sum(fW * sv) / sum(fW) })
set.seed(3); N <- rbinom(length(p_out), n, p_out)
cat(sprintf("predictive: mean %.1f, 90%% (%.0f, %.0f)\n", n * mean(p_out), quantile(N, .05), quantile(N, .95)))
## marginal posterior of the shapes from long chains, priors as in the cut posterior
lpri <- function(lsig, xi) dnorm(lsig, 0, 10, log = TRUE) + dnorm(xi, 0, 0.5, log = TRUE)
mlong <- function(j, iter = 200000, burn = 10000) { g <- gpd_fit(Y[, j], S[[j]]$u); cur <- c(g$lsig, g$xi); lc <- gpd_ll(S[[j]], cur[1], cur[2]) + lpri(cur[1], cur[2]); out <- numeric(iter)
  for (it in 1:iter) { pr <- cur + 1.2 * g$se * rnorm(2); lp <- gpd_ll(S[[j]], pr[1], pr[2]) + lpri(pr[1], pr[2]); if (is.finite(lp) && log(runif(1)) < lp - lc) { cur <- pr; lc <- lp }; out[it] <- cur[2] }
  out[(burn + 1):iter] }
set.seed(11); x1 <- mlong(1); x2 <- mlong(2)
cat(sprintf("shape posteriors: PM2.5 %.2f (%.2f, %.2f); O3 %.2f (%.2f, %.2f), P(xi_O3 < 0) = %.2f\n", median(x1), quantile(x1, .05), quantile(x1, .95), median(x2), quantile(x2, .05), quantile(x2, .95), mean(x2 < 0)))
## ---- panel (f) ----
fb <- readLines("fig-beijing.tex"); k0 <- max(grep("^\\\\begin\\{tikzpicture\\}", fb))
out_tex <- fb[1:(k0 - 1)]; add <- function(...) out_tex <<- c(out_tex, sprintf(...))
coords <- function(x, y, d = 3) paste(sprintf(paste0("(%.", d, "f,%.", d, "f)"), x, y), collapse = " ")
PW <- 5.0; PH <- 4.9; L <- 0.95; B <- 1.1; W <- 3.85; H <- 3.5
panel_begin <- function(x0, y0, sx, sy) { add("\\begin{tikzpicture}"); add("\\path[use as bounding box] (0,0) rectangle (%.2f,%.2f);", PW, PH)
  add("\\begin{scope}[shift={(%.4f,%.4f)}]\\begin{scope}[xscale=%.5f,yscale=%.5f]", L - sx * x0, B - sy * y0, sx, sy) }
dot <- function(x, y, col = "gray!65", size = 0.35) for (i in seq_along(x)) add("\\node[circle,fill=%s,inner sep=%.2fpt] at (%.3f,%.3f) {};", col, size, x[i], y[i])
sxf <- W / 6.3; syf <- H / 5.3
panel_begin(0, 0, sxf, syf)
add("\\begin{scope}\\clip (-0.05,-0.05) rectangle (6.3,5.3);")
add("\\fill[blue!15] plot coordinates {%s} -- plot coordinates {%s} -- cycle;", coords(l5[, 1] / 100, l5[, 2] / 100), coords(rev(h5[, 1]) / 100, rev(h5[, 2]) / 100))
dot(a$PM2.5 / 100, a$O3 / 100)
add("\\draw[red!75!black,densely dotted,semithick] plot coordinates {%s};", coords(mth[, 1] / 100, mth[, 2] / 100))
add("\\draw[blue!70!black,semithick] plot coordinates {%s};", coords(m1[, 1] / 100, m1[, 2] / 100))
add("\\draw[blue!70!black,semithick,dashed] plot coordinates {%s};", coords(m5[, 1] / 100, m5[, 2] / 100))
add("\\end{scope}")
add("\\node[blue!70!black,anchor=north west,inner sep=1pt] at (%.3f,%.3f) {\\scriptsize 1};", m1[jd, 1] / 100 + 0.02, m1[jd, 2] / 100 - 0.02)
add("\\node[blue!70!black,anchor=south west,inner sep=1pt] at (%.3f,%.3f) {\\scriptsize 5};", m5[jd, 1] / 100 + 0.05, m5[jd, 2] / 100 + 0.02)
add("\\draw[->,thick] (0,0) -- (6.3,0) node[below left=-1pt] {\\scriptsize PM$_{2.5}$};")
add("\\draw[->,thick] (0,0) -- (0,5.3) node[left=-1pt] {\\scriptsize O$_3$};")
add("\\foreach \\t/\\l in {2/200,4/400} \\draw (\\t,%.4f) -- (\\t,%.4f) node[below] {\\scriptsize \\l};", 0.05 / syf, -0.05 / syf)
add("\\foreach \\t/\\l in {1/100,2/200,3/300,4/400} \\draw (%.4f,\\t) -- (%.4f,\\t) node[left] {\\scriptsize \\l};", 0.05 / sxf, -0.05 / sxf)
add("\\end{scope}\\end{scope}")
add("\\node[anchor=base] at (%.2f,0.15) {\\scriptsize (f) return level sets, $\\mu$g/m$^3$};", PW / 2)
add("\\end{tikzpicture}")
writeLines(out_tex, "fig-beijing-new.tex"); cat("wrote fig-beijing-new.tex,", length(out_tex), "lines\n")
