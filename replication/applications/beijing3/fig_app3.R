## fig_app3.R -- Figure 6 of the main paper: the trivariate Wanliu fit, drawn at the printed width of the Bayesian Analysis template.
source(system.file("trivariate", "star3.R", package = "star")); K <- 4
fit <- readRDS("results/fit_app3.rds"); S <- readRDS("results/app3_summary.rds"); E <- S$E; n <- nrow(E); R <- rowSums(E); W <- E / R
TH <- S$TH; nd <- length(TH)
labs <- c(expression(PM[2.5]), expression(PM[10]), expression(O[3]))
tern <- function(W) cbind(W[, 2] + W[, 3] / 2, W[, 3] * sqrt(3) / 2)
G <- simplex_grid(120); TG <- tern(G$W)
fW <- rowMeans(sapply(TH, function(th) fW_grid(th, K, G))) / 2      # relative to the uniform density (2 in Lebesgue units)
sig <- apply(sapply(TH, function(th) sigma3(G$W, th$chib, th$nu)), 1, median)
rho <- apply(sapply(TH, function(th) rho3(G$W, th$u, th$p)), 1, median)
pm <- list(u = apply(sapply(TH, `[[`, "u"), 1, median), p = median(sapply(TH, `[[`, "p")))
tail <- R > rho3(W, pm$u, pm$p)
pal <- hcl.colors(64, "YlGnBu", rev = TRUE)
tern_panel <- function(z, main, pts = NULL, keylab, digits = 1) {
  plot(NA, xlim = c(-0.1, 1.1), ylim = c(-0.42, 0.95), asp = 1, axes = FALSE, xlab = "", ylab = "")
  title(main, cex.main = 1.05, font.main = 1, line = 0.2)
  br <- seq(min(z), max(z), length.out = 65); col <- pal[cut(z, br, include.lowest = TRUE)]
  points(TG, pch = 15, col = col, cex = 0.42)
  polygon(c(0, 1, 0.5), c(0, 0, sqrt(3) / 2), border = "grey20")
  if (!is.null(pts)) points(tern(pts), pch = 16, cex = 0.3)
  text(0, 0, labs[1], pos = 1, offset = 0.4, cex = 1); text(1, 0, labs[2], pos = 1, offset = 0.4, cex = 1)
  text(0.5, sqrt(3) / 2, labs[3], pos = 3, offset = 0.3, cex = 1, xpd = NA)
  ## horizontal colour key under the triangle
  xs <- seq(0.1, 0.9, length.out = 65); rect(xs[-65], -0.25, xs[-1], -0.19, col = pal, border = NA); rect(0.1, -0.25, 0.9, -0.19, border = "grey40")
  text(c(0.1, 0.9), -0.25, formatC(range(z), format = "f", digits = digits), pos = 1, offset = 0.3, cex = 0.9)
  text(0.5, -0.25, keylab, pos = 1, offset = 0.3, cex = 0.9)
}
dd <- if (file.exists("results/fig_app3_data.rds")) readRDS("results/fig_app3_data.rds") else stop("run the projection step first")
PR <- dd$PR; thg <- dd$thg; biv <- dd$biv
pdf("fig-app3.pdf", width = 5.6, height = 4.45, pointsize = 9)
layout(rbind(1:3, 4:6, c(7, 7, 7)), heights = c(1.05, 1.1, 0.13))
par(mar = c(0.2, 0.3, 1.5, 0.3))
tern_panel(fW, "(a) angular density", pts = W[tail, ], keylab = "ratio to uniform")
tern_panel(sig, expression("(b) tail scale" ~ sigma(bold(w))), keylab = "radius", digits = 2)
tern_panel(rho, expression("(c) threshold" ~ rho(bold(w))), pts = W[tail, ], keylab = "radius")
par(mar = c(3.3, 3.5, 1.5, 0.6), mgp = c(2.1, 0.6, 0), tcl = -0.3)
pairs <- list(c(1, 2), c(1, 3), c(2, 3)); pl <- c("(d)", "(e)", "(f)")
for (q in 1:3) { jk <- pairs[[q]]; X <- E[, jk] / log(n)
  plot(X, pch = 16, cex = 0.3, col = "grey60", xlim = c(0, 1.25), ylim = c(0, 1.25), xlab = labs[jk[1]], ylab = labs[jk[2]],
       asp = 1, xaxt = "n", yaxt = "n", cex.lab = 1.1)
  axis(1, at = c(0, 0.5, 1), cex.axis = 0.95); axis(2, at = c(0, 0.5, 1), cex.axis = 0.95, las = 1)
  title(pl[q], cex.main = 1.05, font.main = 1, line = 0.3)
  polygon(c(0, 1, 1, 0), c(0, 0, 1, 1), border = "grey40", lty = 2)
  lines(c(0, 1), c(1, 0), col = "grey75")
  lo <- apply(PR[[q]], 1, quantile, .05); hi <- apply(PR[[q]], 1, quantile, .95); md <- apply(PR[[q]], 1, median)
  polygon(c(lo * cos(thg), rev(hi * cos(thg))), c(lo * sin(thg), rev(hi * sin(thg))), col = adjustcolor("steelblue", 0.28), border = NA)
  lines(md * cos(thg), md * sin(thg), col = "steelblue4", lwd = 2)
  if (q == 2 && !is.null(biv)) lines(biv * cos(thg), biv * sin(thg), col = "firebrick", lwd = 2, lty = 2)
  e <- median(S$ETA[, paste(jk, collapse = "")])                     # posterior median of eta_C, as in Table 5
  text(0.02, 1.17, bquote(eta == ~ .(sprintf("%.2f", e))), adj = 0, cex = 1, col = "steelblue4")
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", horiz = TRUE, bty = "n", cex = 1, seg.len = 2.5, x.intersp = 0.6,
       lwd = c(2, 8, 2), lty = c(1, 1, 2), col = c("steelblue4", adjustcolor("steelblue", 0.28), "firebrick"),
       legend = c("projected limit set, median", "90% band", "bivariate fit of Section 6.1"), text.width = c(0.25, 0.1, 0.27))
dev.off()
cat("done\n")
