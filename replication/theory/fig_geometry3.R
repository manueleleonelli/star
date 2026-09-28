## fig_geometry3.R -- Figure 4 of the main paper: trivariate limit sets of the additive-pairwise scale family.
## Drawn at the printed width (6.5 in), so fonts are at their final size.
source(system.file("trivariate", "star3.R", package = "star"))
g3 <- function(x1, x2, x3, chib, nu) { s <- x1 + x2 + x3; W <- cbind(x1, x2, x3) / s; s * gauge3(W, chib, nu) }
height <- function(x1, x2, chib, nu) {           # largest x3 with g <= 1 over (x1, x2)
  G <- expand.grid(x1 = x1, x2 = x2); lo <- rep(0, nrow(G)); hi <- rep(1.2, nrow(G))
  inside0 <- g3(G$x1, G$x2, 0, chib, nu) <= 1 | (G$x1 + G$x2 == 0)
  for (i in 1:40) { mid <- (lo + hi) / 2; ins <- g3(G$x1, G$x2, mid, chib, nu) <= 1; lo[ins] <- mid[ins]; hi[!ins] <- mid[!ins] }
  z <- lo; z[!inside0] <- NA; z[G$x1 + G$x2 == 0] <- 1
  matrix(z, length(x1), length(x2))
}
## projection of the limit set onto the (j, k) plane. G is star-shaped about the origin, so its projection is the
## union of the projected rays to the boundary points sigma(w) w, and the radius in direction theta is the largest
## projected radius among boundary points whose projection has that direction. Computed on a fine simplex grid
## with the barycentre and the pair's edge included, so that the conical maxima of sigma (Proposition S3) are hit.
SIMPLEX_FINE <- local({ m <- 600; G <- expand.grid(i = 0:m, j = 0:m); G <- G[G$i + G$j <= m, ]; W <- cbind(G$i, G$j, m - G$i - G$j) / m
  rbind(W, c(1, 1, 1) / 3, c(0.5, 0.5, 0), c(0.5, 0, 0.5), c(0, 0.5, 0.5)) })
proj_curve <- function(theta, chib, nu, j = 1, k = 2, face = FALSE) {
  W <- SIMPLEX_FINE; if (face) W <- W[W[, setdiff(1:3, c(j, k))] == 0, , drop = FALSE]
  sig <- 1 / gauge3(W, chib, nu); X <- sig * W; a <- atan2(X[, k], X[, j]); r <- sqrt(X[, j]^2 + X[, k]^2)
  bins <- findInterval(a, c(-Inf, (theta[-1] + theta[-length(theta)]) / 2, Inf))
  rr <- tapply(r, factor(bins, levels = seq_along(theta)), max); rr <- as.numeric(rr)
  ## the boundary of a star-shaped projection: the running envelope also has to respect points assigned to neighbouring bins
  ifelse(is.na(rr), 0, rr)
}
biv_radius <- function(theta, chib, nu = 0.5) { u <- c(cos(theta), sin(theta)); w <- u[1] / sum(u); (1 + chib * tau2(w)^(1 + nu)) / sum(u) }

pdf("fig-geometry3.pdf", width = 5.6, height = 4.3, pointsize = 9)
layout(rbind(1:3, 4:6, c(7, 7, 7)), heights = c(1.05, 1, 0.14)); par(oma = c(0, 0.6, 0, 0))
x <- seq(0, 1, length.out = 51); nu <- c(0.5, 0.5, 0.5)
sets <- list(list(chib = c(0, 0, 0), lab = "(a) independence"),
             list(chib = c(0.8, 0, 0), lab = "(b) one dependent pair"),
             list(chib = c(0.5, 0.5, 0.5), lab = "(c) three pairs at 1/2"))
cube <- rbind(c(0,0,0),c(1,0,0),c(1,1,0),c(0,1,0),c(0,0,0),c(0,0,1),c(1,0,1),c(1,1,1),c(0,1,1),c(0,0,1))
par(mar = c(0.2, 0.2, 1.6, 0.2))
for (s in sets) {
  z <- height(x, x, s$chib, nu)
  pm <- persp(x, x, z, theta = 35, phi = 22, expand = 0.85, col = "grey88", border = "grey45", lwd = 0.25,
              xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), box = FALSE, axes = FALSE, scale = FALSE, r = 4)
  lines(trans3d(cube[, 1], cube[, 2], cube[, 3], pm), lty = 2, col = "grey30", lwd = 0.8)
  for (i in 2:4) { e <- rbind(cube[i, ], cube[i, ] + c(0, 0, 1)); lines(trans3d(e[, 1], e[, 2], e[, 3], pm), lty = 2, col = "grey30", lwd = 0.8) }
  points(trans3d(1, 1, 1, pm), pch = 4, cex = 1.3, lwd = 1.5)
  text(trans3d(0.5, -0.12, 0, pm), expression(x[1]), cex = 1.1)
  text(trans3d(1.12, 0.5, 0, pm), expression(x[2]), cex = 1.1)
  corners <- rbind(c(0,0), c(1,0), c(1,1), c(0,1)); px <- sapply(1:4, function(i) trans3d(corners[i,1], corners[i,2], 0.5, pm)$x)
  cL <- corners[which.min(px), ]; pL <- trans3d(cL[1], cL[2], 0.5, pm)
  text(pL$x, pL$y, expression(x[3]), cex = 1.1, pos = 2, offset = 0.3, xpd = NA)
  title(s$lab, cex.main = 1.05, font.main = 1, line = 0.3)
}
par(mar = c(3.2, 3.4, 1.6, 0.6), mgp = c(2, 0.6, 0), tcl = -0.3)
th <- seq(0.005, pi / 2 - 0.005, length.out = 181)
cases <- list(list(chib = c(0.6, 0.2, 0.2), j = 1, k = 2, lab = "(d) the strong pair of (0.6, 0.2, 0.2)"),
              list(chib = c(0.6, 0.2, 0.2), j = 1, k = 3, lab = "(e) a weak pair of (0.6, 0.2, 0.2)"),
              list(chib = c(0.5, 0.4, 0.4), j = 1, k = 2, lab = "(f) a pair of (0.5, 0.4, 0.4)"))
for (cs in cases) { q <- which(PAIRS[, 1] == cs$j & PAIRS[, 2] == cs$k)
  rp <- proj_curve(th, cs$chib, nu, cs$j, cs$k)
  rf <- proj_curve(th, cs$chib, nu, cs$j, cs$k, face = TRUE)
  rb <- sapply(th, biv_radius, chib = cs$chib[q])
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = bquote(x[.(cs$j)]), ylab = bquote(x[.(cs$k)]), asp = 1,
       xaxt = "n", yaxt = "n", cex.lab = 1.1)
  axis(1, at = c(0, 0.5, 1), cex.axis = 0.95); axis(2, at = c(0, 0.5, 1), cex.axis = 0.95, las = 1)
  title(cs$lab, cex.main = 1.05, font.main = 1, line = 0.3)
  polygon(c(0, 1, 1, 0), c(0, 0, 1, 1), border = "grey60", lty = 2)
  lines(rb * cos(th), rb * sin(th), col = "grey60", lwd = 3.5)
  lines(rf * cos(th), rf * sin(th), lty = 3, lwd = 1.4)
  lines(rp * cos(th), rp * sin(th), lwd = 1.6)
  Wf <- SIMPLEX_FINE; sf <- 1 / gauge3(Wf, cs$chib, nu); e <- max(pmin(sf * Wf[, cs$j], sf * Wf[, cs$k])); eb <- (1 + cs$chib[q]) / 2
  text(0.05, 0.30, bquote(eta == ~ .(sprintf("%.2f", e))), adj = 0, cex = 1)
  text(0.05, 0.16, bquote(eta[bivariate] == ~ .(sprintf("%.2f", eb))), adj = 0, cex = 1, col = "grey40")
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", horiz = TRUE, bty = "n", cex = 1, lwd = c(1.6, 1.4, 3.5), lty = c(1, 3, 1), col = c("black", "black", "grey60"),
       legend = c("projection", "face", "bivariate family"), seg.len = 2.5, x.intersp = 0.6, text.width = c(0.15, 0.1, 0.2))
dev.off()
cat("done\n")
