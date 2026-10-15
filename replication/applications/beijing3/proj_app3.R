## proj_app3.R -- projections of the posterior limit sets of the trivariate Wanliu fit onto the coordinate planes,
## computed on the simplex (see fig_geometry3.R), one curve per posterior draw; and eta_C on the fine grid
source(system.file("trivariate", "star3.R", package = "star")); S <- readRDS("results/app3_summary.rds"); TH <- S$TH
SIMPLEX_FINE <- local({ m <- 400; G <- expand.grid(i = 0:m, j = 0:m); G <- G[G$i + G$j <= m, ]; W <- cbind(G$i, G$j, m - G$i - G$j) / m
  rbind(W, c(1, 1, 1) / 3, c(0.5, 0.5, 0), c(0.5, 0, 0.5), c(0, 0.5, 0.5)) })
thg <- seq(0.01, pi / 2 - 0.01, length.out = 181); cuts <- c(-Inf, (thg[-1] + thg[-length(thg)]) / 2, Inf)
pairs <- list(c(1, 2), c(1, 3), c(2, 3)); idx <- round(seq(1, length(TH), length.out = 100))
PR <- lapply(1:3, function(q) matrix(NA, length(thg), length(idx))); ETAf <- matrix(NA, length(idx), 4)
for (s in seq_along(idx)) { th <- TH[[idx[s]]]; W <- SIMPLEX_FINE; sig <- 1 / gauge3(W, th$chib, th$nu); X <- sig * W
  for (q in 1:3) { j <- pairs[[q]][1]; k <- pairs[[q]][2]; a <- atan2(X[, k], X[, j]); r <- sqrt(X[, j]^2 + X[, k]^2)
    bins <- findInterval(a, cuts); rr <- as.numeric(tapply(r, factor(bins, levels = seq_along(thg)), max)); PR[[q]][, s] <- ifelse(is.na(rr), 0, rr)
    ETAf[s, q] <- max(pmin(X[, j], X[, k])) }
  ETAf[s, 4] <- max(pmin(X[, 1], pmin(X[, 2], X[, 3]))) }
## bivariate fit of Section 6.1 for the PM2.5-ozone pair: posterior median limit set boundary on the same angles (400 draws)
x2 <- readRDS("../beijing/results/fit_Wanliu.rds")$fit; D2 <- x2$draws; zc <- paste0("z", 0:x2$K["K_nu"])
w2 <- cos(thg) / (cos(thg) + sin(thg)); id2 <- round(seq(1, nrow(D2), length.out = 400))
biv <- apply(sapply(id2, function(s) star::sigma_star(w2, D2[s, "chibar"], D2[s, zc])), 1, median) / (cos(thg) + sin(thg))
saveRDS(list(PR = PR, thg = thg, biv = biv, ETAf = ETAf), "results/fig_app3_data.rds")
cat("eta medians on the fine grid (12, 13, 23, 123):", round(apply(ETAf, 2, median), 3), "\n")
cat("eta medians of Table 5 (400 draws):", round(apply(S$ETA, 2, median), 3), "\n")
cat("eta of the median projected curve:", round(sapply(1:3, function(q) { md <- apply(PR[[q]], 1, median); max(pmin(md * cos(thg), md * sin(thg))) }), 3), "\n")
