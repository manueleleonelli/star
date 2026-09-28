## app_common.R -- shared machinery for the two applications (source first)
library(star)

## semiparametric margins: empirical below the (1-p_u) quantile, GPD above
gpd_fit <- function(x, u) {
  y <- x[x > u] - u
  nll <- function(p) { s <- exp(p[1]); xi <- p[2]
    if (any(1 + xi*y/s <= 0)) return(1e10)
    length(y)*log(s) + (1 + 1/xi)*sum(log1p(xi*y/s)) }
  o <- optim(c(log(sd(y)), 0.1), nll)
  list(sigma = exp(o$par[1]), xi = o$par[2], u = u, p_u = mean(x > u))
}
to_exp <- function(x, g) {
  F <- ecdf(x)(x) * length(x)/(length(x)+1); ab <- x > g$u
  F[ab] <- 1 - g$p_u * pmax(1 + g$xi*(x[ab]-g$u)/g$sigma, 1e-12)^(-1/g$xi)
  -log(1 - F)
}
standardize <- function(x, y, q_u = 0.9) {
  gx <- gpd_fit(x, quantile(x, q_u)); gy <- gpd_fit(y, quantile(y, q_u))
  list(E = cbind(to_exp(x, gx), to_exp(y, gy)), gx = gx, gy = gy)
}

## everything the paper reports, from one fit
report <- function(fit, E, label, q_chi = c(0.9, 0.95, 0.99), omegas = c(0.2, 0.35, 0.5, 0.65, 0.8),
                   regions = list(c(.9,.9), c(.95,.95), c(.99,.9), c(.9,.99), c(.99,.99)), days = NA) {
  D <- fit$draws; r <- rowSums(E); w <- E[,1]/r
  idx <- round(seq(1, nrow(D), length.out = 300)); q <- function(x) quantile(x, c(.05,.5,.95))
  gl <- seq(0.5/201, 1-0.5/201, length.out = 201); glw <- rep(1/201, 201)
  cat("\n=====", label, "  n =", nrow(E), "=====\n")
  cp <- class_probs(fit); cat("class probabilities:", paste(names(cp), round(cp,3)), "\n")
  cat("chibar (5/50/95%):", round(q(D[,"chibar"]),3), " | eta:", round(q((1+D[,"chibar"])/2),3), "\n")
  sl <- D[,"chibar"]; sl <- sl[sl > 0 & sl < 1]; if (length(sl) > 10) cat("slab-conditional chibar:", round(q(sl),3), "\n")
  cat("p:", round(q(D[,"p"]),2), " | u1:", round(q(D[,"u1"]),2), " | u2:", round(q(D[,"u2"]),2), " | alpha:", round(q(D[,"alpha"]),2), "\n")
  ex <- sapply(idx, function(s){ th <- star:::theta_row(fit,s); mean(r > rho_star(w,th$u1,th$u2,th$p)) })
  rh <- sapply(idx, function(s){ th <- star:::theta_row(fit,s); rho_star(0.5,th$u1,th$u2,th$p) })
  cat("exceedance fraction:", round(q(ex),3), " | rho(1/2):", round(q(rh),2), "\n")
  cq <- chi_q(fit, q_chi, ndraw = 200)
  emp <- sapply(q_chi, function(qq) mean(E[,1] > -log(1-qq) & E[,2] > -log(1-qq))/(1-qq))
  cat("chi(q) model mean:", round(cq[1,],3), "\n   90% lo:", round(cq[2,],3), " hi:", round(cq[3,],3), "\n   empirical:", round(emp,3), "\n")
  ad <- adf_lambda(fit, omega = omegas, ndraw = 200)
  cat("Lambda(omega) mean:", round(ad$summary[1,],3), " lo:", round(ad$summary[2,],3), " hi:", round(ad$summary[3,],3), "\n")
  cat("joint exceedance probabilities (levels are marginal quantiles):\n")
  tab <- t(sapply(regions, function(rg) { t1 <- -log(1-rg[1]); t2 <- -log(1-rg[2])
    v <- sapply(idx, function(s){ th <- star:::theta_row(fit,s)
      star:::joint_surv_cpp(t1,t2,gl,glw,th$phi,th$alpha,th$b,th$u1,th$u2,th$p,th$chibar,th$z) })
    c(lev1 = rg[1], lev2 = rg[2], mean = mean(v), lo = quantile(v,.05), hi = quantile(v,.95),
      empirical = mean(E[,1] > t1 & E[,2] > t2), days = if (is.na(days)) NA else days*mean(v)) }))
  print(round(tab, 4))
  ## figure coordinates: scaled cloud + 40 posterior sigma curves, raw data + 40 rho curves
  ln <- log(nrow(E)); wg <- seq(0.003, 0.997, length.out = 61); sub <- round(seq(1, nrow(D), length.out = 40))
  fig <- list(pts_scaled = round(E/ln, 3), pts_raw = round(E, 3),
    sig = lapply(sub, function(s){ th <- star:::theta_row(fit,s); sg <- sigma_star(wg,th$chibar,th$z); round(cbind(sg*wg, sg*(1-wg)),3) }),
    rho = lapply(sub, function(s){ th <- star:::theta_row(fit,s); rr <- rho_star(wg,th$u1,th$u2,th$p); round(cbind(rr*wg, rr*(1-wg)),3) }))
  invisible(list(cp = cp, chibar = q(D[,"chibar"]), p = q(D[,"p"]), exc = q(ex), rho_half = q(rh),
                 chi = cq, chi_emp = emp, adf = ad, regions = tab, fig = fig, acc = fit$acc_rates))
}
