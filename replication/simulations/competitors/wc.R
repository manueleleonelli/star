## Wadsworth-Campbell (2024) geometric approach, through the authors' package geometricMVE:
## threshold r0(w) = rolling-window 0.95 quantile of R | W (QR.2d, method "empirical"),
## truncated-Gamma radial model above it with a parametric gauge, chosen by AIC among the
## logistic-type, Gaussian, inverted-logistic and square gauges and the additive mixture of
## the first two (the example script of the package); eta from the fitted gauge, chi(q) by
## simulation from the fitted exceedance model (sim.2d).
suppressMessages(library(geometricMVE))
wc_fit <- function(E, tau = 0.95) {
  r <- rowSums(E); w <- E[, 1] / r
  qr <- QR.2d(r = r, w = w, method = "empirical", tau = tau)
  exc <- which(r > qr$r0w & !is.na(qr$r0w))
  rexc <- r[exc]; wexc <- w[exc]; r0w <- qr$r0w[exc]
  cand <- list(
    rvad    = list(gfun = gauge_rvad,        init = c(2, 0.5), lo = c(0, 0.01), hi = c(Inf, 1)),
    gauss   = list(gfun = gauge_gaussian,    init = c(2, 0.5), lo = c(0, 0), hi = c(Inf, 0.999)),
    invlog  = list(gfun = gauge_invlogistic, init = c(2, 2),   lo = c(0, 1), hi = c(Inf, Inf)),
    square  = list(gfun = gauge_square,      init = c(2, 0.5), lo = c(0, 0.01), hi = c(Inf, 0.99)))
  fits <- lapply(cand, function(cd) {
    f <- try(fit.geometricMVE.2d(r = rexc, w = wexc, r0w = r0w, init.val = cd$init, gfun = cd$gfun,
                                 lower.limit = cd$lo, upper.limit = cd$hi), silent = TRUE)
    if (inherits(f, "try-error") || f$convergence != 0) return(NULL)
    list(mle = f$mle, aic = 2 * f$nllh + 2 * length(f$mle), gfun = cd$gfun, add = FALSE)
  })
  fm <- try(fit.geometricMVE.2d(r = rexc, w = wexc, r0w = r0w, init.val = c(2, 0.5, 0.5, 1),
                                add.gauge = TRUE, gauge1 = gauge_rvad, gauge2 = gauge_gaussian,
                                par1ind = 2, par2ind = 3, lower.limit = c(0, 0.01, 0, 0), upper.limit = c(Inf, 1, 0.999, Inf)),
            silent = TRUE)
  fits$mix <- if (inherits(fm, "try-error") || fm$convergence != 0) NULL else
    list(mle = fm$mle, aic = 2 * fm$nllh + 2 * length(fm$mle), add = TRUE)
  fits <- fits[!sapply(fits, is.null)]
  best <- names(fits)[which.min(sapply(fits, `[[`, "aic"))]
  list(fit = fits[[best]], which = best, wexc = wexc, r0w = r0w, pexc = length(exc) / nrow(E))
}
## gauge value at (w, 1 - w) for the selected fit
wc_gauge <- function(f, wg) {
  if (!f$add) return(sapply(wg, function(w) f$gfun(c(w, 1 - w), f$mle[-1])))
  ms <- additivegauge.scaling(gauge1 = gauge_rvad, gauge2 = gauge_gaussian, weight = f$mle[4], par1 = f$mle[2], par2 = f$mle[3])
  sapply(wg, function(w) additivegauge.rescale(c(w, 1 - w), gauge1 = gauge_rvad, gauge2 = gauge_gaussian,
                                               par1 = f$mle[2], par2 = f$mle[3], weight = f$mle[4], ms = ms))
}
wc_eta <- function(wf) { wg <- seq(0.001, 0.999, length.out = 999); max(pmin(wg, 1 - wg) / wc_gauge(wf$fit, wg)) }
## chi(q) = P(E1 > t_q, E2 > t_q) / (1 - q): exceedance fraction times the conditional probability by simulation
wc_chi <- function(wf, qs, nsim = 2e5) {
  f <- wf$fit
  xs <- if (f$add) sim.2d(w = wf$wexc, r0w = wf$r0w, nsim = nsim, par = f$mle, add.gauge = TRUE,
                          gauge1 = gauge_rvad, gauge2 = gauge_gaussian, par1ind = 2, par2ind = 3)
        else sim.2d(w = wf$wexc, r0w = wf$r0w, nsim = nsim, par = f$mle, gfun = f$gfun)
  sapply(qs, function(q) { tq <- -log(1 - q); wf$pexc * mean(xs[, 1] > tq & xs[, 2] > tq) / (1 - q) })
}
