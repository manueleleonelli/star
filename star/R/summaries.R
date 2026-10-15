# Derived posterior quantities ----------------------------------------------

theta_row <- function(fit, s) {
  d <- fit$draws[s, ]; K <- fit$K
  list(phi = unname(d[1:K["K_W"]]),
       alpha = unname(d["alpha"]),
       b = unname(d[paste0("b", 0:K["K_lam"])]),
       u1 = unname(d["u1"]), u2 = unname(d["u2"]), p = unname(d["p"]),
       chibar = unname(d["chibar"]),
       z = unname(d[paste0("z", 0:K["K_nu"])]))
}

#' Posterior draws of the residual tail dependence coefficient eta = (1+chibar)/2
#' @export
eta_draws <- function(fit) (1 + fit$draws[, "chibar"]) / 2

#' Posterior probabilities of the dependence classes (three-component prior)
#' @export
class_probs <- function(fit) {
  if (fit$prior$route != "threecomp")
    stop("class probabilities require the three-component prior")
  c(P_AI_exact = mean(fit$component == "zero"),
    P_slab = mean(fit$component == "slab"),
    P_AD = mean(fit$component == "one"))
}

gl_nodes <- function(M = 101) {          # Gauss-Legendre on (0,1)
  g <- seq(0.5 / M, 1 - 0.5 / M, length.out = M)   # midpoint rule, robust and simple
  list(nodes = g, wts = rep(1 / M, M))
}

#' Posterior of the sub-asymptotic dependence curve chi(q)
#' @export
chi_q <- function(fit, q, ndraw = 200, M = 201) {
  gl <- gl_nodes(M); tq <- -log(1 - q)
  idx <- round(seq(1, nrow(fit$draws), length.out = min(ndraw, nrow(fit$draws))))
  sapply(q, function(qq) {
    tq <- -log(1 - qq)
    vals <- sapply(idx, function(s) {
      th <- theta_row(fit, s)
      joint_surv_cpp(tq, tq, gl$nodes, gl$wts, th$phi, th$alpha, th$b,
                     th$u1, th$u2, th$p, th$chibar, th$z) / (1 - qq)
    })
    c(mean = mean(vals), lo = quantile(vals, 0.025), hi = quantile(vals, 0.975))
  })
}

#' Posterior of the angular dependence function Lambda(omega) of Proposition 4
#' @export
adf_lambda <- function(fit, omega = seq(0.05, 0.95, by = 0.05),
                       ndraw = 200, grid = seq(0.001, 0.999, length.out = 999)) {
  idx <- round(seq(1, nrow(fit$draws), length.out = min(ndraw, nrow(fit$draws))))
  out <- sapply(omega, function(om) {
    vals <- sapply(idx, function(s) {
      th <- theta_row(fit, s)
      sig <- sigma_star(grid, th$chibar, th$z)
      min(pmax(om / grid, (1 - om) / (1 - grid)) / sig)
    })
    c(mean = mean(vals), lo = quantile(vals, 0.025), hi = quantile(vals, 0.975))
  })
  list(omega = omega, summary = out)
}

#' Posterior mean joint survival probability P(E1 > t1, E2 > t2)
#' @export
joint_surv <- function(fit, t1, t2, ndraw = 200, M = 201) {
  gl <- gl_nodes(M)
  idx <- round(seq(1, nrow(fit$draws), length.out = min(ndraw, nrow(fit$draws))))
  vals <- sapply(idx, function(s) {
    th <- theta_row(fit, s)
    joint_surv_cpp(t1, t2, gl$nodes, gl$wts, th$phi, th$alpha, th$b,
                   th$u1, th$u2, th$p, th$chibar, th$z)
  })
  c(mean = mean(vals), lo = quantile(vals, 0.025), hi = quantile(vals, 0.975))
}

# Plots ----------------------------------------------------------------------

#' Scaled sample with the posterior of the limit set boundary (check of fit of Section 4.5)
#' @export
plot_silhouette <- function(fit, ndraw = 100, wgrid = seq(0.002, 0.998, length.out = 199)) {
  n <- length(fit$r); ln <- log(n)
  plot(fit$r * fit$w / ln, fit$r * (1 - fit$w) / ln, pch = 16, cex = 0.3,
       col = "grey60", xlab = expression(x[1]), ylab = expression(x[2]),
       xlim = c(0, 1.3), ylim = c(0, 1.3), asp = 1)
  idx <- round(seq(1, nrow(fit$draws), length.out = min(ndraw, nrow(fit$draws))))
  for (s in idx) {
    th <- theta_row(fit, s)
    sig <- sigma_star(wgrid, th$chibar, th$z)
    lines(sig * wgrid, sig * (1 - wgrid), col = rgb(0, 0, 1, 0.06))
  }
  abline(h = 1, v = 1, lty = 2, col = "grey40")
}

#' Data with the posterior of the threshold curve
#' @export
plot_threshold <- function(fit, ndraw = 100, wgrid = seq(0.002, 0.998, length.out = 199)) {
  plot(fit$r * fit$w, fit$r * (1 - fit$w), pch = 16, cex = 0.3, col = "grey60",
       xlab = expression(x[1]), ylab = expression(x[2]))
  idx <- round(seq(1, nrow(fit$draws), length.out = min(ndraw, nrow(fit$draws))))
  for (s in idx) {
    th <- theta_row(fit, s)
    rho <- rho_star(wgrid, th$u1, th$u2, th$p)
    lines(rho * wgrid, rho * (1 - wgrid), col = rgb(1, 0, 0, 0.06))
  }
}
