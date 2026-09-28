## joint.R -- cut (two-step) versus joint posterior over margins and dependence (Section 4 of the main paper).
## Margins: empirical below the sample 0.9 quantile, GPD(sigma_j, xi_j) above. The cut fixes (sigma_j, xi_j) at
## their maximum likelihood estimates and runs star_mcmc on the standardized data; the joint arm adds one
## random-walk block per margin to the same sampler, with the Jacobian of the transformation in the target.
suppressMessages(library(star))
q_u <- 0.9
gpd_fit <- function(x, u) { y <- x[x > u] - u
  nll <- function(p) { s <- exp(p[1]); xi <- p[2]; if (any(1 + xi * y / s <= 0)) return(1e10); length(y) * log(s) + (1 + 1/xi) * sum(log1p(xi * y / s)) }
  o <- optim(c(log(sd(y)), 0.1), nll, hessian = TRUE)
  list(lsig = o$par[1], xi = o$par[2], u = u, p_u = mean(x > u), se = sqrt(diag(solve(o$hessian)))) }
## standardization to exponential margins for given (log sigma, xi); the bulk part is fixed
std_setup <- function(y) { u <- quantile(y, q_u, names = FALSE); F0 <- ecdf(y)(y) * length(y) / (length(y) + 1)
  list(y = y, u = u, ab = y > u, e0 = -log(1 - F0), p_u = mean(y > u), x = pmax(y - u, 0)) }
to_exp <- function(S, lsig, xi) { e <- S$e0; s <- exp(lsig); z <- 1 + xi * S$x[S$ab] / s
  e[S$ab] <- -log(S$p_u) + log(z) / xi; e }
## log of the marginal factor of the joint density over the exceedances: sum log h - log S = -n_u log sigma - sum log(1 + xi x / sigma)
marg_term <- function(S, lsig, xi) { s <- exp(lsig); z <- 1 + xi * S$x[S$ab] / s; if (any(z <= 0) || xi < -0.5) return(-Inf)
  -sum(S$ab) * lsig - sum(log(z)) }
marg_prior <- function(lsig, xi) dnorm(lsig, 0, 3, log = TRUE) + dnorm(xi, 0, 0.5, log = TRUE)

joint_mcmc <- function(Y, prior, K_W = 5, K_lam = 2, K_nu = 2, n_iter = 20000, burnin = 8000, thin = 5, beta_move = c(1, 1, 1) / 3, verbose = FALSE) {
  S1 <- std_setup(Y[, 1]); S2 <- std_setup(Y[, 2]); g1 <- gpd_fit(Y[, 1], S1$u); g2 <- gpd_fit(Y[, 2], S2$u)
  m1 <- c(g1$lsig, g1$xi); m2 <- c(g2$lsig, g2$xi)
  E <- cbind(to_exp(S1, m1[1], m1[2]), to_exp(S2, m2[1], m2[2])); r <- rowSums(E); w <- E[, 1] / r
  q0 <- prior$q_kappa0
  state <- list(psi = rep(0, K_W - 1), la = log(2), b = rep(0, K_lam + 1), lo1 = log(max(quantile(r, 0.9) / 2 - q0, 0.2)),
                lo2 = log(max(quantile(r, 0.9) / 2 - q0, 0.2)), lp = log(2), lchi = qlogis(0.5), z = rep(0, K_nu + 1), m = "slab", q0 = q0)
  s2t <- star:::state_to_theta; lpr <- star:::log_prior
  ll_cur <- loglik_star(r, w, s2t(state)); lp_cur <- lpr(state, prior); lm_cur <- marg_term(S1, m1[1], m1[2]) + marg_term(S2, m2[1], m2[2]) + marg_prior(m1[1], m1[2]) + marg_prior(m2[1], m2[2])
  blocks <- list(ang = "psi", body = c("la", "b"), thr = c("lo1", "lo2", "lp"), tail = c("lchi", "z"))
  ls <- c(ang = -1, body = -1, thr = -2, tail = -1, m1 = log(0.3), m2 = log(0.3)); acc <- num <- c(ang = 0, body = 0, thr = 0, tail = 0, jump = 0, m1 = 0, m2 = 0)
  keep <- seq(burnin + 1, n_iter, by = thin); P <- K_W + 1 + (K_lam + 1) + 3 + 1 + (K_nu + 1)
  draws <- matrix(NA_real_, length(keep), P); mdraws <- matrix(NA_real_, length(keep), 4); comp <- character(length(keep)); s <- 0
  colnames(draws) <- c(paste0("phi", 1:K_W), "alpha", paste0("b", 0:K_lam), "u1", "u2", "p", "chibar", paste0("z", 0:K_nu)); colnames(mdraws) <- c("lsig1", "xi1", "lsig2", "xi2")
  perturb <- function(state, nm, sc) { for (v in nm) state[[v]] <- state[[v]] + rnorm(length(state[[v]]), 0, sc); state }
  Lse <- list(g1$se, g2$se)
  for (t in seq_len(n_iter)) {
    gam <- if (t <= burnin) t^(-0.6) else 0
    for (bl in names(blocks)) { nm <- blocks[[bl]]
      if (bl == "tail" && state$m != "slab") nm <- "z"
      if (bl == "tail" && state$m == "zero") { state$z <- rnorm(K_nu + 1, 0, prior$v_nu); lp_cur <- lpr(state, prior); next }
      prop <- perturb(state, nm, exp(ls[bl])); ll_p <- loglik_star(r, w, s2t(prop)); lp_p <- lpr(prop, prior)
      a <- min(1, exp(ll_p + lp_p - ll_cur - lp_cur))
      if (runif(1) < a) { state <- prop; ll_cur <- ll_p; lp_cur <- lp_p; acc[bl] <- acc[bl] + 1 }
      num[bl] <- num[bl] + 1; ls[bl] <- ls[bl] + gam * (a - 0.234) }
    ## between-component move
    m_new <- sample(c("zero", "slab", "one"), 1, prob = beta_move); prop <- state; prop$m <- m_new
    if (m_new == "slab") prop$lchi <- qlogis(rbeta(1, prior$a_chi, prior$b_chi))
    ll_p <- loglik_star(r, w, s2t(prop)); lq_fw <- log(beta_move[match(m_new, c("zero", "slab", "one"))]); lq_bw <- log(beta_move[match(state$m, c("zero", "slab", "one"))])
    lpi <- function(m) switch(m, zero = log(prior$pi0), one = log(prior$pi1), slab = log(1 - prior$pi0 - prior$pi1))
    a <- min(1, exp(ll_p + lpi(m_new) + lq_bw - ll_cur - lpi(state$m) - lq_fw))
    if (runif(1) < a) { state <- prop; ll_cur <- ll_p; lp_cur <- lpr(state, prior); acc["jump"] <- acc["jump"] + 1 }
    num["jump"] <- num["jump"] + 1
    ## marginal blocks: (log sigma_j, xi_j) with scaled random walk; the dependence likelihood is re-evaluated on the re-standardized data
    for (j in 1:2) { nm <- c("m1", "m2")[j]; S <- list(S1, S2)[[j]]; cur <- list(m1, m2)[[j]]
      prop_m <- cur + exp(ls[nm]) * Lse[[j]] * rnorm(2)
      lm_j_cur <- marg_term(S, cur[1], cur[2]) + marg_prior(cur[1], cur[2]); lm_j_p <- marg_term(S, prop_m[1], prop_m[2]) + marg_prior(prop_m[1], prop_m[2])
      if (is.finite(lm_j_p)) { Ep <- E; Ep[, j] <- to_exp(S, prop_m[1], prop_m[2]); rp <- rowSums(Ep); wp <- Ep[, 1] / rp
        ll_p <- loglik_star(rp, wp, s2t(state)); a <- min(1, exp(ll_p + lm_j_p - ll_cur - lm_j_cur)) } else a <- 0
      if (runif(1) < a) { if (j == 1) m1 <- prop_m else m2 <- prop_m; E <- Ep; r <- rp; w <- wp; ll_cur <- ll_p; acc[nm] <- acc[nm] + 1 }
      num[nm] <- num[nm] + 1; ls[nm] <- ls[nm] + gam * (a - 0.234) }
    if (t %in% keep) { s <- s + 1; th <- s2t(state); draws[s, ] <- c(th$phi, th$alpha, th$b, th$u1, th$u2, th$p, th$chibar, th$z); mdraws[s, ] <- c(m1, m2); comp[s] <- state$m }
    if (verbose && t %% 5000 == 0) cat(sprintf("iter %d  loglik %.1f  chibar %.3f  xi %.2f %.2f\n", t, ll_cur, s2t(state)$chibar, m1[2], m2[2]))
  }
  out <- list(draws = draws, component = comp, mdraws = mdraws, acc_rates = acc / pmax(num, 1), prior = prior, r = r, w = w, K = c(K_W = K_W, K_lam = K_lam, K_nu = K_nu),
              mle = rbind(m1 = c(g1$lsig, g1$xi), m2 = c(g2$lsig, g2$xi)), se = rbind(g1$se, g2$se), S = list(S1, S2))
  class(out) <- "star_fit"; out
}
## the cut arm: package sampler on the data standardized at the MLE
cut_mcmc <- function(Y, prior, ...) { S1 <- std_setup(Y[, 1]); S2 <- std_setup(Y[, 2]); g1 <- gpd_fit(Y[, 1], S1$u); g2 <- gpd_fit(Y[, 2], S2$u)
  E <- cbind(to_exp(S1, g1$lsig, g1$xi), to_exp(S2, g2$lsig, g2$xi))
  fit <- star_mcmc(E, prior = prior, ...); fit$mle <- rbind(m1 = c(g1$lsig, g1$xi), m2 = c(g2$lsig, g2$xi)); fit$se <- rbind(g1$se, g2$se); fit$S <- list(S1, S2); fit }

## ---- exact cut posterior: psi from the posterior of the marginal model alone, theta | psi by the dependence sampler ----
## Segmented sampler: one adaptive burn-in, then for each of M marginal draws a warm-started segment whose first
## seg_burn iterations are discarded; the kept draws of all segments are pooled with equal weight.
gpd_ll <- function(S, lsig, xi) { s <- exp(lsig); z <- 1 + xi * S$x[S$ab] / s; if (any(z <= 0) || xi < -0.5) return(-Inf); sum(-lsig - (1 + 1 / xi) * log(z)) }
marg_post_draws <- function(y, S, nd, iter = 20000, burn = 5000) { g <- gpd_fit(y, S$u); cur <- c(g$lsig, g$xi)
  lc <- gpd_ll(S, cur[1], cur[2]) + marg_prior(cur[1], cur[2]); out <- matrix(NA, iter, 2)
  for (it in 1:iter) { pr <- cur + 1.2 * g$se * rnorm(2); lp <- gpd_ll(S, pr[1], pr[2]) + marg_prior(pr[1], pr[2])
    if (is.finite(lp) && log(runif(1)) < lp - lc) { cur <- pr; lc <- lp }; out[it, ] <- cur }
  out[round(seq(burn + 1, iter, length.out = nd)), , drop = FALSE] }
cutx_mcmc <- function(Y, prior, M = 12, burnin = 8000, seg_burn = 1500, seg_keep = 2000, thin = 5, K_W = 5, K_lam = 2, K_nu = 2, beta_move = c(1, 1, 1) / 3) {
  S1 <- std_setup(Y[, 1]); S2 <- std_setup(Y[, 2])
  md1 <- marg_post_draws(Y[, 1], S1, M); md2 <- marg_post_draws(Y[, 2], S2, M); PSI <- cbind(md1, md2)
  stdz <- function(m) { E <- cbind(to_exp(S1, PSI[m, 1], PSI[m, 2]), to_exp(S2, PSI[m, 3], PSI[m, 4])); r <- rowSums(E); list(r = r, w = E[, 1] / r) }
  cur_d <- stdz(1); r <- cur_d$r; w <- cur_d$w; q0 <- prior$q_kappa0
  state <- list(psi = rep(0, K_W - 1), la = log(2), b = rep(0, K_lam + 1), lo1 = log(max(quantile(r, 0.9) / 2 - q0, 0.2)),
                lo2 = log(max(quantile(r, 0.9) / 2 - q0, 0.2)), lp = log(2), lchi = qlogis(0.5), z = rep(0, K_nu + 1), m = "slab", q0 = q0)
  s2t <- star:::state_to_theta; lpr <- star:::log_prior
  ll_cur <- loglik_star(r, w, s2t(state)); lp_cur <- lpr(state, prior)
  blocks <- list(ang = "psi", body = c("la", "b"), thr = c("lo1", "lo2", "lp"), tail = c("lchi", "z"))
  ls <- c(ang = -1, body = -1, thr = -2, tail = -1); acc <- num <- c(ang = 0, body = 0, thr = 0, tail = 0, jump = 0)
  perturb <- function(state, nm, sc) { for (v in nm) state[[v]] <- state[[v]] + rnorm(length(state[[v]]), 0, sc); state }
  lpi <- function(m) switch(m, zero = log(prior$pi0), one = log(prior$pi1), slab = log(1 - prior$pi0 - prior$pi1))
  step <- function(adapt_gam) {
    for (bl in names(blocks)) { nm <- blocks[[bl]]
      if (bl == "tail" && state$m != "slab") nm <- "z"
      if (bl == "tail" && state$m == "zero") { state$z <<- rnorm(K_nu + 1, 0, prior$v_nu); lp_cur <<- lpr(state, prior); next }
      prop <- perturb(state, nm, exp(ls[bl])); ll_p <- loglik_star(r, w, s2t(prop)); lp_p <- lpr(prop, prior)
      a <- min(1, exp(ll_p + lp_p - ll_cur - lp_cur))
      if (runif(1) < a) { state <<- prop; ll_cur <<- ll_p; lp_cur <<- lp_p; acc[bl] <<- acc[bl] + 1 }
      num[bl] <<- num[bl] + 1; ls[bl] <<- ls[bl] + adapt_gam * (a - 0.234) }
    m_new <- sample(c("zero", "slab", "one"), 1, prob = beta_move); prop <- state; prop$m <- m_new
    if (m_new == "slab") prop$lchi <- qlogis(rbeta(1, prior$a_chi, prior$b_chi))
    ll_p <- loglik_star(r, w, s2t(prop)); lq_fw <- log(beta_move[match(m_new, c("zero", "slab", "one"))]); lq_bw <- log(beta_move[match(state$m, c("zero", "slab", "one"))])
    a <- min(1, exp(ll_p + lpi(m_new) + lq_bw - ll_cur - lpi(state$m) - lq_fw))
    if (runif(1) < a) { state <<- prop; ll_cur <<- ll_p; lp_cur <<- lpr(state, prior); acc["jump"] <<- acc["jump"] + 1 }
    num["jump"] <<- num["jump"] + 1 }
  for (t in seq_len(burnin)) step(t^(-0.6))
  nk <- floor(seg_keep / thin); P <- K_W + 1 + (K_lam + 1) + 3 + 1 + (K_nu + 1)
  draws <- matrix(NA_real_, M * nk, P); mdraws <- matrix(NA_real_, M * nk, 4); comp <- character(M * nk); s <- 0
  colnames(draws) <- c(paste0("phi", 1:K_W), "alpha", paste0("b", 0:K_lam), "u1", "u2", "p", "chibar", paste0("z", 0:K_nu)); colnames(mdraws) <- c("lsig1", "xi1", "lsig2", "xi2")
  for (m in seq_len(M)) { d <- stdz(m); r <- d$r; w <- d$w; ll_cur <- loglik_star(r, w, s2t(state))
    for (t in seq_len(seg_burn)) step(0)
    for (t in seq_len(seg_keep)) { step(0); if (t %% thin == 0) { s <- s + 1; th <- s2t(state)
      draws[s, ] <- c(th$phi, th$alpha, th$b, th$u1, th$u2, th$p, th$chibar, th$z); mdraws[s, ] <- PSI[m, ]; comp[s] <- state$m } } }
  out <- list(draws = draws, component = comp, mdraws = mdraws, seg = rep(seq_len(M), each = nk), acc_rates = acc / pmax(num, 1), prior = prior,
              K = c(K_W = K_W, K_lam = K_lam, K_nu = K_nu), S = list(S1, S2), PSI = PSI)
  class(out) <- "star_fit"; out
}
