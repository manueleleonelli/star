# Adaptive Metropolis-within-Gibbs sampler ----------------------------------

#' Fit the self-thresholding angular-radial model
#'
#' @param E n x 2 matrix on standard exponential margins
#' @param prior output of \code{star_prior()}
#' @param K_W,K_lam,K_nu basis sizes (angular density; log body rate; log shape)
#' @param n_iter,burnin,thin chain settings; adaptation runs during burn-in
#' @param beta_move proposal probabilities of (zero, slab, one) for the
#'   between-component move under the three-component prior
#' @param fix_thr NULL, or list(u1, u2, p): fixes the threshold curve at the
#'   given values and removes the threshold block from the sampler (used by
#'   the fixed-threshold ablation of the simulation study)
#' @param verbose print progress
#' @return list with draws (matrix), component (character, threecomp only),
#'   loglik trace, acceptance rates, prior, data, basis sizes
#' @export
star_mcmc <- function(E, prior = star_prior(), K_W = 5, K_lam = 2, K_nu = 2,
                       n_iter = 30000, burnin = 10000, thin = 5,
                       beta_move = c(1, 1, 1) / 3, fix_thr = NULL, verbose = TRUE) {
  r <- rowSums(E); w <- E[, 1] / r
  n <- length(r)

  # --- initialization (data driven) ---
  q0 <- prior$q_kappa0
  state <- list(
    psi = rep(0, K_W - 1), la = log(2), b = rep(0, K_lam + 1),
    lo1 = log(max(quantile(r, 0.9) / 2 - q0, 0.2)),
    lo2 = log(max(quantile(r, 0.9) / 2 - q0, 0.2)),
    lp = log(2), lchi = qlogis(0.5), z = rep(0, K_nu + 1),
    m = if (prior$route == "threecomp") "slab" else "slab", q0 = q0
  )
  if (!is.null(fix_thr)) {
    state$q0 <- 0
    state$lo1 <- log(fix_thr$u1); state$lo2 <- log(fix_thr$u2)
    state$lp <- log(fix_thr$p)
  }
  ll_cur <- loglik_star(r, w, state_to_theta(state))
  lp_cur <- log_prior(state, prior)

  # --- blocks: names of state components updated jointly ---
  blocks <- list(ang = "psi", body = c("la", "b"), thr = c("lo1", "lo2", "lp"),
                 tail = c("lchi", "z"))
  if (!is.null(fix_thr)) blocks$thr <- NULL
  ls <- c(ang = -1, body = -1, thr = -2, tail = -1)   # log proposal scales
  acc <- num <- c(ang = 0, body = 0, thr = 0, tail = 0, jump = 0)
  targ <- c(ang = 0.234, body = 0.234, thr = 0.234, tail = 0.234)

  keep <- seq(burnin + 1, n_iter, by = thin)
  P <- K_W + 1 + (K_lam + 1) + 3 + 1 + (K_nu + 1)
  draws <- matrix(NA_real_, length(keep), P)
  colnames(draws) <- c(paste0("phi", 1:K_W), "alpha", paste0("b", 0:K_lam),
                       "u1", "u2", "p", "chibar", paste0("z", 0:K_nu))
  comp <- character(length(keep)); llt <- numeric(length(keep)); s <- 0

  perturb <- function(state, nm, sc) {
    for (v in nm) state[[v]] <- state[[v]] + rnorm(length(state[[v]]), 0, sc)
    state
  }

  for (t in seq_len(n_iter)) {
    gam <- if (t <= burnin) t^(-0.6) else 0
    # ---- random-walk blocks ----
    for (bl in names(blocks)) {
      nm <- blocks[[bl]]
      if (bl == "tail" && prior$route == "threecomp" && state$m != "slab")
        nm <- "z"                                   # chibar pinned at an atom
      if (bl == "tail" && prior$route == "threecomp" && state$m == "zero") {
        # likelihood free of z: Gibbs refresh from the prior
        state$z <- rnorm(K_nu + 1, 0, prior$v_nu)
        lp_cur <- log_prior(state, prior); next
      }
      prop <- perturb(state, nm, exp(ls[bl]))
      ll_p <- loglik_star(r, w, state_to_theta(prop))
      lp_p <- log_prior(prop, prior)
      a <- min(1, exp(ll_p + lp_p - ll_cur - lp_cur))
      if (runif(1) < a) { state <- prop; ll_cur <- ll_p; lp_cur <- lp_p; acc[bl] <- acc[bl] + 1 }
      num[bl] <- num[bl] + 1
      ls[bl] <- ls[bl] + gam * (a - targ[bl])
    }
    # ---- between-component move (three-component prior) ----
    if (prior$route == "threecomp") {
      m_new <- sample(c("zero", "slab", "one"), 1, prob = beta_move)
      prop <- state; prop$m <- m_new
      if (m_new == "slab") prop$lchi <- qlogis(rbeta(1, prior$a_chi, prior$b_chi))
      ll_p <- loglik_star(r, w, state_to_theta(prop))
      # proposal densities: atoms have mass beta_j; slab proposes from the slab
      # prior, so the Beta densities cancel and only the component masses and
      # prior odds remain in the ratio (Appendix B of the paper).
      lq_fw <- log(beta_move[match(m_new, c("zero", "slab", "one"))])
      lq_bw <- log(beta_move[match(state$m, c("zero", "slab", "one"))])
      lpi <- function(m) switch(m, zero = log(prior$pi0), one = log(prior$pi1),
                                slab = log(1 - prior$pi0 - prior$pi1))
      a <- min(1, exp(ll_p + lpi(m_new) + lq_bw - ll_cur - lpi(state$m) - lq_fw))
      if (runif(1) < a) {
        state <- prop; ll_cur <- ll_p; lp_cur <- log_prior(state, prior)
        acc["jump"] <- acc["jump"] + 1
      }
      num["jump"] <- num["jump"] + 1
    }
    # ---- store ----
    if (t %in% keep) {
      s <- s + 1
      th <- state_to_theta(state)
      draws[s, ] <- c(th$phi, th$alpha, th$b, th$u1, th$u2, th$p, th$chibar, th$z)
      comp[s] <- state$m; llt[s] <- ll_cur
    }
    if (verbose && t %% 5000 == 0)
      cat(sprintf("iter %d  loglik %.1f  chibar %.3f\n", t, ll_cur,
                  state_to_theta(state)$chibar))
  }
  out <- list(draws = draws, component = comp, loglik = llt,
              acc_rates = acc / pmax(num, 1), prior = prior,
              r = r, w = w, K = c(K_W = K_W, K_lam = K_lam, K_nu = K_nu))
  class(out) <- "star_fit"
  out
}
