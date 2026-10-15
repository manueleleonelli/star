# Model functions (R wrappers around the C++ core) --------------------------

#' Tail scale function sigma(w) = 1 + chibar * tau(w)^(1+nu(w))
#' @param w angles in (0,1); @param chibar in [0,1]; @param z Bernstein coefficients of log nu
#' @export
sigma_star <- function(w, chibar, z) sigma_cpp(w, chibar, z)

#' Power-mean threshold curve rho(w)
#' @export
rho_star <- function(w, u1, u2, p) rho_cpp(w, u1, u2, p)

#' Body rate lambda(w) = exp(Bernstein(b))
#' @export
lambda_star <- function(w, b) exp(bern_eval_vec(w, b))

#' Relaxation function nu(w) = exp(Bernstein(z))
#' @export
nu_star <- function(w, z) exp(bern_eval_vec(w, z))

#' Bernstein angular density f_W(w)
#' @export
fW_star <- function(w, phi) fW_cpp(w, phi)

#' Exact spliced log-likelihood of the model (d = 2)
#' @param r radii; @param w angles; @param theta named list with components
#'   phi, alpha, b, u1, u2, p, chibar, z
#' @export
loglik_star <- function(r, w, theta) {
  loglik_cpp(r, w, theta$phi, theta$alpha, theta$b,
             theta$u1, theta$u2, theta$p, theta$chibar, theta$z)
}

# Simulation -----------------------------------------------------------------

#' Simulate from the model (bivariate)
#' @return a list with E (n x 2 matrix), r, w
#' @export
rstar <- function(n, theta) {
  K <- length(theta$phi)
  comp <- sample.int(K, n, replace = TRUE, prob = theta$phi)
  w <- rbeta(n, comp, K - comp + 1)
  lam <- lambda_star(w, theta$b)
  rho <- rho_star(w, theta$u1, theta$u2, theta$p)
  sig <- sigma_star(w, theta$chibar, theta$z)
  pi_w <- pgamma(rho, shape = theta$alpha, rate = lam)
  body <- runif(n) < pi_w
  r <- numeric(n)
  # body: Gamma truncated to (0, rho(w)) by inverse cdf
  ub <- runif(sum(body)) * pi_w[body]
  r[body] <- qgamma(ub, shape = theta$alpha, rate = lam[body])
  # tail: rho + Exp(sigma)
  r[!body] <- rho[!body] + rexp(sum(!body)) * sig[!body]
  list(E = cbind(r * w, r * (1 - w)), r = r, w = w)
}

# Priors ---------------------------------------------------------------------

#' Prior specification
#' @param route "beta" (continuous prior on chibar) or "threecomp" (atoms at 0 and 1 and a Beta slab)
#' @param v_p sd of the N(log 2, v_p^2) prior on log p, supporting p in (0, Inf)
#' @param p_min lower bound on the exponent p (0 = unrestricted; 1 restricts to convex thresholds)
#' @export
star_prior <- function(route = c("beta", "threecomp"),
                        a_dir = 1, v_alpha = 0.5, v_lambda = 1, v_nu = 1,
                        kappa0 = 0.8, s_u = 1, v_p = 0.75, p_min = 0,
                        a_chi = 1, b_chi = 1, pi0 = 0.25, pi1 = 0.25) {
  route <- match.arg(route)
  list(route = route, a_dir = a_dir, v_alpha = v_alpha, v_lambda = v_lambda,
       v_nu = v_nu, kappa0 = kappa0, q_kappa0 = -log(1 - kappa0), s_u = s_u,
       v_p = v_p, p_min = p_min, a_chi = a_chi, b_chi = b_chi, pi0 = pi0, pi1 = pi1)
}

# Log prior on the sampling scale (unconstrained parameters), Jacobians included.
# psi: softmax coordinates of phi (length K-1, reference = last category)
# state: list(psi, la = log alpha, b, lo1 = log o1, lo2 = log o2, lz = qlogis(zeta), lchi = qlogis(chibar), z, m)
log_prior <- function(state, prior) {
  phi <- softmax_phi(state$psi)
  lp <- prior$a_dir * sum(log(phi))                               # Dirichlet(a) + softmax Jacobian
  lp <- lp + dnorm(state$la, log(2), prior$v_alpha, log = TRUE)   # log alpha
  lp <- lp + sum(dnorm(state$b, 0, prior$v_lambda, log = TRUE))   # body coefficients
  for (lo in c(state$lo1, state$lo2)) {                           # half-normal offsets, log scale
    o <- exp(lo)
    lp <- lp - o^2 / (2 * prior$s_u^2) + lo
  }
  if (!is.null(prior$p_min) && exp(state$lp) < prior$p_min) return(-Inf)  # hard floor on p
  lp <- lp + dnorm(state$lp, log(2), prior$v_p, log = TRUE)      # log p ~ N(log 2, v_p^2)
  lp <- lp + sum(dnorm(state$z, 0, prior$v_nu, log = TRUE))       # shape coefficients
  # chibar contribution depends on the route / component
  if (prior$route == "beta" || state$m == "slab") {
    chib <- plogis(state$lchi)
    lp <- lp + dbeta(chib, prior$a_chi, prior$b_chi, log = TRUE) +
      log(chib) + log(1 - chib)                                   # slab density + logit Jacobian
    if (prior$route == "threecomp") lp <- lp + log(1 - prior$pi0 - prior$pi1)
  } else if (state$m == "zero") {
    lp <- lp + log(prior$pi0)
  } else {                                                        # m == "one"
    lp <- lp + log(prior$pi1)
  }
  lp
}

softmax_phi <- function(psi) {
  e <- exp(c(psi, 0))
  e / sum(e)
}

state_to_theta <- function(state) {
  chib <- switch(state$m, zero = 0, one = 1, slab = plogis(state$lchi))
  list(phi = softmax_phi(state$psi), alpha = exp(state$la), b = state$b,
       u1 = state$q0 + exp(state$lo1), u2 = state$q0 + exp(state$lo2),
       p = exp(state$lp), chibar = chib, z = state$z)
}
