## star3.R -- the trivariate implementation of STAR used in Section 6.2 of the paper (details in Supplement S3.7).
## Angular density: Dirichlet (Bernstein) basis of degree K on the 2-simplex.
## Body: Gamma(alpha, lambda(w)) with log lambda linear in w.  Threshold: p-ball, rho(w) = {sum (w_j/u_j)^p}^{-1/p}.
## Tail: exponential (xi = 0) with scale sigma(w) = 1/g(w), the additive-pairwise gauge
##   g(w) = 1 - sum_{j<k} (w_j + w_k) {1 - 1/sigma_jk(v)},  sigma_jk(v) = 1 + chibar_jk tau(v)^(1+nu_jk),  v = w_j/(w_j+w_k),
## subject to the envelope g(w) >= max_j w_j on the simplex (checked on a grid at every tail proposal).
## Inference: adaptive random-walk Metropolis within Gibbs on unconstrained scales, slab prior on each chibar_jk.

PAIRS <- rbind(c(1, 2), c(1, 3), c(2, 3))

## ---- bases and functional parameters ----------------------------------------------------------
basis_index <- function(K) { ix <- expand.grid(k1 = 0:K, k2 = 0:K); ix <- ix[ix$k1 + ix$k2 <= K, ]; cbind(ix$k1, ix$k2, K - ix$k1 - ix$k2) }
basis3 <- function(W, K) {                        # n x M matrix of Dirichlet densities with parameters k+1
  ix <- basis_index(K); lw <- log(pmax(W, 1e-300))
  sapply(seq_len(nrow(ix)), function(m) { a <- ix[m, ] + 1
    exp(lgamma(sum(a)) - sum(lgamma(a)) + lw %*% (a - 1)) })
}
lambda3 <- function(W, b) exp(W %*% b)[, 1]
rho3 <- function(W, u, p) rowSums(sweep(W, 2, u, "/")^p)^(-1/p)
tau2 <- function(v) pmin(v, 1 - v) / pmax(v, 1 - v)
gauge3 <- function(W, chib, nu) {                 # g(w) for rows of W on the simplex (sum w = 1)
  g <- rep(1, nrow(W))
  for (q in 1:3) { j <- PAIRS[q, 1]; k <- PAIRS[q, 2]; s <- W[, j] + W[, k]; v <- ifelse(s > 0, W[, j] / s, 0.5)
    sjk <- 1 + chib[q] * tau2(v)^(1 + nu[q]); g <- g - s * (1 - 1 / sjk) }
  g
}
sigma3 <- function(W, chib, nu) 1 / gauge3(W, chib, nu)
simplex_grid <- function(m) { ix <- expand.grid(i = 1:(m - 1), j = 1:(m - 1)); ix <- ix[ix$i + ix$j < m, ]
  W <- cbind(ix$i, ix$j, m - ix$i - ix$j) / m; list(W = W, h2 = 1 / m^2) }
ADM_GRID <- simplex_grid(80)$W; ADM_MAX <- pmax(ADM_GRID[, 1], ADM_GRID[, 2], ADM_GRID[, 3])
admissible <- function(chib, nu, W = ADM_GRID) all(gauge3(W, chib, nu) >= ADM_MAX - 1e-9)

## ---- parameters --------------------------------------------------------------------------------
unpack <- function(x, K) {
  M <- nrow(basis_index(K)); i <- 0
  psi <- x[i + seq_len(M - 1)]; i <- i + M - 1
  e <- exp(c(psi, 0)); phi <- e / sum(e)
  a <- x[i + 1]; b <- x[i + 2:4]; o <- x[i + 5:7]; lp <- x[i + 8]; cc <- x[i + 9:11]; z <- x[i + 12]
  list(phi = phi, psi = psi, alpha = exp(a), a = a, b = b, u = -log(1 - 0.8) + exp(o), o = o, p = exp(lp), lp = lp,
       chib = plogis(cc), cc = cc, nu = rep(exp(z), 3), z = z)
}
pack <- function(th) c(log(th$phi[-length(th$phi)]) - log(th$phi[length(th$phi)]), log(th$alpha), th$b,
                        log(pmax(th$u + log(1 - 0.8), 1e-6)), log(th$p), qlogis(th$chib), log(th$nu[1]))
blocks <- function(K) { M <- nrow(basis_index(K)); m1 <- M - 1
  list(ang1 = seq_len(ceiling(m1 / 2)), ang2 = (ceiling(m1 / 2) + 1):m1, body = m1 + 1:4, thr = m1 + 5:8, tail = m1 + 9:12) }

## ---- likelihood and prior ----------------------------------------------------------------------
loglik3 <- function(th, dat) {
  fW <- as.numeric(dat$B %*% th$phi)
  lam <- lambda3(dat$W, th$b); rho <- rho3(dat$W, th$u, th$p); sig <- sigma3(dat$W, th$chib, th$nu)
  tail <- dat$R > rho
  ll <- sum(log(fW)) - 2 * sum(dat$logR)
  ll <- ll + sum(dgamma(dat$R[!tail], shape = th$alpha, rate = lam[!tail], log = TRUE))
  ll + sum(pgamma(rho[tail], shape = th$alpha, rate = lam[tail], lower.tail = FALSE, log.p = TRUE) - log(sig[tail]) - (dat$R[tail] - rho[tail]) / sig[tail])
}
logprior3 <- function(th, hyper) {
  with(hyper, sum(a_dir * log(th$phi)) + dnorm(th$a, log(3), v_alpha, log = TRUE) + sum(dnorm(th$b, 0, v_lambda, log = TRUE)) +
    sum(-exp(2 * th$o) / (2 * s_u^2) + th$o) + dnorm(th$lp, log(2), v_p, log = TRUE) +
    sum(log(th$chib * (1 - th$chib))) + sum(dnorm(th$z, 0, v_nu, log = TRUE)))
}
default_hyper <- list(a_dir = 1, v_alpha = 0.5, v_lambda = 1, v_nu = 1, v_p = 0.75, s_u = 1, p_min = 1)

## ---- sampler -------------------------------------------------------------------------------------
init3 <- function(dat, K) {
  M <- nrow(basis_index(K)); R <- dat$R
  u0 <- rep(quantile(R, 0.8) / 3^(1 - 1/2), 3)      # so that rho at the barycentre sits near the 0.8 quantile of R
  th <- list(phi = rep(1 / M, M), alpha = 3, b = c(0, 0, 0), u = pmax(u0, -log(0.2) + 0.05), p = 2, chib = c(0.3, 0.3, 0.3), nu = c(1, 1, 1))
  pack(th)
}
mcmc3 <- function(E, n_iter = 30000, burnin = 10000, thin = 5, K = 4, hyper = default_hyper, seed = 1, verbose = TRUE, init = NULL) {
  set.seed(seed)
  R <- rowSums(E); W <- E / R; dat <- list(R = R, logR = log(R), W = W, B = basis3(W, K))
  x <- if (is.null(init)) init3(dat, K) else init; bl <- blocks(K)
  post <- function(x, check = TRUE) { th <- unpack(x, K); if (th$p < hyper$p_min || (check && !admissible(th$chib, th$nu))) return(-Inf)
    loglik3(th, dat) + logprior3(th, hyper) }
  lp <- post(x); if (!is.finite(lp)) stop("bad initial state")
  sc <- c(ang1 = 0.3, ang2 = 0.3, body = 0.05, thr = 0.05, tail = 0.2); acc <- att <- sapply(bl, function(b) 0)
  L <- lapply(bl, function(b) diag(length(b)))                       # Cholesky factors of the proposal covariances
  H <- matrix(NA, burnin, length(x))
  keep <- seq(burnin + 1, n_iter, by = thin); out <- matrix(NA, length(keep), length(x)); llout <- numeric(length(keep)); r <- 0
  for (it in seq_len(n_iter)) {
    for (nb in names(bl)) { b <- bl[[nb]]; xp <- x; xp[b] <- x[b] + sc[nb] * as.numeric(L[[nb]] %*% rnorm(length(b))); lpp <- post(xp, check = nb == "tail"); att[nb] <- att[nb] + 1
      if (is.finite(lpp) && log(runif(1)) < lpp - lp) { x <- xp; lp <- lpp; acc[nb] <- acc[nb] + 1 } }
    if (it <= burnin) { H[it, ] <- x
      if (it %% 50 == 0) { for (nb in names(bl)) { rate <- acc[nb] / att[nb]; sc[nb] <- sc[nb] * exp(2 * (rate - 0.25) / sqrt(it / 50)) }; acc[] <- 0; att[] <- 0 }
      if (it >= 1000 && it %% 500 == 0) for (nb in names(bl)) { b <- bl[[nb]]; S <- cov(H[(it %/% 2):it, b, drop = FALSE]) + diag(1e-6, length(b))
        Lnew <- t(chol(S)); s0 <- sqrt(mean(diag(L[[nb]] %*% t(L[[nb]])))); s1 <- sqrt(mean(diag(S))); L[[nb]] <- Lnew; sc[nb] <- sc[nb] * s0 / s1 } }
    if (it %in% keep) { r <- r + 1; out[r, ] <- x; llout[r] <- lp }
    if (verbose && it %% 5000 == 0) cat(sprintf("iter %d  logpost %.1f  scales %s\n", it, lp, paste(signif(sc, 2), collapse = " ")))
  }
  list(draws = out, logpost = llout, K = K, acc = acc / pmax(att, 1), scales = sc, hyper = hyper, E = E)
}
theta_at <- function(fit, s) unpack(fit$draws[s, ], fit$K)

## ---- derived quantities --------------------------------------------------------------------------
QGRID <- simplex_grid(240)                         # quadrature grid on the simplex (Lebesgue in (w1, w2))
QB <- NULL
fW_grid <- function(th, K, G = QGRID) { if (is.null(QB) || ncol(QB) != length(th$phi)) QB <<- basis3(G$W, K); as.numeric(QB %*% th$phi) }
radial_surv <- function(s, W, th) {                # P(R > s | W = w), the spliced survivor, vectorized over rows
  lam <- lambda3(W, th$b); rho <- rho3(W, th$u, th$p); sig <- sigma3(W, th$chib, th$nu)
  body <- s <= rho
  out <- numeric(length(s))
  out[body] <- pgamma(s[body], shape = th$alpha, rate = lam[body], lower.tail = FALSE)
  out[!body] <- pgamma(rho[!body], shape = th$alpha, rate = lam[!body], lower.tail = FALSE) * exp(-(s[!body] - rho[!body]) / sig[!body])
  out
}
joint_surv3 <- function(th, K, t, G = QGRID) {     # P(E_j > t_j for j with t_j > 0), t a length-3 vector with zeros allowed
  W <- G$W; s <- rep(0, nrow(W)); for (j in 1:3) if (t[j] > 0) s <- pmax(s, t[j] / W[, j])
  sum(fW_grid(th, K, G) * radial_surv(s, W, th)) * G$h2
}
chi3 <- function(th, K, q, C = 1:3) { t <- rep(0, 3); t[C] <- -log(1 - q); joint_surv3(th, K, t) / (1 - q) }
ETA_GRID <- { m <- 300; ix <- expand.grid(i = 0:m, j = 0:m); ix <- ix[ix$i + ix$j <= m, ]; cbind(ix$i, ix$j, m - ix$i - ix$j) / m }
eta3 <- function(th, C = 1:3, W = ETA_GRID) { sig <- sigma3(W, th$chib, th$nu); min(1, max(do.call(pmin, lapply(C, function(j) W[, j])) * sig)) }
exc_frac3 <- function(th, K, G = QGRID) { W <- G$W; rho <- rho3(W, th$u, th$p); lam <- lambda3(W, th$b)
  sum(fW_grid(th, K, G) * pgamma(rho, shape = th$alpha, rate = lam, lower.tail = FALSE)) * G$h2 }

## ---- simulation from the model -------------------------------------------------------------------
rdirichlet <- function(n, a) { x <- matrix(rgamma(n * length(a), rep(a, each = n)), n); x / rowSums(x) }
rstar3 <- function(n, th, K) {
  ix <- basis_index(K); comp <- sample.int(nrow(ix), n, TRUE, th$phi)
  W <- t(sapply(comp, function(m) rdirichlet(1, ix[m, ] + 1)))
  lam <- lambda3(W, th$b); rho <- rho3(W, th$u, th$p); sig <- sigma3(W, th$chib, th$nu)
  pb <- pgamma(rho, shape = th$alpha, rate = lam); tail <- runif(n) > pb
  R <- numeric(n); R[!tail] <- qgamma(runif(sum(!tail)) * pb[!tail], shape = th$alpha, rate = lam[!tail])
  R[tail] <- rho[tail] + rexp(sum(tail)) * sig[tail]
  list(E = R * W, tail = tail)
}
