\name{star_mcmc}
\alias{star_mcmc}
\title{Posterior sampling for the STAR model}
\description{Adaptive Metropolis-within-Gibbs sampler for the bivariate STAR model, with the blocks, adaptation and between-component moves described in Section S3 of the supplementary material.}
\usage{
star_mcmc(E, prior = star_prior(), K_W = 5, K_lam = 2, K_nu = 2,
          n_iter = 30000, burnin = 10000, thin = 5,
          beta_move = c(1, 1, 1)/3, fix_thr = NULL, verbose = TRUE)
}
\arguments{
\item{E}{an \eqn{n \times 2} matrix of observations on standard exponential margins.}
\item{prior}{output of \code{\link{star_prior}}.}
\item{K_W}{number of Beta densities in the angular density.}
\item{K_lam, K_nu}{degrees of the Bernstein bases of \eqn{\log\lambda(w)} and \eqn{\log\nu(w)}.}
\item{n_iter, burnin, thin}{chain length, burn-in (during which the proposal scales are adapted) and thinning.}
\item{beta_move}{proposal probabilities of the components (zero, slab, one) in the between-component move of the three-component prior.}
\item{fix_thr}{\code{NULL}, or a list with components \code{u1}, \code{u2} and \code{p} that fixes the threshold curve and removes its block from the sampler.}
\item{verbose}{logical; print progress every 5000 iterations.}
}
\value{An object of class \code{"star_fit"}: a list with \code{draws} (matrix of retained draws of the angular weights, \code{alpha}, the body coefficients \code{b0, b1, ...}, \code{u1}, \code{u2}, \code{p}, \code{chibar} and the coefficients \code{z0, z1, ...}), \code{component} (the component label of each draw under the three-component prior), \code{loglik} (log-likelihood trace), \code{acc_rates}, \code{prior}, the radii \code{r} and angles \code{w}, and the basis sizes \code{K}.}
\examples{
\donttest{
set.seed(1)
th <- list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0),
           u1 = 2.2, u2 = 2.2, p = 2, chibar = 0.6, z = log(c(0.5, 0.5, 0.5)))
dat <- rstar(1000, th)
fit <- star_mcmc(dat$E, star_prior("threecomp"), n_iter = 6000, burnin = 2000,
                 thin = 4, verbose = FALSE)
class_probs(fit)
chi_q(fit, c(0.95, 0.99))
}
}
