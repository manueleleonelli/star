\name{summaries}
\alias{eta_draws}
\alias{class_probs}
\alias{chi_q}
\alias{adf_lambda}
\alias{joint_surv}
\title{Posterior summaries of a STAR fit}
\description{Derived quantities of Sections 3 and 4.5 of the paper, computed from the retained draws of \code{\link{star_mcmc}}.}
\usage{
eta_draws(fit)
class_probs(fit)
chi_q(fit, q, ndraw = 200, M = 201)
adf_lambda(fit, omega = seq(0.05, 0.95, by = 0.05), ndraw = 200,
           grid = seq(0.001, 0.999, length.out = 999))
joint_surv(fit, t1, t2, ndraw = 200, M = 201)
}
\arguments{
\item{fit}{an object returned by \code{\link{star_mcmc}}.}
\item{q}{vector of levels in (0,1).}
\item{omega}{vector of rays in (0,1) at which the angular dependence function is evaluated.}
\item{t1, t2}{levels on the standard exponential scale.}
\item{ndraw}{number of equally spaced posterior draws used.}
\item{M}{number of quadrature nodes for the joint survivor (12).}
\item{grid}{grid of angles for the minimization in Proposition 4.}
}
\value{\code{eta_draws}: draws of \eqn{\eta = (1+\bar\chi)/2}. \code{class_probs}: posterior probabilities of the atom at zero, the slab and the atom at one (three-component prior only). \code{chi_q}: a matrix with the posterior mean and 95\% interval of \eqn{\chi(q)} for each \code{q}. \code{adf_lambda}: a list with \code{omega} and the posterior mean and 95\% interval of \eqn{\Lambda(\omega)}. \code{joint_surv}: the posterior mean and 95\% interval of \eqn{P(E_1 > t_1, E_2 > t_2)}.}
