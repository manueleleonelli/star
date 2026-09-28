\name{star-package}
\alias{star-package}
\alias{star}
\docType{package}
\title{Self-Thresholding Angular-Radial Models of Extremal Dependence}
\description{
Bayesian fitting of the self-thresholding angular-radial (STAR) model of Leonelli and Gamerman. After marginal standardization to standard exponential margins, a bivariate observation is written as a radius \eqn{R = E_1 + E_2} and an angle \eqn{W = E_1/R}. The angle follows a Bernstein polynomial density, and the radius given the angle follows a spliced distribution: a Gamma body with angle-dependent rate below the threshold curve \eqn{\rho(w)}, a power-mean contour with anchors \eqn{u_1, u_2} and exponent \eqn{p}, and an exponential tail with scale \eqn{\sigma(w) = 1 + \bar\chi \tau(w)^{1+\nu(w)}} above it. The threshold curve is estimated with all other parameters, \eqn{\bar\chi} is the coefficient of extremal association, and \eqn{\sigma(\cdot)} is the reciprocal gauge of the implied limit set.

The main functions are \code{\link{star_prior}} and \code{\link{star_mcmc}}; \code{\link{rstar}} simulates from the model, and \code{\link{chi_q}}, \code{\link{adf_lambda}}, \code{\link{joint_surv}}, \code{\link{eta_draws}} and \code{\link{class_probs}} summarize a fit. The trivariate implementation of Section 6.2 of the paper is in \code{system.file("trivariate", "star3.R", package = "star")}.
}
\references{
Leonelli, M. and Gamerman, D. Self-thresholding angular-radial models of extremal dependence. Submitted to \emph{Bayesian Analysis}.
}
\keyword{package}
