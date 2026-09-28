\name{model_functions}
\alias{sigma_star}
\alias{rho_star}
\alias{lambda_star}
\alias{nu_star}
\alias{fW_star}
\title{Functional components of the STAR model}
\description{The tail scale, threshold curve, body rate, relaxation function and angular density, evaluated at angles \code{w} in (0,1).}
\usage{
sigma_star(w, chibar, z)
rho_star(w, u1, u2, p)
lambda_star(w, b)
nu_star(w, z)
fW_star(w, phi)
}
\arguments{
\item{w}{vector of angles in (0,1).}
\item{chibar}{coefficient of extremal association, in [0,1].}
\item{z}{Bernstein coefficients of \eqn{\log\nu(w)}.}
\item{u1, u2, p}{anchors and exponent of the power-mean threshold curve.}
\item{b}{Bernstein coefficients of \eqn{\log\lambda(w)}.}
\item{phi}{weights of the Beta mixture defining the angular density.}
}
\value{A vector of the same length as \code{w}: \eqn{\sigma(w) = 1 + \bar\chi\tau(w)^{1+\nu(w)}}, \eqn{\rho(w) = \{(w/u_1)^p + ((1-w)/u_2)^p\}^{-1/p}}, \eqn{\lambda(w)}, \eqn{\nu(w)} or \eqn{f_W(w)}.}
\examples{
w <- seq(0.01, 0.99, by = 0.01)
plot(w, sigma_star(w, 0.8, log(c(0.5, 0.5, 0.5))), type = "l")
}
