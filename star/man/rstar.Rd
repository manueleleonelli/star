\name{rstar}
\alias{rstar}
\title{Simulate from the bivariate STAR model}
\description{Draws an angle from the Bernstein density, a body or tail indicator, and a truncated Gamma or shifted exponential radius.}
\usage{rstar(n, theta)}
\arguments{
\item{n}{sample size.}
\item{theta}{a list with components \code{phi} (angular weights), \code{alpha} (body shape), \code{b} (Bernstein coefficients of \eqn{\log\lambda}), \code{u1}, \code{u2}, \code{p} (threshold curve), \code{chibar} and \code{z} (Bernstein coefficients of \eqn{\log\nu}).}
}
\value{A list with the \eqn{n \times 2} matrix \code{E} on the exponential scale, the radii \code{r} and the angles \code{w}.}
\examples{
th <- list(phi = rep(0.2, 5), alpha = 2, b = c(0, 0, 0), u1 = 2.2, u2 = 2.2,
           p = 2, chibar = 0.6, z = log(c(0.5, 0.5, 0.5)))
x <- rstar(500, th)
}
