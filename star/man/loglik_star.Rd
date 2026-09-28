\name{loglik_star}
\alias{loglik_star}
\title{Log-likelihood of the bivariate STAR model}
\description{Exact spliced log-likelihood (11) of the paper, evaluated in C++.}
\usage{loglik_star(r, w, theta)}
\arguments{
\item{r, w}{radii and angles of the standardized observations.}
\item{theta}{parameter list as in \code{\link{rstar}}.}
}
\value{The log-likelihood, a scalar.}
