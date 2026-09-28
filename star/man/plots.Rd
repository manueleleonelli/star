\name{plots}
\alias{plot_silhouette}
\alias{plot_threshold}
\title{Diagnostic plots for a STAR fit}
\description{\code{plot_silhouette} draws the scaled sample \eqn{E/\log n} with posterior draws of the limit set boundary; \code{plot_threshold} draws the standardized data with posterior draws of the threshold curve.}
\usage{
plot_silhouette(fit, ndraw = 100, wgrid = seq(0.002, 0.998, length.out = 199))
plot_threshold(fit, ndraw = 100, wgrid = seq(0.002, 0.998, length.out = 199))
}
\arguments{
\item{fit}{an object returned by \code{\link{star_mcmc}}.}
\item{ndraw}{number of posterior draws drawn.}
\item{wgrid}{grid of angles.}
}
\value{Called for their side effect.}
