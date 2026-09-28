\name{star_prior}
\alias{star_prior}
\title{Prior specification for the STAR model}
\description{Collects the hyperparameters of the priors of Sections 4.2 and 4.3 of the paper.}
\usage{
star_prior(route = c("beta", "threecomp"), a_dir = 1, v_alpha = 0.5,
           v_lambda = 1, v_nu = 1, kappa0 = 0.8, s_u = 1, v_p = 0.75,
           p_min = 0, a_chi = 1, b_chi = 1, pi0 = 0.25, pi1 = 0.25)
}
\arguments{
\item{route}{\code{"beta"} for a continuous Beta prior on \eqn{\bar\chi}, or \code{"threecomp"} for the three-component prior with atoms at 0 and 1 and a Beta slab.}
\item{a_dir}{parameter of the symmetric Dirichlet prior on the angular weights.}
\item{v_alpha}{standard deviation of the normal prior on \eqn{\log\alpha}, the log body shape.}
\item{v_lambda}{standard deviation of the normal priors on the Bernstein coefficients of \eqn{\log\lambda(w)}.}
\item{v_nu}{standard deviation of the normal priors on the Bernstein coefficients of \eqn{\log\nu(w)}.}
\item{kappa0}{reference level: the anchors are \eqn{u_j = -\log(1-\kappa_0) + o_j}.}
\item{s_u}{scale of the half-normal priors on the anchor offsets \eqn{o_j}.}
\item{v_p}{standard deviation of the \eqn{N(\log 2, v_p^2)} prior on \eqn{\log p}.}
\item{p_min}{lower bound on the exponent \eqn{p}; 0 leaves it unrestricted and 1 restricts the threshold curve to convex shapes.}
\item{a_chi, b_chi}{parameters of the Beta prior (or slab) on \eqn{\bar\chi}.}
\item{pi0, pi1}{prior probabilities of the atoms \eqn{\bar\chi = 0} and \eqn{\bar\chi = 1} under \code{route = "threecomp"}.}
}
\value{A list of hyperparameters, passed to \code{\link{star_mcmc}}.}
\examples{
pr <- star_prior("threecomp", p_min = 1)
}
