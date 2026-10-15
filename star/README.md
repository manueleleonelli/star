# star

Self-thresholding angular-radial (STAR) models of extremal dependence (Leonelli and Gamerman).

```r
remotes::install_github("manueleleonelli/star", subdir = "star")   # needs Rcpp and a C++ compiler
library(star)
th  <- list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0),
            u1 = 2.2, u2 = 2.2, p = 2, chibar = 0.6, z = log(c(0.5, 0.5, 0.5)))
dat <- rstar(1500, th)
fit <- star_mcmc(dat$E, prior = star_prior("threecomp", p_min = 1))
class_probs(fit); chi_q(fit, c(0.9, 0.95, 0.99))
```

The trivariate implementation used in Section 6.2 of the paper is `inst/trivariate/star3.R`:
`source(system.file("trivariate", "star3.R", package = "star"))`.
