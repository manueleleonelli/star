## wanliu_cut60.R -- cut posterior for Section 6 with 60 marginal draws: (sigma_j, xi_j) from the posterior of the
## marginal model alone (priors N(0, 10^2) on log sigma_j and N(0, 0.5^2) on xi_j), theta refitted on each standardized sample
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
source("joint.R"); library(parallel)
x <- readRDS("results/fit_Wanliu.rds"); a <- x$data; Y <- cbind(a$PM2.5, a$O3); n <- nrow(Y)
prior <- star_prior("threecomp", p_min = 1); M <- 60
lpri <- function(lsig, xi) dnorm(lsig, 0, 10, log = TRUE) + dnorm(xi, 0, 0.5, log = TRUE)
S <- list(std_setup(Y[, 1]), std_setup(Y[, 2]))
mdraw <- function(j, iter = 70000, burn = 10000) { g <- gpd_fit(Y[, j], S[[j]]$u); cur <- c(g$lsig, g$xi); lc <- gpd_ll(S[[j]], cur[1], cur[2]) + lpri(cur[1], cur[2])
  out <- matrix(NA, iter, 2); for (it in 1:iter) { pr <- cur + 1.2 * g$se * rnorm(2); lp <- gpd_ll(S[[j]], pr[1], pr[2]) + lpri(pr[1], pr[2])
    if (is.finite(lp) && log(runif(1)) < lp - lc) { cur <- pr; lc <- lp }; out[it, ] <- cur }
  out[round(seq(burn + 1, iter, length.out = M)), ] }
## resumable: the fits are skipped when results/wanliu_cut60_fits.rds exists
if (file.exists("results/wanliu_cut60_fits.rds")) cat("results/wanliu_cut60_fits.rds exists, nothing to do\n") else {
set.seed(8); md1 <- mdraw(1); md2 <- mdraw(2)
fits <- mclapply(1:M, function(m) { E <- cbind(to_exp(S[[1]], md1[m, 1], md1[m, 2]), to_exp(S[[2]], md2[m, 1], md2[m, 2]))
  set.seed(200 + m); f <- star_mcmc(E, prior = prior, n_iter = 16000, burnin = 6000, thin = 5, verbose = FALSE)
  idx <- round(seq(1, nrow(f$draws), length.out = 40)); list(draws = f$draws[idx, ], component = f$component[idx], K = f$K, m = c(md1[m, ], md2[m, ])) }, mc.cores = 2, mc.preschedule = FALSE)
saveRDS(list(fits = fits, md1 = md1, md2 = md2), "results/wanliu_cut60_fits.rds"); cat("done", format(Sys.time()), "\n") }
