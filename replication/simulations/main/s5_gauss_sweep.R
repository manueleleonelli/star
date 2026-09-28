# s5_gauss_sweep.R -- the classification "money figure": posterior probability
# of asymptotic dependence for Gaussian data across rho, the case where
# classical procedures misclassify AI as AD. Saves results for plotting.
source("00_common.R"); source("00_extra.R")
RHOS <- c(0.3, 0.5, 0.7, 0.9); NS <- c(2000, 8000)
if (SCALE == "test") { RHOS <- c(0.3, 0.9); NS <- c(700) }
one_rep <- function(seed, rho, n) {
  set.seed(seed)
  fit <- star_mcmc(sim_gauss(n, rho), prior = star_prior("threecomp"),
                    n_iter = SET_NEW$ITER, burnin = SET_NEW$BURN, thin = SET_NEW$THIN, verbose = FALSE)
  cp <- class_probs(fit)
  data.frame(rho = rho, n = n, seed = seed, P1 = cp[3],
             chibar_mean = mean(fit$draws[, "chibar"]))
}
res <- list(); i <- 0
for (rho in RHOS) for (n in NS) { i <- i+1
  res[[i]] <- run_cell("s5", sprintf("rho%02.0f_n%d", 100*rho, n), SET$R,
                       function(r) one_rep(5e4*i + r, rho, n)) }
res <- do.call(rbind, res); saveRDS(res, "results/s5_gauss_sweep.rds")
cat("\n===== Study 5: Gaussian sweep, P(AD | data) =====\n")
print(aggregate(cbind(P1, chibar_mean) ~ rho + n, res, function(x) round(mean(x),3)), row.names = FALSE)
cat("\ntrue eta = (1+rho)/2, all AI (chi = 0); the curve is the money figure.\n")
