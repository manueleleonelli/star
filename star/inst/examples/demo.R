library(star)
set.seed(11)
theta0 <- list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0),
               u1 = 2.2, u2 = 2.2, p = 2, chibar = 0.6, z = c(log(0.5), log(0.5), log(0.5)))
dat <- rstar(1500, theta0)
fit <- star_mcmc(dat$E, prior = star_prior("beta"), n_iter = 12000, burnin = 4000, thin = 4)
print(fit$acc_rates)
print(quantile(fit$draws[, "chibar"], c(0.05, 0.5, 0.95)))
