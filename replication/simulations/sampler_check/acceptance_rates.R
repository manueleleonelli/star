## acceptance_rates.R -- the verification runs quoted in Section S3.6 of the supplement: realized acceptance rates of the
## four random-walk blocks and of the between-component move, ten datasets of n = 1500 from the model under the
## three-component prior. Output in acceptance_rates.log. About ten minutes on two cores.
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
suppressMessages(library(star)); library(parallel)
th <- function(cb) list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0), u1 = 2.2, u2 = 2.2, p = 2, chibar = cb, z = log(c(0.5, 0.5, 0.5)))
jobs <- expand.grid(cb = c(0, 0.3, 0.6, 0.9, 1), seed = 1:2)
res <- mclapply(seq_len(nrow(jobs)), function(i) { set.seed(jobs$seed[i] * 100 + i)
  d <- rstar(1500, th(jobs$cb[i]))
  f <- star_mcmc(d$E, star_prior("threecomp"), n_iter = 12000, burnin = 4000, thin = 4, verbose = FALSE)
  c(cb = jobs$cb[i], f$acc_rates) }, mc.cores = 2)
print(round(do.call(rbind, res), 3))
