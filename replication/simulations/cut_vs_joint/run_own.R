## run_own.R -- P(AD | data) for the same copula samples as run_all.R, fitted on the model's own margins
## (rstar output, no marginal transformation), three-component prior, 20 replicates at n = 2000 and 400
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
source("joint.R"); library(parallel)
s <- readLines("run_joint.R"); i <- grep("^prior <- star_prior", s); eval(parse(text = s[4:(i - 1)]))
prior <- star_prior("threecomp", p_min = 1); dir.create("res_own", showWarnings = FALSE)
jobs <- c(lapply(1:20, function(r) list(n = 400, rep = r)), lapply(1:20, function(r) list(n = 2000, rep = r)))
invisible(mclapply(jobs, function(j) { f <- sprintf("res_own/%d_%02d.rds", j$n, j$rep); if (file.exists(f)) return(NULL)
  set.seed(5000 + 100 * (j$n == 400) + j$rep); E <- rstar(j$n, th0)$E
  fit <- star_mcmc(E, prior = prior, n_iter = 20000, burnin = 8000, thin = 5, verbose = FALSE)
  idx <- round(seq(1, nrow(fit$draws), length.out = 300)); th <- lapply(idx, function(k) star:::theta_row(fit, k))
  chi <- sapply(th, function(t) c(js(t, -log(.05), -log(.05)) / .05, js(t, -log(.01), -log(.01)) / .01))
  saveRDS(list(n = j$n, rep = j$rep, pad = mean(fit$component == "one"), pzero = mean(fit$component == "zero"), chibar = median(fit$draws[, "chibar"]),
               chi95 = median(chi[1, ]), chi99 = median(chi[2, ])), f); NULL }, mc.cores = 2, mc.preschedule = FALSE))
res <- lapply(list.files("res_own", full.names = TRUE), readRDS)
for (n in c(400, 2000)) { R <- Filter(function(z) z$n == n, res)
  cat(sprintf("n = %d (%d reps): P(AD) %.2f  P(zero) %.2f  median chibar %.2f  chi95 %.3f  chi99 %.3f\n", n, length(R), mean(sapply(R, `[[`, "pad")), mean(sapply(R, `[[`, "pzero")),
      mean(sapply(R, `[[`, "chibar")), mean(sapply(R, `[[`, "chi95")), mean(sapply(R, `[[`, "chi99")))) }
