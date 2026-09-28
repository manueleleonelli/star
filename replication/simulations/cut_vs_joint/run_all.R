## run_all.R -- the three arms (plug-in cut, exact cut, joint) on the same datasets, saving thinned draws so that
## any functional can be computed afterwards. Two data-generating processes share the copula samples (same seeds):
##   ln: lognormal bulk with a GPD(1, 0.3) tail above the 0.9 quantile (the design of run_joint.R; plug-in and joint
##       arms reproduce res/, the exact cut reproduces res_cutx/);
##   gp: GPD(1, 0.3) over the whole range, so that the marginal model is exact above any threshold.
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
source("joint.R"); library(parallel)
s <- readLines("run_joint.R"); i <- grep("^prior <- star_prior", s); eval(parse(text = s[4:(i - 1)]))   # th0, F1, F2, truth pieces
prior <- star_prior("threecomp", p_min = 1)
q_ln <- function(U) ifelse(U <= 0.9, qlnorm(U), u0 + sig0 / xi0 * ((1 - (U - 0.9) / 0.1)^(-xi0) - 1))
q_gp <- function(U) sig0 / xi0 * ((1 - U)^(-xi0) - 1)
M <- as.integer(Sys.getenv("MCUT", "10"))
dir.create("res_all", showWarnings = FALSE)
keep_draws <- function(fit, joint) { idx <- round(seq(1, nrow(fit$draws), length.out = 300))
  list(th = fit$draws[idx, ], comp = fit$component[idx], psi = if (joint) fit$mdraws[idx, ] else matrix(c(fit$mle[1, ], fit$mle[2, ]), 300, 4, byrow = TRUE), acc = fit$acc_rates) }
marg_long <- function(y, S, iter = 40000, burn = 5000, nkeep = 2000) { g <- gpd_fit(y, S$u); cur <- c(g$lsig, g$xi)
  lc <- gpd_ll(S, cur[1], cur[2]) + marg_prior(cur[1], cur[2]); out <- matrix(NA, iter, 2)
  for (it in 1:iter) { pr <- cur + 1.2 * g$se * rnorm(2); lp <- gpd_ll(S, pr[1], pr[2]) + marg_prior(pr[1], pr[2])
    if (is.finite(lp) && log(runif(1)) < lp - lc) { cur <- pr; lc <- lp }; out[it, ] <- cur }
  out[round(seq(burn + 1, iter, length.out = nkeep)), ] }
one_all <- function(job) { dgp <- job$dgp; n <- job$n; rep <- job$rep
  f <- sprintf("res_all/%s_%d_%02d.rds", dgp, n, rep); if (file.exists(f)) return(invisible(NULL))
  qm <- if (dgp == "ln") q_ln else q_gp
  set.seed(5000 + 100 * (n == 400) + rep); E <- rstar(n, th0)$E; Y <- cbind(qm(F1(E[, 1])), qm(F2(E[, 2])))
  t0 <- Sys.time(); fc <- cut_mcmc(Y, prior, n_iter = 20000, burnin = 8000, thin = 5, verbose = FALSE); tc <- as.numeric(Sys.time() - t0, units = "mins")
  t0 <- Sys.time(); fj <- joint_mcmc(Y, prior, n_iter = 20000, burnin = 8000, thin = 5); tj <- as.numeric(Sys.time() - t0, units = "mins")
  set.seed(9000 + 100 * (n == 400) + rep); t0 <- Sys.time(); fx <- cutx_mcmc(Y, prior, M = M, burnin = 8000, seg_burn = 1500, seg_keep = 2000, thin = 5); tx <- as.numeric(Sys.time() - t0, units = "mins")
  set.seed(77 + rep + 100 * (n == 400)); ml <- cbind(marg_long(Y[, 1], fc$S[[1]]), marg_long(Y[, 2], fc$S[[2]])); colnames(ml) <- c("lsig1", "xi1", "lsig2", "xi2")
  out <- list(dgp = dgp, n = n, rep = rep, S = fc$S, mle = fc$mle, se = fc$se, cut = keep_draws(fc, FALSE), joint = keep_draws(fj, TRUE), cutx = keep_draws(fx, TRUE),
              cutx_psi = fx$PSI, marg = ml, mins = c(cut = tc, joint = tj, cutx = tx))
  saveRDS(out, f); invisible(NULL) }
jobs <- c(lapply(1:20, function(r) list(dgp = "gp", n = 400, rep = r)), lapply(1:20, function(r) list(dgp = "ln", n = 400, rep = r)),
          lapply(1:20, function(r) list(dgp = "gp", n = 2000, rep = r)), lapply(1:20, function(r) list(dgp = "ln", n = 2000, rep = r)))
if (nzchar(Sys.getenv("ONLY"))) jobs <- Filter(function(j) paste(j$dgp, j$n) == Sys.getenv("ONLY"), jobs)
invisible(mclapply(jobs, function(j) tryCatch(one_all(j), error = function(e) cat("ERR", j$dgp, j$n, j$rep, conditionMessage(e), "\n")), mc.cores = 2, mc.preschedule = FALSE))
cat("done", format(Sys.time()), "\n")
