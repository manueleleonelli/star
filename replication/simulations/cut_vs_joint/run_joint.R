## run_joint.R -- cut versus joint posterior: two settings, NREP replicates each
source("joint.R"); library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
NREP <- as.integer(Sys.getenv("NREP", "20")); ITER <- as.integer(Sys.getenv("ITER", "20000")); BURN <- as.integer(Sys.getenv("BURN", "8000"))
th0 <- list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0), u1 = 2.2, u2 = 2.2, p = 2, chibar = 0.6, z = rep(log(0.5), 3))
sig0 <- 1; xi0 <- 0.3; u0 <- qlnorm(0.9)                     # true margins: lognormal bulk, GPD(1, 0.3) tail above the 0.9 quantile
qmarg <- function(U) ifelse(U <= 0.9, qlnorm(U), u0 + sig0 / xi0 * ((1 - (U - 0.9) / 0.1)^(-xi0) - 1))
gl <- seq(0.5 / 201, 1 - 0.5 / 201, length.out = 201); gw <- rep(1 / 201, 201)
js <- function(th, t1, t2) star:::joint_surv_cpp(t1, t2, gl, gw, th$phi, th$alpha, th$b, th$u1, th$u2, th$p, th$chibar, th$z)
## exact marginal distribution functions of the model, so that Y_j = qmarg{F_j(E_j)} has margin G exactly and the copula of the model
tg <- seq(0, 30, by = 0.01); lS1 <- log(pmax(sapply(tg, function(t) js(th0, t, 0)), 1e-300)); lS2 <- log(pmax(sapply(tg, function(t) js(th0, 0, t)), 1e-300))
F1 <- function(x) 1 - exp(approxfun(tg, lS1, rule = 2)(x)); F2 <- function(x) 1 - exp(approxfun(tg, lS2, rule = 2)(x))
tq <- function(q, j) uniroot(function(t) (if (j == 1) js(th0, t, 0) else js(th0, 0, t)) - (1 - q), c(0.01, 30))$root   # marginal quantiles of the model
chi_cop <- function(q) js(th0, tq(q, 1), tq(q, 2)) / (1 - q)                                                          # chi(q) of the copula of the model
ystar <- qmarg(0.995); pstar0 <- js(th0, tq(0.995, 1), tq(0.995, 2))
truth <- c(chibar = 0.6, chi95 = chi_cop(0.95), chi99 = chi_cop(0.99), pstar = pstar0, xi1 = xi0, xi2 = xi0)
prior <- star_prior("threecomp", p_min = 1)
summarize <- function(fit, joint) {
  D <- fit$draws; idx <- round(seq(1, nrow(D), length.out = 300)); S <- fit$S
  Tlev <- function(j, lsig, xi) { s <- exp(lsig); z <- 1 + xi * (ystar - S[[j]]$u) / s; if (z <= 0) return(60); if (abs(xi) < 1e-8) return(-log(S[[j]]$p_u) + (ystar - S[[j]]$u) / s); -log(S[[j]]$p_u) + log(z) / xi }
  out <- t(sapply(idx, function(k) { th <- star:::theta_row(fit, k)
    m <- if (joint) fit$mdraws[k, ] else c(fit$mle[1, ], fit$mle[2, ])
    c(chibar = th$chibar, chi95 = js(th, -log(.05), -log(.05)) / 0.05, chi99 = js(th, -log(.01), -log(.01)) / 0.01,
      pstar = js(th, Tlev(1, m[1], m[2]), Tlev(2, m[3], m[4])), xi1 = m[2], xi2 = m[4]) }))
  q <- apply(out, 2, quantile, c(.05, .5, .95))
  if (!joint) { for (j in 1:2) { nm <- paste0("xi", j); q[, nm] <- fit$mle[j, 2] + c(-1, 0, 1) * qnorm(.95) * c(fit$se[j, 2], 0, fit$se[j, 2]) } }
  list(lo = q[1, ], med = q[2, ], hi = q[3, ], pad = mean(fit$component == "one"), acc = fit$acc_rates)
}
one <- function(job) { case <- job$case; rep <- job$rep; n <- job$n
  f <- sprintf("res/%s_%02d.rds", case, rep); if (file.exists(f)) return(readRDS(f))
  set.seed(5000 + 100 * (case == "B") + rep); E <- rstar(n, th0)$E; Y <- cbind(qmarg(F1(E[, 1])), qmarg(F2(E[, 2])))
  t0 <- Sys.time(); fc <- cut_mcmc(Y, prior, n_iter = ITER, burnin = BURN, thin = 5, verbose = FALSE); tc <- as.numeric(Sys.time() - t0, units = "mins")
  t0 <- Sys.time(); fj <- joint_mcmc(Y, prior, n_iter = ITER, burnin = BURN, thin = 5); tj <- as.numeric(Sys.time() - t0, units = "mins")
  out <- list(case = case, rep = rep, n = n, cut = summarize(fc, FALSE), joint = summarize(fj, TRUE), mins = c(cut = tc, joint = tj))
  saveRDS(out, f); out }
dir.create("res", showWarnings = FALSE)
jobs <- c(lapply(seq_len(NREP), function(r) list(case = "A", rep = r, n = 2000)), lapply(seq_len(NREP), function(r) list(case = "B", rep = r, n = 400)))
res <- mclapply(jobs, one, mc.cores = 2, mc.preschedule = FALSE)
tab <- function(case, arm) { R <- Filter(function(x) x$case == case, res); med <- t(sapply(R, function(x) x[[arm]]$med)); lo <- t(sapply(R, function(x) x[[arm]]$lo)); hi <- t(sapply(R, function(x) x[[arm]]$hi))
  tr <- truth[colnames(med)]
  data.frame(truth = tr, bias = colMeans(med) - tr, rmse = sqrt(colMeans(sweep(med, 2, tr)^2)), cover = colMeans(sweep(lo, 2, tr, "<=") & sweep(hi, 2, tr, ">=")), width = colMeans(hi - lo)) }
for (cs in c("A", "B")) for (arm in c("cut", "joint")) { cat("\n==== case", cs, arm, "====\n"); print(round(tab(cs, arm), 4)) }
cat("\nminutes:", round(colMeans(t(sapply(res, `[[`, "mins"))), 1), "\n")
saveRDS(list(res = res, truth = truth, th0 = th0), "joint_results.rds")
