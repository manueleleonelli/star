## tab_all.R -- functionals of the saved draws of run_all.R, and bias / RMSE / coverage / width over replicates.
## Every target is a property of the copula of the model at th0 (or of the true margins). Estimators:
##   std  : on the standardized scale of each arm, at exponential quantiles t_q (as in the paper);
##   cop  : the copula of the fitted model, at its own marginal quantiles (free of any standardization);
##   dir  : on the original scale, true levels y mapped by T(y; psi) and read off the fitted joint survivor
##          (the paper's route to original-scale summaries);
##   skl  : on the original scale, true levels y mapped by F(y; psi) and read off the fitted copula (Sklar route:
##          the fitted copula joined to the marginal model).
suppressMessages(library(star))
s <- readLines("run_joint.R"); i <- grep("^prior <- star_prior", s); eval(parse(text = s[4:(i - 1)]))   # th0, F1, F2, js, tq
q_ln <- function(U) ifelse(U <= 0.9, qlnorm(U), u0 + sig0 / xi0 * ((1 - (U - 0.9) / 0.1)^(-xi0) - 1))
q_gp <- function(U) sig0 / xi0 * ((1 - U)^(-xi0) - 1)
th_of <- function(d) list(phi = unname(d[1:5]), alpha = unname(d[6]), b = unname(d[7:9]), u1 = unname(d[10]), u2 = unname(d[11]), p = unname(d[12]), chibar = unname(d[13]), z = unname(d[14:16]))
Tj <- function(y, S, lsig, xi) { z <- 1 + xi * (y - S$u) / exp(lsig); if (z <= 0) return(60); -log(S$p_u) + if (abs(xi) < 1e-8) (y - S$u) / exp(lsig) else log(z) / xi }
## the model's own marginal quantile at exceedance probability pe (margin j)
mqp <- function(th, pe, j) uniroot(function(t) log(max(if (j == 1) js(th, t, 0) else js(th, 0, t), 1e-300)) - log(max(pe, 1e-250)), c(0.001, 80), tol = 1e-8)$root
tq95 <- c(tq(0.95, 1), tq(0.95, 2)); tq99 <- c(tq(0.99, 1), tq(0.99, 2)); tq995 <- c(tq(0.995, 1), tq(0.995, 2))
chi95_0 <- js(th0, tq95[1], tq95[2]) / 0.05; chi99_0 <- js(th0, tq99[1], tq99[2]) / 0.01
truth <- c(chibar = 0.6, chi95 = chi95_0, chi99 = chi99_0, pstar = js(th0, tq995[1], tq995[2]), pmix = js(th0, tq995[1], tq95[2]), mstar = 0.005, xi1 = 0.3, xi2 = 0.3)
fun_draws <- function(A, S, qm) { yq <- qm(c(0.95, 0.99, 0.995))
  t(sapply(seq_len(nrow(A$th)), function(k) { th <- th_of(A$th[k, ]); ps <- unname(A$psi[k, ])
    T1 <- sapply(yq, Tj, S = S[[1]], lsig = ps[1], xi = ps[2]); T2 <- sapply(yq, Tj, S = S[[2]], lsig = ps[3], xi = ps[4])
    C1 <- sapply(exp(-T1), mqp, th = th, j = 1); C2 <- sapply(exp(-T2), mqp, th = th, j = 2)       # Sklar: model quantiles at the marginal-model levels
    c(chibar = th$chibar,
      chi95_std = js(th, -log(.05), -log(.05)) / .05, chi99_std = js(th, -log(.01), -log(.01)) / .01,
      chi95_cop = js(th, mqp(th, .05, 1), mqp(th, .05, 2)) / .05, chi99_cop = js(th, mqp(th, .01, 1), mqp(th, .01, 2)) / .01,
      chi95_dir = js(th, T1[1], T2[1]) / .05, chi99_dir = js(th, T1[2], T2[2]) / .01,
      chi95_skl = js(th, C1[1], C2[1]) / .05, chi99_skl = js(th, C1[2], C2[2]) / .01,
      pstar_dir = js(th, T1[3], T2[3]), pstar_skl = js(th, C1[3], C2[3]), pmix_dir = js(th, T1[3], T2[1]), pmix_skl = js(th, C1[3], C2[1]),
      mstar_dir = js(th, T1[3], 0), mstar_skl = exp(-T1[3]), xi1 = ps[2], xi2 = ps[4]) })) }
summ_one <- function(x) { qm <- if (x$dgp == "ln") q_ln else q_gp
  lapply(c(cut = "cut", cutx = "cutx", joint = "joint"), function(arm) { A <- x[[arm]]; F <- fun_draws(A, x$S, qm)
    Q <- apply(F, 2, quantile, c(.05, .5, .95))
    if (arm == "cut") for (j in 1:2) Q[, paste0("xi", j)] <- x$mle[j, 2] + c(-1, 0, 1) * qnorm(.95) * x$se[j, 2]
    if (arm == "cutx") for (j in 1:2) Q[, paste0("xi", j)] <- quantile(x$marg[, paste0("xi", j)], c(.05, .5, .95))
    list(Q = Q, pad = mean(A$comp == "one"), pzero = mean(A$comp == "zero")) }) }
cache <- if (file.exists("tab_all_cache.rds")) readRDS("tab_all_cache.rds") else list()
for (f in list.files("res_all", full.names = TRUE)) if (is.null(cache[[basename(f)]])) { x <- readRDS(f)
  cache[[basename(f)]] <- c(list(dgp = x$dgp, n = x$n, rep = x$rep, mins = x$mins), summ_one(x)) }
saveRDS(cache, "tab_all_cache.rds")
agg <- function(dgp, n, arm) { R <- Filter(function(z) z$dgp == dgp && z$n == n, cache); if (!length(R)) return(NULL)
  lo <- t(sapply(R, function(z) z[[arm]]$Q[1, ])); med <- t(sapply(R, function(z) z[[arm]]$Q[2, ])); hi <- t(sapply(R, function(z) z[[arm]]$Q[3, ]))
  tr <- truth[sub("_.*", "", colnames(med))]
  list(n_rep = length(R), bias = colMeans(med) - tr, rmse = sqrt(colMeans(sweep(med, 2, tr)^2)), cover = colMeans(sweep(lo, 2, tr, "<=") & sweep(hi, 2, tr, ">=")),
       width = colMeans(hi - lo), pad = mean(sapply(R, function(z) z[[arm]]$pad)), pzero = mean(sapply(R, function(z) z[[arm]]$pzero)),
       mins = mean(sapply(R, function(z) z$mins[arm]))) }
show <- function(dgp, n) { A <- lapply(c(cut = "cut", cutx = "cutx", joint = "joint"), function(a) agg(dgp, n, a)); if (is.null(A$cut)) return(invisible(NULL))
  cat(sprintf("\n===== %s margins, n = %d, %d replicates =====\n", dgp, n, A$cut$n_rep))
  cat(sprintf("P(AD): plug-in %.2f  exact cut %.2f  joint %.2f   P(chibar = 0): %.2f %.2f %.2f   minutes: %.1f %.1f %.1f\n", A$cut$pad, A$cutx$pad, A$joint$pad,
      A$cut$pzero, A$cutx$pzero, A$joint$pzero, A$cut$mins, A$cutx$mins, A$joint$mins))
  nm <- names(A$cut$bias); sc <- ifelse(grepl("^(pstar|pmix|mstar)", nm), 1e4, 1)
  tab <- data.frame(row.names = nm, truth = truth[sub("_.*", "", nm)] * sc)
  for (a in names(A)) { tab[[paste0(a, ".bias")]] <- A[[a]]$bias * sc; tab[[paste0(a, ".rmse")]] <- A[[a]]$rmse * sc; tab[[paste0(a, ".cov")]] <- A[[a]]$cover; tab[[paste0(a, ".wid")]] <- A[[a]]$width * sc }
  print(round(tab, 3)); cat("(pstar, pmix, mstar in units of 1e-4)\n"); invisible(tab) }
if (sys.nframe() == 0) for (d in c("gp", "ln")) for (n in c(400, 2000)) show(d, n)
