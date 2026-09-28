# s3_misspecification.R -- Study 3: estimation of true tail quantities when
# the data come from copulas outside the model class. Reports error in eta,
# chi(q) at q in {0.95, 0.99, 0.999} and the ADF at omega = 0.3 against the
# analytic or Monte Carlo truths.

source("00_common.R")
set.seed(33)

gens <- list(
  gauss05  = list(f = function(n) sim_gauss(n, 0.5),       gen = "gauss",       par = 0.5),
  invlog05 = list(f = function(n) sim_invlogistic(n, 0.5), gen = "invlogistic", par = 0.5),
  frank5   = list(f = function(n) sim_frank(n, 5),         gen = "frank",       par = 5),
  t05_4    = list(f = function(n) sim_t(n, 0.5, 4),        gen = "t",           par = c(0.5, 4)))

QS <- c(0.95, 0.99, 0.999)
truths <- lapply(gens, function(g) list(
  eta = true_eta(g$gen, g$par[1]),
  chi = sapply(QS, function(q) true_chi_q(g$gen, g$par, q,
               nmc = if (SCALE == "test") 2e5 else 2e7)),
  lam03 = true_lambda_adf(g$gen, g$par[1], 0.3)))
saveRDS(truths, "results/s3_truths.rds")

one_rep <- function(seed, gname, n) {
  set.seed(seed); g <- gens[[gname]]; tr <- truths[[gname]]
  fit <- star_mcmc(g$f(n), prior = star_prior("beta"), n_iter = SET$ITER,
                    burnin = SET$BURN, thin = SET$THIN, verbose = FALSE)
  et <- eta_draws(fit)
  cq <- chi_q(fit, QS, ndraw = 120)
  ad <- adf_lambda(fit, omega = 0.3, ndraw = 120)
  data.frame(gen = gname, n = n, seed = seed,
             eta_err = median(et) - tr$eta,
             eta_cover = as.numeric(tr$eta >= quantile(et, .05) &
                                    tr$eta <= quantile(et, .95)),
             chi95_err  = cq[1, 1] - tr$chi[1],
             chi99_err  = cq[1, 2] - tr$chi[2],
             chi999_err = cq[1, 3] - tr$chi[3],
             lam03_err  = ad$summary[1, 1] - tr$lam03)
}

res <- list(); i <- 0
for (gname in names(gens)) {
  i <- i + 1; n <- SET$N2
  res[[i]] <- run_cell("s3", sprintf("%s_n%d", gname, n), SET$R,
                       function(r) one_rep(3e4*i + r, gname, n))
}
res <- do.call(rbind, res); saveRDS(res, "results/s3_all.rds")

cat("\n===== Study 3: misspecification (errors vs truth) =====\n")
print(aggregate(cbind(eta_err, eta_cover, chi95_err, chi99_err, chi999_err,
                      lam03_err) ~ gen, res, mean), digits = 2, row.names = FALSE)
cat("\nreading: invlog05 tests the exact-class gap of the paper; frank5 tests\n",
    "the body/tail separation (eta_err should be near zero despite body\n",
    "association); t05_4 shows the behaviour at the AD boundary, eta -> 1.\n")
