# s2_classification.R -- Study 2: classification operating characteristics
# Three-component prior fitted to model truths (both atoms and two slab
# values) and to real copulas; reports posterior class probabilities, the
# 0.5-rule error and the log score. True class: AD iff eta = 1.

source("00_common.R")
set.seed(22)

gens <- list(
  m_chi0    = list(f = function(n) rstar(n, theta_model(chibar = 0))$E,   AD = FALSE),
  m_chi05   = list(f = function(n) rstar(n, theta_model(chibar = 0.5))$E, AD = FALSE),
  m_chi095  = list(f = function(n) rstar(n, theta_model(chibar = 0.95))$E,AD = FALSE),
  m_chi1    = list(f = function(n) rstar(n, theta_model(chibar = 1))$E,   AD = TRUE),
  gauss05   = list(f = function(n) sim_gauss(n, 0.5),        AD = FALSE),
  gauss09   = list(f = function(n) sim_gauss(n, 0.9),        AD = FALSE),
  invlog05  = list(f = function(n) sim_invlogistic(n, 0.5),  AD = FALSE),
  logist05  = list(f = function(n) sim_logistic(n, 0.5),     AD = TRUE),
  t05_4     = list(f = function(n) sim_t(n, 0.5, 4),         AD = TRUE))

NS <- c(SET$N1 * 2, SET$N2 * 2)   # classification needs more data: 2000, 8000 at beach scale

one_rep <- function(seed, gname, n) {
  set.seed(seed); g <- gens[[gname]]
  fit <- star_mcmc(g$f(n), prior = star_prior("threecomp"),
                    n_iter = round(SET$ITER * 1.25), burnin = SET$BURN,
                    thin = SET$THIN, verbose = FALSE)
  cp <- class_probs(fit)
  pAD <- cp["P_AD"]; pAI_side <- cp["P_AI_exact"] + cp["P_slab"]
  data.frame(gen = gname, n = n, seed = seed,
             P0 = cp[1], Pslab = cp[2], P1 = cp[3],
             chibar_mean = mean(fit$draws[, "chibar"]),
             correct_05 = as.numeric((pAD > 0.5) == g$AD),
             logscore = log(pmax(if (g$AD) pAD else pAI_side, 1e-4)))
}

res <- list(); i <- 0
for (gname in names(gens)) for (n in NS) {
  i <- i + 1
  res[[i]] <- run_cell("s2", sprintf("%s_n%d", gname, n), SET$R,
                       function(r) one_rep(2e4*i + r, gname, n))
}
res <- do.call(rbind, res); saveRDS(res, "results/s2_all.rds")

cat("\n===== Study 2: classification =====\n")
print(aggregate(cbind(P0, Pslab, P1, chibar_mean, correct_05, logscore) ~ gen + n,
                res, mean), digits = 2, row.names = FALSE)
cat("\nreading: P1 should be high for m_chi1, logist05, t05_4 and near zero\n",
    "elsewhere; gauss09 is the classical hard case, report its P1 as it is.\n")
