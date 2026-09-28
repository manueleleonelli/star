# s6_hill.R -- frequentist Hill-type benchmark for eta (mean excess of the
# structure variable), for context alongside the model estimates of Study 3.
source("00_common.R"); source("00_extra.R")
gens <- list(
  gauss05  = list(f=function(n) sim_gauss(n,0.5),       eta=0.75),
  invlog05 = list(f=function(n) sim_invlogistic(n,0.5), eta=2^-0.5),
  frank5   = list(f=function(n) sim_frank(n,5),         eta=0.50),
  t05_4    = list(f=function(n) sim_t(n,0.5,4),         eta=1.00),
  model06  = list(f=function(n) rstar(n, theta_model(chibar=0.6))$E, eta=0.80))
N <- if (SCALE=="test") 700 else 4000
one_rep <- function(seed, gname) {
  set.seed(seed); g <- gens[[gname]]; h <- hill_eta(g$f(N), kappa = 0.9)
  data.frame(gen=gname, seed=seed, eta_true=g$eta, eta_hat=h["eta"],
             err=h["eta"]-g$eta, cover=as.numeric(g$eta>=h["lo"] & g$eta<=h["hi"]), k=h["k"])
}
res <- list()
for (gname in names(gens))
  res[[gname]] <- run_cell("s6", gname, SET$R, function(r) one_rep(6e4+37*r+nchar(gname), gname))
res <- do.call(rbind, res); saveRDS(res, "results/s6_hill.rds")
cat("\n===== Study 6: Hill-type eta benchmark (kappa = 0.9) =====\n")
print(aggregate(cbind(eta_true, eta_hat, err, cover) ~ gen, res, mean), digits=2, row.names=FALSE)
cat("\ncompare err/coverage with the model's eta in Study 3: the Hill estimator\n",
    "has no threshold-curve uncertainty and no body/tail separation.\n")
