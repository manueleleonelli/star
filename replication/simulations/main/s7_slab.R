# s7_slab.R -- slab sensitivity at the boundary: does concentrating the slab
# near one sharpen the AD/AI contest for strongly dependent AI generators?
source("00_common.R"); source("00_extra.R")
gens <- list(gauss09  = function(n) sim_gauss(n,0.9),
             m_chi095 = function(n) rstar(n, theta_model(chibar=0.95))$E,
             m_chi1   = function(n) rstar(n, theta_model(chibar=1))$E)
NS <- c(2000, 8000); if (SCALE=="test") NS <- 700
slabs <- list(uniform=c(1,1), concentrated=c(2,0.8))
one_rep <- function(seed, gname, n, slab) {
  set.seed(seed); ab <- slabs[[slab]]
  fit <- star_mcmc(gens[[gname]](n),
                    prior = star_prior("threecomp", a_chi=ab[1], b_chi=ab[2]),
                    n_iter=SET_NEW$ITER, burnin=SET_NEW$BURN, thin=SET_NEW$THIN, verbose=FALSE)
  data.frame(gen=gname, n=n, slab=slab, seed=seed, P1=class_probs(fit)[3])
}
res <- list(); i <- 0
for (gname in names(gens)) for (n in NS) for (slab in names(slabs)) { i <- i+1
  res[[i]] <- run_cell("s7", sprintf("%s_n%d_%s", gname, n, slab), SET$R,
                       function(r) one_rep(7e4*i+r, gname, n, slab)) }
res <- do.call(rbind, res); saveRDS(res, "results/s7_slab.rds")
cat("\n===== Study 7: slab sensitivity, P(AD | data) =====\n")
print(reshape(aggregate(P1 ~ gen+n+slab, res, mean), idvar=c("gen","n"),
              timevar="slab", direction="wide"), digits=2, row.names=FALSE)
