# s8_misspec_ablation.R -- the ablation of Study 4 under MISSPECIFICATION:
# does conditioning away threshold uncertainty undercover when the data are
# not from the model? Data from Gaussian(0.5) and Frank(5); truth from Study 3.
source("00_common.R"); source("00_extra.R")
truths <- readRDS("results/s3_truths.rds")
gens <- list(gauss05 = function(n) sim_gauss(n,0.5),
             frank5  = function(n) sim_frank(n,5))
N <- if (SCALE=="test") 700 else 4000
QS <- c(0.95, 0.99, 0.999); gl <- seq(0.5/201,1-0.5/201,length.out=201)
calib_u <- function(r,w,k) uniroot(function(u) mean(r>rho_star(w,u,u,2))-(1-k), c(0.3,12))$root
one_rep <- function(seed, gname, variant) {
  set.seed(seed); E <- gens[[gname]](N); r <- rowSums(E); w <- E[,1]/r
  fx <- if (variant=="full") NULL else { u <- calib_u(r,w,as.numeric(variant)); list(u1=u,u2=u,p=2) }
  fit <- star_mcmc(E, prior=star_prior("beta"), n_iter=SET_NEW$ITER, burnin=SET_NEW$BURN,
                    thin=SET_NEW$THIN, fix_thr=fx, verbose=FALSE)
  tr <- truths[[gname]]$chi
  cq <- chi_q(fit, QS, ndraw=120)
  data.frame(gen=gname, variant=variant, seed=seed,
             chi95_cover =as.numeric(tr[1]>=cq[2,1] & tr[1]<=cq[3,1]),
             chi99_cover =as.numeric(tr[2]>=cq[2,2] & tr[2]<=cq[3,2]),
             chi999_cover=as.numeric(tr[3]>=cq[2,3] & tr[3]<=cq[3,3]),
             chi999_err  =cq[1,3]-tr[3], chi999_w=cq[3,3]-cq[2,3])
}
res <- list(); i <- 0
for (gname in names(gens)) for (v in c("full","0.8","0.9","0.95")) { i <- i+1
  res[[i]] <- run_cell("s8", sprintf("%s_v%s", gname, gsub("[.]","",v)), SET$R,
                       function(r) one_rep(8e4*i+r, gname, v)) }
res <- do.call(rbind, res); saveRDS(res, "results/s8_misspec_ablation.rds")
cat("\n===== Study 8: misspecified ablation (coverage of true chi(q)) =====\n")
print(aggregate(cbind(chi95_cover, chi99_cover, chi999_cover, chi999_err, chi999_w) ~ gen+variant,
                res, mean), digits=2, row.names=FALSE)
cat("\nexpectation: fixed-kappa variants undercover chi(q) and have narrower\n",
    "intervals; the full model, propagating threshold uncertainty, does not.\n")
