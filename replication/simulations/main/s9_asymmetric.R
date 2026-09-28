# s9_asymmetric.R -- scenario (iv): asymmetric structure. Two truths: (a) a
# within-model asymmetric truth (angle-varying nu, unequal anchors), testing
# recovery of asymmetric sigma and rho; (b) an out-of-model non-exchangeable
# asymptotically dependent copula (Khoudraji-Gumbel), testing the classifier
# and the asymmetric fit.
source("00_common.R"); source("00_extra.R")
N <- if (SCALE=="test") 700 else 4000
one_rep_within <- function(seed) {
  set.seed(seed); th0 <- theta_asym(chibar=0.6); dat <- rstar(N, th0)
  fit <- star_mcmc(dat$E, prior=star_prior("beta"), n_iter=SET_NEW$ITER,
                    burnin=SET_NEW$BURN, thin=SET_NEW$THIN, verbose=FALSE)
  sig <- function(wq) sapply(round(seq(1,nrow(fit$draws),length.out=120)), function(s){
    th <- star:::theta_row(fit,s); sigma_star(wq, th$chibar, th$z)})
  rho <- function(wq) sapply(round(seq(1,nrow(fit$draws),length.out=120)), function(s){
    th <- star:::theta_row(fit,s); rho_star(wq, th$u1, th$u2, th$p)})
  s03 <- sig(0.3); s07 <- sig(0.7); r03 <- rho(0.3); r07 <- rho(0.7)
  data.frame(seed=seed,
             sig_asym_true = sigma_star(0.3,0.6,th0$z) - sigma_star(0.7,0.6,th0$z),
             sig_asym_hat  = median(s03) - median(s07),
             rho_asym_true = rho_star(0.3,th0$u1,th0$u2,2) - rho_star(0.7,th0$u1,th0$u2,2),
             rho_asym_hat  = median(r03) - median(r07))
}
one_rep_khoud <- function(seed) {
  set.seed(seed); E <- sim_khoudraji(N, theta_gumbel=2, a=0.9, b=0.5)
  fit <- star_mcmc(E, prior=star_prior("threecomp"), n_iter=SET_NEW$ITER,
                    burnin=SET_NEW$BURN, thin=SET_NEW$THIN, verbose=FALSE)
  r03 <- sapply(round(seq(1,nrow(fit$draws),length.out=120)), function(s){
    th <- star:::theta_row(fit,s); rho_star(0.3,th$u1,th$u2,th$p)})
  r07 <- sapply(round(seq(1,nrow(fit$draws),length.out=120)), function(s){
    th <- star:::theta_row(fit,s); rho_star(0.7,th$u1,th$u2,th$p)})
  data.frame(seed=seed, P1=class_probs(fit)[3],
             rho_asym_hat=median(r03)-median(r07))
}
within <- run_cell("s9", "within", SET$R, function(r) one_rep_within(9e4+r))
khoud  <- run_cell("s9", "khoud",  SET$R, function(r) one_rep_khoud(95e3+r))
saveRDS(list(within=within, khoud=khoud), "results/s9_asymmetric.rds")
cat("\n===== Study 9a: within-model asymmetric recovery =====\n")
print(colMeans(within[,-1]), digits=3)
cat("(sig_asym and rho_asym: sign and magnitude of hat should track true)\n")
cat("\n===== Study 9b: Khoudraji-Gumbel (asymmetric AD) =====\n")
cat(sprintf("P(AD | data) = %.2f ; fitted rho asymmetry = %.2f (nonzero = asymmetric fit)\n",
    mean(khoud$P1), mean(khoud$rho_asym_hat)))
