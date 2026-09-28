# s11_theory_checks.R -- confirm Proposition 3 in the fitted object: at the
# atom, chi(q) should decay like c/log{1/(1-q)} with c the closed-form
# constant. Fit to model data at chibar = 1 and compare fitted chi(q) with the
# theoretical slowly varying form across a grid of q.
source("00_common.R"); source("00_extra.R")
N <- if (SCALE=="test") 700 else 4000
QS <- c(0.9, 0.99, 0.999, 0.9999)
c_prop3 <- function(th) {                     # constant of Proposition 3
  eta <- (1+th$chibar)/2; nub <- exp(star:::bern_eval_vec(0.5, th$z))
  s <- 1 + 3*th$chibar + 2*th$chibar*nub
  rho_half <- rho_star(0.5, th$u1, th$u2, th$p); lam_half <- lambda_star(0.5, th$b)
  A <- pgamma(rho_half, th$alpha, rate=lam_half, lower.tail=FALSE)
  fW_half <- fW_star(0.5, th$phi)
  2*eta^2 * fW_half * A * exp(rho_half/(1+th$chibar)) / s
}
one_rep <- function(seed) {
  set.seed(seed); th0 <- theta_model(chibar=1)
  fit <- star_mcmc(rstar(N, th0)$E, prior=star_prior("beta"),
                    n_iter=SET_NEW$ITER, burnin=SET_NEW$BURN, thin=SET_NEW$THIN, verbose=FALSE)
  idx <- round(seq(1,nrow(fit$draws),length.out=100))
  chi_num <- sapply(QS, function(q){ tq <- -log(1-q); gl <- seq(0.5/201,1-0.5/201,length.out=201)
    mean(sapply(idx, function(s){ th <- star:::theta_row(fit,s)
      star:::joint_surv_cpp(tq,tq,gl,rep(1/201,201),th$phi,th$alpha,th$b,th$u1,th$u2,th$p,th$chibar,th$z)/(1-q)}))})
  c_hat <- mean(sapply(idx, function(s) c_prop3(star:::theta_row(fit,s))))
  chi_theory <- c_hat / log(1/(1-QS))
  data.frame(seed=seed, q=QS, chi_num=chi_num, chi_theory=chi_theory,
             ratio=chi_num/chi_theory)
}
NR <- if (SCALE=="test") 3 else 15
res <- run_cell("s11", "atom", NR, function(r) one_rep(11e4+r), chunk = 5)
saveRDS(res, "results/s11_theory_checks.rds")
cat("\n===== Study 11: Proposition 3 in the fitted object (chibar = 1) =====\n")
print(aggregate(cbind(chi_num, chi_theory, ratio) ~ q, res, mean), digits=3, row.names=FALSE)
cat("\nas q -> 1 the ratio chi_num / (c/log) should approach 1, confirming the\n",
    "slowly varying form; deep-q agreement is the check that matters.\n")
