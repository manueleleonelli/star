## ht_sens_run.R -- sensitivity of HT to its conditioning threshold and margins (Table S8); writes results/htsens_*.rds
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
source("ht_sens_fun.R")
## Study 3: copulas at n = 4000, exponential margins
s3 <- readRDS("../main/results/s3_all.rds"); tr <- readRDS("../main/results/s3_truths.rds"); QS <- c(0.95, 0.99, 0.999)
gens <- list(gauss05 = function(n) sim_gauss(n, 0.5), invlog05 = function(n) sim_invlogistic(n, 0.5), frank5 = function(n) sim_frank(n, 5), t05_4 = function(n) sim_t(n, 0.5, 4))
r3 <- do.call(rbind, mclapply(seq_len(nrow(s3)), function(i) {
  g <- as.character(s3$gen[i]); sd <- s3$seed[i]; set.seed(sd); E <- gens[[g]](s3$n[i])
  cbind(gen = g, seed = sd, truth = tr[[g]]$chi, one(E, QS, exp_margins, sd)) }, mc.cores = 2))
saveRDS(r3, "results/htsens_s3.rds"); cat("s3 done\n")
## Study 1: within the model, q = 0.99
s1 <- readRDS("../main/results/s1_all.rds"); x <- s1[s1$quantity == "chi_q", ]
cells <- expand.grid(chibar = c(0.15, 0.45, 0.75, 0.95), n = c(1000, 4000)); cells$u1 <- 2.2; cells$u2 <- 2.2
cells <- rbind(cells, data.frame(chibar = 0.45, n = 4000, u1 = 2.8, u2 = 1.7))
TM <- lapply(seq_len(nrow(cells)), function(i) true_margins(theta_model(cells$chibar[i], u1 = cells$u1[i], u2 = cells$u2[i])))
r1 <- do.call(rbind, mclapply(seq_len(nrow(x)), function(k) {
  sd <- x$seed[k]; i <- floor((sd - 1) / 1e4); cl <- cells[i, ]
  set.seed(sd); E <- rstar(cl$n, theta_model(cl$chibar, u1 = cl$u1, u2 = cl$u2))$E
  cbind(chibar0 = cl$chibar, n = cl$n, u1 = cl$u1, seed = sd, truth = x$truth[k], one(E, 0.99, TM[[i]], sd)) }, mc.cores = 2))
saveRDS(r1, "results/htsens_s1.rds"); cat("s1 done\n")
## Study 4: chibar = 0.6, n = 4000, q = 0.999
s4 <- readRDS("../main/results/s4_all.rds"); f4 <- s4[s4$variant == "full", ]; th0 <- theta_model(0.6); TM0 <- true_margins(th0)
t3 <- -log(1 - 0.999); tr4 <- star:::joint_surv_cpp(t3, t3, gl, gw, th0$phi, th0$alpha, th0$b, th0$u1, th0$u2, th0$p, th0$chibar, th0$z) / 0.001
r4 <- do.call(rbind, mclapply(f4$seed, function(sd) { set.seed(sd); E <- rstar(4000, th0)$E
  cbind(seed = sd, truth = tr4, one(E, 0.999, TM0, sd)) }, mc.cores = 2))
saveRDS(r4, "results/htsens_s4.rds"); cat("s4 done\n")
