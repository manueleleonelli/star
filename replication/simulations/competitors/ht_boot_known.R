## bootstrap coverage of HT with the true margins, conditioning threshold 0.9, within the model (Sections S3.2 and S3.6 datasets)
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
source("ht_sens_fun.R")
B <- 100
bootF <- function(E, qs, Ft, seed, nsim = 10000) {
  set.seed(seed); n <- nrow(E); st <- ht_chi_F(E, qs, Ft$F1, Ft$F2, 0.9, fits = TRUE)
  b <- replicate(B, ht_chi_F(E[sample.int(n, n, replace = TRUE), , drop = FALSE], qs, Ft$F1, Ft$F2, 0.9, nsim = nsim, starts = st))
  b <- matrix(b, nrow = length(qs)); cbind(lo = apply(b, 1, quantile, 0.05, na.rm = TRUE), hi = apply(b, 1, quantile, 0.95, na.rm = TRUE))
}
s1 <- readRDS("../main/results/s1_all.rds"); x <- s1[s1$quantity == "chi_q", ]
cells <- expand.grid(chibar = c(0.15, 0.45, 0.75, 0.95), n = c(1000, 4000)); cells$u1 <- 2.2; cells$u2 <- 2.2
cells <- rbind(cells, data.frame(chibar = 0.45, n = 4000, u1 = 2.8, u2 = 1.7))
TM <- lapply(seq_len(nrow(cells)), function(i) true_margins(theta_model(cells$chibar[i], u1 = cells$u1[i], u2 = cells$u2[i])))
r1 <- do.call(rbind, mclapply(seq_len(nrow(x)), function(k) {
  sd <- x$seed[k]; i <- floor((sd - 1) / 1e4); cl <- cells[i, ]
  set.seed(sd); E <- rstar(cl$n, theta_model(cl$chibar, u1 = cl$u1, u2 = cl$u2))$E
  ci <- bootF(E, 0.99, TM[[i]], sd + 4e6); tr <- x$truth[k]
  data.frame(chibar0 = cl$chibar, n = cl$n, u1 = cl$u1, seed = sd, cover = as.numeric(tr >= ci[1] & tr <= ci[2])) }, mc.cores = 2))
saveRDS(r1, "results/htkboot_s1.rds"); cat("boot s1 done\n")
s4 <- readRDS("../main/results/s4_all.rds"); f4 <- s4[s4$variant == "full", ]; th0 <- theta_model(0.6); TM0 <- true_margins(th0)
t3 <- -log(1 - 0.999); tr4 <- star:::joint_surv_cpp(t3, t3, gl, gw, th0$phi, th0$alpha, th0$b, th0$u1, th0$u2, th0$p, th0$chibar, th0$z) / 0.001
r4 <- do.call(rbind, mclapply(f4$seed, function(sd) { set.seed(sd); E <- rstar(4000, th0)$E
  ci <- bootF(E, 0.999, TM0, sd + 4e6); data.frame(seed = sd, cover = as.numeric(tr4 >= ci[1] & tr4 <= ci[2])) }, mc.cores = 2))
saveRDS(r4, "results/htkboot_s4.rds"); cat("boot s4 done\n")
