## run_ht.R -- conditional extremes model (HT) on the datasets of Studies 1, 3 and 4, regenerated from their seeds:
## chi(q) with estimated and with true margins, and bootstrap intervals (Tables 2, 3 and S7). Writes chunks to htres/;
## then run assemble.R. Run from this folder.
library(parallel); if (.Platform$OS.type == "windows") mclapply <- function(X, FUN, ..., mc.cores = 1, mc.preschedule = TRUE) lapply(X, FUN, ...)  # no forking on Windows
suppressMessages(library(star)); library(parallel)
source("ht.R"); source("gens.R")
theta_model <- function(chibar = 0.6, p = 2, u1 = 2.2, u2 = 2.2, nubar = 0.5)
  list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0), u1 = u1, u2 = u2, p = p, chibar = chibar, z = rep(log(nubar), 3))
ht_chi_known <- function(E, qs, kappa = 0.9, nsim = 20000) {
  y1 <- lap_q(1 - exp(-E[, 1])); y2 <- lap_q(1 - exp(-E[, 2])); u <- lap_q(kappa)
  f12 <- ht_fit(y1, y2, u); f21 <- ht_fit(y2, y1, u)
  sapply(qs, function(q) { v <- lap_q(q); 0.5 * (ht_joint(f12, v, v, nsim) + ht_joint(f21, v, v, nsim)) / (1 - q) })
}
B <- 100; QS <- c(0.95, 0.99, 0.999); dir.create("htres", showWarnings = FALSE)
run_chunks <- function(tag, N, FUN, chunk = 10) {
  for (k in seq_len(ceiling(N / chunk))) {
    f <- sprintf("htres/%s_%03d.rds", tag, k); if (file.exists(f)) next
    idx <- ((k - 1) * chunk + 1):min(k * chunk, N)
    res <- mclapply(idx, function(i) tryCatch(FUN(i), error = function(e) paste("ERR:", conditionMessage(e))), mc.cores = 2)
    saveRDS(res, f); cat(tag, k, format(Sys.time(), "%H:%M:%S"), "\n")
  }
}
## Study 3
s3 <- readRDS("../main/results/s3_all.rds"); tr <- readRDS("../main/results/s3_truths.rds")
gens <- list(gauss05 = function(n) sim_gauss(n, 0.5), invlog05 = function(n) sim_invlogistic(n, 0.5), frank5 = function(n) sim_frank(n, 5), t05_4 = function(n) sim_t(n, 0.5, 4))
run_chunks("s3", nrow(s3), function(i) {
  g <- as.character(s3$gen[i]); sd <- s3$seed[i]; set.seed(sd); E <- gens[[g]](s3$n[i])
  set.seed(sd + 2e6); est <- ht_chi(E, QS); set.seed(sd + 3e6); kn <- ht_chi_known(E, QS)
  ci <- ht_boot(E, QS, B, seed = sd + 4e6); truth <- tr[[g]]$chi
  data.frame(gen = g, seed = sd, q = QS, truth = truth, ht_err = est - truth, ht_cover = as.numeric(truth >= ci[, 1] & truth <= ci[, 2]), htk_err = kn - truth)
})
## Study 4
s4 <- readRDS("../main/results/s4_all.rds"); f4 <- s4[s4$variant == "full", ]; th0 <- theta_model(0.6)
gl <- seq(0.5/201, 1 - 0.5/201, length.out = 201); t3 <- -log(1 - 0.999)
tr4 <- star:::joint_surv_cpp(t3, t3, gl, rep(1/201, 201), th0$phi, th0$alpha, th0$b, th0$u1, th0$u2, th0$p, th0$chibar, th0$z) / 0.001
run_chunks("s4", nrow(f4), function(i) {
  sd <- f4$seed[i]; set.seed(sd); E <- rstar(4000, th0)$E; set.seed(sd + 2e6); est <- ht_chi(E, 0.999); ci <- ht_boot(E, 0.999, B, seed = sd + 4e6)
  data.frame(seed = sd, truth = tr4, ht_err = est - tr4, ht_cover = as.numeric(tr4 >= ci[1] & tr4 <= ci[2]))
})
## Study 1
s1 <- readRDS("../main/results/s1_all.rds"); x <- s1[s1$quantity == "chi_q", ]
cells <- expand.grid(chibar = c(0.15, 0.45, 0.75, 0.95), n = c(1000, 4000)); cells$u1 <- 2.2; cells$u2 <- 2.2
cells <- rbind(cells, data.frame(chibar = 0.45, n = 4000, u1 = 2.8, u2 = 1.7))
run_chunks("s1", nrow(x), function(k) {
  sd <- x$seed[k]; i <- floor((sd - 1) / 1e4); cl <- cells[i, ]
  set.seed(sd); E <- rstar(cl$n, theta_model(cl$chibar, u1 = cl$u1, u2 = cl$u2))$E
  set.seed(sd + 2e6); est <- ht_chi(E, 0.99); ci <- ht_boot(E, 0.99, B, seed = sd + 4e6); tr <- x$truth[k]
  data.frame(chibar0 = cl$chibar, n = cl$n, u1 = cl$u1, seed = sd, truth = tr, ht_err = est - tr, ht_cover = as.numeric(tr >= ci[1] & tr <= ci[2]))
})
cat("ALL DONE\n")
