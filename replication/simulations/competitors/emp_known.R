## Empirical chi(q) on the data scale, n^{-1} sum 1{E1 > t_q, E2 > t_q} / (1 - q), t_q = -log(1 - q),
## with a 90% percentile bootstrap interval (B = 300), on the exact datasets of Studies 1, 3 and 4.
suppressMessages(library(star))
theta_model <- function(chibar = 0.6, p = 2, u1 = 2.2, u2 = 2.2, nubar = 0.5)
  list(phi = c(0.1, 0.2, 0.4, 0.2, 0.1), alpha = 2, b = c(0, 0, 0),
       u1 = u1, u2 = u2, p = p, chibar = chibar, z = rep(log(nubar), 3))
emp_ci <- function(E, q, B = 300) {
  x <- E[, 1] > -log(1 - q) & E[, 2] > -log(1 - q); n <- length(x)
  b <- replicate(B, mean(x[sample.int(n, n, replace = TRUE)])) / (1 - q)
  c(est = mean(x) / (1 - q), quantile(b, c(0.05, 0.95), names = FALSE))
}
gl <- seq(0.5/201, 1 - 0.5/201, length.out = 201)
## Study 1: within the model, q = 0.99
s1 <- readRDS("../main/results/s1_all.rds"); x <- s1[s1$quantity == "chi_q", ]
cells <- expand.grid(chibar = c(0.15, 0.45, 0.75, 0.95), n = c(1000, 4000)); cells$u1 <- 2.2; cells$u2 <- 2.2
cells <- rbind(cells, data.frame(chibar = 0.45, n = 4000, u1 = 2.8, u2 = 1.7))
r1 <- do.call(rbind, lapply(seq_len(nrow(x)), function(k) {
  sd <- x$seed[k]; i <- floor((sd - 1) / 1e4); cl <- cells[i, ]
  stopifnot(abs(cl$chibar - x$chibar0[k]) < 1e-9, cl$n == x$n[k], cl$u1 == x$u1[k])
  set.seed(sd); E <- rstar(cl$n, theta_model(cl$chibar, u1 = cl$u1, u2 = cl$u2))$E
  set.seed(sd + 1e6); e <- emp_ci(E, 0.99); tr <- x$truth[k]
  data.frame(chibar0 = cl$chibar, n = cl$n, u1 = cl$u1, seed = sd, truth = tr,
             model_err = x$err[k], model_cover = x$cover[k],
             emp_err = e[1] - tr, emp_cover = as.numeric(tr >= e[2] & tr <= e[3]))
}))
## Study 4: within the model, chibar = 0.6, n = 4000, q = 0.999 (estimated-threshold fits)
s4 <- readRDS("../main/results/s4_all.rds"); f4 <- s4[s4$variant == "full", ]
th0 <- theta_model(0.6); t3 <- -log(1 - 0.999)
tr4 <- star:::joint_surv_cpp(t3, t3, gl, rep(1/201, 201), th0$phi, th0$alpha, th0$b,
                              th0$u1, th0$u2, th0$p, th0$chibar, th0$z) / 0.001
r4 <- do.call(rbind, lapply(f4$seed, function(sd) {
  set.seed(sd); E <- rstar(4000, th0)$E; set.seed(sd + 1e6); e <- emp_ci(E, 0.999)
  data.frame(seed = sd, emp_err = e[1] - tr4, emp_cover = as.numeric(tr4 >= e[2] & tr4 <= e[3]))
}))
r4 <- merge(r4, data.frame(seed = f4$seed, model_err = f4$chi999_med - tr4, model_cover = f4$chi999_cover))
saveRDS(list(s1 = r1, s4 = r4, truth4 = tr4), "results/emp_known.rds")
rmse <- function(v) sqrt(mean(v^2))
cat("== Study 1, q = 0.99 ==\n")
print(merge(aggregate(cbind(model_err, emp_err) ~ chibar0 + n + u1, r1, rmse),
            aggregate(cbind(model_cover, emp_cover) ~ chibar0 + n + u1, r1, mean)), digits = 4)
cat("mean errors\n"); print(aggregate(cbind(model_err, emp_err) ~ chibar0 + n + u1, r1, mean), digits = 3)
cat("== Study 4, q = 0.999, truth", tr4, "==\n")
print(c(model_rmse = rmse(r4$model_err), emp_rmse = rmse(r4$emp_err), model_bias = mean(r4$model_err),
        emp_bias = mean(r4$emp_err), model_cover = mean(r4$model_cover), emp_cover = mean(r4$emp_cover)), digits = 3)
