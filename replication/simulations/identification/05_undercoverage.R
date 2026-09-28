# 05_undercoverage.R -- Study 0b: resolving the diagonal undercoverage at
# n = 8000. Two candidate causes, crossed: chain length (mixing) and the
# prior scale v_p on log p (regularization of the right-skewed exponent).
# Also reports batch-means effective sample sizes for the threshold block.
# Runtime at defaults: roughly 3-4 hours on 8 cores.

source("00_utils.R")
library(parallel)
set.seed(505)
CORES <- as.integer(Sys.getenv("STAR_CORES", max(1, detectCores() - 1)))
SCALE <- Sys.getenv("STAR_SCALE", "beach")
R <- if (SCALE == "test") 2 else 30
N <- if (SCALE == "test") 600 else 8000
ITERS <- if (SCALE == "test") c(short = 1200, long = 2400) else c(short = 8000, long = 30000)

th0 <- theta_true(chibar = 0.6, p = 2)
rho_half_true <- rho_star(0.5, th0$u1, th0$u2, th0$p)

ess_bm <- function(x) {                        # batch-means ESS
  S <- length(x); B <- floor(sqrt(S)); nb <- floor(S / B)
  bm <- colMeans(matrix(x[1:(B * nb)], nrow = B))
  tau <- max(1, B * var(bm) / var(x)); S / tau
}

run_par2 <- function(X, FUN) {
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(CORES); on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, { library(star); source("00_utils.R") })
    clusterExport(cl, ls(globalenv()), envir = globalenv())
    parLapply(cl, X, FUN)
  } else mclapply(X, FUN, mc.cores = CORES)
}

cells <- expand.grid(iter = names(ITERS), v_p = c(1.5, 0.75),
                     stringsAsFactors = FALSE)
dir.create("results05", showWarnings = FALSE)
out <- list()
for (i in seq_len(nrow(cells))) {
  it <- ITERS[cells$iter[i]]; vp <- cells$v_p[i]
  f <- sprintf("results05/cell_%s_vp%.2f.rds", cells$iter[i], vp)
  if (file.exists(f)) { out[[i]] <- readRDS(f); next }
  cat(sprintf("cell %d/4: iter = %d, v_p = %.2f (%s)\n", i, it, vp,
              format(Sys.time(), "%H:%M")))
  res <- run_par2(seq_len(R), function(rr) {
    set.seed(5e4 * i + rr)
    dat <- rstar(N, th0)
    fit <- star_mcmc(dat$E, prior = star_prior("beta", v_p = vp),
                      n_iter = it, burnin = round(it * 0.4), thin = 5,
                      verbose = FALSE)
    D <- fit$draws
    rho_half <- sapply(seq_len(nrow(D)), function(s) {
      th <- star:::theta_row(fit, s); rho_star(0.5, th$u1, th$u2, th$p) })
    ci_r <- quantile(rho_half, c(.05, .95)); ci_p <- quantile(D[, "p"], c(.05, .95))
    data.frame(iter = it, v_p = vp, seed = rr,
               rho_cover = as.numeric(rho_half_true >= ci_r[1] & rho_half_true <= ci_r[2]),
               rho_w = diff(ci_r), p_cover = as.numeric(2 >= ci_p[1] & 2 <= ci_p[2]),
               p_med = median(D[, "p"]),
               ess_logp = ess_bm(log(D[, "p"])), ess_u1 = ess_bm(D[, "u1"]))
  })
  bad <- !vapply(res, is.data.frame, logical(1))
  if (any(bad)) stop("worker error: ",
                     conditionMessage(attr(res[[which(bad)[1]]], "condition")))
  out[[i]] <- do.call(rbind, res); saveRDS(out[[i]], f)
}
res <- do.call(rbind, out); saveRDS(res, "results05/undercoverage.rds")

cat("\n===== Study 0b: undercoverage resolution (truth rho(1/2) =",
    round(rho_half_true, 3), ") =====\n")
print(aggregate(cbind(rho_cover, rho_w, p_cover, p_med, ess_logp, ess_u1)
                ~ iter + v_p, res, mean), digits = 3, row.names = FALSE)
cat("\nreading: if longer chains alone restore rho_cover to ~0.9 the cause was\n",
    "mixing (check ess_logp: values below ~200 at iter = 8000 confirm it);\n",
    "if only v_p = 0.75 restores it, the cause is the right skew of p and the\n",
    "tighter default is the one to adopt; if both help, both matter.\n")
