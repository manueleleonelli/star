# 02_contraction.R -- Component B: posterior contraction rates
# CI width and posterior sd of parameters and functionals versus n; the
# log-log slope against n identifies the contraction rate (-1/2 = regular).
# The settings below are those of the run reported in the supplement (Section S4.1, Table S3).

source("00_utils.R")
library(parallel)
set.seed(202)

REPS <- 30
NS <- c(500, 2000, 8000)
ITER <- 8000; BURN <- 3000; THIN <- 5
CORES <- getOption("mc.cores", 8)

th0 <- theta_true(chibar = 0.6, p = 2)
truth <- functionals(th0)

one_rep <- function(seed, n) {
  set.seed(seed)
  dat <- rstar(n, th0)
  fit <- star_mcmc(dat$E, prior = star_prior("beta"), n_iter = ITER,
                    burnin = BURN, thin = THIN, verbose = FALSE)
  cbind(summarize_rep(fit, truth), n = n, seed = seed)
}

grid <- expand.grid(seed = seq_len(REPS), n = NS)
cat(sprintf("Component B: %d fits (%d replicates x %d sample sizes); progress by chunk...\n", nrow(grid), REPS, length(NS)))
res <- run_par(seq_len(nrow(grid)),
               function(i) one_rep(grid$seed[i] + 1000 * grid$n[i], grid$n[i]),
               cores = CORES)
res <- do.call(rbind, res)
saveRDS(res, "contraction.rds")

## contraction slopes: log median width ~ log n, per quantity
agg <- aggregate(width ~ quantity + n, res, median)
cat("\ncontraction slopes (regular identification: about -0.5):\n")
for (q in unique(agg$quantity)) {
  a <- agg[agg$quantity == q, ]
  s <- coef(lm(log(width) ~ log(n), a))[2]
  cat(sprintf("  %-10s slope %5.2f   widths:", q, s))
  cat(sprintf(" %.3g", a$width), "\n")
}
cat("\ncoverage of 90% intervals by n:\n")
cv <- aggregate(cover ~ quantity + n, res, mean)
print(round(xtabs(cover ~ quantity + n, cv), 2))
