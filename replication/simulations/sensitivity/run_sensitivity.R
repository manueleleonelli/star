## run_sensitivity.R -- all sensitivity studies in one file, no arguments needed.
##   Needs the star package. Run from this folder: source("run_sensitivity.R").
## Checkpointed per study (results/*.rds); re-sourcing skips finished studies.
## Knobs: NSEED (paired seeds per setting) and ITER. Defaults are supplement-grade
## and take roughly 1.5-2 hours sequentially at NSEED = 5.

suppressMessages(library(star))
NSEED <- 5; ITER <- 10000; BURN <- 3000; THIN <- 4; N <- 1000
dir.create("results", showWarnings = FALSE)

theta0 <- function(chibar = 0.6) list(phi = c(0.1,0.2,0.4,0.2,0.1), alpha = 2, b = c(0,0,0),
                                      u1 = 2.2, u2 = 2.2, p = 2, chibar = chibar, z = rep(log(0.5),3))
summ <- function(fit, E) {
  D <- fit$draws; r <- rowSums(E); w <- E[,1]/r; idx <- round(seq(1, nrow(D), length.out = 150))
  ex <- mean(sapply(idx, function(s){ th <- star:::theta_row(fit,s); mean(r > rho_star(w,th$u1,th$u2,th$p)) }))
  rh <- median(sapply(idx, function(s){ th <- star:::theta_row(fit,s); rho_star(0.5,th$u1,th$u2,th$p) }))
  tq <- -log(0.01); gl <- seq(0.5/201, 1-0.5/201, length.out = 201)
  cq <- mean(sapply(idx, function(s){ th <- star:::theta_row(fit,s)
    star:::joint_surv_cpp(tq,tq,gl,rep(1/201,201),th$phi,th$alpha,th$b,th$u1,th$u2,th$p,th$chibar,th$z)/0.01 }))
  data.frame(chibar = median(D[,"chibar"]), chibar_w = unname(diff(quantile(D[,"chibar"], c(.05,.95)))),
             exc = ex, rho_half = rh, p = median(D[,"p"]), chi99 = cq, alpha = median(D[,"alpha"]))
}
fit_beta <- function(E, ...) star_mcmc(E, prior = star_prior("beta", ...), n_iter = ITER, burnin = BURN, thin = THIN, verbose = FALSE)

## ---------------- 1. kappa0 and half-normal scale (paired) ----------------
if (!file.exists("results/kappa.rds")) {
  cells <- list(k07 = c(0.7,1), k08 = c(0.8,1), k09 = c(0.9,1), s05 = c(0.8,0.5), s2 = c(0.8,2)); out <- NULL
  for (s in seq_len(NSEED)) { set.seed(100 + s); dat <- rstar(N, theta0())
    for (nm in names(cells)) { k <- cells[[nm]]
      out <- rbind(out, data.frame(seed = s, cell = nm, summ(fit_beta(dat$E, kappa0 = k[1], s_u = k[2]), dat$E)))
      cat(sprintf("[1 kappa] seed %d %s done (%s)\n", s, nm, format(Sys.time(), "%H:%M"))) } }
  saveRDS(out, "results/kappa.rds")
}
## ---------------- 2. class prior: spike weights and slab (paired) ----------------
if (!file.exists("results/prior.rds")) {
  priors <- list(spike10 = c(0.1,0.1,1,1), default = c(0.25,0.25,1,1), spike40 = c(0.4,0.4,1,1), slab22 = c(0.25,0.25,2,2)); out <- NULL
  for (cb in c(0.5, 0.95)) for (s in seq_len(NSEED)) { set.seed(200 + 10*s + round(100*cb)); dat <- rstar(N, theta0(cb))
    for (nm in names(priors)) { pr <- priors[[nm]]
      fit <- star_mcmc(dat$E, prior = star_prior("threecomp", pi0 = pr[1], pi1 = pr[2], a_chi = pr[3], b_chi = pr[4]),
                        n_iter = round(ITER*1.25), burnin = BURN, thin = THIN, verbose = FALSE)
      cp <- class_probs(fit)
      out <- rbind(out, data.frame(chibar0 = cb, seed = s, prior = nm, P0 = cp[[1]], Pslab = cp[[2]], P1 = cp[[3]], chibar_mean = mean(fit$draws[,"chibar"])))
      cat(sprintf("[2 prior] chibar %.2f seed %d %s: P1 = %.2f (%s)\n", cb, s, nm, cp[[3]], format(Sys.time(), "%H:%M"))) } }
  saveRDS(out, "results/prior.rds")
}
## ---------------- 3. basis orders (paired) ----------------
if (!file.exists("results/basis.rds")) {
  cfgs <- list(small = c(3,1,1), default = c(5,2,2), large = c(8,4,4)); out <- NULL
  for (s in seq_len(NSEED)) { set.seed(400 + s); dat <- rstar(N, theta0())
    for (nm in names(cfgs)) { K <- cfgs[[nm]]
      fit <- star_mcmc(dat$E, prior = star_prior("beta"), K_W = K[1], K_lam = K[2], K_nu = K[3], n_iter = ITER, burnin = BURN, thin = THIN, verbose = FALSE)
      out <- rbind(out, data.frame(seed = s, basis = nm, summ(fit, dat$E))); cat(sprintf("[3 basis] seed %d %s done\n", s, nm)) } }
  saveRDS(out, "results/basis.rds")
}
## ---------------- 4. angle-dependent body shape: is a constant alpha enough? ----------------
if (!file.exists("results/alphaw.rds")) {
  rb_alphaw <- function(n, th, af) {
    K <- length(th$phi); comp <- sample.int(K, n, TRUE, th$phi); w <- rbeta(n, comp, K-comp+1)
    lam <- lambda_star(w, th$b); rho <- rho_star(w, th$u1, th$u2, th$p); sig <- sigma_star(w, th$chibar, th$z); al <- af(w)
    pi_w <- pgamma(rho, shape = al, rate = lam); body <- runif(n) < pi_w; r <- numeric(n)
    r[body] <- qgamma(runif(sum(body))*pi_w[body], shape = al[body], rate = lam[body])
    r[!body] <- rho[!body] + rexp(sum(!body))*sig[!body]; cbind(r*w, r*(1-w)) }
  af <- function(w) 1.5 + 1.5*w; th <- theta0(0.6)
  set.seed(7); Eb <- rb_alphaw(2e6, th, af); tq <- -log(0.01); rb <- rowSums(Eb)
  truth <- c(chi99 = mean(Eb[,1] > tq & Eb[,2] > tq)/0.01, exc = mean(rb > rho_star(Eb[,1]/rb, th$u1, th$u2, th$p)), rho_half = rho_star(0.5, th$u1, th$u2, th$p))
  out <- NULL
  for (s in seq_len(NSEED)) { set.seed(300 + s); E <- rb_alphaw(N, th, af)
    out <- rbind(out, data.frame(seed = s, summ(fit_beta(E), E))); cat(sprintf("[4 alpha(w)] seed %d done\n", s)) }
  saveRDS(list(fits = out, truth = truth), "results/alphaw.rds")
}
## ---------------- 5. small-n ridge variant ----------------
if (!file.exists("results/ridge300.rds")) {
  out <- NULL
  for (s in seq_len(NSEED)) { set.seed(303 + s); dat <- rstar(300, theta0(0.6))
    for (nm in c("anchor","weak")) { pri <- if (nm == "anchor") star_prior("beta") else star_prior("beta", kappa0 = 0.5, s_u = 5)
      fit <- star_mcmc(dat$E, prior = pri, n_iter = ITER, burnin = BURN, thin = THIN, verbose = FALSE)
      D <- with(as.data.frame(fit$draws), cbind(log_u1 = log(u1), log_u2 = log(u2), log_p = log(p),
                logit_chibar = qlogis(pmin(pmax(chibar,1e-6),1-1e-6)), log_alpha = log(alpha), b0 = b0))
      C <- cor(D); C[upper.tri(C, diag = TRUE)] <- 0; wm <- which(abs(C) == max(abs(C)), arr.ind = TRUE)[1,]
      r <- rowSums(dat$E); ex <- mean(sapply(round(seq(1,nrow(fit$draws),length.out=100)), function(k){ th2 <- star:::theta_row(fit,k); mean(r > rho_star(dat$E[,1]/r, th2$u1, th2$u2, th2$p)) }))
      out <- rbind(out, data.frame(seed = s, prior = nm, max_corr = max(abs(C)), pair = paste(colnames(D)[wm], collapse = "-"),
                   sd_log_p = sd(D[,"log_p"]), sd_log_u1 = sd(D[,"log_u1"]), exc = ex))
      cat(sprintf("[5 ridge n=300] seed %d %s done\n", s, nm)) } }
  saveRDS(out, "results/ridge300.rds")
}

## ======================= consolidated summary =======================
cat("\n\n================= SENSITIVITY STUDIES: SUMMARY (n = 1000, NSEED =", NSEED, ") =================\n")
k <- readRDS("results/kappa.rds"); k$cell <- factor(k$cell, levels = c("k07","k08","k09","s05","s2"))
cat("\n--- 1. kappa0 / half-normal scale (truth: chibar 0.6, rho(1/2) 3.11) ---\n")
a <- aggregate(cbind(chibar, chibar_w, exc, rho_half, p, chi99) ~ cell, k, mean)
a$setting <- c("kappa0=0.7","kappa0=0.8 default","kappa0=0.9","s_u=0.5","s_u=2")[as.integer(a$cell)]
print(a[, c("setting","chibar","chibar_w","exc","rho_half","p","chi99")], digits = 3, row.names = FALSE)
cat("paired shifts vs default (mean [sd] of within-seed differences):\n")
ref <- k[k$cell == "k08",]
for (cc in c("k07","k09","s05","s2")) { m <- merge(k[k$cell == cc,], ref, by = "seed")
  cat(sprintf("  %-4s d.chibar %+.3f [%.3f]  d.exc %+.3f [%.3f]  d.rho_half %+.2f [%.2f]  d.chi99 %+.4f [%.4f]\n", cc,
    mean(m$chibar.x-m$chibar.y), sd(m$chibar.x-m$chibar.y), mean(m$exc.x-m$exc.y), sd(m$exc.x-m$exc.y),
    mean(m$rho_half.x-m$rho_half.y), sd(m$rho_half.x-m$rho_half.y), mean(m$chi99.x-m$chi99.y), sd(m$chi99.x-m$chi99.y))) }

p <- readRDS("results/prior.rds"); p$prior <- factor(p$prior, levels = c("spike10","default","spike40","slab22"))
cat("\n--- 2. class prior: P(AD | data) by spike weights and slab (truths chibar 0.5 interior, 0.95 near boundary) ---\n")
a <- aggregate(cbind(P0, Pslab, P1, chibar_mean) ~ chibar0 + prior, p, mean); a$PAI <- a$P0 + a$Pslab
print(a[order(a$chibar0, a$prior), c("chibar0","prior","P0","Pslab","P1","PAI","chibar_mean")], digits = 2, row.names = FALSE)

b <- readRDS("results/basis.rds"); b$basis <- factor(b$basis, levels = c("small","default","large"))
cat("\n--- 3. basis orders (K_W, K_lam, K_nu) = (3,1,1) / (5,2,2) / (8,4,4) ---\n")
print(aggregate(cbind(chibar, chibar_w, exc, rho_half, p, chi99) ~ basis, b, mean), digits = 3, row.names = FALSE)

al <- readRDS("results/alphaw.rds")
cat("\n--- 4. angle-dependent body shape alpha(w) = 1.5 + 1.5 w, fitted with constant alpha ---\n")
cat(sprintf("truth: chibar 0.6, chi(0.99) %.4f, exceedance %.3f, rho(1/2) %.2f\n", al$truth["chi99"], al$truth["exc"], al$truth["rho_half"]))
print(round(colMeans(al$fits[, -1]), 3))

rg <- readRDS("results/ridge300.rds")
cat("\n--- 5. ridge variant at n = 300: anchoring vs weak prior ---\n")
print(aggregate(cbind(max_corr, sd_log_p, sd_log_u1, exc) ~ prior, rg, mean), digits = 3, row.names = FALSE)
cat("most correlated pair (anchor):", names(sort(table(rg$pair[rg$prior=="anchor"]), decreasing = TRUE))[1],
    "| (weak):", names(sort(table(rg$pair[rg$prior=="weak"]), decreasing = TRUE))[1], "\n")
cat("\nAll five studies complete; results in results/*.rds\n")
