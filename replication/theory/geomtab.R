## geomtab.R -- Table S1 of the supplement: eta_C of the ratio scale function and of the additive-pairwise gauge (11),
## with the envelope check on a grid of the simplex. Writes geomtab.csv.
source(system.file("trivariate", "star3.R", package = "star"))
W <- ETA_GRID; mn <- function(C) do.call(pmin, lapply(C, function(j) W[, j]))
eta_of <- function(sig) c(eta12 = max(mn(1:2) * sig), eta13 = max(mn(c(1,3)) * sig), eta23 = max(mn(2:3) * sig), eta123 = max(mn(1:3) * sig))
adm <- function(sig) all(1/sig >= apply(W, 1, max) - 1e-9)
naive <- function(chib, nu = 0.5) 1 + chib * (apply(W, 1, min) / apply(W, 1, max))^(1 + nu)
rows <- list()
for (cb in c(0.5, 0.9)) { s <- naive(cb); rows[[length(rows)+1]] <- c(family = "naive", chib = sprintf("%.1f, all pairs", cb), round(eta_of(s), 3), adm = adm(s)) }
for (ch in list(c(0.6,0,0), c(0.6,0.2,0.2), c(0.5,0.4,0.4), c(0.5,0.5,0.5), c(0.55,0.55,0.55), c(1,0,0))) {
  s <- sigma3(W, ch, rep(0.5,3)); rows[[length(rows)+1]] <- c(family = "additive", chib = paste(ch, collapse = ", "), round(eta_of(s), 3), adm = adm(s)) }
tab <- do.call(rbind, rows); print(tab, quote = FALSE)
cat("bivariate reference (1+chibar)/2: 0.5->0.75, 0.6->0.80, 0.9->0.95, 1->1\n")
cat("sigma at barycentre, additive 0.5 all pairs:", sigma3(matrix(1/3,1,3), rep(.5,3), rep(.5,3)), "\n")
write.csv(tab, "geomtab.csv", row.names = FALSE)
