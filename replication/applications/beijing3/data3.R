## Trivariate Beijing data: warm-season daily maxima of all six pollutants at Wanliu, with the
## cross-station sensor-fault cleaning of the main paper applied to each pollutant.
dir <- "../../data/beijing"
fs <- list.files(dir, pattern = "\\.csv\\.xz$"); L <- lapply(file.path(dir, fs), read.csv); names(L) <- sub("\\.csv\\.xz$", "", fs)
key <- function(d) paste(d$year, d$month, d$day, d$hour); k0 <- key(L[[1]])
grab <- function(v) { M <- sapply(L, function(d) d[[v]][match(k0, key(d))]); colnames(M) <- names(L); M }
clean <- function(M, lev, fac = 3) {
  for (j in seq_len(ncol(M))) { om <- suppressWarnings(apply(M[, -j, drop = FALSE], 1, max, na.rm = TRUE))
    bad <- !is.na(M[, j]) & M[, j] > lev & is.finite(om) & M[, j] > fac * om
    M[bad, j] <- NA }
  M }
lev <- c(PM2.5 = 300, PM10 = 300, SO2 = 300, NO2 = 300, CO = 300 * 10, O3 = 300)   # CO is reported in larger units
date <- as.Date(sprintf("%d-%02d-%02d", L[[1]]$year, L[[1]]$month, L[[1]]$day))
st <- "Wanliu"
daily <- NULL
for (v in names(lev)) { M <- clean(grab(v), lev[v]); a <- aggregate(list(x = M[, st]), list(date = date), max, na.rm = TRUE)
  names(a)[2] <- v; daily <- if (is.null(daily)) a else merge(daily, a, by = "date") }
daily <- daily[as.integer(format(daily$date, "%m")) %in% 4:9, ]
ok <- apply(daily[, -1], 1, function(r) all(is.finite(r)))
cat("warm-season days:", nrow(daily), " complete on all six:", sum(ok), " complete on PM2.5,O3:", sum(is.finite(daily$PM2.5) & is.finite(daily$O3)), "\n")
d <- daily[ok, ]
print(round(cor(d[, -1], method = "spearman"), 2))
chi_emp <- function(x, y, q) { u <- quantile(x, q); v <- quantile(y, q); mean(x > u & y > v) / (1 - q) }
cat("\nempirical chi(0.9), chi(0.95) by pair:\n")
nm <- names(lev)
for (i in 1:5) for (j in (i+1):6) cat(sprintf("%-6s %-6s  %.2f  %.2f\n", nm[i], nm[j], chi_emp(d[[nm[i]]], d[[nm[j]]], .9), chi_emp(d[[nm[i]]], d[[nm[j]]], .95)))
saveRDS(daily, "results/daily6_wanliu.rds")
