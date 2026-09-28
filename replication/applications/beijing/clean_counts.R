## clean_counts.R -- hours removed by the cross-station cleaning rule, per station and pollutant (Section S5.1)
dir <- "../../data/beijing"
fs <- list.files(dir, pattern = "\\.csv\\.xz$"); L <- lapply(file.path(dir, fs), read.csv); names(L) <- sub("\\.csv\\.xz$", "", fs)
key <- function(d) paste(d$year, d$month, d$day, d$hour); k0 <- key(L[[1]])
grab <- function(v) { M <- sapply(L, function(d) d[[v]][match(k0, key(d))]); colnames(M) <- names(L); M }
date <- as.Date(sprintf("%d-%02d-%02d", L[[1]]$year, L[[1]]$month, L[[1]]$day)); warm <- as.integer(format(date, "%m")) %in% 4:9
for (v in c("O3", "PM2.5")) { M <- grab(v)
  for (j in seq_len(ncol(M))) { om <- suppressWarnings(apply(M[, -j, drop = FALSE], 1, max, na.rm = TRUE))
    bad <- !is.na(M[, j]) & M[, j] > 300 & is.finite(om) & M[, j] > 3 * om
    if (colnames(M)[j] %in% c("Wanliu", "Dongsi") || sum(bad) > 0) cat(sprintf("%-6s %-14s flagged %3d (warm %3d)  max flagged %s  dates %s\n", v, colnames(M)[j], sum(bad), sum(bad & warm),
      if (any(bad)) round(max(M[bad, j])) else "-", if (any(bad)) paste(unique(format(date[bad], "%Y-%m")), collapse = ",") else "-")) } }
## raw daily maxima at Wanliu without cleaning, for comparison
x <- readRDS("results/daily12.rds")[["Wanliu"]]; cat("Wanliu daily n", nrow(x), " max PM", max(x$PM2.5), " max O3", max(x$O3), "\n")
