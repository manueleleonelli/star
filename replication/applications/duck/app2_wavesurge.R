## app2_wavesurge.R -- wave height and storm surge at Duck, North Carolina, from
## public NOAA sources: NDBC buoy 44100 (waves) and CO-OPS gauge 8651370 (water level).
## NEEDS INTERNET. Surge = observed hourly water level minus astronomical prediction.
## Output: results/wavesurge.rds. The pair is expected near the asymptotic dependence
## boundary; watch the exceedance fraction and rho(1/2) (Section 4.4 diagnostic).
source("app_common.R"); dir.create("results", showWarnings = FALSE)
ITER <- as.integer(Sys.getenv("STAR_ITER", "30000")); BURN <- round(ITER/3)
years <- 2010:2023; buoy <- "44100"; gauge <- "8651370"
## resumable: download and fit are skipped when results/wavesurge.rds exists (tables_duck.R reads it)
if (file.exists("results/wavesurge.rds")) cat("results/wavesurge.rds exists, nothing to do; run tables_duck.R\n") else {

## --- waves: NDBC standard meteorological files, one per year ---
wave <- do.call(rbind, lapply(years, function(y) {
  u <- sprintf("https://www.ndbc.noaa.gov/data/historical/stdmet/%sh%d.txt.gz", buoy, y)
  f <- tempfile(); ok <- tryCatch({ download.file(u, f, quiet = TRUE); TRUE }, error = function(e) FALSE)
  if (!ok) { cat("no wave file for", y, "\n"); return(NULL) }
  d <- read.table(gzfile(f), header = FALSE, skip = 2, fill = TRUE)
  names(d)[1:9] <- c("YY","MM","DD","hh","mm","WDIR","WSPD","GST","WVHT")
  d <- d[!is.na(d$WVHT) & d$WVHT < 90, ]                                    # 99 = missing
  data.frame(date = as.Date(sprintf("%d-%02d-%02d", d$YY, d$MM, d$DD)), wvht = d$WVHT) }))
wave <- aggregate(wvht ~ date, wave, max)

## --- surge: CO-OPS API, hourly observed and predicted water level, one year per request ---
coops <- function(product, y) {
  u <- sprintf(paste0("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=%s",
       "&application=star&begin_date=%d0101&end_date=%d1231&datum=MSL&station=%s",
       "&time_zone=GMT&units=metric&interval=h&format=csv"), product, y, y, gauge)
  d <- tryCatch(read.csv(u), error = function(e) NULL); if (is.null(d) || ncol(d) < 2) return(NULL)
  names(d)[1:2] <- c("t", "v"); data.frame(t = as.POSIXct(d$t, tz = "GMT"), v = as.numeric(d$v)) }
obs <- do.call(rbind, lapply(years, function(y) coops("hourly_height", y)))
prd <- do.call(rbind, lapply(years, function(y) coops("predictions", y)))
m <- merge(obs, prd, by = "t", suffixes = c(".obs", ".prd")); m$surge <- m$v.obs - m$v.prd
m$date <- as.Date(m$t); surge <- aggregate(surge ~ date, m[!is.na(m$surge), ], max)

## --- pair, standardize, fit ---
a <- merge(wave, surge, by = "date"); cat("Duck NC paired days:", nrow(a), "\n")
cat("lag-1 autocorrelation wave", round(acf(a$wvht, plot=FALSE)$acf[2],2), "surge", round(acf(a$surge, plot=FALSE)$acf[2],2), "\n")
st <- standardize(a$wvht, a$surge); E <- st$E
cat("GPD shapes: wave", round(st$gx$xi,2), " surge", round(st$gy$xi,2), "\n")
set.seed(2026)
fit <- star_mcmc(E, prior = star_prior("threecomp", p_min = 1),
                  n_iter = ITER, burnin = BURN, thin = 5, verbose = TRUE)
rep <- report(fit, E, "DUCK NC WAVE-SURGE", days = 365)
saveRDS(list(data = a, E = E, margins = st[c("gx","gy")], fit = fit, report = rep), "results/wavesurge.rds")
cat("\nsaved results/wavesurge.rds\n") }
