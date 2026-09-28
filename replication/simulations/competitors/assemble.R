# assemble.R -- rebuilds results/ht_all.rds and results/wc_all.rds from the
# per-chunk outputs of run_ht.R (htres/) and run_wc.R (wcres/).
# Run from the competitors folder:  Rscript assemble.R
# Each chunk file <dir>/<tag>_<chunk>.rds is a list of data frames, one per
# replicate; the assembled object is list(s1 =, s3 =, s4 =), each element the
# rbind of all chunks of that tag in chunk order.
# Optional argument: a suffix for the output names (e.g. "_check"), used to
# verify the rebuild without overwriting the saved files.

args <- commandArgs(trailingOnly = TRUE)
suffix <- if (length(args)) args[1] else ""
TAGS <- c("s1", "s3", "s4")

assemble <- function(dir) {
  out <- lapply(TAGS, function(tag) {
    files <- sort(list.files(dir, pattern = sprintf("^%s_[0-9]+\\.rds$", tag),
                             full.names = TRUE))
    if (!length(files)) stop("no chunk files for tag ", tag, " in ", dir)
    rows <- unlist(lapply(files, readRDS), recursive = FALSE)
    bad <- !vapply(rows, is.data.frame, logical(1))
    if (any(bad)) stop(sum(bad), " failed replicate(s) in ", dir, "/", tag, ": ",
                       paste(unique(unlist(rows[bad])), collapse = "; "))
    do.call(rbind, rows)
  })
  setNames(out, TAGS)
}

ht_all <- assemble("htres")
wc_all <- assemble("wcres")
saveRDS(ht_all, sprintf("results/ht_all%s.rds", suffix))
saveRDS(wc_all, sprintf("results/wc_all%s.rds", suffix))
cat(sprintf("ht_all: %s rows; wc_all: %s rows\n",
            paste(sapply(ht_all, nrow), collapse = "/"),
            paste(sapply(wc_all, nrow), collapse = "/")))
