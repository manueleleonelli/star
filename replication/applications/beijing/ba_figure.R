## ba_figure.R -- Bayesian Analysis layout of Figure 5: takes fig-beijing-new.tex (written by sec6_cut.R) and
## wraps each row of three panels in \resizebox{\textwidth}{!}{...}. Writes fig-beijing-ba.tex, which the
## paper inputs as fig-beijing.tex. Run from this folder after app1_showcase.R and sec6_cut.R.
fb <- readLines("fig-beijing-new.tex")
k1 <- grep("^\\\\begin\\{tikzpicture\\}", fb)[1]                  # first panel
kr <- which(fb == "\\end{tikzpicture}\\\\[4pt]%")                  # end of the first row
kl <- max(which(fb == "\\end{tikzpicture}"))                        # end of the last panel
stopifnot(length(kr) == 1, k1 < kr, kr < kl)
fb[1] <- "% BA version of fig-beijing.tex: each row scaled to the text width"
fb[kl] <- "\\end{tikzpicture}}"
fb <- c(fb[1:(k1 - 1)], "\\resizebox{\\textwidth}{!}{%", fb[k1:(kr - 1)], "\\end{tikzpicture}}\\\\[4pt]", "\\resizebox{\\textwidth}{!}{%", fb[(kr + 1):length(fb)])
## no newline after the last line, as in the file of the paper
con <- file("fig-beijing-ba.tex", "wb"); writeChar(paste(fb, collapse = "\n"), con, eos = NULL); close(con)
cat("wrote fig-beijing-ba.tex,", length(fb), "lines\n")
