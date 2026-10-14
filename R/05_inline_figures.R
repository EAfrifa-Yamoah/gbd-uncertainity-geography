#!/usr/bin/env Rscript
# Move the figure environments from the end of main.tex to the paragraph where
# each is first cited, and write main_inline.tex. Placement is [htbp] so LaTeX
# sets each figure on or after the page of first mention.
if (!dir.exists("manuscript")) stop("run from the repository root: Rscript R/05_inline_figures.R")
x <- readLines("manuscript/main.tex")
start <- grep("^\\\\section\\*\\{Figures\\}", x)
clear <- start - 1  # the \clearpage before the Figures section
end_doc <- grep("^\\\\end\\{document\\}", x)
fig_block <- x[(start + 1):(end_doc - 1)]
body <- x[1:(clear - 1)]
# split the figure block into environments
b <- grep("^\\\\begin\\{figure\\}", fig_block); e <- grep("^\\\\end\\{figure\\}", fig_block)
figs <- lapply(seq_along(b), function(i) sub("\\[p\\]", "[htbp]", fig_block[b[i]:e[i]]))
labels <- vapply(figs, function(f) sub(".*\\\\label\\{(fig:[a-z]+)\\}.*", "\\1", f[grep("\\\\label", f)]), "")
# for each figure, find the first body line citing it and insert after that paragraph
out <- body
# insert in reverse order of first citation, so that two figures cited in the
# same paragraph keep their citation order after insertion
first_cite <- vapply(labels, function(l) grep(paste0("\\\\ref\\{", l, "\\}"), body)[1], 1)
for (i in order(first_cite, seq_along(figs), decreasing = TRUE)) {
  pat <- paste0("\\\\ref\\{", labels[i], "\\}")
  hit <- grep(pat, out)[1]
  # end of the paragraph: next blank line at or after hit
  k <- hit; while (k <= length(out) && nzchar(out[k])) k <- k + 1
  out <- append(out, c("", figs[[i]]), after = k - 1)
}
out <- c(out, "", "\\end{document}")
writeLines(out, "manuscript/main_inline.tex")
cat("figures inlined:", paste(labels, collapse = ", "), "\n")
