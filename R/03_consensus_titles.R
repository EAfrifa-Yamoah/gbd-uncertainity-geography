# Print titles and first sentence of abstracts from a saved Consensus search result
suppressPackageStartupMessages(library(jsonlite))
f <- commandArgs(trailingOnly = TRUE)[1]
t <- fromJSON(f)$text[1]
m <- regmatches(t, gregexpr("\\[\\d+\\] \\[[^\\]]+\\]\\([^)]+\\) \\([^)]+\\)", t))[[1]]
cat(m, sep = "\n")
