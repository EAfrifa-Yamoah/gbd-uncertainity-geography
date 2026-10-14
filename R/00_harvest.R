#!/usr/bin/env Rscript
# Harvest saved GBD connector tool results (JSON {result: "<csv>"}) into
# data/raw/, naming each file from its own content so that parallel tool
# calls cannot be mislabelled. Files already harvested are renamed *.harvested.
suppressPackageStartupMessages({library(jsonlite)})
dir_in  <- "/root/.claude/projects/-home-claude/4a32cf05-c4a9-5296-84b9-e2caa28798ef/tool-results"
dir_out <- "/home/claude/gbd_uncertainty/data/raw"
files <- list.files(dir_in, pattern = "^mcp-GBD_2023_Collab-.*\\.txt$", full.names = TRUE)
slug <- function(x) gsub("[^A-Za-z0-9]+", "_", tolower(x))
for (f in files) {
  txt <- tryCatch(fromJSON(f)$result, error = function(e) NULL)
  if (is.null(txt) || !nzchar(txt)) { cat("skip (no result):", basename(f), "\n"); next }
  d <- tryCatch(read.csv(text = txt, check.names = FALSE, stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) { cat("skip (not csv):", basename(f), "\n"); next }
  tool <- sub("^mcp-GBD_2023_Collab-(.*)-[0-9]+\\.txt$", "\\1", basename(f))
  parts <- c(tool)
  if ("risk_name" %in% names(d))       parts <- c(parts, slug(d$risk_name[1]))
  if ("rei_name"  %in% names(d))       parts <- c(parts, slug(d$rei_name[1]))
  if ("cause_name" %in% names(d))      parts <- c(parts, slug(d$cause_name[1]))
  if ("impairment_name" %in% names(d)) parts <- c(parts, slug(d$impairment_name[1]))
  if ("measure" %in% names(d))         parts <- c(parts, d$measure[1])
  if ("age_group" %in% names(d))       parts <- c(parts, d$age_group[1])
  if ("year" %in% names(d))            parts <- c(parts, paste0("y", d$year[1]))
  if ("covariate_name" %in% names(d))  parts <- c(parts, slug(d$covariate_name[1]))
  # location scope: admin1 harvests carry a parent column, country harvests do not
  if ("parent_name" %in% names(d))     parts <- c(parts, "within", slug(d$parent_name[1]))
  if ("super_region_name" %in% names(d) && length(unique(d$super_region_name)) == 1 &&
      tool == "get_top_locations")     parts <- c(parts, "within", slug(d$super_region_name[1]))
  name <- paste(parts, collapse = "__")
  out <- file.path(dir_out, paste0(name, ".csv"))
  k <- 1; base <- out
  while (file.exists(out)) { k <- k + 1; out <- sub("\\.csv$", paste0("_", k, ".csv"), base) }
  writeLines(txt, out)
  file.rename(f, paste0(f, ".harvested"))
  cat(sprintf("%-90s %5d rows\n", basename(out), nrow(d)))
}
