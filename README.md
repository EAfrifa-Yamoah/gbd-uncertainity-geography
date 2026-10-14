# gbd-uncertainty-geography

Code, data extracts and manuscript for

> Afrifa-Yamoah E. The geography of uncertainty in the Global Burden of
> Disease Study 2023: two geographies in one interval. Manuscript, 2026.

The study analyses the relative width of the 95% uncertainty intervals
published with GBD 2023 estimates across 204 countries, four measures, the 30
leading causes, 1990 to 2023, 508 subnational units and 16 risk quantities, and
shows that a DALY interval superimposes a data driven mortality component on a
globally shared disability weight component.

## The practical tool

`R/decompose_ui.R` provides `decompose_daly_ui()`, which reads any published
DALY uncertainty interval in terms of its YLL and YLD components from the three
published intervals alone:

```r
source("R/decompose_ui.R")
# Ghana 2023, all-cause, age standardised rates per 100,000: c(point, lower, upper)
decompose_daly_ui(daly = c(49775, 46069, 54329),
                  yll  = c(37964, 36285, 39789),
                  yld  = c(11811,  8991, 15592))
```

It returns the three relative widths, the YLD share of DALYs, the share of the
DALY interval variance that comes from YLDs under independent components, the
ratio of the observed DALY width to the width expected under independence, the
implied YLL to YLD draw correlation, and a one line reading rule. Section 4.2
of the manuscript (Box 1) describes the procedure.

## Layout

```
R/00_harvest.R         names and copies connector extracts into data/raw (only needed to add extracts)
R/01_build_panel.R     builds the five analysis panels in data/panel
R/02_analysis.R        every statistic in the paper (output/stats.txt) and the figures
R/04_supplement.R      every table in the main text and supplement, written to manuscript/tables
R/05_inline_figures.R  derives manuscript/main_inline.tex from main.tex
R/decompose_ui.R       the decomposition function
data/raw/              134 verbatim CSV extracts from the IHME GBD 2023 Collab connector, 14 October 2026
data/panel/            derived panels (CSV)
figures/               eight figures, 300 dpi, colour blind safe palettes
output/                stats.txt, stats_supplement.txt, summary CSVs
manuscript/            main.tex (figures at end), main_inline.tex (figures in text),
                       supplement.tex, refs.bib, tables/, compiled PDFs
```

## Reproduce

R >= 4.3 with jsonlite, dplyr, tidyr, readr, purrr, stringr, ggplot2, scales,
mgcv, patchwork, sf, quantreg, broom, rnaturalearthdata.

```
Rscript R/01_build_panel.R
LANG=C.UTF-8 Rscript R/02_analysis.R
Rscript R/04_supplement.R
cd manuscript && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Data source and terms

All estimates are from the Institute for Health Metrics and Evaluation (IHME),
Global Burden of Disease Study 2023, obtained through the GBD 2023 Collab
connector. Use of GBD results is subject to the IHME free of charge
non commercial user agreement; check the current terms before redistributing
the extracts in `data/` and remove them from the public repository if the
terms require it (the scripts regenerate them from the connector).

## Licence

Code: MIT (see LICENSE). Manuscript text and figures: all rights reserved
pending publication.
