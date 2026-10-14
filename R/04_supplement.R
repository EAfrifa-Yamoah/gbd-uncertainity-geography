#!/usr/bin/env Rscript
# Supplementary analyses and LaTeX tables for the manuscript. Every table in
# the main text and supplement is written here from the data; nothing is typed.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(mgcv)})
if (!dir.exists("data/raw")) stop("run from the repository root: Rscript R/04_supplement.R")
P <- readRDS("data/panel/panels.rds"); cc <- P$country_cause; cv <- P$covars
dir.create("manuscript/tables", showWarnings = FALSE, recursive = TRUE)
sink("output/stats_supplement.txt", split = TRUE)
sr_levels <- c("High-income", "Central Europe, Eastern Europe, and Central Asia", "Latin America and Caribbean",
               "North Africa and Middle East", "South Asia", "Southeast Asia, East Asia, and Oceania", "Sub-Saharan Africa")
sr_tex <- c("High income", "Central Europe, Eastern Europe, Central Asia", "Latin America and Caribbean",
            "North Africa and Middle East", "South Asia", "Southeast Asia, East Asia, Oceania", "Sub-Saharan Africa")
cc <- cc %>% mutate(sr = factor(super_region_name, sr_levels, sr_tex)) %>% left_join(cv, by = "location_id")
f3 <- function(x) formatC(x, format = "f", digits = 3)
f2 <- function(x) formatC(x, format = "f", digits = 2)
f1 <- function(x) formatC(x, format = "f", digits = 1)
tex_table <- function(df, file, align) {
  hdr <- paste(names(df), collapse = " & ")
  body <- apply(df, 1, function(r) paste(r, collapse = " & "))
  writeLines(c(sprintf("\\begin{tabular}{%s}", align), "\\toprule", paste0(hdr, " \\\\"), "\\midrule",
               paste0(body, " \\\\"), "\\bottomrule", "\\end{tabular}"), file.path("manuscript/tables", file))
}
ac <- cc %>% filter(cause_id == 294, year == 2023)

# ---- Table 1: median rw by super region and measure, 2023 -------------------
t1 <- ac %>% group_by(sr, measure) %>% summarise(m = median(rw), .groups = "drop") %>%
  pivot_wider(names_from = measure, values_from = m) %>%
  select(sr, deaths, ylls, ylds, dalys)
n_sr <- ac %>% filter(measure == "dalys") %>% count(sr)
t1 <- t1 %>% left_join(n_sr, by = "sr") %>%
  transmute(`Super region` = as.character(sr), `$n$` = n, Deaths = f3(deaths), YLLs = f3(ylls), YLDs = f3(ylds), DALYs = f3(dalys))
glob <- ac %>% group_by(measure) %>% summarise(m = median(rw)) %>% pivot_wider(names_from = measure, values_from = m)
t1 <- bind_rows(t1, tibble(`Super region` = "All countries", `$n$` = 204L, Deaths = f3(glob$deaths),
                           YLLs = f3(glob$ylls), YLDs = f3(glob$ylds), DALYs = f3(glob$dalys)))
tex_table(t1, "tab1_region_measure.tex", "lrrrrr")

# ---- Table 2: cause ranking ---------------------------------------------------
cl <- cc %>% filter(cause_id != 294, year == 2023, measure == "dalys", lower > 0, point > 0)
cs <- cl %>% group_by(cause_name, cause_type) %>%
  summarise(n = n(), med = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75),
            rho = cor(rw, sdi, method = "spearman"), .groups = "drop") %>% arrange(desc(med))
dd <- cc %>% filter(cause_id != 294, year == 2023, measure %in% c("deaths", "dalys"), lower > 0) %>%
  select(location_name, cause_name, measure, rw) %>% pivot_wider(names_from = measure, values_from = rw) %>%
  filter(!is.na(deaths)) %>% group_by(cause_name) %>% summarise(ratio = median(deaths / dalys))
t2 <- cs %>% left_join(dd, by = "cause_name") %>%
  transmute(Cause = cause_name, Type = cause_type, `$n$` = n,
            `Median (IQR)` = sprintf("%s (%s to %s)", f2(med), f2(q25), f2(q75)),
            `$\\rho$ with SDI` = f2(rho), `Deaths/DALYs` = ifelse(is.na(ratio), "--", f2(ratio)))
tex_table(t2, "tab2_causes.tex", "llrlrr")

# ---- Table 3: composition test ----------------------------------------------
comp <- read_csv("output/composition_allcause_2023.csv", show_col_types = FALSE) %>%
  mutate(sr = factor(super_region, sr_tex))
t3 <- comp %>% group_by(sr) %>%
  summarise(n = n(), share = median(yld_share), rwd = median(rw_dalys), contrib = median(yld_contrib),
            ratio = median(rw_dalys / rw_pred_indep), rho = median(rho_implied), .groups = "drop")
t3 <- bind_rows(t3, comp %>% summarise(sr = "All countries", n = n(), share = median(yld_share), rwd = median(rw_dalys),
                                     contrib = median(yld_contrib), ratio = median(rw_dalys / rw_pred_indep), rho = median(rho_implied)))
t3 <- t3 %>% transmute(`Super region` = as.character(sr), `$n$` = n, `YLD share` = f3(share), `rw(DALYs)` = f3(rwd),
                       `YLD variance share` = f3(contrib), `Observed/independent` = f3(ratio), `Implied $\\rho$` = f3(rho))
tex_table(t3, "tab3_composition.tex", "lrrrrrr")

# ---- S1: robustness of the ranking to the width metric ---------------------
cat("== S1. Metric robustness ==\n")
for (m in c("deaths", "ylds", "dalys")) {
  d <- ac %>% filter(measure == m) %>% mutate(abw = upper - lower)
  cat(sprintf("%s: Spearman rank(rw) ~ rank(lrw) = %.3f; rw ~ absolute width = %.3f; rw ~ SDI = %.3f; lrw ~ SDI = %.3f; abs ~ SDI = %.3f\n",
              m, cor(d$rw, d$lrw, method = "spearman"), cor(d$rw, d$abw, method = "spearman"),
              cor(d$rw, d$sdi, method = "spearman"), cor(d$lrw, d$sdi, method = "spearman"), cor(d$abw, d$sdi, method = "spearman")))
}
s1 <- ac %>% mutate(abw = upper - lower) %>% group_by(sr, measure) %>%
  summarise(rw = median(rw), lrw = median(lrw), abw = median(abw), .groups = "drop") %>%
  filter(measure %in% c("deaths", "dalys")) %>%
  mutate(measure = toupper(measure)) %>%
  pivot_wider(names_from = measure, values_from = c(rw, lrw, abw)) %>%
  transmute(`Super region` = as.character(sr), `rw deaths` = f3(rw_DEATHS), `log ratio deaths` = f3(lrw_DEATHS),
            `absolute deaths` = f1(abw_DEATHS), `rw DALYs` = f3(rw_DALYS), `log ratio DALYs` = f3(lrw_DALYS), `absolute DALYs` = f0 <- formatC(abw_DALYS, format = "f", digits = 0))
tex_table(s1, "tabS1_metric_robustness.tex", "lrrrrrr")

# ---- S2: independence test with asymmetric arms -----------------------------
cat("\n== S2. Independence test with separate arms ==\n")
arms <- ac %>% select(location_name, sr, measure, point, lower, upper) %>%
  pivot_wider(names_from = measure, values_from = c(point, lower, upper)) %>%
  mutate(up_yll = upper_ylls - point_ylls, up_yld = upper_ylds - point_ylds, up_daly = upper_dalys - point_dalys,
         lo_yll = point_ylls - lower_ylls, lo_yld = point_ylds - lower_ylds, lo_daly = point_dalys - lower_dalys,
         r_up = up_daly / sqrt(up_yll^2 + up_yld^2), r_lo = lo_daly / sqrt(lo_yll^2 + lo_yld^2),
         rho_up = (up_daly^2 - up_yll^2 - up_yld^2) / (2 * up_yll * up_yld),
         rho_lo = (lo_daly^2 - lo_yll^2 - lo_yld^2) / (2 * lo_yll * lo_yld))
cat(sprintf("upper arm: observed/independent median %.3f (IQR %.3f to %.3f), implied rho %.3f\n",
            median(arms$r_up), quantile(arms$r_up, .25), quantile(arms$r_up, .75), median(arms$rho_up)))
cat(sprintf("lower arm: observed/independent median %.3f (IQR %.3f to %.3f), implied rho %.3f\n",
            median(arms$r_lo), quantile(arms$r_lo, .25), quantile(arms$r_lo, .75), median(arms$rho_lo)))
s2 <- arms %>% group_by(sr) %>% summarise(n = n(), r_up = median(r_up), r_lo = median(r_lo), rho_up = median(rho_up), rho_lo = median(rho_lo)) %>%
  transmute(`Super region` = as.character(sr), `$n$` = n, `Upper arm ratio` = f3(r_up), `Lower arm ratio` = f3(r_lo),
            `Implied $\\rho$ upper` = f3(rho_up), `Implied $\\rho$ lower` = f3(rho_lo))
tex_table(s2, "tabS2_arms.tex", "lrrrrr")

# ---- S3: GAM deviance table -------------------------------------------------
cat("\n== S3. GAM deviance explained ==\n")
rows <- list()
for (m in c("deaths", "ylls", "ylds", "dalys")) {
  d <- ac %>% filter(measure == m) %>% mutate(y = log(rw), lpop = log10(population))
  dev <- function(f) 100 * summary(gam(f, data = d, method = "REML"))$dev.expl
  rows[[m]] <- tibble(Measure = toupper(m),
    `SDI` = f1(dev(y ~ s(sdi, k = 6))), `Population` = f1(dev(y ~ s(lpop, k = 6))),
    `Super region` = f1(dev(y ~ sr)), `SDI + population` = f1(dev(y ~ s(sdi, k = 6) + s(lpop, k = 6))),
    `SDI + population + region` = f1(dev(y ~ s(sdi, k = 6) + s(lpop, k = 6) + sr)),
    `HAQ + population + region` = f1(dev(y ~ s(haq_value, k = 6) + s(lpop, k = 6) + sr)))
}
s3 <- bind_rows(rows); print(as.data.frame(s3))
tex_table(s3, "tabS3_gam.tex", "lrrrrrr")

# ---- S4: time table (median rw by measure and year) --------------------------
tm <- cc %>% filter(cause_id == 294) %>% group_by(measure, year) %>% summarise(m = median(rw), .groups = "drop") %>%
  pivot_wider(names_from = measure, values_from = m) %>% arrange(year) %>%
  transmute(Year = year, Deaths = ifelse(is.na(deaths), "--", f3(deaths)), YLLs = ifelse(is.na(ylls), "--", f3(ylls)),
            YLDs = ifelse(is.na(ylds), "--", f3(ylds)), DALYs = f3(dalys))
tex_table(tm, "tabS4_time.tex", "lrrrr")

# ---- S5: admin1 summary -------------------------------------------------------
ads <- read_csv("output/admin1_summary.csv", show_col_types = FALSE)
nat_key <- c(brazil = "Brazil", ethiopia = "Ethiopia", india = "India", indonesia = "Indonesia", iran = "Iran", italy = "Italy",
             japan = "Japan", kenya = "Kenya", mexico = "Mexico", nigeria = "Nigeria", norway = "Norway", pakistan = "Pakistan",
             philippines = "Philippines", poland = "Poland", south_africa = "South Africa", uk = "United Kingdom", usa = "United States")
s5 <- ads %>% arrange(desc(med_over_nat)) %>%
  transmute(Country = nat_key[country], Units = n_units, `National rw` = f3(nat_rw), `Admin1 median` = f3(med_rw),
            `Admin1 range` = sprintf("%s to %s", f3(min_rw), f3(max_rw)), `Median/national` = f2(med_over_nat),
            `$\\rho$ with size` = f2(rho_rw_count))
tex_table(s5, "tabS5_admin1.tex", "lrrrlrr")

# ---- S6: risk tables ---------------------------------------------------------
sev <- P$country_sev %>% left_join(cv, by = "location_id")
rd  <- P$country_riskdaly %>% left_join(cv, by = "location_id") %>% filter(lower > 0)
s6a <- sev %>% group_by(risk_name) %>% summarise(med = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75),
                                                 rho_sdi = cor(rw, sdi, method = "spearman"), rho_lev = cor(rw, point, method = "spearman")) %>%
  arrange(desc(med)) %>% transmute(Risk = risk_name, Quantity = "Exposure (SEV)", `Median (IQR)` = sprintf("%s (%s to %s)", f2(med), f2(q25), f2(q75)),
                                   `$\\rho$ with SDI` = f2(rho_sdi), `$\\rho$ with level` = f2(rho_lev))
s6b <- rd %>% group_by(risk_name) %>% summarise(med = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75),
                                                rho_sdi = cor(rw, sdi, method = "spearman"), rho_lev = cor(rw, point, method = "spearman")) %>%
  arrange(desc(med)) %>% transmute(Risk = risk_name, Quantity = "Attributable DALYs", `Median (IQR)` = sprintf("%s (%s to %s)", f2(med), f2(q25), f2(q75)),
                                   `$\\rho$ with SDI` = f2(rho_sdi), `$\\rho$ with level` = f2(rho_lev))
tex_table(bind_rows(s6a, s6b), "tabS6_risks.tex", "llllr")

# ---- Practical helper: decompose a published DALY interval -------------------
cat("\n== Helper function test ==\n")
source("R/decompose_ui.R")
ex <- ac %>% filter(location_name == "Ghana") %>% select(measure, point, lower, upper)
g <- function(m, col) ex[[col]][ex$measure == m]
res <- decompose_daly_ui(daly = c(g("dalys","point"), g("dalys","lower"), g("dalys","upper")),
                         yll  = c(g("ylls","point"),  g("ylls","lower"),  g("ylls","upper")),
                         yld  = c(g("ylds","point"),  g("ylds","lower"),  g("ylds","upper")))
print(res, digits = 3)
sink()
cat("tables written to manuscript/tables\n")
