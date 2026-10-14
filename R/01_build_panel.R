#!/usr/bin/env Rscript
# Build tidy analysis panels from the raw GBD 2023 connector extracts.
# Every row keeps the connector's own point estimate and 95% uncertainty
# interval (UI). Derived quantities:
#   rw    = (upper - lower) / point          relative UI width (dimensionless)
#   lrw   = log(upper / lower)               log interval ratio (scale free)
#   asym  = (upper - point) / (point - lower) asymmetry (> 1: right skewed)
# Only the combined sex, age standardised series is used for the main panels.
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(purrr); library(stringr)})
setwd("/home/claude/gbd_uncertainty")
raw <- "data/raw"
files <- list.files(raw, full.names = TRUE)

derive <- function(d, point, lo, hi) {
  d$point <- d[[point]]; d$lower <- d[[lo]]; d$upper <- d[[hi]]
  d %>% mutate(rw   = (upper - lower) / point,
               lrw  = log(upper / lower),
               asym = (upper - point) / (point - lower))
}

# ---- 1. Country level cause panel (get_top_locations, cause mode) ----------
cause_files <- grep("get_top_locations__.*__(dalys|deaths|ylls|ylds|prevalence)__", files, value = TRUE)
country_cause <- map_dfr(cause_files, function(f) {
  d <- read_csv(f, show_col_types = FALSE, progress = FALSE)
  d %>% select(location_id, location_name, super_region_name, cause_id, cause_name,
               cause_type_id, year, measure, age_group,
               rate_combined, rate_lower_combined, rate_upper_combined,
               count_combined, rate_male, rate_lower_male, rate_upper_male,
               rate_female, rate_lower_female, rate_upper_female) %>%
    derive("rate_combined", "rate_lower_combined", "rate_upper_combined")
})
country_cause <- country_cause %>%
  mutate(cause_type = recode(as.character(cause_type_id), "294" = "All causes",
                             "295" = "CMNN", "409" = "NCD", "687" = "Injuries"))
cat("country_cause rows:", nrow(country_cause), " causes:", n_distinct(country_cause$cause_id),
    " measures:", paste(unique(country_cause$measure), collapse = ","), "\n")

# ---- 2. Country level SEV panel --------------------------------------------
sev_files <- grep("get_top_locations__.*__sev__", files, value = TRUE)
country_sev <- map_dfr(sev_files, function(f) read_csv(f, show_col_types = FALSE, progress = FALSE)) %>%
  select(location_id, location_name, super_region_name, rei_id, risk_name, year,
         sev_combined, sev_lower_combined, sev_upper_combined) %>%
  derive("sev_combined", "sev_lower_combined", "sev_upper_combined")
cat("country_sev rows:", nrow(country_sev), " risks:", n_distinct(country_sev$rei_id), "\n")

# ---- 3. Country level risk attributable DALYs -------------------------------
risk_files <- grep("riskdalys__", files, value = TRUE)
country_riskdaly <- map_dfr(risk_files, function(f) read_csv(f, show_col_types = FALSE, progress = FALSE)) %>%
  select(location_id, location_name, super_region_name, rei_id, risk_name, year, measure,
         rate_combined, rate_lower_combined, rate_upper_combined, percent_combined,
         percent_lower_combined, percent_upper_combined) %>%
  derive("rate_combined", "rate_lower_combined", "rate_upper_combined")
cat("country_riskdaly rows:", nrow(country_riskdaly), "\n")

# ---- 4. Covariates ----------------------------------------------------------
sdi <- read_csv(file.path(raw, "covariate__sdi_2023.csv"), show_col_types = FALSE) %>%
  select(location_id, sdi, sdi_quintile)
cov_read <- function(f, nm) read_csv(file.path(raw, f), show_col_types = FALSE) %>%
  select(location_id, value_combined, lower_combined, upper_combined) %>%
  rename_with(~ paste0(nm, "_", sub("_combined", "", .x)), -location_id)
haq <- cov_read("covariate__haqi_2023.csv", "haq")
uhc <- cov_read("covariate__uhc_2023.csv", "uhc")
ldi <- cov_read("covariate__ldi_2023.csv", "ldi")
pop <- read_csv(file.path(raw, "covariate__population_2023.csv"), show_col_types = FALSE)
pop <- pop %>% select(location_id, population = value_combined)
covars <- sdi %>% left_join(haq, by = "location_id") %>% left_join(uhc, by = "location_id") %>%
  left_join(ldi, by = "location_id") %>% left_join(pop, by = "location_id") %>%
  mutate(sdi_quintile = factor(sdi_quintile,
         levels = c("Low SDI", "Low-middle SDI", "Middle SDI", "High-middle SDI", "High SDI")))
cat("covars rows:", nrow(covars), " missing HAQ:", sum(is.na(covars$haq_value)),
    " missing pop:", sum(is.na(covars$population)), "\n")
print(table(covars$sdi_quintile, useNA = "ifany"))

# ---- 5. Admin1 panels -------------------------------------------------------
adm_daly_files <- grep("admin1__.*__dalys\\.csv", files, value = TRUE)
admin1_daly <- map_dfr(adm_daly_files, function(f) {
  read_csv(f, show_col_types = FALSE, progress = FALSE) %>%
    mutate(country = str_match(basename(f), "admin1__(.*)__dalys")[, 2]) %>%
    select(country, location_id, location_name, year, measure,
           rate_combined, rate_lower_combined, rate_upper_combined, count_combined) %>%
    derive("rate_combined", "rate_lower_combined", "rate_upper_combined")
})
adm_sev_files <- grep("admin1__.*__sbp_sev\\.csv", files, value = TRUE)
admin1_sev <- map_dfr(adm_sev_files, function(f) {
  read_csv(f, show_col_types = FALSE, progress = FALSE) %>%
    mutate(country = str_match(basename(f), "admin1__(.*)__sbp_sev")[, 2]) %>%
    select(country, location_id, location_name, year,
           sev_combined, sev_lower_combined, sev_upper_combined) %>%
    derive("sev_combined", "sev_lower_combined", "sev_upper_combined")
})
cat("admin1_daly rows:", nrow(admin1_daly), " countries:", n_distinct(admin1_daly$country), "\n")
cat("admin1_sev rows:", nrow(admin1_sev), "\n")

dir.create("data/panel", showWarnings = FALSE)
write_csv(country_cause,    "data/panel/country_cause.csv")
write_csv(country_sev,      "data/panel/country_sev.csv")
write_csv(country_riskdaly, "data/panel/country_riskdaly.csv")
write_csv(covars,           "data/panel/covariates.csv")
write_csv(admin1_daly,      "data/panel/admin1_daly.csv")
write_csv(admin1_sev,       "data/panel/admin1_sev.csv")
saveRDS(list(country_cause = country_cause, country_sev = country_sev,
             country_riskdaly = country_riskdaly, covars = covars,
             admin1_daly = admin1_daly, admin1_sev = admin1_sev), "data/panel/panels.rds")
cat("done\n")
