#!/usr/bin/env Rscript
# The geography of epistemic uncertainty in GBD 2023: analysis script.
# Outcome throughout: relative 95% UI width rw = (upper - lower)/point of the
# combined sex, age standardised estimate. Secondary: lrw = log(upper/lower).
# All numbers quoted in the manuscript are written to output/stats.txt.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(scales)
  library(mgcv); library(patchwork); library(sf); library(quantreg); library(broom)
})
if (!dir.exists("data/raw")) stop("run from the repository root: Rscript R/02_analysis.R")
P <- readRDS("data/panel/panels.rds")
cc <- P$country_cause; cv <- P$covars; sev <- P$country_sev; rd <- P$country_riskdaly
ad <- P$admin1_daly; as <- P$admin1_sev
sink("output/stats.txt", split = TRUE)
cat("GBD 2023 uncertainty width analysis, run", format(Sys.time()), "\n\n")

# ---- palettes (fixed meanings across all figures) ---------------------------
# super regions: Okabe-Ito (colour blind safe); never reused for anything else
sr_levels <- c("High-income", "Central Europe, Eastern Europe, and Central Asia",
               "Latin America and Caribbean", "North Africa and Middle East",
               "South Asia", "Southeast Asia, East Asia, and Oceania", "Sub-Saharan Africa")
sr_short  <- c("High income", "Central Europe, Eastern Europe, Central Asia",
               "Latin America and Caribbean", "North Africa and Middle East",
               "South Asia", "Southeast Asia, East Asia, Oceania", "Sub-Saharan Africa")
sr_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
names(sr_pal) <- sr_short
# measures: greys/shapes only, so measure never competes with super region colour
meas_levels <- c("deaths", "ylls", "ylds", "dalys", "prevalence")
meas_labels <- c("Deaths", "YLLs", "YLDs", "DALYs", "Prevalence")
cc <- cc %>% mutate(super_region = factor(super_region_name, sr_levels, sr_short),
                    measure_f = factor(measure, meas_levels, meas_labels)) %>%
  left_join(cv, by = "location_id")
sev <- sev %>% mutate(super_region = factor(super_region_name, sr_levels, sr_short)) %>% left_join(cv, by = "location_id")
rd  <- rd  %>% mutate(super_region = factor(super_region_name, sr_levels, sr_short)) %>% left_join(cv, by = "location_id")
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
          legend.position = "bottom", strip.text = element_text(face = "bold")))
fig <- function(p, name, w = 8, h = 5) ggsave(file.path("figures", name), p, width = w, height = h, dpi = 300, bg = "white")

# =============================================================================
# 1. All cause burden, 2023: how wide are the intervals, and where?
# =============================================================================
ac <- cc %>% filter(cause_id == 294, year == 2023)
cat("== 1. All-cause 2023, relative UI width by measure (204 countries) ==\n")
s1 <- ac %>% group_by(measure_f) %>%
  summarise(n = n(), median = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75),
            min = min(rw), max = max(rw), p_gt_0.25 = mean(rw > 0.25), p_gt_0.5 = mean(rw > 0.5),
            med_asym = median(asym), .groups = "drop")
print(as.data.frame(s1), digits = 3)
cat("\nSpearman correlation of rw across measures within country:\n")
w1 <- ac %>% select(location_name, measure, rw) %>% pivot_wider(names_from = measure, values_from = rw)
print(round(cor(w1[, meas_levels[1:4]], method = "spearman"), 2))

cat("\nWidest and narrowest all-cause DALY intervals, 2023:\n")
d23 <- ac %>% filter(measure == "dalys") %>% arrange(desc(rw))
print(as.data.frame(d23 %>% select(location_name, super_region, sdi, point, lower, upper, rw) %>% head(12)), digits = 4)
print(as.data.frame(d23 %>% select(location_name, super_region, sdi, point, lower, upper, rw) %>% tail(8)), digits = 4)
cat("\nBy super region, all-cause DALYs 2023:\n")
print(as.data.frame(d23 %>% group_by(super_region) %>%
  summarise(n = n(), median_rw = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75), .groups = "drop")), digits = 3)
cat("\nBy SDI quintile (GBD assignment), all-cause DALYs 2023:\n")
print(as.data.frame(d23 %>% group_by(sdi_quintile) %>%
  summarise(n = n(), median_rw = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75), .groups = "drop")), digits = 3)
cat("\nBy super region, all-cause deaths 2023:\n")
print(as.data.frame(ac %>% filter(measure == "deaths") %>% group_by(super_region) %>%
  summarise(n = n(), median_rw = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75), .groups = "drop")), digits = 3)

# ---- 1b. Composition: is DALY uncertainty the YLD share in disguise? --------
cat("\n== 1b. Composition test: all-cause DALYs 2023 ==\n")
comp <- ac %>% select(location_name, super_region, sdi, measure, point, lower, upper, rw) %>%
  pivot_wider(names_from = measure, values_from = c(point, lower, upper, rw)) %>%
  mutate(yld_share = point_ylds / point_dalys,
         hw_yll = (upper_ylls - lower_ylls) / 2, hw_yld = (upper_ylds - lower_ylds) / 2,
         hw_daly = (upper_dalys - lower_dalys) / 2,
         # predicted DALY half width if YLL and YLD draws were independent, or perfectly correlated
         hw_indep = sqrt(hw_yll^2 + hw_yld^2), hw_perf = hw_yll + hw_yld,
         rw_pred_indep = 2 * hw_indep / point_dalys, rw_pred_perf = 2 * hw_perf / point_dalys,
         # implied correlation between YLL and YLD draws from the observed DALY width
         rho_implied = (hw_daly^2 - hw_yll^2 - hw_yld^2) / (2 * hw_yll * hw_yld),
         yld_contrib = hw_yld^2 / (hw_yll^2 + hw_yld^2))
cat("YLD share of DALYs: median", round(median(comp$yld_share), 3), " range", round(range(comp$yld_share), 3), "\n")
cat("Spearman rw(DALYs) ~ YLD share:", round(cor(comp$rw_dalys, comp$yld_share, method = "spearman"), 3), "\n")
cat("Spearman YLD share ~ SDI:", round(cor(comp$yld_share, comp$sdi, method = "spearman"), 3), "\n")
cat("Share of DALY interval variance (half width squared) contributed by YLDs under independence: median",
    round(median(comp$yld_contrib), 3), " IQR", round(quantile(comp$yld_contrib, c(.25, .75)), 3), "\n")
cat("Observed DALY rw / predicted under independence: median", round(median(comp$rw_dalys / comp$rw_pred_indep), 3),
    " IQR", round(quantile(comp$rw_dalys / comp$rw_pred_indep, c(.25, .75)), 3), "\n")
cat("Observed DALY rw / predicted under perfect correlation: median", round(median(comp$rw_dalys / comp$rw_pred_perf), 3), "\n")
cat("Implied YLL-YLD draw correlation: median", round(median(comp$rho_implied), 3),
    " IQR", round(quantile(comp$rho_implied, c(.25, .75)), 3), "\n")
lm_c <- lm(log(rw_dalys) ~ log(yld_share), data = comp); lm_s <- lm(log(rw_dalys) ~ sdi, data = comp)
lm_cs <- lm(log(rw_dalys) ~ log(yld_share) + sdi, data = comp)
cat(sprintf("R2 log rw(DALYs): YLD share alone %.3f | SDI alone %.3f | both %.3f (SDI coefficient given share: %.3f, p = %.3g)\n",
            summary(lm_c)$r.squared, summary(lm_s)$r.squared, summary(lm_cs)$r.squared,
            coef(lm_cs)["sdi"], summary(lm_cs)$coefficients["sdi", 4]))
# does the YLD rw itself vary with SDI once share is accounted for?
cat("Spearman rw(YLDs) ~ SDI:", round(cor(comp$rw_ylds, comp$sdi, method = "spearman"), 3),
    "; rw(YLLs) ~ SDI:", round(cor(comp$rw_ylls, comp$sdi, method = "spearman"), 3), "\n")
pC <- ggplot(comp, aes(yld_share, rw_dalys)) +
  geom_point(aes(colour = super_region), size = 1.4, alpha = 0.85) +
  scale_colour_manual(values = sr_pal, name = NULL) +
  labs(x = "YLD share of all-cause DALYs, 2023", y = "Relative 95% UI width, all-cause DALYs") +
  guides(colour = guide_legend(nrow = 4))
fig(pC, "fig1b_daly_rw_vs_yld_share.png", 8, 6)
write_csv(comp, "output/composition_allcause_2023.csv")

# Figure 1: distribution by super region and measure
p1 <- ggplot(ac, aes(x = super_region, y = rw, fill = super_region)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
  facet_wrap(~ measure_f, nrow = 1) +
  scale_fill_manual(values = sr_pal, name = NULL) +
  scale_y_continuous(trans = "log10", breaks = c(0.003, 0.01, 0.03, 0.1, 0.3), labels = c("0.003", "0.01", "0.03", "0.1", "0.3")) +
  labs(x = NULL, y = "Relative 95% UI width, (upper - lower) / estimate") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 3))
fig(p1, "fig1_allcause_rw_by_region_measure.png", 10, 5)

# Figure 2: world maps of rw, DALYs and deaths
data(countries50, package = "rnaturalearthdata")
world <- st_as_sf(countries50) %>% filter(admin != "Antarctica") %>%
  select(name_long, admin, iso_a3, geometry)
xw <- c("Lao People's Democratic Republic" = "Lao PDR", "Viet Nam" = "Vietnam",
        "Sao Tome and Principe" = "São Tomé and Principe",
        "Micronesia (Federated States of)" = "Federated States of Micronesia",
        "Democratic People's Republic of Korea" = "Dem. Rep. Korea",
        "Republic of Moldova" = "Moldova", "Iran (Islamic Republic of)" = "Iran",
        "Bolivia (Plurinational State of)" = "Bolivia",
        "Venezuela (Bolivarian Republic of)" = "Venezuela", "Türkiye" = "Turkey",
        "Syrian Arab Republic" = "Syria", "Côte d'Ivoire" = "Côte d'Ivoire",
        "Russian Federation" = "Russian Federation", "United Republic of Tanzania" = "Tanzania",
        "Brunei Darussalam" = "Brunei Darussalam", "Czechia" = "Czech Republic",
        "Eswatini" = "eSwatini",
        "Congo" = "Republic of the Congo", "Timor-Leste" = "Timor-Leste",
        "Gambia" = "The Gambia", "Bahamas" = "Bahamas")
mapd <- ac %>% filter(measure %in% c("dalys", "deaths")) %>%
  mutate(ne_name = ifelse(location_name %in% names(xw), xw[location_name], location_name))
ascii <- function(x) { x <- enc2utf8(x); x[grepl("Ivoire", x)] <- "Cote d'Ivoire"; x[grepl("rkiye", x)] <- "Turkey"
  iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "") }
mapd$ne_name <- ascii(mapd$ne_name); world$name_long <- ascii(world$name_long); world$admin <- ascii(world$admin)
unmatched <- setdiff(unique(mapd$ne_name), c(world$name_long, world$admin))
cat("\nCountries not matched to map polygons (small islands expected):", length(unmatched), "\n")
cat(paste(unmatched, collapse = "; "), "\n")
wm <- world %>% mutate(key = ifelse(name_long %in% mapd$ne_name, name_long, admin))
mapj <- mapd %>% select(ne_name, measure_f, rw) %>%
  inner_join(wm, by = c("ne_name" = "key")) %>% st_as_sf()
p2 <- ggplot() +
  geom_sf(data = world, fill = "grey92", colour = "white", linewidth = 0.1) +
  geom_sf(data = mapj, aes(fill = rw), colour = "white", linewidth = 0.1) +
  facet_wrap(~ measure_f, ncol = 1) +
  scale_fill_viridis_c(option = "viridis", trans = "log10", name = "Relative 95% UI width",
                       breaks = c(0.05, 0.1, 0.2, 0.4), labels = c("0.05", "0.1", "0.2", "0.4")) +
  coord_sf(crs = "+proj=robin") + theme_void(base_size = 11) +
  theme(legend.position = "bottom", legend.key.width = unit(1.6, "cm"), strip.text = element_text(face = "bold"))
fig(p2, "fig2_map_allcause_rw.png", 8, 8)

# =============================================================================
# 2. What predicts interval width? SDI, HAQ, population size
# =============================================================================
cat("\n== 2. Predictors of log relative width, all-cause 2023 ==\n")
for (m in c("dalys", "deaths", "ylds")) {
  d <- ac %>% filter(measure == m) %>% mutate(lrw_ = log(rw), lpop = log10(population))
  g0 <- gam(lrw_ ~ s(sdi, k = 6) + s(lpop, k = 6) + super_region, data = d, method = "REML")
  g1 <- gam(lrw_ ~ s(haq_value, k = 6) + s(lpop, k = 6) + super_region, data = d, method = "REML")
  g2 <- gam(lrw_ ~ s(sdi, k = 6) + s(lpop, k = 6), data = d, method = "REML")
  g3 <- gam(lrw_ ~ s(sdi, k = 6), data = d, method = "REML")
  g4 <- gam(lrw_ ~ super_region, data = d, method = "REML")
  cat(sprintf("\n%s: deviance explained  SDI only %.1f%% | region only %.1f%% | SDI+pop %.1f%% | SDI+pop+region %.1f%% | HAQ+pop+region %.1f%%\n",
              m, 100*summary(g3)$dev.expl, 100*summary(g4)$dev.expl, 100*summary(g2)$dev.expl,
              100*summary(g0)$dev.expl, 100*summary(g1)$dev.expl))
  cat("  Spearman rw~SDI:", round(cor(d$rw, d$sdi, method = "spearman"), 3),
      " rw~HAQ:", round(cor(d$rw, d$haq_value, method = "spearman"), 3),
      " rw~log pop:", round(cor(d$rw, d$lpop, method = "spearman"), 3), "\n")
  # linear summary for interpretability: ratio of rw between SDI 0.3 and 0.8 at median pop
  pr <- predict(g2, newdata = data.frame(sdi = c(0.3, 0.8), lpop = median(d$lpop)))
  cat(sprintf("  GAM (SDI+pop): rw at SDI 0.3 is %.2f times rw at SDI 0.8 (median population)\n", exp(pr[1] - pr[2])))
  pr2 <- predict(g2, newdata = data.frame(sdi = median(d$sdi), lpop = c(6, 8)))
  cat(sprintf("  GAM (SDI+pop): rw at population 1e6 is %.2f times rw at 1e8 (median SDI)\n", exp(pr2[1] - pr2[2])))
  print(summary(g0)$s.table)
}
# Figure 3: rw vs SDI and vs population
d3 <- ac %>% filter(measure %in% c("dalys", "deaths", "ylds"))
p3a <- ggplot(d3, aes(sdi, rw)) + geom_point(aes(colour = super_region), size = 1.3, alpha = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), colour = "black", linewidth = 0.6, se = TRUE) +
  facet_wrap(~ measure_f, nrow = 1, scales = "free_y") + scale_colour_manual(values = sr_pal, name = NULL) +
  scale_y_continuous(trans = "log10") + labs(x = "Socio-demographic Index, 2023", y = "Relative 95% UI width")
p3b <- ggplot(d3, aes(population, rw)) + geom_point(aes(colour = super_region), size = 1.3, alpha = 0.8) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), colour = "black", linewidth = 0.6, se = TRUE) +
  facet_wrap(~ measure_f, nrow = 1, scales = "free_y") + scale_colour_manual(values = sr_pal, name = NULL) +
  scale_x_continuous(trans = "log10", labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(trans = "log10") + labs(x = "Population, 2023", y = "Relative 95% UI width")
p3 <- (p3a / p3b) + plot_layout(guides = "collect") & theme(legend.position = "bottom") & guides(colour = guide_legend(nrow = 3))
fig(p3, "fig3_rw_vs_sdi_population.png", 10, 7.5)

# =============================================================================
# 3. Cause structure: is uncertainty a property of places or of causes?
# =============================================================================
cat("\n== 3. Cause level DALYs 2023 (30 level 3 causes x 204 countries) ==\n")
cl <- cc %>% filter(cause_id != 294, year == 2023, measure == "dalys", lower > 0, point > 0) %>%
  mutate(lrw_ = log(rw))
cat("cells:", nrow(cl), " dropped (lower<=0):", sum(cc$cause_id != 294 & cc$year == 2023 & cc$measure == "dalys" & !(cc$lower > 0 & cc$point > 0)), "\n")
# two way additive decomposition of log rw
a2 <- anova(lm(lrw_ ~ factor(cause_id) + location_name, data = cl))
ss <- a2$`Sum Sq`; names(ss) <- rownames(a2)
cat("Variance (SS) share of log rw: cause", round(100*ss[1]/sum(ss), 1), "% | country",
    round(100*ss[2]/sum(ss), 1), "% | residual", round(100*ss[3]/sum(ss), 1), "%\n")
a2b <- anova(lm(lrw_ ~ location_name + factor(cause_id), data = cl))
ssb <- a2b$`Sum Sq`
cat("(order reversed) country", round(100*ssb[1]/sum(ssb), 1), "% | cause", round(100*ssb[2]/sum(ssb), 1), "%\n")
cat("\nMedian rw by cause (DALYs 2023), with Spearman rw~SDI across countries:\n")
cs <- cl %>% group_by(cause_name, cause_type) %>%
  summarise(n = n(), median_rw = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75),
            rho_sdi = cor(rw, sdi, method = "spearman"), rho_rate = cor(rw, point, method = "spearman"), .groups = "drop") %>%
  arrange(desc(median_rw))
print(as.data.frame(cs), digits = 3)
cat("\nMedian rw by cause type:\n")
print(as.data.frame(cl %>% group_by(cause_type) %>% summarise(median_rw = median(rw), n = n())), digits = 3)

# deaths vs DALYs vs prevalence for the same cause and country
cat("\n== 3b. Fatal versus nonfatal: deaths, DALYs, prevalence for the same cause and country ==\n")
cmp <- cc %>% filter(cause_id != 294, year == 2023, lower > 0) %>% select(location_name, cause_name, cause_type, sdi, measure, rw) %>%
  pivot_wider(names_from = measure, values_from = rw)
cmp_dd <- cmp %>% filter(!is.na(deaths), !is.na(dalys))
cat("causes with both deaths and DALYs:", n_distinct(cmp_dd$cause_name), " cells:", nrow(cmp_dd), "\n")
cat("median ratio rw(deaths)/rw(DALYs):", round(median(cmp_dd$deaths / cmp_dd$dalys), 2),
    " share with deaths wider:", round(mean(cmp_dd$deaths > cmp_dd$dalys), 2), "\n")
print(as.data.frame(cmp_dd %>% group_by(cause_name) %>% summarise(med_ratio = median(deaths/dalys), n = n()) %>% arrange(med_ratio)), digits = 2)
cmp_pd <- cmp %>% filter(!is.na(prevalence), !is.na(dalys))
cat("\ncauses with prevalence and DALYs:", n_distinct(cmp_pd$cause_name), "\n")
print(as.data.frame(cmp_pd %>% group_by(cause_name) %>% summarise(med_rw_prev = median(prevalence), med_rw_daly = median(dalys), n = n())), digits = 2)

# Figure 4: heatmap of log rw, countries ordered by SDI, causes ordered by median rw
ord_c <- cs$cause_name
ord_l <- cl %>% distinct(location_name, sdi) %>% arrange(sdi) %>% pull(location_name)
hm <- cl %>% mutate(cause_name = factor(cause_name, rev(ord_c)), location_name = factor(location_name, ord_l))
p4 <- ggplot(hm, aes(location_name, cause_name, fill = rw)) + geom_tile() +
  scale_fill_viridis_c(trans = "log10", name = "Relative 95% UI width",
                       breaks = c(0.1, 0.3, 1, 3), labels = c("0.1", "0.3", "1", "3"), oob = squish, limits = c(0.08, 4)) +
  labs(x = "204 countries, ordered by Socio-demographic Index (low to high)", y = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.key.width = unit(1.6, "cm"))
fig(p4, "fig4_heatmap_cause_country_rw.png", 11, 7)

# =============================================================================
# 4. Time: has uncertainty narrowed or widened across the 1990 to 2023 series?
# =============================================================================
cat("\n== 4. All-cause DALYs, relative width by year ==\n")
tt <- cc %>% filter(cause_id == 294, measure == "dalys")
print(as.data.frame(tt %>% group_by(year) %>% summarise(median_rw = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75))), digits = 3)
cat("\nBy SDI quintile and year (median rw):\n")
print(as.data.frame(tt %>% group_by(sdi_quintile, year) %>% summarise(m = median(rw), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = m)), digits = 3)
ttw <- tt %>% select(location_name, super_region, sdi, year, rw) %>% pivot_wider(names_from = year, values_from = rw, names_prefix = "y") %>%
  mutate(ratio_1990_2023 = y1990 / y2023)
cat("\nratio rw(1990)/rw(2023): median", round(median(ttw$ratio_1990_2023), 2),
    " share >1:", round(mean(ttw$ratio_1990_2023 > 1), 2), "\n")
cat("Spearman ratio_1990_2023 ~ SDI:", round(cor(ttw$ratio_1990_2023, ttw$sdi, method = "spearman"), 3), "\n")
cat("rank stability: Spearman rw1990 ~ rw2023:", round(cor(ttw$y1990, ttw$y2023, method = "spearman"), 3), "\n")
cat("\nAll-cause by measure and year (median rw, IQR):\n")
tm <- cc %>% filter(cause_id == 294)
print(as.data.frame(tm %>% group_by(measure_f, year) %>%
  summarise(median_rw = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75), .groups = "drop")), digits = 3)
cat("\nDeaths: median rw by super region and year:\n")
print(as.data.frame(tm %>% filter(measure == "deaths") %>% group_by(super_region, year) %>%
  summarise(m = median(rw), .groups = "drop") %>% pivot_wider(names_from = year, values_from = m)), digits = 3)
tdw <- tm %>% filter(measure == "deaths") %>% select(location_name, sdi, super_region, year, rw) %>%
  pivot_wider(names_from = year, values_from = rw, names_prefix = "y")
cat("\nDeaths recency: ratio rw(2023)/rw(2019): median", round(median(tdw$y2023/tdw$y2019), 2),
    " IQR", round(quantile(tdw$y2023/tdw$y2019, c(.25, .75)), 2),
    " share > 1:", round(mean(tdw$y2023 > tdw$y2019), 2), "\n")
cat("Deaths recency: ratio rw(2023)/rw(2022): median", round(median(tdw$y2023/tdw$y2022), 2), "\n")
cat("Deaths: ratio rw(2019)/rw(1990): median", round(median(tdw$y2019/tdw$y1990), 2), "\n")
cat("Spearman ratio(2023/2019) ~ SDI:", round(cor(tdw$y2023/tdw$y2019, tdw$sdi, method = "spearman"), 3), "\n")
cat("Deaths 2023/2019 ratio by super region (median):\n")
print(as.data.frame(tdw %>% group_by(super_region) %>% summarise(med_ratio = median(y2023/y2019), n = n())), digits = 3)
tdd <- tm %>% filter(measure == "dalys") %>% select(location_name, year, rw) %>% pivot_wider(names_from = year, values_from = rw, names_prefix = "y")
cat("DALYs recency: ratio rw(2023)/rw(2019): median", round(median(tdd$y2023/tdd$y2019), 2), "\n")
tyd <- tm %>% filter(measure == "ylds") %>% select(location_name, year, rw) %>% pivot_wider(names_from = year, values_from = rw, names_prefix = "y")
cat("YLDs: rw 1990 median", round(median(tyd$y1990), 3), " 2010", round(median(tyd$y2010), 3), " 2023", round(median(tyd$y2023), 3),
    "; IQR 2023", round(quantile(tyd$y2023, c(.25,.75)), 3), "\n")
p5 <- ggplot(tm %>% filter(measure != "prevalence"), aes(year, rw, group = location_name, colour = super_region)) +
  geom_line(alpha = 0.3, linewidth = 0.25) +
  stat_summary(aes(group = super_region), fun = median, geom = "line", linewidth = 1.1) +
  facet_wrap(~ measure_f, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = sr_pal, name = NULL) + scale_y_continuous(trans = "log10") +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2020)) +
  labs(x = NULL, y = "Relative 95% UI width, all-cause") + guides(colour = guide_legend(nrow = 3))
fig(p5, "fig5_rw_time_allcause.png", 11, 5)

# =============================================================================
# 5. Subnational: does uncertainty inherit the national level or vary within?
# =============================================================================
cat("\n== 5. Admin1 all-cause DALYs 2023 versus national ==\n")
nat_key <- c(brazil = "Brazil", ethiopia = "Ethiopia", india = "India", indonesia = "Indonesia",
             iran = "Iran (Islamic Republic of)", italy = "Italy", japan = "Japan", kenya = "Kenya",
             mexico = "Mexico", nigeria = "Nigeria", norway = "Norway", pakistan = "Pakistan",
             philippines = "Philippines", poland = "Poland", south_africa = "South Africa",
             uk = "United Kingdom", usa = "United States of America")
nat <- ac %>% filter(measure == "dalys", location_name %in% nat_key) %>%
  mutate(country = names(nat_key)[match(location_name, nat_key)]) %>% select(country, nat_rw = rw, nat_point = point, sdi)
ads <- ad %>% group_by(country) %>%
  summarise(n_units = n(), med_rw = median(rw), min_rw = min(rw), max_rw = max(rw),
            cv_rw = sd(rw)/mean(rw), max_over_min = max(rw)/min(rw),
            rho_rw_count = cor(rw, count_combined, method = "spearman"), .groups = "drop") %>%
  left_join(nat, by = "country") %>% mutate(med_over_nat = med_rw / nat_rw) %>% arrange(desc(med_over_nat))
print(as.data.frame(ads), digits = 3)
cat("\nPooled: median admin1 rw / national rw =", round(median(ads$med_over_nat), 2),
    "; countries where every admin1 unit is wider than national:", sum(ads$min_rw > ads$nat_rw), "of 17\n")
cat("Spearman admin1 rw ~ DALY count (unit size), pooled within country (median rho):", round(median(ads$rho_rw_count), 2), "\n")
# ICC like decomposition
a5 <- anova(lm(log(rw) ~ country, data = ad)); s5 <- a5$`Sum Sq`
cat("Share of admin1 log rw variance between countries:", round(100*s5[1]/sum(s5), 1), "%\n")
cat("\nAdmin1 high SBP SEV 2023 (relative width):\n")
ass <- as %>% group_by(country) %>% summarise(n = n(), med_rw = median(rw), min_rw = min(rw), max_rw = max(rw), max_over_min = max(rw)/min(rw)) %>% arrange(desc(max_over_min))
print(as.data.frame(ass), digits = 3)
cat("\nUnits with widest SEV intervals:\n")
print(as.data.frame(as %>% arrange(desc(rw)) %>% select(country, location_name, point, lower, upper, rw) %>% head(10)), digits = 3)
ad2 <- ad %>% mutate(country_lab = nat_key[country], country_lab = factor(country_lab, levels = ads %>% arrange(med_rw) %>% pull(country) %>% nat_key[.]))
nat2 <- nat %>% mutate(country_lab = factor(nat_key[country], levels = levels(ad2$country_lab)))
p6 <- ggplot(ad2, aes(rw, country_lab)) +
  geom_point(size = 1.1, alpha = 0.6, colour = "grey40") +
  geom_point(data = nat2, aes(nat_rw, country_lab), shape = 23, size = 2.8, fill = "#CC79A7", colour = "black") +
  scale_x_continuous(trans = "log10") +
  labs(x = "Relative 95% UI width, all-cause DALYs 2023", y = NULL)
fig(p6, "fig6_admin1_vs_national_rw.png", 8, 5.5)

# =============================================================================
# 6. Risk factors: exposure (SEV) and attributable burden
# =============================================================================
cat("\n== 6. Risk exposure (SEV) and risk attributable DALYs, 2023 ==\n")
print(as.data.frame(sev %>% group_by(risk_name) %>%
  summarise(median_rw = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75),
            rho_sdi = cor(rw, sdi, method = "spearman"), rho_sev = cor(rw, point, method = "spearman")) %>% arrange(desc(median_rw))), digits = 3)
cat("\nRisk attributable DALYs (rows with lower <= 0 excluded from rw summaries):", sum(rd$lower <= 0), "\n")
print(as.data.frame(rd %>% filter(lower > 0) %>% group_by(risk_name) %>%
  summarise(n = n(), median_rw = median(rw), q25 = quantile(rw, .25), q75 = quantile(rw, .75),
            rho_sdi = cor(rw, sdi, method = "spearman")) %>% arrange(desc(median_rw))), digits = 3)
cat("\nAll-cause DALY rw for comparison (median):", round(median(d23$rw), 3), "\n")
cat("\nSEV rw by super region, high systolic blood pressure:\n")
print(as.data.frame(sev %>% filter(rei_id == 107) %>% group_by(super_region) %>% summarise(median_rw = median(rw), n = n())), digits = 3)
rr <- bind_rows(sev %>% mutate(kind = "Exposure (SEV)") %>% select(risk_name, super_region, rw, kind),
                rd %>% filter(lower > 0) %>% mutate(kind = "Attributable DALYs") %>% select(risk_name, super_region, rw, kind))
rr <- rr %>% mutate(risk_name = factor(risk_name, levels = rr %>% group_by(risk_name) %>% summarise(m = median(rw)) %>% arrange(m) %>% pull(risk_name)))
p7 <- ggplot(rr, aes(rw, risk_name, fill = super_region)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.25, position = position_dodge(preserve = "single")) +
  facet_wrap(~ kind, scales = "free_y") + scale_fill_manual(values = sr_pal, name = NULL) +
  scale_x_continuous(trans = "log10") + labs(x = "Relative 95% UI width", y = NULL) + guides(fill = guide_legend(nrow = 3))
fig(p7, "fig7_risk_rw.png", 11, 7)

# =============================================================================
# 7. Asymmetry of intervals (shape of the draw distribution)
# =============================================================================
cat("\n== 7. Interval asymmetry (upper - point)/(point - lower), all-cause and cause level 2023 ==\n")
print(as.data.frame(ac %>% group_by(measure_f) %>% summarise(median_asym = median(asym), q25 = quantile(asym, .25), q75 = quantile(asym, .75), p_gt1 = mean(asym > 1))), digits = 3)
print(as.data.frame(cl %>% summarise(median_asym = median(asym), q25 = quantile(asym, .25), q75 = quantile(asym, .75), p_gt1 = mean(asym > 1))), digits = 3)
cat("Spearman asym ~ rw (cause level):", round(cor(cl$asym, cl$rw, method = "spearman"), 3), "\n")

# =============================================================================
# 8. A single country scorecard: Ghana (home country interest) and Australia
# =============================================================================
cat("\n== 8. Scorecards ==\n")
for (ctry in c("Ghana", "Australia", "Nigeria")) {
  x <- cc %>% filter(location_name == ctry, year == 2023, measure == "dalys") %>% arrange(desc(rw)) %>%
    select(cause_name, point, lower, upper, rw)
  cat("\n", ctry, " all-cause DALY rw rank (1 = widest of 204):", which(d23$location_name == ctry), " rw =", round(d23$rw[d23$location_name == ctry], 3), "\n")
  print(as.data.frame(x %>% filter(cause_name != "All causes") %>% head(6)), digits = 3)
}
sink()
write_csv(cs, "output/cause_summary_dalys2023.csv")
write_csv(ads, "output/admin1_summary.csv")
write_csv(d23 %>% select(location_name, super_region, sdi, haq_value, population, point, lower, upper, rw, asym), "output/allcause_dalys_2023_rw.csv")
cat("figures written:", paste(list.files("figures"), collapse = ", "), "\n")
