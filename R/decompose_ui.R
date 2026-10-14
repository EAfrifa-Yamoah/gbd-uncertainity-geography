# decompose_daly_ui: read a published GBD DALY uncertainty interval in terms of
# its fatal (YLL) and nonfatal (YLD) components, using only the published point
# estimates and 95% uncertainty intervals for the three quantities.
#
# Arguments: numeric vectors c(point, lower, upper) for DALYs, YLLs and YLDs,
# all on the same scale (rate per 100,000 or count) and for the same location,
# year, sex and age group.
#
# Returns a one row data frame with
#   rw_daly, rw_yll, rw_yld   relative widths (upper - lower) / point
#   yld_share                 YLD point estimate / DALY point estimate
#   yld_var_share             share of DALY interval variance from YLDs, assuming
#                             independent YLL and YLD draws
#   obs_over_indep            observed DALY half width / half width predicted
#                             under independence (1 = consistent with independence)
#   rho_implied               YLL to YLD draw correlation implied by the observed width
#   note                      a reading rule
decompose_daly_ui <- function(daly, yll, yld) {
  stopifnot(length(daly) == 3, length(yll) == 3, length(yld) == 3)
  rw <- function(v) (v[3] - v[2]) / v[1]
  hw <- function(v) (v[3] - v[2]) / 2
  hd <- hw(daly); hl <- hw(yll); hn <- hw(yld)
  indep <- sqrt(hl^2 + hn^2)
  out <- data.frame(rw_daly = rw(daly), rw_yll = rw(yll), rw_yld = rw(yld),
                    yld_share = yld[1] / daly[1],
                    yld_var_share = hn^2 / (hl^2 + hn^2),
                    obs_over_indep = hd / indep,
                    rho_implied = (hd^2 - hl^2 - hn^2) / (2 * hl * hn))
  out$note <- if (out$yld_var_share > 0.8) {
    paste("DALY interval is dominated by the nonfatal component;",
          "use the YLL interval to judge local mortality data.")
  } else if (out$yld_var_share < 0.2) {
    paste("DALY interval is dominated by the fatal component",
          "and reflects mortality data availability.")
  } else paste("DALY interval mixes fatal and nonfatal uncertainty",
              "in comparable amounts.")
  out
}
