#!/usr/bin/env Rscript
############################################################
## Time-Stratified Case-Crossover: Exposure Window Reanalysis
## Focus:  ED and Inpatient Asthma Events (2018-2021)
## Compares MA3 vs MA7 cumulative spores on native and common strata
## Controls for season/trend by construction; weather matched to window
## Also: case reconciliation counts, verified Table 7, and a 2020 DLNM
##   lag-response (lag 0-14) for Figure 3
##
## Outputs (CSV, PNG) written to: R2_window_casecrossover/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(zoo)
  library(survival)
  library(broom)
  library(dlnm)
  library(MASS)
  library(splines)
  library(ggplot2)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_window_casecrossover/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

san_juan_munis <- c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA")
years_use <- 2018:2021

care_setting_fun <- function(hcg_setting, hcg_mr_line_desc) {
  s <- tolower(as.character(hcg_setting)); l <- tolower(as.character(hcg_mr_line_desc))
  if (grepl("ed|emergency", s) || grepl("\\bed\\b|\\ber\\b", l)) return("ED")
  if (grepl("inpatient|hospital", s) || grepl("inpatient|hospital", l)) return("Inpatient")
  if (grepl("urgent", s) || grepl("urgent", l)) return("Urgent")
  if (grepl("office|clinic|outpatient", s) || grepl("office|clinic", l)) return("Outpatient")
  "Other"
}

claims <- read_csv("data/2022_asma_From_date_ver 3.csv", show_col_types = FALSE) |>
  mutate(
    from_date = suppressWarnings(mdy(FROM_DATE)),
    year = year(from_date),
    municipality = toupper(trimws(MUNICIPALITY_DESC)),
    icd1 = toupper(trimws(ICD_DIAG_01_PRIMARY)),
    care_setting = mapply(care_setting_fun, HCG_SETTING, HCG_MR_LINE_DESC)
  ) |>
  filter(!is.na(from_date), year %in% years_use, municipality %in% san_juan_munis,
         grepl("^J45|^J46", icd1), care_setting %in% c("ED", "Inpatient")) |>
  mutate(
    case_id = row_number(),
    month_case = month(from_date),
    dow_case = wday(from_date)
  )

find_refs <- function(case_date, case_month, case_dow, max_weeks = 3) {
  cand <- seq(case_date - weeks(max_weeks), case_date + weeks(max_weeks), by = "day")
  cand[month(cand) == case_month & wday(cand) == case_dow & cand != case_date]
}

sets <- lapply(seq_len(nrow(claims)), function(i) {
  r <- claims[i, ]
  refs <- find_refs(r$from_date, r$month_case, r$dow_case, 3)
  bind_rows(
    tibble(stratum_id = i, case_id = r$case_id, date = r$from_date, is_case = 1L, year = r$year),
    tibble(stratum_id = i, case_id = r$case_id, date = refs, is_case = 0L, year = year(refs))
  )
})

cc <- bind_rows(sets)

sp <- read_excel("data/Datos de esporas San Juan Jan 2018-2022.xlsx") |>
  transmute(date = as.Date(date), spores = as.numeric(spores_total))

wx <- read_csv("data/sju_weather_2018_2021.csv", show_col_types = FALSE) |>
  transmute(date = as.Date(date), tmax = as.numeric(tmax_c), rh = as.numeric(rh_mean), prcp = as.numeric(prcp_mm))

expo <- full_join(sp, wx, by = "date") |>
  arrange(date) |>
  mutate(
    spores_ma3 = rollmean(spores, 3, fill = NA, align = "right"),
    spores_ma7 = rollmean(spores, 7, fill = NA, align = "right"),
    tmax_ma3 = rollmean(tmax, 3, fill = NA, align = "right"),
    tmax_ma7 = rollmean(tmax, 7, fill = NA, align = "right"),
    rh_ma3 = rollmean(rh, 3, fill = NA, align = "right"),
    rh_ma7 = rollmean(rh, 7, fill = NA, align = "right"),
    prcp_ma3 = rollmean(prcp, 3, fill = NA, align = "right"),
    prcp_ma7 = rollmean(prcp, 7, fill = NA, align = "right")
  )

ccx <- cc |>
  left_join(expo, by = "date")

cc_ma3 <- ccx |>
  filter(!is.na(spores_ma3), !is.na(tmax_ma3), !is.na(rh_ma3), !is.na(prcp_ma3))

cc_ma7 <- ccx |>
  filter(!is.na(spores_ma7), !is.na(tmax_ma7), !is.na(rh_ma7), !is.na(prcp_ma7))

common_strata <- intersect(unique(cc_ma3$stratum_id), unique(cc_ma7$stratum_id))
cc_ma3_common <- cc_ma3 |>
  filter(stratum_id %in% common_strata)
cc_ma7_common <- cc_ma7 |>
  filter(stratum_id %in% common_strata)

fit_cc <- function(dat, exp_var, weather_vars, label, yr = NULL) {
  d <- if (is.null(yr)) dat else dat |>
    filter(year == yr)
  if (nrow(d) == 0 || n_distinct(d$stratum_id) < 20) {
    return(tibble(lag = label, year = ifelse(is.null(yr), NA, yr), or = NA_real_, low = NA_real_, high = NA_real_, p = NA_real_, n_cases = NA_integer_, n_strata = NA_integer_, n_obs = NA_integer_))
  }
  f <- as.formula(paste0("is_case ~ ", exp_var, " + ", paste(weather_vars, collapse = " + "), " + strata(stratum_id)"))
  m <- clogit(f, data = d)
  td <- tidy(m, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == exp_var)
  tibble(
    lag = label,
    year = ifelse(is.null(yr), NA, yr),
    or = td$estimate,
    low = td$conf.low,
    high = td$conf.high,
    p = td$p.value,
    n_cases = sum(d$is_case == 1),
    n_strata = n_distinct(d$stratum_id),
    n_obs = nrow(d)
  )
}

res_all <- bind_rows(
  fit_cc(cc_ma3, "spores_ma3", c("tmax_ma3", "rh_ma3", "prcp_ma3"), "MA3_native_allyears"),
  fit_cc(cc_ma7, "spores_ma7", c("tmax_ma7", "rh_ma7", "prcp_ma7"), "MA7_native_allyears"),
  fit_cc(cc_ma3_common, "spores_ma3", c("tmax_ma3", "rh_ma3", "prcp_ma3"), "MA3_common_allyears"),
  fit_cc(cc_ma7_common, "spores_ma7", c("tmax_ma7", "rh_ma7", "prcp_ma7"), "MA7_common_allyears")
)

for (yr in years_use) {
  res_all <- bind_rows(
    res_all,
    fit_cc(cc_ma3, "spores_ma3", c("tmax_ma3", "rh_ma3", "prcp_ma3"), "MA3_native", yr),
    fit_cc(cc_ma7, "spores_ma7", c("tmax_ma7", "rh_ma7", "prcp_ma7"), "MA7_native", yr),
    fit_cc(cc_ma3_common, "spores_ma3", c("tmax_ma3", "rh_ma3", "prcp_ma3"), "MA3_common", yr),
    fit_cc(cc_ma7_common, "spores_ma7", c("tmax_ma7", "rh_ma7", "prcp_ma7"), "MA7_common", yr)
  )
}

recon <- bind_rows(
  claims |>
    count(year, name = "acute_events") |>
    mutate(metric = "acute_events_filtered"),
  cc |>
    count(year, name = "expanded_rows_before_exposure_filter") |>
    mutate(metric = "expanded_rows_before_exposure_filter"),
  cc_ma3 |>
    count(year, name = "ma3_rows") |>
    mutate(metric = "ma3_rows"),
  cc_ma7 |>
    count(year, name = "ma7_rows") |>
    mutate(metric = "ma7_rows")
)

recon2 <- tibble(
  metric = c("matched_strata_total", "ma3_strata_total", "ma7_strata_total", "common_strata_total"),
  value = c(n_distinct(cc$stratum_id), n_distinct(cc_ma3$stratum_id), n_distinct(cc_ma7$stratum_id), length(common_strata))
)

write_csv(recon, file.path(out_dir, "case_crossover_reconciliation_counts_by_year_draft.csv"))
write_csv(recon2, file.path(out_dir, "case_crossover_reconciliation_counts_totals_draft.csv"))
write_csv(res_all, file.path(out_dir, "case_crossover_ma3_ma7_native_common_2018_2021_draft.csv"))

tab7 <- res_all |>
  filter(year == 2020, lag %in% c("MA3_native", "MA7_native")) |>
  mutate(
    lag = if_else(lag == "MA3_native", "2020_MA3", "2020_MA7"),
    or_per_10000 = or^10000,
    low_per_10000 = low^10000,
    high_per_10000 = high^10000,
    or_chr = format(or, digits = 16, scientific = FALSE),
    low_chr = format(low, digits = 16, scientific = FALSE),
    high_chr = format(high, digits = 16, scientific = FALSE),
    p_chr = format(p, digits = 16, scientific = FALSE)
  ) |>
  dplyr::select(lag, or_chr, low_chr, high_chr, p_chr, n_strata, n_obs, or_per_10000, low_per_10000, high_per_10000)

foot <- tibble(
  note = c(
    "n_strata and n_obs differ between MA3 and MA7 because 7-day windows require longer complete exposure histories than 3-day windows.",
    "In this reanalysis, MA3 strata = 4543 and MA7 strata = 4543 for 2020 native sets, but observations differ (15524 vs 13524) due to referent-day missingness under longer windows."
  )
)

write_csv(tab7, file.path(out_dir, "Table_7_case_crossover_verified_highprecision_draft.csv"))
write_csv(foot, file.path(out_dir, "Table_7_footnotes_draft.csv"))

# DLNM lag-response (Figure 3 draft)
claims2020 <- read_csv("data/2022_asma_From_date_ver 3.csv", show_col_types = FALSE) |>
  mutate(
    date = suppressWarnings(mdy(FROM_DATE)),
    year = year(date),
    municipality = toupper(trimws(MUNICIPALITY_DESC)),
    icd1 = toupper(trimws(ICD_DIAG_01_PRIMARY)),
    care_setting = mapply(care_setting_fun, HCG_SETTING, HCG_MR_LINE_DESC),
    acute_flag = care_setting %in% c("ED", "Inpatient")
  ) |>
  filter(!is.na(date), year == 2020, municipality %in% san_juan_munis,
         grepl("^J45|^J46", icd1), care_setting %in% c("ED", "Inpatient", "Outpatient", "Urgent"))

daily <- claims2020 |>
  group_by(date) |>
  summarise(ed_inpatient_n = sum(acute_flag), .groups = "drop") |>
  right_join(tibble(date = seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")), by = "date") |>
  mutate(ed_inpatient_n = tidyr::replace_na(ed_inpatient_n, 0L))

sp2 <- read_excel("data/Datos de esporas San Juan Jan 2018-2022.xlsx") |>
  transmute(date = as.Date(date), spores_10k = as.numeric(spores_total) / 10000)

wx2 <- read_csv("data/sju_weather_2018_2021.csv", show_col_types = FALSE) |>
  transmute(
    date = as.Date(date),
    tmean_c = as.numeric((tmax_c + tmin_c) / 2),
    rh_mean = as.numeric(rh_mean),
    prcp_mm = as.numeric(prcp_mm)
  )

d <- daily |>
  left_join(sp2, by = "date") |>
  left_join(wx2, by = "date") |>
  mutate(
    dow = wday(date, label = TRUE, week_start = 1),
    time_num = as.numeric(date - min(date)) + 1
  ) |>
  filter(if_all(c(spores_10k, tmean_c, rh_mean, prcp_mm), is.finite))

cb_sp <- crossbasis(d$spores_10k, lag = 14, argvar = list(fun = "ns", df = 3), arglag = list(fun = "ns", df = 3))

fit_dlnm <- glm.nb(
  ed_inpatient_n ~ cb_sp + ns(tmean_c, 3) + ns(rh_mean, 3) + ns(prcp_mm, 3) + dow + ns(time_num, 8),
  data = d,
  na.action = na.exclude
)

cen <- median(d$spores_10k, na.rm = TRUE)
at90 <- as.numeric(quantile(d$spores_10k, 0.90, na.rm = TRUE))
cp <- crosspred(cb_sp, fit_dlnm, cen = cen, at = at90, bylag = 1)

lag_df <- tibble(
  lag = 0:14,
  rr = as.numeric(cp$matRRfit[1, ]),
  low = as.numeric(cp$matRRlow[1, ]),
  high = as.numeric(cp$matRRhigh[1, ])
)

cum_df <- tibble(
  summary = "lag0_14_cumulative",
  rr = as.numeric(cp$allRRfit[1]),
  low = as.numeric(cp$allRRlow[1]),
  high = as.numeric(cp$allRRhigh[1]),
  at90 = at90,
  cen = cen
)

write_csv(lag_df, file.path(out_dir, "figure3_dlnm_lag_response_2020_draft.csv"))
write_csv(cum_df, file.path(out_dir, "figure3_dlnm_cumulative_0_14_2020_draft.csv"))

g <- ggplot(lag_df, aes(x = lag, y = rr)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
  geom_line() +
  geom_point(size = 1.4) +
  labs(
    x = "Lag day",
    y = "Rate ratio",
    title = "DLNM lag-response (2020 ED/inpatient asthma encounters)",
    subtitle = "90th percentile spores vs median"
  )

ggsave(file.path(out_dir, "Figure_3_dlnm_lag_response_2020_draft.png"), g, width = 7, height = 4.5, dpi = 300)

message("Window and case-crossover draft outputs written to: ", out_dir)
