#!/usr/bin/env Rscript
############################################################
## ED/Inpatient Fraction: Severity vs Utilization Composition
## Focus:  2020 encounters by care setting (pre-lockdown, lockdown, reopening)
## Tests whether the spore-fraction association survives a 2020 period indicator
##   - Daily binomial GAM (ED/inpatient fraction) + NB GAMs (ED/inpatient and
##     other counts modeled separately), with and without period2020
##   - Controls: cyclic day-of-year spline + day-of-week + weather
## Also: per-day encounters by period/setting, telehealth ascertainment,
##   2020 age-stratified acute ORs with FDR status, 2020 ED/inpatient case-crossover
##
## Outputs (CSV) written to: R2_ed_inpatient_fraction/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(mgcv)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_ed_inpatient_fraction/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

claims <- read_csv("data/2022_asma_From_date_ver 3.csv", show_col_types = FALSE) |>
  mutate(
    from_date = suppressWarnings(mdy(FROM_DATE)),
    year = year(from_date),
    muni = toupper(MUNICIPALITY_DESC),
    icd1 = toupper(trimws(ICD_DIAG_01_PRIMARY)),
    hcg_setting = as.character(HCG_SETTING),
    hcg_mr_line_desc = as.character(HCG_MR_LINE_DESC),
    line_group_desc = as.character(LINE_GROUP_DESC),
    form_type = as.character(FORM_TYPE),
    care_setting = case_when(
      grepl("ed|emergency", tolower(hcg_setting)) | grepl("\\bed\\b|\\ber\\b", tolower(hcg_mr_line_desc)) ~ "ED",
      grepl("inpatient|hospital", tolower(hcg_setting)) | grepl("inpatient|hospital", tolower(hcg_mr_line_desc)) ~ "Inpatient",
      grepl("urgent", tolower(hcg_setting)) | grepl("urgent", tolower(hcg_mr_line_desc)) ~ "Urgent",
      grepl("office|clinic|outpatient", tolower(hcg_setting)) | grepl("office|clinic", tolower(hcg_mr_line_desc)) ~ "Outpatient",
      grepl("dme|supplies", tolower(hcg_mr_line_desc)) ~ "DME",
      TRUE ~ "Other"
    ),
    asthma_primary = grepl("^J45|^J46", icd1),
    telehealth_flag = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone",
                            paste(hcg_setting, hcg_mr_line_desc, line_group_desc, form_type),
                            ignore.case = TRUE),
    period = case_when(
      year < 2020 ~ "pre2020",
      year == 2020 ~ "y2020",
      year >= 2021 ~ "y2021plus"
    )
  ) |>
  filter(!is.na(from_date),
         muni %in% c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA"))

clinical_settings <- c("ED", "Inpatient", "Outpatient", "Urgent")

asthma <- claims |>
  filter(asthma_primary, care_setting %in% clinical_settings) |>
  mutate(
    acute_flag = care_setting %in% c("ED", "Inpatient"),
    other_flag = care_setting %in% c("Outpatient", "Urgent")
  )

period_days <- tibble(
  period = c("pre2020", "y2020", "y2021plus"),
  n_days = c(
    length(seq.Date(as.Date("2018-01-01"), as.Date("2019-12-31"), by = "day")),
    length(seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")),
    length(seq.Date(as.Date("2021-01-01"), as.Date("2022-12-31"), by = "day"))
  )
)

enc <- asthma |>
  group_by(period) |>
  summarise(
    ed_inpatient_n = sum(acute_flag),
    other_n = sum(other_flag),
    total_n = n(),
    telehealth_n = sum(telehealth_flag, na.rm = TRUE),
    telehealth_other_n = sum(telehealth_flag & other_flag, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(period_days, by = "period") |>
  mutate(
    ed_inpatient_per_day = ed_inpatient_n / n_days,
    other_per_day = other_n / n_days,
    total_per_day = total_n / n_days,
    ed_inpatient_encounter_fraction = ed_inpatient_n / (ed_inpatient_n + other_n),
    telehealth_share_total = telehealth_n / total_n,
    telehealth_share_other = telehealth_other_n / pmax(other_n, 1)
  )

write_csv(enc, file.path(out_dir, "table1_period_setting_per_day_draft.csv"))

tele_detail <- asthma |>
  mutate(
    tele_source_hcg_setting = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone", hcg_setting, ignore.case = TRUE),
    tele_source_line_desc = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone", hcg_mr_line_desc, ignore.case = TRUE),
    tele_source_line_group = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone", line_group_desc, ignore.case = TRUE),
    tele_source_form_type = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone", form_type, ignore.case = TRUE)
  ) |>
  group_by(period) |>
  summarise(
    n_total = n(),
    tele_any = sum(telehealth_flag),
    tele_from_hcg_setting = sum(tele_source_hcg_setting),
    tele_from_line_desc = sum(tele_source_line_desc),
    tele_from_line_group = sum(tele_source_line_group),
    tele_from_form_type = sum(tele_source_form_type),
    .groups = "drop"
  )

write_csv(tele_detail, file.path(out_dir, "telehealth_ascertainment_draft.csv"))

sp <- read_excel("data/Datos de esporas San Juan Jan 2018-2022.xlsx") |>
  transmute(date = as.Date(date), spores_10k = as.numeric(spores_total) / 10000)

wx <- read_csv("data/sju_weather_2018_2021.csv", show_col_types = FALSE) |>
  transmute(
    date = as.Date(date),
    tmean_c = as.numeric((tmax_c + tmin_c) / 2),
    rh_mean = as.numeric(rh_mean),
    wind_ms = as.numeric(wind_ms),
    prcp_mm = as.numeric(prcp_mm),
    dewpoint_c = as.numeric(dewpoint_c)
  )

daily2020 <- asthma |>
  filter(year(from_date) == 2020) |>
  group_by(date = from_date) |>
  summarise(
    ed_inpatient_n = sum(acute_flag),
    other_n = sum(other_flag),
    total_n = n(),
    .groups = "drop"
  ) |>
  right_join(tibble(date = seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")), by = "date") |>
  mutate(across(c(ed_inpatient_n, other_n, total_n), ~ tidyr::replace_na(.x, 0L))) |>
  left_join(sp, by = "date") |>
  left_join(wx, by = "date") |>
  mutate(
    dow = wday(date, label = TRUE, week_start = 1),
    doy = yday(date),
    period2020 = case_when(
      date < as.Date("2020-03-15") ~ "pre_lockdown",
      date <= as.Date("2020-06-30") ~ "lockdown",
      TRUE ~ "reopening"
    ),
    period2020 = factor(period2020, levels = c("pre_lockdown", "lockdown", "reopening"))
  ) |>
  filter(if_all(c(spores_10k, tmean_c, rh_mean, wind_ms, prcp_mm, dewpoint_c), is.finite))

fit_frac_base <- gam(
  cbind(ed_inpatient_n, other_n) ~ spores_10k + s(doy, bs = "cc", k = 6) + dow +
    tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
  family = binomial(),
  data = daily2020,
  method = "REML",
  knots = list(doy = c(0.5, 366.5))
)

fit_frac_period <- gam(
  cbind(ed_inpatient_n, other_n) ~ spores_10k + period2020 + s(doy, bs = "cc", k = 6) + dow +
    tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
  family = binomial(),
  data = daily2020,
  method = "REML",
  knots = list(doy = c(0.5, 366.5))
)

fit_ed_base <- gam(
  ed_inpatient_n ~ spores_10k + s(doy, bs = "cc", k = 6) + dow +
    tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
  family = nb(),
  data = daily2020,
  method = "REML",
  knots = list(doy = c(0.5, 366.5))
)

fit_ed_period <- gam(
  ed_inpatient_n ~ spores_10k + period2020 + s(doy, bs = "cc", k = 6) + dow +
    tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
  family = nb(),
  data = daily2020,
  method = "REML",
  knots = list(doy = c(0.5, 366.5))
)

fit_other_base <- gam(
  other_n ~ spores_10k + s(doy, bs = "cc", k = 6) + dow +
    tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
  family = nb(),
  data = daily2020,
  method = "REML",
  knots = list(doy = c(0.5, 366.5))
)

fit_other_period <- gam(
  other_n ~ spores_10k + period2020 + s(doy, bs = "cc", k = 6) + dow +
    tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
  family = nb(),
  data = daily2020,
  method = "REML",
  knots = list(doy = c(0.5, 366.5))
)

extract_spore <- function(fit, model_name, outcome_type) {
  sm <- summary(fit)$p.table
  b <- sm["spores_10k", "Estimate"]
  se <- sm["spores_10k", "Std. Error"]
  tibble(
    model = model_name,
    outcome = outcome_type,
    estimate = exp(b),
    ci_low = exp(b - 1.96 * se),
    ci_high = exp(b + 1.96 * se),
    p_value = sm["spores_10k", "Pr(>|z|)"]
  )
}

modsum <- bind_rows(
  extract_spore(fit_frac_base, "base", "ED/inpatient encounter fraction"),
  extract_spore(fit_frac_period, "with_2020_period_indicator", "ED/inpatient encounter fraction"),
  extract_spore(fit_ed_base, "base", "ED/inpatient count"),
  extract_spore(fit_ed_period, "with_2020_period_indicator", "ED/inpatient count"),
  extract_spore(fit_other_base, "base", "Other encounter count"),
  extract_spore(fit_other_period, "with_2020_period_indicator", "Other encounter count")
)

write_csv(modsum, file.path(out_dir, "models_2020_fraction_and_counts_period_indicator_draft.csv"))

t6_or <- read_csv("statistical_results/pooled_daily_LOGIT_spores_ma7_by_pandemic_phase_age_pphase_age.csv", show_col_types = FALSE)
t6_rr <- read_csv("statistical_results/pooled_daily_NB_spores_ma7_by_pandemic_phase_age_pphase_age.csv", show_col_types = FALSE)

wald_p <- function(est, se) 2 * pnorm(-abs(est / se))

joint <- bind_rows(
  t6_rr |> mutate(type = "RR", est = RR_pphase_age, p_nom = wald_p(log_est_pphase_age, log_se_pphase_age)),
  t6_or |> mutate(type = "OR", est = OR_pphase_age, p_nom = wald_p(log_est_pphase_age, log_se_pphase_age))
) |>
  mutate(p_bh_48 = p.adjust(p_nom, method = "BH"))

a2020 <- joint |>
  filter(type == "OR", pandemic_phase_pphase_age == "y2020") |>
  transmute(
    age_group = age_group_pphase_age,
    or = est,
    p_nominal = p_nom,
    p_bh_48,
    sig_nominal = p_nominal < 0.05,
    sig_bh = p_bh_48 < 0.05
  )

write_csv(a2020, file.path(out_dir, "y2020_age_stratified_acute_or_fdr_status_draft.csv"))

cc <- read_csv("statistical_results/case_crossover_ed_inpatient_results_cases.csv", show_col_types = FALSE)
write_csv(cc, file.path(out_dir, "case_crossover_ed_inpatient_2020_draft.csv"))

message("ED/inpatient encounter fraction drafts written to: ", out_dir)
