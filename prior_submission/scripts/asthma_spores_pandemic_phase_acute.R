############################################################
## Spores ↔ Asthma models: Baseline vs Pandemic Phase (pphase_acute)
## San Juan health region (claims) + San Juan station spores + SJU weather
##
## A) Baseline (2018–2019) modeled jointly, with year-specific spores effect:
##    - Weekly NB (total claims) + weekly binomial logit (acute fraction)
##    - Exposure: weekly cumulative spores (sum of daily spores_10k)
##    - Controls: cyclic week-of-year spline + year indicators + weekly mean weather
##    - Key term: spores_cum_wk_pphase_acute × year_f_pphase_acute
##
## B) Pandemic/post-pandemic pooled daily models:
##    - pandemic_phase_pphase_acute: pre2020, y2020, y2021plus
##    - Exposure: spores_10k_ma7_pphase_acute × pandemic_phase_pphase_acute
##    - Controls: cyclic day-of-year spline + day-of-week + smooth time trend + weather
##
## Exposure QC (applied separately to spores and pollen BEFORE modeling):
##   - If a week has <= 3 missing days: impute missing days with that week’s mean
##   - If a week has  > 3 missing days OR all missing: exclude week (set all days in week to NA)
##
## Outputs (CSV) written to: statistical_results/
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(readxl)
  library(janitor)
  library(rnoaa)
  library(mgcv)
  library(MASS)
  library(splines)
  library(zoo)
})

# ===========================================================
# 1) SETTINGS -------
# ===========================================================

claims_path_pphase_acute  <- "data/2022_asma_From_date_ver 3.xlsx"
claims_sheet_pphase_acute <- 1

spore_files_sj_pphase_acute <- "data/Datos esporas y polen 24-08-19 a 31-12-22 Corregidas 12-09-25.xlsx"

# NOAA token (recommended: environment variable)
# Sys.setenv(NOAATOKEN="YOUR_TOKEN"); options(noaakey = Sys.getenv("NOAATOKEN"))
options(noaakey = Sys.getenv("NOAATOKEN"))

sju_station_id_pphase_acute <- "RQW00011641"

years_all_pphase_acute      <- 2018:2022
years_baseline_pphase_acute <- 2018:2019

# Health region mapping (EDIT to match your definition)
san_juan_munis_pphase_acute <- c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA")
caguas_munis_pphase_acute   <- c("CAGUAS", "GURABO", "JUNCOS", "SAN LORENZO")

clinical_settings_pphase_acute <- c("ED", "Inpatient", "Outpatient", "Urgent")

max_missing_days_week_pphase_acute <- 3

spore_scale_divisor_pphase_acute <- 10000  # spores_10k = spores_total / 10000

weather_covars_pphase_acute <- c("tmean_c_pphase_acute", "rh_mean_pphase_acute", "wind_ms_pphase_acute", "prcp_mm_pphase_acute", "dewpoint_c_pphase_acute")

out_dir_pphase_acute <- "statistical_results"
dir.create(out_dir_pphase_acute, showWarnings = FALSE, recursive = TRUE)

# ===========================================================
# 2) HELPERS  -----
# ===========================================================

parse_date_any_pphase_acute <- function(x_pphase_acute) {
  if (inherits(x_pphase_acute, "Date")) return(x_pphase_acute)
  if (is.numeric(x_pphase_acute)) return(as.Date(x_pphase_acute, origin = "1899-12-30")) # Excel numeric dates
  x_pphase_acute <- as.character(x_pphase_acute)
  d_pphase_acute <- suppressWarnings(ymd(x_pphase_acute))
  d_pphase_acute <- dplyr::coalesce(d_pphase_acute, suppressWarnings(mdy(x_pphase_acute)))
  d_pphase_acute <- dplyr::coalesce(d_pphase_acute, suppressWarnings(dmy(x_pphase_acute)))
  dplyr::coalesce(d_pphase_acute, suppressWarnings(as.Date(x_pphase_acute)))
}

safe_mean_pphase_acute <- function(x_pphase_acute) {
  m_pphase_acute <- mean(x_pphase_acute, na.rm = TRUE)
  if (is.nan(m_pphase_acute)) NA_real_ else m_pphase_acute
}

# ---- claims: care setting and asthma primary ----
care_setting_fun_pphase_acute <- function(hcg_setting_pphase_acute, hcg_mr_line_desc_pphase_acute) {
  s_pphase_acute <- tolower(as.character(hcg_setting_pphase_acute))
  l_pphase_acute <- tolower(as.character(hcg_mr_line_desc_pphase_acute))
  if (grepl("ed|emergency", s_pphase_acute) || grepl("\\bed\\b|\\ber\\b", l_pphase_acute)) return("ED")
  if (grepl("inpatient|hospital", s_pphase_acute) || grepl("inpatient|hospital", l_pphase_acute)) return("Inpatient")
  if (grepl("urgent", s_pphase_acute) || grepl("urgent", l_pphase_acute)) return("Urgent")
  if (grepl("office|clinic|outpatient", s_pphase_acute) || grepl("office|clinic", l_pphase_acute)) return("Outpatient")
  if (grepl("dme|supplies", l_pphase_acute)) return("DME")
  "Other"
}

ensure_care_setting_pphase_acute <- function(df_pphase_acute) {
  if ("care_setting_pphase_acute" %in% names(df_pphase_acute)) return(df_pphase_acute)

  if ("care_setting" %in% names(df_pphase_acute)) {
    return(df_pphase_acute %>% mutate(care_setting_pphase_acute = care_setting))
  }

  needed_pphase_acute <- c("hcg_setting", "hcg_mr_line_desc")
  if (!all(needed_pphase_acute %in% names(df_pphase_acute))) stop("No care_setting and missing hcg_setting/hcg_mr_line_desc to create it.")

  df_pphase_acute %>% mutate(care_setting_pphase_acute = mapply(care_setting_fun_pphase_acute, hcg_setting, hcg_mr_line_desc))
}

ensure_asthma_primary_pphase_acute <- function(df_pphase_acute) {
  if ("asthma_primary_pphase_acute" %in% names(df_pphase_acute)) return(df_pphase_acute)

  if ("asthma_primary" %in% names(df_pphase_acute)) {
    return(df_pphase_acute %>% mutate(asthma_primary_pphase_acute = asthma_primary))
  }

  if (!"icd_diag_01_primary" %in% names(df_pphase_acute)) stop("Missing asthma_primary and icd_diag_01_primary.")

  df_pphase_acute %>%
    mutate(
      icd_diag_01_primary = stringr::str_to_upper(stringr::str_trim(as.character(icd_diag_01_primary))),
      asthma_primary_pphase_acute = if_else(stringr::str_detect(icd_diag_01_primary, "^J45|^J46"), "asthma", "other")
    )
}

# ---- spores/pollen ----
read_spores_sj_onefile_pphase_acute <- function(path_pphase_acute) {
  x_pphase_acute <- read_excel(path_pphase_acute) %>% clean_names()
  date_col_pphase_acute <- if ("date" %in% names(x_pphase_acute)) "date" else names(x_pphase_acute)[1]
  spore_col_pphase_acute  <- grep("esporas|spore", names(x_pphase_acute), ignore.case = TRUE, value = TRUE)[1]
  pollen_col_pphase_acute <- grep("polen|pollen", names(x_pphase_acute), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(spore_col_pphase_acute)) stop("Could not identify a spores column in: ", path_pphase_acute)

  x_pphase_acute %>%
    transmute(
      date = parse_date_any_pphase_acute(.data[[date_col_pphase_acute]]),
      spores_total_raw_pphase_acute = suppressWarnings(as.numeric(.data[[spore_col_pphase_acute]])),
      pollen_total_raw_pphase_acute = if (!is.na(pollen_col_pphase_acute)) suppressWarnings(as.numeric(.data[[pollen_col_pphase_acute]])) else NA_real_
    ) %>%
    filter(!is.na(date))
}

weekly_qc_impute_or_exclude_pphase_acute <- function(df_pphase_acute, value_col_pphase_acute, base_name_pphase_acute, max_missing_pphase_acute = 3) {
  value_sym_pphase_acute <- rlang::ensym(value_col_pphase_acute)

  out_col_pphase_acute        <- paste0(base_name_pphase_acute, "_pphase_acute")
  imputed_col_pphase_acute    <- paste0(base_name_pphase_acute, "_imputed_pphase_acute")
  week_excl_col_pphase_acute  <- paste0(base_name_pphase_acute, "_week_excluded_pphase_acute")
  nmiss_col_pphase_acute      <- paste0(base_name_pphase_acute, "_n_missing_week_pphase_acute")
  weekmean_col_pphase_acute   <- paste0(base_name_pphase_acute, "_week_mean_pphase_acute")

  df_pphase_acute %>%
    mutate(week_start_pphase_acute = floor_date(date, unit = "week", week_start = 1)) %>%
    group_by(week_start_pphase_acute) %>%
    mutate(
      n_missing_week_tmp_pphase_acute = sum(is.na(!!value_sym_pphase_acute)),
      week_mean_tmp_pphase_acute      = safe_mean_pphase_acute(!!value_sym_pphase_acute),
      week_ok_tmp_pphase_acute        = (n_missing_week_tmp_pphase_acute <= max_missing_pphase_acute) & !is.na(week_mean_tmp_pphase_acute),
      week_excluded_tmp_pphase_acute  = !week_ok_tmp_pphase_acute,
      imputed_day_tmp_pphase_acute    = week_ok_tmp_pphase_acute & is.na(!!value_sym_pphase_acute),
      value_clean_tmp_pphase_acute    = case_when(
        week_ok_tmp_pphase_acute & is.na(!!value_sym_pphase_acute) ~ week_mean_tmp_pphase_acute,
        week_ok_tmp_pphase_acute ~ as.numeric(!!value_sym_pphase_acute),
        TRUE ~ NA_real_
      ),
      !!out_col_pphase_acute       := value_clean_tmp_pphase_acute,
      !!imputed_col_pphase_acute   := imputed_day_tmp_pphase_acute,
      !!week_excl_col_pphase_acute := week_excluded_tmp_pphase_acute,
      !!nmiss_col_pphase_acute     := n_missing_week_tmp_pphase_acute,
      !!weekmean_col_pphase_acute  := week_mean_tmp_pphase_acute
    ) %>%
    ungroup() %>%
    dplyr::select(-n_missing_week_tmp_pphase_acute, -week_mean_tmp_pphase_acute, -week_ok_tmp_pphase_acute,
           -week_excluded_tmp_pphase_acute, -imputed_day_tmp_pphase_acute, -value_clean_tmp_pphase_acute)
}

# ---- weather (unit handling) ----
needs_div10_pphase_acute <- function(x_pphase_acute, thr_pphase_acute) {
  x_pphase_acute <- suppressWarnings(as.numeric(x_pphase_acute))
  q_pphase_acute <- suppressWarnings(quantile(abs(x_pphase_acute), probs = 0.99, na.rm = TRUE, names = FALSE))
  if (is.na(q_pphase_acute)) return(FALSE)
  q_pphase_acute > thr_pphase_acute
}

pull_weather_sju_pphase_acute <- function(year_value_pphase_acute, station_id_pphase_acute) {
  vars_pphase_acute <- c("PRCP", "TMAX", "TMIN", "AWND", "ADPT", "RHAV")

  w_pphase_acute <- meteo_pull_monitors(
    monitors   = station_id_pphase_acute,
    var        = vars_pphase_acute,
    date_min   = sprintf("%d-01-01", year_value_pphase_acute),
    date_max   = sprintf("%d-12-31", year_value_pphase_acute),
    keep_flags = FALSE
  ) %>%
    as_tibble() %>%
    mutate(date = as.Date(date))

  for (nm_pphase_acute in c("prcp","tmax","tmin","awnd","adpt","rhav")) {
    if (!nm_pphase_acute %in% names(w_pphase_acute)) w_pphase_acute[[nm_pphase_acute]] <- NA_real_
  }

  prcp_mm_pphase_acute     <- if (needs_div10_pphase_acute(w_pphase_acute$prcp, thr_pphase_acute = 300)) w_pphase_acute$prcp / 10 else w_pphase_acute$prcp
  tmax_c_pphase_acute      <- if (needs_div10_pphase_acute(w_pphase_acute$tmax, thr_pphase_acute = 80))  w_pphase_acute$tmax / 10 else w_pphase_acute$tmax
  tmin_c_pphase_acute      <- if (needs_div10_pphase_acute(w_pphase_acute$tmin, thr_pphase_acute = 80))  w_pphase_acute$tmin / 10 else w_pphase_acute$tmin
  dewpoint_c_pphase_acute  <- if (needs_div10_pphase_acute(w_pphase_acute$adpt, thr_pphase_acute = 80))  w_pphase_acute$adpt / 10 else w_pphase_acute$adpt
  wind_ms_pphase_acute     <- if (needs_div10_pphase_acute(w_pphase_acute$awnd, thr_pphase_acute = 30))  w_pphase_acute$awnd / 10 else w_pphase_acute$awnd

  tibble(
    date = w_pphase_acute$date,
    prcp_mm_pphase_acute = prcp_mm_pphase_acute,
    tmax_c_pphase_acute  = tmax_c_pphase_acute,
    tmin_c_pphase_acute  = tmin_c_pphase_acute,
    tmean_c_pphase_acute = (tmax_c_pphase_acute + tmin_c_pphase_acute) / 2,
    dewpoint_c_pphase_acute = dewpoint_c_pphase_acute,
    rh_mean_pphase_acute = w_pphase_acute$rhav,
    wind_ms_pphase_acute = wind_ms_pphase_acute
  )
}

# ---- linear-combination CI for exp(beta) ----
lincomb_exp_pphase_acute <- function(beta_pphase_acute, V_pphase_acute, idx_pphase_acute, label_pphase_acute = "effect_pphase_acute") {
  b_pphase_acute <- sum(beta_pphase_acute[names(idx_pphase_acute)] * idx_pphase_acute, na.rm = TRUE)
  v_pphase_acute <- as.numeric(t(idx_pphase_acute) %*% V_pphase_acute[names(idx_pphase_acute), names(idx_pphase_acute), drop = FALSE] %*% idx_pphase_acute)
  se_pphase_acute <- sqrt(v_pphase_acute)
  ci_pphase_acute <- b_pphase_acute + c(-1, 1) * 1.96 * se_pphase_acute
  tibble(
    log_est_pphase_acute = b_pphase_acute,
    log_se_pphase_acute  = se_pphase_acute,
    !!label_pphase_acute := exp(b_pphase_acute),
    low_pphase_acute     = exp(ci_pphase_acute[1]),
    high_pphase_acute    = exp(ci_pphase_acute[2])
  )
}

get_group_effects_pphase_acute <- function(fit_pphase_acute, base_term_pphase_acute, group_var_pphase_acute, label_pphase_acute) {
  beta_pphase_acute <- coef(fit_pphase_acute)
  V_pphase_acute <- vcov(fit_pphase_acute)

  lvls_pphase_acute <- levels(model.frame(fit_pphase_acute)[[group_var_pphase_acute]])
  ref_pphase_acute  <- lvls_pphase_acute[1]

  map_dfr(lvls_pphase_acute, function(lv_pphase_acute) {
    idx_pphase_acute <- rep(0, length(beta_pphase_acute)); names(idx_pphase_acute) <- names(beta_pphase_acute)
    idx_pphase_acute[base_term_pphase_acute] <- 1

    if (lv_pphase_acute != ref_pphase_acute) {
      t1_pphase_acute <- paste0(base_term_pphase_acute, ":", group_var_pphase_acute, lv_pphase_acute)
      t2_pphase_acute <- paste0(group_var_pphase_acute, lv_pphase_acute, ":", base_term_pphase_acute)
      if (t1_pphase_acute %in% names(beta_pphase_acute)) idx_pphase_acute[t1_pphase_acute] <- 1
      if (t2_pphase_acute %in% names(beta_pphase_acute)) idx_pphase_acute[t2_pphase_acute] <- 1
    }

    est_pphase_acute <- lincomb_exp_pphase_acute(beta_pphase_acute, V_pphase_acute,
                                     idx_pphase_acute[names(idx_pphase_acute) %in% names(beta_pphase_acute)],
                                     label_pphase_acute = label_pphase_acute)
    tibble(level_pphase_acute = lv_pphase_acute) %>% bind_cols(est_pphase_acute)
  })
}

# ===========================================================
# 3) BUILD DAILY ASTHMA OUTCOMES (San Juan health region)-----
# ===========================================================

claims_raw_pphase_acute <- read_excel(claims_path_pphase_acute, sheet = claims_sheet_pphase_acute) %>%
  clean_names() %>%
  mutate(from_date = parse_date_any_pphase_acute(from_date)) %>%
  ensure_care_setting_pphase_acute() %>%
  ensure_asthma_primary_pphase_acute()

claims_pphase_acute <- claims_raw_pphase_acute %>%
  mutate(
    municipality_desc_pphase_acute = toupper(as.character(municipality_desc)),
    year_pphase_acute = year(from_date),
    health_region_pphase_acute = case_when(
      municipality_desc_pphase_acute %in% san_juan_munis_pphase_acute ~ "San Juan",
      municipality_desc_pphase_acute %in% caguas_munis_pphase_acute   ~ "Caguas",
      TRUE                                                ~ NA_character_
    ),
    health_region_pphase_acute = factor(health_region_pphase_acute, levels = c("San Juan", "Caguas"))
  ) %>%
  filter(year_pphase_acute %in% years_all_pphase_acute) %>%
  filter(!is.na(health_region_pphase_acute))

claims_sj_pphase_acute <- claims_pphase_acute %>%
  filter(
    health_region_pphase_acute == "San Juan",
    asthma_primary_pphase_acute == "asthma",
    care_setting_pphase_acute %in% clinical_settings_pphase_acute
  ) %>%
  mutate(acute_flag_pphase_acute = as.integer(care_setting_pphase_acute %in% c("ED", "Inpatient")))

calendar_pphase_acute <- tibble(
  date = seq.Date(as.Date(sprintf("%d-01-01", min(years_all_pphase_acute))),
                  as.Date(sprintf("%d-12-31", max(years_all_pphase_acute))),
                  by = "day")
) %>%
  mutate(
    year_pphase_acute = year(date),
    doy_pphase_acute  = yday(date),
    dow_pphase_acute  = wday(date, label = TRUE, week_start = 1),
    week_start_pphase_acute = floor_date(date, unit = "week", week_start = 1),
    week_of_year_pphase_acute = isoweek(date),
    date_num_pphase_acute = as.numeric(date)
  ) %>%
  filter(year_pphase_acute %in% years_all_pphase_acute)

daily_asthma_pphase_acute <- claims_sj_pphase_acute %>%
  group_by(date = from_date) %>%
  summarise(
    total_claims_pphase_acute = n(),
    acute_claims_pphase_acute = sum(acute_flag_pphase_acute == 1, na.rm = TRUE),
    nonacute_claims_pphase_acute = total_claims_pphase_acute - acute_claims_pphase_acute,
    .groups = "drop"
  )

daily_outcome_pphase_acute <- calendar_pphase_acute %>%
  left_join(daily_asthma_pphase_acute, by = "date") %>%
  mutate(
    across(c(total_claims_pphase_acute, acute_claims_pphase_acute, nonacute_claims_pphase_acute),
           ~ replace_na(.x, 0L))
  )

# ===========================================================
# 4) BUILD DAILY SPORES/POLLEN with weekly QC (San Juan station)-----
# ===========================================================

spores_raw_pphase_acute <- map_dfr(spore_files_sj_pphase_acute, read_spores_sj_onefile_pphase_acute) %>%
  group_by(date) %>%
  summarise(
    spores_total_raw_pphase_acute  = safe_mean_pphase_acute(spores_total_raw_pphase_acute),
    pollen_total_raw_pphase_acute  = safe_mean_pphase_acute(pollen_total_raw_pphase_acute),
    .groups = "drop"
  ) %>%
  arrange(date)

spores_daily_pphase_acute <- tibble(date = calendar_pphase_acute$date) %>%
  left_join(spores_raw_pphase_acute, by = "date") %>%
  arrange(date)

spores_daily_pphase_acute <- weekly_qc_impute_or_exclude_pphase_acute(
  spores_daily_pphase_acute,
  value_col_pphase_acute = spores_total_raw_pphase_acute,
  base_name_pphase_acute = "spores_total",
  max_missing_pphase_acute = max_missing_days_week_pphase_acute
)

spores_daily_pphase_acute <- weekly_qc_impute_or_exclude_pphase_acute(
  spores_daily_pphase_acute,
  value_col_pphase_acute = pollen_total_raw_pphase_acute,
  base_name_pphase_acute = "pollen_total",
  max_missing_pphase_acute = max_missing_days_week_pphase_acute
)

spores_daily_pphase_acute <- spores_daily_pphase_acute %>%
  mutate(
    spores_10k_pphase_acute = spores_total_pphase_acute / spore_scale_divisor_pphase_acute
  )

# ===========================================================
# 5) WEATHER + MERGE into DAILY dataset + MA7 + pandemic phase -----
# ===========================================================

weather_daily_pphase_acute <- map_dfr(years_all_pphase_acute, ~ pull_weather_sju_pphase_acute(.x, sju_station_id_pphase_acute))

daily_pphase_acute <- daily_outcome_pphase_acute %>%
  left_join(spores_daily_pphase_acute,  by = "date") %>%
  left_join(weather_daily_pphase_acute, by = "date") %>%
  mutate(
    pandemic_phase_pphase_acute = case_when(
      year_pphase_acute < 2020  ~ "pre2020",
      year_pphase_acute == 2020 ~ "y2020",
      year_pphase_acute >= 2021 ~ "y2021plus"
    ),
    pandemic_phase_pphase_acute = factor(pandemic_phase_pphase_acute, levels = c("pre2020", "y2020", "y2021plus"))
  ) %>%
  arrange(date) %>%
  mutate(
    spores_10k_ma7_pphase_acute = zoo::rollapply(spores_10k_pphase_acute, width = 7, FUN = mean,
                                          align = "right", fill = NA_real_, partial = FALSE)
  )

weather_covars_use_pphase_acute <- intersect(weather_covars_pphase_acute, names(daily_pphase_acute))
if (length(weather_covars_use_pphase_acute) == 0) stop("No weather covariates found after weather pull.")

# ===========================================================
# A) BASELINE (2018–2019): weekly NB + weekly logit with spores × year -----
# ===========================================================

weekly_baseline_pphase_acute <- daily_pphase_acute %>%
  filter(year_pphase_acute %in% years_baseline_pphase_acute) %>%
  filter(week_of_year_pphase_acute <= 52) %>%
  mutate(week_start_pphase_acute = floor_date(date, unit = "week", week_start = 1)) %>%
  group_by(year_pphase_acute, week_start_pphase_acute, week_of_year_pphase_acute) %>%
  summarise(
    total_claims_wk_pphase_acute    = sum(total_claims_pphase_acute),
    acute_claims_wk_pphase_acute    = sum(acute_claims_pphase_acute),
    nonacute_claims_wk_pphase_acute = sum(nonacute_claims_pphase_acute),
    spores_cum_wk_pphase_acute = sum(spores_10k_pphase_acute, na.rm = FALSE),
    across(all_of(weather_covars_use_pphase_acute), ~ safe_mean_pphase_acute(.x), .names = "{.col}_wkmean_pphase_acute"),
    spores_week_excluded_pphase_acute = any(spores_total_week_excluded_pphase_acute, na.rm = TRUE),
    pollen_week_excluded_pphase_acute = any(pollen_total_week_excluded_pphase_acute, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!spores_week_excluded_pphase_acute) %>%
  filter(!is.na(spores_cum_wk_pphase_acute)) %>%
  mutate(
    year_f_pphase_acute = factor(year_pphase_acute)
  )

# Set reference level to the first available year
if (nrow(weekly_baseline_pphase_acute) > 0) {
  available_years_pphase_acute <- sort(unique(weekly_baseline_pphase_acute$year_pphase_acute))
  
  # Check if we have at least 2 years for year-specific models
  if (length(available_years_pphase_acute) >= 2) {
    weekly_baseline_pphase_acute <- weekly_baseline_pphase_acute %>%
      mutate(year_f_pphase_acute = relevel(factor(year_pphase_acute), ref = as.character(available_years_pphase_acute[1])))
  } else {
    cat("\nWarning: Only", length(available_years_pphase_acute), "year(s) available after filtering. Fitting pooled baseline model without year interaction.\n")
    weekly_baseline_pphase_acute <- weekly_baseline_pphase_acute %>%
      mutate(year_f_pphase_acute = factor(year_pphase_acute))
  }
}

knots_week_pphase_acute <- list(week_of_year_pphase_acute = c(0.5, 52.5))
wk_weather_cols_pphase_acute <- names(weekly_baseline_pphase_acute)[stringr::str_ends(names(weekly_baseline_pphase_acute), "_wkmean_pphase_acute")]

# Dynamically adjust k based on available data
n_weeks_pphase_acute <- nrow(weekly_baseline_pphase_acute)
k_week_pphase_acute <- min(20, max(3, floor(n_weeks_pphase_acute * 0.8)))

# Build formula based on number of years available
if (length(available_years_pphase_acute) >= 2) {
  f_wk_nb_pphase_acute <- as.formula(
    paste0(
      "total_claims_wk_pphase_acute ~ spores_cum_wk_pphase_acute * year_f_pphase_acute + s(week_of_year_pphase_acute, bs='cc', k=", k_week_pphase_acute, ")",
      if (length(wk_weather_cols_pphase_acute) > 0) paste0(" + ", paste(wk_weather_cols_pphase_acute, collapse = " + ")) else ""
    )
  )
} else {
  f_wk_nb_pphase_acute <- as.formula(
    paste0(
      "total_claims_wk_pphase_acute ~ spores_cum_wk_pphase_acute + s(week_of_year_pphase_acute, bs='cc', k=", k_week_pphase_acute, ")",
      if (length(wk_weather_cols_pphase_acute) > 0) paste0(" + ", paste(wk_weather_cols_pphase_acute, collapse = " + ")) else ""
    )
  )
}

fit_wk_nb_pphase_acute <- gam(
  f_wk_nb_pphase_acute,
  data   = weekly_baseline_pphase_acute,
  family = nb(),
  method = "REML",
  knots  = knots_week_pphase_acute
)

weekly_baseline_logit_pphase_acute <- weekly_baseline_pphase_acute %>% filter(total_claims_wk_pphase_acute > 0)

if (length(available_years_pphase_acute) >= 2) {
  f_wk_logit_pphase_acute <- as.formula(
    paste0(
      "cbind(acute_claims_wk_pphase_acute, nonacute_claims_wk_pphase_acute) ~ spores_cum_wk_pphase_acute * year_f_pphase_acute + s(week_of_year_pphase_acute, bs='cc', k=", k_week_pphase_acute, ")",
      if (length(wk_weather_cols_pphase_acute) > 0) paste0(" + ", paste(wk_weather_cols_pphase_acute, collapse = " + ")) else ""
    )
  )
} else {
  f_wk_logit_pphase_acute <- as.formula(
    paste0(
      "cbind(acute_claims_wk_pphase_acute, nonacute_claims_wk_pphase_acute) ~ spores_cum_wk_pphase_acute + s(week_of_year_pphase_acute, bs='cc', k=", k_week_pphase_acute, ")",
      if (length(wk_weather_cols_pphase_acute) > 0) paste0(" + ", paste(wk_weather_cols_pphase_acute, collapse = " + ")) else ""
    )
  )
}

fit_wk_logit_pphase_acute <- gam(
  f_wk_logit_pphase_acute,
  data   = weekly_baseline_logit_pphase_acute,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = knots_week_pphase_acute
)

# Only extract year-specific effects if we have multiple years
if (length(available_years_pphase_acute) >= 2) {
  baseline_nb_effects_by_year_pphase_acute <- get_group_effects_pphase_acute(
    fit_pphase_acute = fit_wk_nb_pphase_acute,
    base_term_pphase_acute = "spores_cum_wk_pphase_acute",
    group_var_pphase_acute = "year_f_pphase_acute",
    label_pphase_acute = "RR_pphase_acute"
  )

  baseline_logit_effects_by_year_pphase_acute <- get_group_effects_pphase_acute(
    fit_pphase_acute = fit_wk_logit_pphase_acute,
    base_term_pphase_acute = "spores_cum_wk_pphase_acute",
    group_var_pphase_acute = "year_f_pphase_acute",
    label_pphase_acute = "OR_pphase_acute"
  )
} else {
  # Extract pooled effect only
  beta_nb <- coef(fit_wk_nb_pphase_acute)
  V_nb <- vcov(fit_wk_nb_pphase_acute)
  baseline_nb_effects_by_year_pphase_acute <- lincomb_exp_pphase_acute(
    beta_nb, V_nb,
    setNames(1, "spores_cum_wk_pphase_acute"),
    label_pphase_acute = "RR_pphase_acute"
  ) %>% mutate(level_pphase_acute = as.character(available_years_pphase_acute[1]))
  
  beta_logit <- coef(fit_wk_logit_pphase_acute)
  V_logit <- vcov(fit_wk_logit_pphase_acute)
  baseline_logit_effects_by_year_pphase_acute <- lincomb_exp_pphase_acute(
    beta_logit, V_logit,
    setNames(1, "spores_cum_wk_pphase_acute"),
    label_pphase_acute = "OR_pphase_acute"
  ) %>% mutate(level_pphase_acute = as.character(available_years_pphase_acute[1]))
}

# ===========================================================
# B) PANDEMIC PHASE: daily NB + daily logit with spores × pandemic_phase -----
# ===========================================================

daily_phase_pphase_acute <- daily_pphase_acute %>%
  filter(!is.na(spores_10k_ma7_pphase_acute)) %>%
  mutate(
    pandemic_phase_pphase_acute = relevel(pandemic_phase_pphase_acute, ref = "pre2020")
  )

knots_doy_pphase_acute <- list(doy_pphase_acute = c(0.5, 365.5))
day_weather_cols_pphase_acute <- intersect(weather_covars_use_pphase_acute, names(daily_phase_pphase_acute))

# Adjust k for day-of-year spline based on available data
n_days_pphase_acute <- nrow(daily_phase_pphase_acute)
k_doy_pphase_acute <- min(30, max(5, floor(n_days_pphase_acute / 50)))

# Daily NB model formula
f_day_nb_pphase_acute <- as.formula(
  paste0(
    "total_claims_pphase_acute ~ spores_10k_ma7_pphase_acute * pandemic_phase_pphase_acute + ",
    "s(doy_pphase_acute, bs='cc', k=", k_doy_pphase_acute, ") + dow_pphase_acute + s(date_num_pphase_acute, bs='cr', k=10)",
    if (length(day_weather_cols_pphase_acute) > 0) paste0(" + ", paste(day_weather_cols_pphase_acute, collapse = " + ")) else ""
  )
)

fit_day_nb_pphase_acute <- gam(
  f_day_nb_pphase_acute,
  data   = daily_phase_pphase_acute,
  family = nb(),
  method = "REML",
  knots  = knots_doy_pphase_acute
)

# Daily logit model (only days with claims > 0)
daily_phase_logit_cc_pphase_acute <- daily_phase_pphase_acute %>%
  filter(total_claims_pphase_acute > 0)

f_day_logit_pphase_acute <- as.formula(
  paste0(
    "cbind(acute_claims_pphase_acute, nonacute_claims_pphase_acute) ~ spores_10k_ma7_pphase_acute * pandemic_phase_pphase_acute + ",
    "s(doy_pphase_acute, bs='cc', k=", k_doy_pphase_acute, ") + dow_pphase_acute + s(date_num_pphase_acute, bs='cr', k=10)",
    if (length(day_weather_cols_pphase_acute) > 0) paste0(" + ", paste(day_weather_cols_pphase_acute, collapse = " + ")) else ""
  )
)

fit_day_logit_pphase_acute <- gam(
  f_day_logit_pphase_acute,
  data   = daily_phase_logit_cc_pphase_acute,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = knots_doy_pphase_acute
)

phase_nb_effects_pphase_acute <- get_group_effects_pphase_acute(
  fit_pphase_acute = fit_day_nb_pphase_acute,
  base_term_pphase_acute = "spores_10k_ma7_pphase_acute",
  group_var_pphase_acute = "pandemic_phase_pphase_acute",
  label_pphase_acute = "RR_pphase_acute"
)

phase_logit_effects_pphase_acute <- get_group_effects_pphase_acute(
  fit_pphase_acute = fit_day_logit_pphase_acute,
  base_term_pphase_acute = "spores_10k_ma7_pphase_acute",
  group_var_pphase_acute = "pandemic_phase_pphase_acute",
  label_pphase_acute = "OR_pphase_acute"
)

write.csv(phase_nb_effects_pphase_acute,
          file.path(out_dir_pphase_acute, "pooled_daily_NB_spores_ma7_by_pandemic_phase_pphase_acute.csv"),
          row.names = FALSE)

write.csv(phase_logit_effects_pphase_acute,
          file.path(out_dir_pphase_acute, "pooled_daily_LOGIT_spores_ma7_by_pandemic_phase_pphase_acute.csv"),
          row.names = FALSE)

cat("\nDone. Outputs written to:", out_dir_pphase_acute, "\n")
