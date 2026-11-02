############################################################
## Spores ↔ Asthma models: Baseline vs Pandemic Phase (pphase)
## San Juan health region (claims) + San Juan station spores + SJU weather
##
## A) Baseline (2018–2019) modeled jointly, with year-specific spores effect:
##    - Weekly NB (total claims) + weekly binomial logit (acute fraction)
##    - Exposure: weekly cumulative spores (sum of daily spores_10k)
##    - Controls: cyclic week-of-year spline + year indicators + weekly mean weather
##    - Key term: spores_cum_wk_pphase × year_f_pphase
##
## B) Pandemic/post-pandemic pooled daily models:
##    - pandemic_phase_pphase: pre2020, y2020, y2021plus
##    - Exposure: spores_10k_ma7_pphase × pandemic_phase_pphase
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

claims_path_pphase  <- "data/2022_asma_From_date_ver 3.xlsx"
claims_sheet_pphase <- 1

spore_files_sj_pphase <- "data/Datos esporas y polen 24-08-19 a 31-12-22 Corregidas 12-09-25.xlsx"

# NOAA token (recommended: environment variable)
# Sys.setenv(NOAATOKEN="YOUR_TOKEN"); options(noaakey = Sys.getenv("NOAATOKEN"))
options(noaakey = Sys.getenv("NOAATOKEN"))

sju_station_id_pphase <- "RQW00011641"

years_all_pphase      <- 2018:2022
years_baseline_pphase <- 2018:2019

# Health region mapping (EDIT to match your definition)
san_juan_munis_pphase <- c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA")
caguas_munis_pphase   <- c("CAGUAS", "GURABO", "JUNCOS", "SAN LORENZO")

clinical_settings_pphase <- c("ED", "Inpatient", "Outpatient", "Urgent")

max_missing_days_week_pphase <- 3

spore_scale_divisor_pphase <- 10000  # spores_10k = spores_total / 10000

weather_covars_pphase <- c("tmean_c_pphase", "rh_mean_pphase", "wind_ms_pphase", "prcp_mm_pphase", "dewpoint_c_pphase")

out_dir_pphase <- "statistical_results"
dir.create(out_dir_pphase, showWarnings = FALSE, recursive = TRUE)

# ===========================================================
# 2) HELPERS  -----
# ===========================================================

parse_date_any_pphase <- function(x_pphase) {
  if (inherits(x_pphase, "Date")) return(x_pphase)
  if (is.numeric(x_pphase)) return(as.Date(x_pphase, origin = "1899-12-30")) # Excel numeric dates
  x_pphase <- as.character(x_pphase)
  d_pphase <- suppressWarnings(ymd(x_pphase))
  d_pphase <- dplyr::coalesce(d_pphase, suppressWarnings(mdy(x_pphase)))
  d_pphase <- dplyr::coalesce(d_pphase, suppressWarnings(dmy(x_pphase)))
  dplyr::coalesce(d_pphase, suppressWarnings(as.Date(x_pphase)))
}

safe_mean_pphase <- function(x_pphase) {
  m_pphase <- mean(x_pphase, na.rm = TRUE)
  if (is.nan(m_pphase)) NA_real_ else m_pphase
}

# ---- claims: care setting and asthma primary ----
care_setting_fun_pphase <- function(hcg_setting_pphase, hcg_mr_line_desc_pphase) {
  s_pphase <- tolower(as.character(hcg_setting_pphase))
  l_pphase <- tolower(as.character(hcg_mr_line_desc_pphase))
  if (grepl("ed|emergency", s_pphase) || grepl("\\bed\\b|\\ber\\b", l_pphase)) return("ED")
  if (grepl("inpatient|hospital", s_pphase) || grepl("inpatient|hospital", l_pphase)) return("Inpatient")
  if (grepl("urgent", s_pphase) || grepl("urgent", l_pphase)) return("Urgent")
  if (grepl("office|clinic|outpatient", s_pphase) || grepl("office|clinic", l_pphase)) return("Outpatient")
  if (grepl("dme|supplies", l_pphase)) return("DME")
  "Other"
}

ensure_care_setting_pphase <- function(df_pphase) {
  if ("care_setting_pphase" %in% names(df_pphase)) return(df_pphase)

  if ("care_setting" %in% names(df_pphase)) {
    return(df_pphase %>% mutate(care_setting_pphase = care_setting))
  }

  needed_pphase <- c("hcg_setting", "hcg_mr_line_desc")
  if (!all(needed_pphase %in% names(df_pphase))) stop("No care_setting and missing hcg_setting/hcg_mr_line_desc to create it.")

  df_pphase %>% mutate(care_setting_pphase = mapply(care_setting_fun_pphase, hcg_setting, hcg_mr_line_desc))
}

ensure_asthma_primary_pphase <- function(df_pphase) {
  if ("asthma_primary_pphase" %in% names(df_pphase)) return(df_pphase)

  if ("asthma_primary" %in% names(df_pphase)) {
    return(df_pphase %>% mutate(asthma_primary_pphase = asthma_primary))
  }

  if (!"icd_diag_01_primary" %in% names(df_pphase)) stop("Missing asthma_primary and icd_diag_01_primary.")

  df_pphase %>%
    mutate(
      icd_diag_01_primary = stringr::str_to_upper(stringr::str_trim(as.character(icd_diag_01_primary))),
      asthma_primary_pphase = if_else(stringr::str_detect(icd_diag_01_primary, "^J45|^J46"), "asthma", "other")
    )
}

# ---- spores/pollen ----
read_spores_sj_onefile_pphase <- function(path_pphase) {
  x_pphase <- read_excel(path_pphase) %>% clean_names()
  date_col_pphase <- if ("date" %in% names(x_pphase)) "date" else names(x_pphase)[1]
  spore_col_pphase  <- grep("esporas|spore", names(x_pphase), ignore.case = TRUE, value = TRUE)[1]
  pollen_col_pphase <- grep("polen|pollen", names(x_pphase), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(spore_col_pphase)) stop("Could not identify a spores column in: ", path_pphase)

  x_pphase %>%
    transmute(
      date = parse_date_any_pphase(.data[[date_col_pphase]]),
      spores_total_raw_pphase = suppressWarnings(as.numeric(.data[[spore_col_pphase]])),
      pollen_total_raw_pphase = if (!is.na(pollen_col_pphase)) suppressWarnings(as.numeric(.data[[pollen_col_pphase]])) else NA_real_
    ) %>%
    filter(!is.na(date))
}

weekly_qc_impute_or_exclude_pphase <- function(df_pphase, value_col_pphase, base_name_pphase, max_missing_pphase = 3) {
  value_sym_pphase <- rlang::ensym(value_col_pphase)

  out_col_pphase        <- paste0(base_name_pphase, "_pphase")
  imputed_col_pphase    <- paste0(base_name_pphase, "_imputed_pphase")
  week_excl_col_pphase  <- paste0(base_name_pphase, "_week_excluded_pphase")
  nmiss_col_pphase      <- paste0(base_name_pphase, "_n_missing_week_pphase")
  weekmean_col_pphase   <- paste0(base_name_pphase, "_week_mean_pphase")

  df_pphase %>%
    mutate(week_start_pphase = floor_date(date, unit = "week", week_start = 1)) %>%
    group_by(week_start_pphase) %>%
    mutate(
      n_missing_week_tmp_pphase = sum(is.na(!!value_sym_pphase)),
      week_mean_tmp_pphase      = safe_mean_pphase(!!value_sym_pphase),
      week_ok_tmp_pphase        = (n_missing_week_tmp_pphase <= max_missing_pphase) & !is.na(week_mean_tmp_pphase),
      week_excluded_tmp_pphase  = !week_ok_tmp_pphase,
      imputed_day_tmp_pphase    = week_ok_tmp_pphase & is.na(!!value_sym_pphase),
      value_clean_tmp_pphase    = case_when(
        week_ok_tmp_pphase & is.na(!!value_sym_pphase) ~ week_mean_tmp_pphase,
        week_ok_tmp_pphase ~ as.numeric(!!value_sym_pphase),
        TRUE ~ NA_real_
      ),
      !!out_col_pphase       := value_clean_tmp_pphase,
      !!imputed_col_pphase   := imputed_day_tmp_pphase,
      !!week_excl_col_pphase := week_excluded_tmp_pphase,
      !!nmiss_col_pphase     := n_missing_week_tmp_pphase,
      !!weekmean_col_pphase  := week_mean_tmp_pphase
    ) %>%
    ungroup() %>%
    dplyr::select(-n_missing_week_tmp_pphase, -week_mean_tmp_pphase, -week_ok_tmp_pphase,
           -week_excluded_tmp_pphase, -imputed_day_tmp_pphase, -value_clean_tmp_pphase)
}

# ---- weather (unit handling) ----
needs_div10_pphase <- function(x_pphase, thr_pphase) {
  x_pphase <- suppressWarnings(as.numeric(x_pphase))
  q_pphase <- suppressWarnings(quantile(abs(x_pphase), probs = 0.99, na.rm = TRUE, names = FALSE))
  if (is.na(q_pphase)) return(FALSE)
  q_pphase > thr_pphase
}

pull_weather_sju_pphase <- function(year_value_pphase, station_id_pphase) {
  vars_pphase <- c("PRCP", "TMAX", "TMIN", "AWND", "ADPT", "RHAV")

  w_pphase <- meteo_pull_monitors(
    monitors   = station_id_pphase,
    var        = vars_pphase,
    date_min   = sprintf("%d-01-01", year_value_pphase),
    date_max   = sprintf("%d-12-31", year_value_pphase),
    keep_flags = FALSE
  ) %>%
    as_tibble() %>%
    mutate(date = as.Date(date))

  for (nm_pphase in c("prcp","tmax","tmin","awnd","adpt","rhav")) {
    if (!nm_pphase %in% names(w_pphase)) w_pphase[[nm_pphase]] <- NA_real_
  }

  prcp_mm_pphase     <- if (needs_div10_pphase(w_pphase$prcp, thr_pphase = 300)) w_pphase$prcp / 10 else w_pphase$prcp
  tmax_c_pphase      <- if (needs_div10_pphase(w_pphase$tmax, thr_pphase = 80))  w_pphase$tmax / 10 else w_pphase$tmax
  tmin_c_pphase      <- if (needs_div10_pphase(w_pphase$tmin, thr_pphase = 80))  w_pphase$tmin / 10 else w_pphase$tmin
  dewpoint_c_pphase  <- if (needs_div10_pphase(w_pphase$adpt, thr_pphase = 80))  w_pphase$adpt / 10 else w_pphase$adpt
  wind_ms_pphase     <- if (needs_div10_pphase(w_pphase$awnd, thr_pphase = 30))  w_pphase$awnd / 10 else w_pphase$awnd

  tibble(
    date = w_pphase$date,
    prcp_mm_pphase = prcp_mm_pphase,
    tmax_c_pphase  = tmax_c_pphase,
    tmin_c_pphase  = tmin_c_pphase,
    tmean_c_pphase = (tmax_c_pphase + tmin_c_pphase) / 2,
    dewpoint_c_pphase = dewpoint_c_pphase,
    rh_mean_pphase = w_pphase$rhav,
    wind_ms_pphase = wind_ms_pphase
  )
}

# ---- linear-combination CI for exp(beta) ----
lincomb_exp_pphase <- function(beta_pphase, V_pphase, idx_pphase, label_pphase = "effect_pphase") {
  b_pphase <- sum(beta_pphase[names(idx_pphase)] * idx_pphase, na.rm = TRUE)
  v_pphase <- as.numeric(t(idx_pphase) %*% V_pphase[names(idx_pphase), names(idx_pphase), drop = FALSE] %*% idx_pphase)
  se_pphase <- sqrt(v_pphase)
  ci_pphase <- b_pphase + c(-1, 1) * 1.96 * se_pphase
  tibble(
    log_est_pphase = b_pphase,
    log_se_pphase  = se_pphase,
    !!label_pphase := exp(b_pphase),
    low_pphase     = exp(ci_pphase[1]),
    high_pphase    = exp(ci_pphase[2])
  )
}

get_group_effects_pphase <- function(fit_pphase, base_term_pphase, group_var_pphase, label_pphase) {
  beta_pphase <- coef(fit_pphase)
  V_pphase <- vcov(fit_pphase)

  lvls_pphase <- levels(model.frame(fit_pphase)[[group_var_pphase]])
  ref_pphase  <- lvls_pphase[1]

  map_dfr(lvls_pphase, function(lv_pphase) {
    idx_pphase <- rep(0, length(beta_pphase)); names(idx_pphase) <- names(beta_pphase)
    idx_pphase[base_term_pphase] <- 1

    if (lv_pphase != ref_pphase) {
      t1_pphase <- paste0(base_term_pphase, ":", group_var_pphase, lv_pphase)
      t2_pphase <- paste0(group_var_pphase, lv_pphase, ":", base_term_pphase)
      if (t1_pphase %in% names(beta_pphase)) idx_pphase[t1_pphase] <- 1
      if (t2_pphase %in% names(beta_pphase)) idx_pphase[t2_pphase] <- 1
    }

    est_pphase <- lincomb_exp_pphase(beta_pphase, V_pphase,
                                     idx_pphase[names(idx_pphase) %in% names(beta_pphase)],
                                     label_pphase = label_pphase)
    tibble(level_pphase = lv_pphase) %>% bind_cols(est_pphase)
  })
}

# ===========================================================
# 3) BUILD DAILY ASTHMA OUTCOMES (San Juan health region)-----
# ===========================================================

claims_raw_pphase <- read_excel(claims_path_pphase, sheet = claims_sheet_pphase) %>%
  clean_names() %>%
  mutate(from_date = parse_date_any_pphase(from_date)) %>%
  ensure_care_setting_pphase() %>%
  ensure_asthma_primary_pphase()

claims_pphase <- claims_raw_pphase %>%
  mutate(
    municipality_desc_pphase = toupper(as.character(municipality_desc)),
    year_pphase = year(from_date),
    health_region_pphase = case_when(
      municipality_desc_pphase %in% san_juan_munis_pphase ~ "San Juan",
      municipality_desc_pphase %in% caguas_munis_pphase   ~ "Caguas",
      TRUE                                                ~ NA_character_
    ),
    health_region_pphase = factor(health_region_pphase, levels = c("San Juan", "Caguas"))
  ) %>%
  filter(year_pphase %in% years_all_pphase) %>%
  filter(!is.na(health_region_pphase))

claims_sj_pphase <- claims_pphase %>%
  filter(
    health_region_pphase == "San Juan",
    asthma_primary_pphase == "asthma",
    care_setting_pphase %in% clinical_settings_pphase
  ) %>%
  mutate(acute_flag_pphase = as.integer(care_setting_pphase %in% c("ED", "Inpatient")))

calendar_pphase <- tibble(
  date = seq.Date(as.Date(sprintf("%d-01-01", min(years_all_pphase))),
                  as.Date(sprintf("%d-12-31", max(years_all_pphase))),
                  by = "day")
) %>%
  mutate(
    year_pphase = year(date),
    doy_pphase  = yday(date),
    dow_pphase  = wday(date, label = TRUE, week_start = 1),
    week_start_pphase = floor_date(date, unit = "week", week_start = 1),
    week_of_year_pphase = isoweek(date),
    date_num_pphase = as.numeric(date)
  ) %>%
  filter(year_pphase %in% years_all_pphase)

daily_asthma_pphase <- claims_sj_pphase %>%
  group_by(date = from_date) %>%
  summarise(
    total_claims_pphase = n(),
    acute_claims_pphase = sum(acute_flag_pphase == 1, na.rm = TRUE),
    nonacute_claims_pphase = total_claims_pphase - acute_claims_pphase,
    .groups = "drop"
  )

daily_outcome_pphase <- calendar_pphase %>%
  left_join(daily_asthma_pphase, by = "date") %>%
  mutate(
    across(c(total_claims_pphase, acute_claims_pphase, nonacute_claims_pphase),
           ~ replace_na(.x, 0L))
  )

# ===========================================================
# 4) BUILD DAILY SPORES/POLLEN with weekly QC (San Juan station)-----
# ===========================================================

spores_raw_pphase <- map_dfr(spore_files_sj_pphase, read_spores_sj_onefile_pphase) %>%
  group_by(date) %>%
  summarise(
    spores_total_raw_pphase  = safe_mean_pphase(spores_total_raw_pphase),
    pollen_total_raw_pphase  = safe_mean_pphase(pollen_total_raw_pphase),
    .groups = "drop"
  ) %>%
  arrange(date)

spores_daily_pphase <- tibble(date = calendar_pphase$date) %>%
  left_join(spores_raw_pphase, by = "date") %>%
  arrange(date)

spores_daily_pphase <- weekly_qc_impute_or_exclude_pphase(
  spores_daily_pphase,
  value_col_pphase = spores_total_raw_pphase,
  base_name_pphase = "spores_total",
  max_missing_pphase = max_missing_days_week_pphase
)

spores_daily_pphase <- weekly_qc_impute_or_exclude_pphase(
  spores_daily_pphase,
  value_col_pphase = pollen_total_raw_pphase,
  base_name_pphase = "pollen_total",
  max_missing_pphase = max_missing_days_week_pphase
)

spores_daily_pphase <- spores_daily_pphase %>%
  mutate(
    spores_10k_pphase = spores_total_pphase / spore_scale_divisor_pphase
  )

# ===========================================================
# 5) WEATHER + MERGE into DAILY dataset + MA7 + pandemic phase -----
# ===========================================================

weather_daily_pphase <- map_dfr(years_all_pphase, ~ pull_weather_sju_pphase(.x, sju_station_id_pphase))

daily_pphase <- daily_outcome_pphase %>%
  left_join(spores_daily_pphase,  by = "date") %>%
  left_join(weather_daily_pphase, by = "date") %>%
  mutate(
    pandemic_phase_pphase = case_when(
      year_pphase < 2020  ~ "pre2020",
      year_pphase == 2020 ~ "y2020",
      year_pphase >= 2021 ~ "y2021plus"
    ),
    pandemic_phase_pphase = factor(pandemic_phase_pphase, levels = c("pre2020", "y2020", "y2021plus"))
  ) %>%
  arrange(date) %>%
  mutate(
    spores_10k_ma7_pphase = zoo::rollapply(spores_10k_pphase, width = 7, FUN = mean,
                                          align = "right", fill = NA_real_, partial = FALSE)
  )

weather_covars_use_pphase <- intersect(weather_covars_pphase, names(daily_pphase))
if (length(weather_covars_use_pphase) == 0) stop("No weather covariates found after weather pull.")

# ===========================================================
# A) BASELINE (2018–2019): weekly NB + weekly logit with spores × year -----
# ===========================================================

weekly_baseline_pphase <- daily_pphase %>%
  filter(year_pphase %in% years_baseline_pphase) %>%
  filter(week_of_year_pphase <= 52) %>%
  mutate(week_start_pphase = floor_date(date, unit = "week", week_start = 1)) %>%
  group_by(year_pphase, week_start_pphase, week_of_year_pphase) %>%
  summarise(
    total_claims_wk_pphase    = sum(total_claims_pphase),
    acute_claims_wk_pphase    = sum(acute_claims_pphase),
    nonacute_claims_wk_pphase = sum(nonacute_claims_pphase),
    spores_cum_wk_pphase = sum(spores_10k_pphase, na.rm = FALSE),
    across(all_of(weather_covars_use_pphase), ~ safe_mean_pphase(.x), .names = "{.col}_wkmean_pphase"),
    spores_week_excluded_pphase = any(spores_total_week_excluded_pphase, na.rm = TRUE),
    pollen_week_excluded_pphase = any(pollen_total_week_excluded_pphase, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!spores_week_excluded_pphase) %>%
  filter(!is.na(spores_cum_wk_pphase)) %>%
  mutate(
    year_f_pphase = factor(year_pphase)
  )

# Set reference level to the first available year
if (nrow(weekly_baseline_pphase) > 0) {
  available_years_pphase <- sort(unique(weekly_baseline_pphase$year_pphase))
  
  # Check if we have at least 2 years for year-specific models
  if (length(available_years_pphase) >= 2) {
    weekly_baseline_pphase <- weekly_baseline_pphase %>%
      mutate(year_f_pphase = relevel(factor(year_pphase), ref = as.character(available_years_pphase[1])))
  } else {
    cat("\nWarning: Only", length(available_years_pphase), "year(s) available after filtering. Fitting pooled baseline model without year interaction.\n")
    weekly_baseline_pphase <- weekly_baseline_pphase %>%
      mutate(year_f_pphase = factor(year_pphase))
  }
}

knots_week_pphase <- list(week_of_year_pphase = c(0.5, 52.5))
wk_weather_cols_pphase <- names(weekly_baseline_pphase)[stringr::str_ends(names(weekly_baseline_pphase), "_wkmean_pphase")]

# Dynamically adjust k based on available data
n_weeks_pphase <- nrow(weekly_baseline_pphase)
k_week_pphase <- min(20, max(3, floor(n_weeks_pphase * 0.8)))

# Build formula based on number of years available
if (length(available_years_pphase) >= 2) {
  f_wk_nb_pphase <- as.formula(
    paste0(
      "total_claims_wk_pphase ~ spores_cum_wk_pphase * year_f_pphase + s(week_of_year_pphase, bs='cc', k=", k_week_pphase, ")",
      if (length(wk_weather_cols_pphase) > 0) paste0(" + ", paste(wk_weather_cols_pphase, collapse = " + ")) else ""
    )
  )
} else {
  f_wk_nb_pphase <- as.formula(
    paste0(
      "total_claims_wk_pphase ~ spores_cum_wk_pphase + s(week_of_year_pphase, bs='cc', k=", k_week_pphase, ")",
      if (length(wk_weather_cols_pphase) > 0) paste0(" + ", paste(wk_weather_cols_pphase, collapse = " + ")) else ""
    )
  )
}

fit_wk_nb_pphase <- gam(
  f_wk_nb_pphase,
  data   = weekly_baseline_pphase,
  family = nb(),
  method = "REML",
  knots  = knots_week_pphase
)

weekly_baseline_logit_pphase <- weekly_baseline_pphase %>% filter(total_claims_wk_pphase > 0)

if (length(available_years_pphase) >= 2) {
  f_wk_logit_pphase <- as.formula(
    paste0(
      "cbind(acute_claims_wk_pphase, nonacute_claims_wk_pphase) ~ spores_cum_wk_pphase * year_f_pphase + s(week_of_year_pphase, bs='cc', k=", k_week_pphase, ")",
      if (length(wk_weather_cols_pphase) > 0) paste0(" + ", paste(wk_weather_cols_pphase, collapse = " + ")) else ""
    )
  )
} else {
  f_wk_logit_pphase <- as.formula(
    paste0(
      "cbind(acute_claims_wk_pphase, nonacute_claims_wk_pphase) ~ spores_cum_wk_pphase + s(week_of_year_pphase, bs='cc', k=", k_week_pphase, ")",
      if (length(wk_weather_cols_pphase) > 0) paste0(" + ", paste(wk_weather_cols_pphase, collapse = " + ")) else ""
    )
  )
}

fit_wk_logit_pphase <- gam(
  f_wk_logit_pphase,
  data   = weekly_baseline_logit_pphase,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = knots_week_pphase
)

# Only extract year-specific effects if we have multiple years
if (length(available_years_pphase) >= 2) {
  baseline_nb_effects_by_year_pphase <- get_group_effects_pphase(
    fit_pphase = fit_wk_nb_pphase,
    base_term_pphase = "spores_cum_wk_pphase",
    group_var_pphase = "year_f_pphase",
    label_pphase = "RR_pphase"
  )

  baseline_logit_effects_by_year_pphase <- get_group_effects_pphase(
    fit_pphase = fit_wk_logit_pphase,
    base_term_pphase = "spores_cum_wk_pphase",
    group_var_pphase = "year_f_pphase",
    label_pphase = "OR_pphase"
  )
} else {
  # Extract pooled effect only
  beta_nb <- coef(fit_wk_nb_pphase)
  V_nb <- vcov(fit_wk_nb_pphase)
  baseline_nb_effects_by_year_pphase <- lincomb_exp_pphase(
    beta_nb, V_nb,
    setNames(1, "spores_cum_wk_pphase"),
    label_pphase = "RR_pphase"
  ) %>% mutate(level_pphase = as.character(available_years_pphase[1]))
  
  beta_logit <- coef(fit_wk_logit_pphase)
  V_logit <- vcov(fit_wk_logit_pphase)
  baseline_logit_effects_by_year_pphase <- lincomb_exp_pphase(
    beta_logit, V_logit,
    setNames(1, "spores_cum_wk_pphase"),
    label_pphase = "OR_pphase"
  ) %>% mutate(level_pphase = as.character(available_years_pphase[1]))
}

# ===========================================================
# B) PANDEMIC PHASE: daily NB + daily logit with spores × pandemic_phase -----
# ===========================================================

daily_phase_pphase <- daily_pphase %>%
  filter(!is.na(spores_10k_ma7_pphase)) %>%
  mutate(
    pandemic_phase_pphase = relevel(pandemic_phase_pphase, ref = "pre2020")
  )

knots_doy_pphase <- list(doy_pphase = c(0.5, 365.5))
day_weather_cols_pphase <- intersect(weather_covars_use_pphase, names(daily_phase_pphase))

# Adjust k for day-of-year spline based on available data
n_days_pphase <- nrow(daily_phase_pphase)
k_doy_pphase <- min(30, max(5, floor(n_days_pphase / 50)))

# Daily NB model formula
f_day_nb_pphase <- as.formula(
  paste0(
    "total_claims_pphase ~ spores_10k_ma7_pphase * pandemic_phase_pphase + ",
    "s(doy_pphase, bs='cc', k=", k_doy_pphase, ") + dow_pphase + s(date_num_pphase, bs='cr', k=10)",
    if (length(day_weather_cols_pphase) > 0) paste0(" + ", paste(day_weather_cols_pphase, collapse = " + ")) else ""
  )
)

fit_day_nb_pphase <- gam(
  f_day_nb_pphase,
  data   = daily_phase_pphase,
  family = nb(),
  method = "REML",
  knots  = knots_doy_pphase
)

# Daily logit model (only days with claims > 0)
daily_phase_logit_cc_pphase <- daily_phase_pphase %>%
  filter(total_claims_pphase > 0)

f_day_logit_pphase <- as.formula(
  paste0(
    "cbind(acute_claims_pphase, nonacute_claims_pphase) ~ spores_10k_ma7_pphase * pandemic_phase_pphase + ",
    "s(doy_pphase, bs='cc', k=", k_doy_pphase, ") + dow_pphase + s(date_num_pphase, bs='cr', k=10)",
    if (length(day_weather_cols_pphase) > 0) paste0(" + ", paste(day_weather_cols_pphase, collapse = " + ")) else ""
  )
)

fit_day_logit_pphase <- gam(
  f_day_logit_pphase,
  data   = daily_phase_logit_cc_pphase,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = knots_doy_pphase
)

phase_nb_effects_pphase <- get_group_effects_pphase(
  fit_pphase = fit_day_nb_pphase,
  base_term_pphase = "spores_10k_ma7_pphase",
  group_var_pphase = "pandemic_phase_pphase",
  label_pphase = "RR_pphase"
)

phase_logit_effects_pphase <- get_group_effects_pphase(
  fit_pphase = fit_day_logit_pphase,
  base_term_pphase = "spores_10k_ma7_pphase",
  group_var_pphase = "pandemic_phase_pphase",
  label_pphase = "OR_pphase"
)

write.csv(phase_nb_effects_pphase,
          file.path(out_dir_pphase, "pooled_daily_NB_spores_ma7_by_pandemic_phase_pphase.csv"),
          row.names = FALSE)

write.csv(phase_logit_effects_pphase,
          file.path(out_dir_pphase, "pooled_daily_LOGIT_spores_ma7_by_pandemic_phase_pphase.csv"),
          row.names = FALSE)

cat("\nDone. Outputs written to:", out_dir_pphase, "\n")
