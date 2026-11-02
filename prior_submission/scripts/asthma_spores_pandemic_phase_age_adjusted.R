############################################################
## Spores ↔ Asthma models: Baseline vs Pandemic Phase BY AGE GROUP
## San Juan health region (claims) + San Juan station spores + SJU weather
##
## A) Baseline (2018–2019) modeled jointly, with year-specific spores effect BY AGE:
##    - Weekly NB (total claims) + weekly binomial logit (acute fraction)
##    - Exposure: weekly cumulative spores (sum of daily spores_10k)
##    - Controls: cyclic week-of-year spline + year indicators + weekly mean weather
##    - Key term: spores_cum_wk_pphase_age × year_f_pphase_age × age_group_pphase_age
##
## B) Pandemic/post-pandemic pooled daily models BY AGE:
##    - pandemic_phase_pphase_age: pre2020, y2020, y2021plus
##    - Exposure: spores_10k_ma7_pphase_age × pandemic_phase_pphase_age × age_group_pphase_age
##    - Controls: cyclic day-of-year spline + day-of-week + smooth time trend + weather
##
## Exposure QC (applied separately to spores and pollen BEFORE modeling):
##   - If a week has <= 3 missing days: impute missing days with that week's mean
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
  library(rlang)
})

# ===========================================================
# 1) SETTINGS
# ===========================================================

claims_path_pphase_age  <- "data/2022_asma_From_date_ver 3.xlsx"
claims_sheet_pphase_age <- 1

# Allow a single file or multiple files
spore_files_sj_pphase_age <- c("data/Datos esporas y polen 24-08-19 a 31-12-22 Corregidas 12-09-25.xlsx")

# NOAA token (recommended: environment variable)
# Sys.setenv(NOAATOKEN="YOUR_TOKEN"); options(noaakey = Sys.getenv("NOAATOKEN"))
options(noaakey = Sys.getenv("NOAATOKEN"))

sju_station_id_pphase_age <- "RQW00011641"

years_all_pphase_age      <- 2018:2022
years_baseline_pphase_age <- 2018:2019

# Health region mapping (EDIT to match your definition)
san_juan_munis_pphase_age <- c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA")
caguas_munis_pphase_age   <- c("CAGUAS", "GURABO", "JUNCOS", "SAN LORENZO")

clinical_settings_pphase_age <- c("ED", "Inpatient", "Outpatient", "Urgent")

max_missing_days_week_pphase_age <- 3
spore_scale_divisor_pphase_age <- 10000  # spores_10k = spores_total / 10000

# Weather covariates to adjust for (created by pull_weather_sju_pphase_age)
weather_covars_pphase_age <- c("tmean_c_pphase_age", "rh_mean_pphase_age", "wind_ms_pphase_age", "prcp_mm_pphase_age", "dewpoint_c_pphase_age")

# Age group definitions
age_breaks_pphase_age <- c(0, 4, 11, 17, 39, 64, 74, 84, Inf)
age_labels_pphase_age <- c("0-4", "5-11", "12-17", "18-39", "40-64", "65-74", "75-84", "85+")

# Alternative simpler age grouping (uncomment if preferred):
# age_breaks_pphase_age <- c(0, 18, 65, Inf)
# age_labels_pphase_age <- c("0-17", "18-64", "65+")

out_dir_pphase_age <- "statistical_results"
dir.create(out_dir_pphase_age, showWarnings = FALSE, recursive = TRUE)

# ===========================================================
# 2) HELPERS
# ===========================================================

parse_date_any_pphase_age <- function(x_pphase_age) {
  if (inherits(x_pphase_age, "Date")) return(x_pphase_age)
  if (is.numeric(x_pphase_age)) return(as.Date(x_pphase_age, origin = "1899-12-30")) # Excel numeric dates
  x_pphase_age <- as.character(x_pphase_age)
  d_pphase_age <- suppressWarnings(ymd(x_pphase_age))
  d_pphase_age <- dplyr::coalesce(d_pphase_age, suppressWarnings(mdy(x_pphase_age)))
  d_pphase_age <- dplyr::coalesce(d_pphase_age, suppressWarnings(dmy(x_pphase_age)))
  dplyr::coalesce(d_pphase_age, suppressWarnings(as.Date(x_pphase_age)))
}

safe_mean_pphase_age <- function(x_pphase_age) {
  m_pphase_age <- mean(x_pphase_age, na.rm = TRUE)
  if (is.nan(m_pphase_age)) NA_real_ else m_pphase_age
}

# ---- claims: care setting and asthma primary ----
care_setting_fun_pphase_age <- function(hcg_setting_pphase_age, hcg_mr_line_desc_pphase_age) {
  s_pphase_age <- tolower(as.character(hcg_setting_pphase_age))
  l_pphase_age <- tolower(as.character(hcg_mr_line_desc_pphase_age))
  if (grepl("ed|emergency", s_pphase_age) || grepl("\\bed\\b|\\ber\\b", l_pphase_age)) return("ED")
  if (grepl("inpatient|hospital", s_pphase_age) || grepl("inpatient|hospital", l_pphase_age)) return("Inpatient")
  if (grepl("urgent", s_pphase_age) || grepl("urgent", l_pphase_age)) return("Urgent")
  if (grepl("office|clinic|outpatient", s_pphase_age) || grepl("office|clinic", l_pphase_age)) return("Outpatient")
  if (grepl("dme|supplies", l_pphase_age)) return("DME")
  "Other"
}

ensure_care_setting_pphase_age <- function(df_pphase_age) {
  if ("care_setting_pphase_age" %in% names(df_pphase_age)) return(df_pphase_age)

  if ("care_setting" %in% names(df_pphase_age)) {
    return(df_pphase_age %>% mutate(care_setting_pphase_age = care_setting))
  }

  needed_pphase_age <- c("hcg_setting", "hcg_mr_line_desc")
  if (!all(needed_pphase_age %in% names(df_pphase_age))) stop("No care_setting and missing hcg_setting/hcg_mr_line_desc to create it.")

  df_pphase_age %>% mutate(care_setting_pphase_age = mapply(care_setting_fun_pphase_age, hcg_setting, hcg_mr_line_desc))
}

ensure_asthma_primary_pphase_age <- function(df_pphase_age) {
  if ("asthma_primary_pphase_age" %in% names(df_pphase_age)) return(df_pphase_age)

  if ("asthma_primary" %in% names(df_pphase_age)) {
    return(df_pphase_age %>% mutate(asthma_primary_pphase_age = asthma_primary))
  }

  if (!"icd_diag_01_primary" %in% names(df_pphase_age)) stop("Missing asthma_primary and icd_diag_01_primary.")

  df_pphase_age %>%
    mutate(
      icd_diag_01_primary = stringr::str_to_upper(stringr::str_trim(as.character(icd_diag_01_primary))),
      asthma_primary_pphase_age = if_else(stringr::str_detect(icd_diag_01_primary, "^J45|^J46"), "asthma", "other")
    )
}

# ---- age group helper ----
ensure_age_group_pphase_age <- function(df_pphase_age) {
  if ("age_group_pphase_age" %in% names(df_pphase_age)) return(df_pphase_age)

  if ("age_group" %in% names(df_pphase_age)) {
    return(df_pphase_age %>% mutate(age_group_pphase_age = as.factor(age_group)))
  }

  if (!"age" %in% names(df_pphase_age)) {
    warning("No age column found; cannot create age_group_pphase_age. Setting to NA.")
    return(df_pphase_age %>% mutate(age_group_pphase_age = NA_character_))
  }

  df_pphase_age %>%
    mutate(
      age_raw_pphase_age = stringr::str_trim(as.character(age)),
      age_num_pphase_age = dplyr::case_when(
        stringr::str_detect(age_raw_pphase_age, "\\+$") ~ suppressWarnings(as.numeric(stringr::str_remove(age_raw_pphase_age, "\\+"))),
        TRUE ~ suppressWarnings(as.numeric(age_raw_pphase_age))
      ),
      age_group_pphase_age = cut(
        age_num_pphase_age,
        breaks = age_breaks_pphase_age,
        labels = age_labels_pphase_age,
        right  = TRUE,
        include.lowest = TRUE
      )
    ) %>%
    dplyr::select(-age_raw_pphase_age, -age_num_pphase_age)
}

# ---- spores/pollen ----
read_spores_sj_onefile_pphase_age <- function(path_pphase_age) {
  x_pphase_age <- read_excel(path_pphase_age) %>% clean_names()
  date_col_pphase_age <- if ("date" %in% names(x_pphase_age)) "date" else names(x_pphase_age)[1]
  spore_col_pphase_age  <- grep("esporas|spore", names(x_pphase_age), ignore.case = TRUE, value = TRUE)[1]
  pollen_col_pphase_age <- grep("polen|pollen", names(x_pphase_age), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(spore_col_pphase_age)) stop("Could not identify a spores column in: ", path_pphase_age)

  x_pphase_age %>%
    transmute(
      date = parse_date_any_pphase_age(.data[[date_col_pphase_age]]),
      spores_total_raw_pphase_age = suppressWarnings(as.numeric(.data[[spore_col_pphase_age]])),
      pollen_total_raw_pphase_age = if (!is.na(pollen_col_pphase_age))
        suppressWarnings(as.numeric(.data[[pollen_col_pphase_age]])) else NA_real_
    ) %>%
    filter(!is.na(date))
}

weekly_qc_impute_or_exclude_pphase_age <- function(df_pphase_age, value_col_pphase_age, base_name_pphase_age, max_missing_pphase_age = 3) {
  value_sym_pphase_age <- rlang::ensym(value_col_pphase_age)

  out_col_pphase_age       <- paste0(base_name_pphase_age, "_pphase_age")
  imputed_col_pphase_age   <- paste0(base_name_pphase_age, "_imputed_pphase_age")
  week_excl_col_pphase_age <- paste0(base_name_pphase_age, "_week_excluded_pphase_age")
  nmiss_col_pphase_age     <- paste0(base_name_pphase_age, "_n_missing_week_pphase_age")
  weekmean_col_pphase_age  <- paste0(base_name_pphase_age, "_week_mean_pphase_age")

  df_pphase_age %>%
    mutate(week_start_pphase_age = floor_date(date, unit = "week", week_start = 1)) %>%
    group_by(week_start_pphase_age) %>%
    mutate(
      n_missing_week_tmp_pphase_age = sum(is.na(!!value_sym_pphase_age)),
      week_mean_tmp_pphase_age      = safe_mean_pphase_age(!!value_sym_pphase_age),
      week_ok_tmp_pphase_age        = (n_missing_week_tmp_pphase_age <= max_missing_pphase_age) & !is.na(week_mean_tmp_pphase_age),
      week_excluded_tmp_pphase_age  = !week_ok_tmp_pphase_age,
      imputed_day_tmp_pphase_age    = week_ok_tmp_pphase_age & is.na(!!value_sym_pphase_age),
      value_clean_tmp_pphase_age    = case_when(
        week_ok_tmp_pphase_age & is.na(!!value_sym_pphase_age) ~ week_mean_tmp_pphase_age,
        week_ok_tmp_pphase_age ~ as.numeric(!!value_sym_pphase_age),
        TRUE ~ NA_real_
      ),
      !!rlang::sym(out_col_pphase_age)       := value_clean_tmp_pphase_age,
      !!rlang::sym(imputed_col_pphase_age)   := imputed_day_tmp_pphase_age,
      !!rlang::sym(week_excl_col_pphase_age) := week_excluded_tmp_pphase_age,
      !!rlang::sym(nmiss_col_pphase_age)     := n_missing_week_tmp_pphase_age,
      !!rlang::sym(weekmean_col_pphase_age)  := week_mean_tmp_pphase_age
    ) %>%
    ungroup() %>%
    dplyr::select(-n_missing_week_tmp_pphase_age, -week_mean_tmp_pphase_age, -week_ok_tmp_pphase_age,
           -week_excluded_tmp_pphase_age, -imputed_day_tmp_pphase_age, -value_clean_tmp_pphase_age)
}

# ---- weather (unit handling) ----
needs_div10_pphase_age <- function(x_pphase_age, thr_pphase_age) {
  x_pphase_age <- suppressWarnings(as.numeric(x_pphase_age))
  q_pphase_age <- suppressWarnings(quantile(abs(x_pphase_age), probs = 0.99, na.rm = TRUE, names = FALSE))
  if (is.na(q_pphase_age)) return(FALSE)
  q_pphase_age > thr_pphase_age
}

pull_weather_sju_pphase_age <- function(year_value_pphase_age, station_id_pphase_age) {
  vars_pphase_age <- c("PRCP", "TMAX", "TMIN", "AWND", "ADPT", "RHAV")

  w_pphase_age <- meteo_pull_monitors(
    monitors   = station_id_pphase_age,
    var        = vars_pphase_age,
    date_min   = sprintf("%d-01-01", year_value_pphase_age),
    date_max   = sprintf("%d-12-31", year_value_pphase_age),
    keep_flags = FALSE
  ) %>%
    as_tibble() %>%
    mutate(date = as.Date(date))

  for (nm_pphase_age in c("prcp","tmax","tmin","awnd","adpt","rhav")) {
    if (!nm_pphase_age %in% names(w_pphase_age)) w_pphase_age[[nm_pphase_age]] <- NA_real_
  }

  prcp_mm_pphase_age     <- if (needs_div10_pphase_age(w_pphase_age$prcp, thr_pphase_age = 300)) w_pphase_age$prcp / 10 else w_pphase_age$prcp
  tmax_c_pphase_age      <- if (needs_div10_pphase_age(w_pphase_age$tmax, thr_pphase_age = 80))  w_pphase_age$tmax / 10 else w_pphase_age$tmax
  tmin_c_pphase_age      <- if (needs_div10_pphase_age(w_pphase_age$tmin, thr_pphase_age = 80))  w_pphase_age$tmin / 10 else w_pphase_age$tmin
  dewpoint_c_pphase_age  <- if (needs_div10_pphase_age(w_pphase_age$adpt, thr_pphase_age = 80))  w_pphase_age$adpt / 10 else w_pphase_age$adpt
  wind_ms_pphase_age     <- if (needs_div10_pphase_age(w_pphase_age$awnd, thr_pphase_age = 30))  w_pphase_age$awnd / 10 else w_pphase_age$awnd

  tibble(
    date = w_pphase_age$date,
    prcp_mm_pphase_age     = prcp_mm_pphase_age,
    tmax_c_pphase_age      = tmax_c_pphase_age,
    tmin_c_pphase_age      = tmin_c_pphase_age,
    tmean_c_pphase_age     = (tmax_c_pphase_age + tmin_c_pphase_age) / 2,
    dewpoint_c_pphase_age  = dewpoint_c_pphase_age,
    rh_mean_pphase_age     = w_pphase_age$rhav,
    wind_ms_pphase_age     = wind_ms_pphase_age
  )
}

# ---- linear-combination CI for exp(beta) ----
lincomb_exp_pphase_age <- function(beta_pphase_age, V_pphase_age, idx_pphase_age, label_pphase_age = "effect_pphase_age") {
  b_pphase_age <- sum(beta_pphase_age[names(idx_pphase_age)] * idx_pphase_age, na.rm = TRUE)
  v_pphase_age <- as.numeric(t(idx_pphase_age) %*% V_pphase_age[names(idx_pphase_age), names(idx_pphase_age), drop = FALSE] %*% idx_pphase_age)
  se_pphase_age <- sqrt(v_pphase_age)
  ci_pphase_age <- b_pphase_age + c(-1, 1) * 1.96 * se_pphase_age

  tibble(
    log_est_pphase_age = b_pphase_age,
    log_se_pphase_age  = se_pphase_age,
    !!rlang::sym(label_pphase_age) := exp(b_pphase_age),
    low_pphase_age  = exp(ci_pphase_age[1]),
    high_pphase_age = exp(ci_pphase_age[2])
  )
}

# Extract slope of base_term across two factor variables in a (base * f1 * f2) model
get_group_effects_2way_pphase_age <- function(fit_pphase_age, base_term_pphase_age, group_var1_pphase_age, group_var2_pphase_age, label_pphase_age) {
  beta_pphase_age <- coef(fit_pphase_age)
  V_pphase_age <- vcov(fit_pphase_age)

  lvls1_pphase_age <- levels(model.frame(fit_pphase_age)[[group_var1_pphase_age]])
  lvls2_pphase_age <- levels(model.frame(fit_pphase_age)[[group_var2_pphase_age]])
  ref1_pphase_age  <- lvls1_pphase_age[1]
  ref2_pphase_age  <- lvls2_pphase_age[1]

  map_dfr(lvls1_pphase_age, function(lv1_pphase_age) {
    map_dfr(lvls2_pphase_age, function(lv2_pphase_age) {

      idx_pphase_age <- rep(0, length(beta_pphase_age)); names(idx_pphase_age) <- names(beta_pphase_age)
      idx_pphase_age[base_term_pphase_age] <- 1

      # base:group1
      if (lv1_pphase_age != ref1_pphase_age) {
        t1 <- paste0(base_term_pphase_age, ":", group_var1_pphase_age, lv1_pphase_age)
        t2 <- paste0(group_var1_pphase_age, lv1_pphase_age, ":", base_term_pphase_age)
        if (t1 %in% names(beta_pphase_age)) idx_pphase_age[t1] <- 1
        if (t2 %in% names(beta_pphase_age)) idx_pphase_age[t2] <- 1
      }

      # base:group2
      if (lv2_pphase_age != ref2_pphase_age) {
        t1 <- paste0(base_term_pphase_age, ":", group_var2_pphase_age, lv2_pphase_age)
        t2 <- paste0(group_var2_pphase_age, lv2_pphase_age, ":", base_term_pphase_age)
        if (t1 %in% names(beta_pphase_age)) idx_pphase_age[t1] <- 1
        if (t2 %in% names(beta_pphase_age)) idx_pphase_age[t2] <- 1
      }

      # base:group1:group2 (3-way)
      if (lv1_pphase_age != ref1_pphase_age && lv2_pphase_age != ref2_pphase_age) {
        patterns <- c(
          paste0(base_term_pphase_age, ":", group_var1_pphase_age, lv1_pphase_age, ":", group_var2_pphase_age, lv2_pphase_age),
          paste0(base_term_pphase_age, ":", group_var2_pphase_age, lv2_pphase_age, ":", group_var1_pphase_age, lv1_pphase_age),
          paste0(group_var1_pphase_age, lv1_pphase_age, ":", base_term_pphase_age, ":", group_var2_pphase_age, lv2_pphase_age),
          paste0(group_var2_pphase_age, lv2_pphase_age, ":", base_term_pphase_age, ":", group_var1_pphase_age, lv1_pphase_age)
        )
        for (pat in patterns) {
          if (pat %in% names(beta_pphase_age)) idx_pphase_age[pat] <- 1
        }
      }

      est <- lincomb_exp_pphase_age(beta_pphase_age, V_pphase_age,
                                idx_pphase_age[names(idx_pphase_age) %in% names(beta_pphase_age)],
                                label_pphase_age = label_pphase_age)

      tibble(
        !!rlang::sym(group_var1_pphase_age) := lv1_pphase_age,
        !!rlang::sym(group_var2_pphase_age) := lv2_pphase_age
      ) %>% bind_cols(est)
    })
  })
}

# Extract slope across ONE factor variable in a (base * f) model
get_group_effects_pphase_age <- function(fit_pphase_age, base_term_pphase_age, group_var_pphase_age, label_pphase_age) {
  beta_pphase_age <- coef(fit_pphase_age)
  V_pphase_age <- vcov(fit_pphase_age)

  lvls_pphase_age <- levels(model.frame(fit_pphase_age)[[group_var_pphase_age]])
  ref_pphase_age  <- lvls_pphase_age[1]

  map_dfr(lvls_pphase_age, function(lv_pphase_age) {
    idx_pphase_age <- rep(0, length(beta_pphase_age)); names(idx_pphase_age) <- names(beta_pphase_age)
    idx_pphase_age[base_term_pphase_age] <- 1

    if (lv_pphase_age != ref_pphase_age) {
      t1 <- paste0(base_term_pphase_age, ":", group_var_pphase_age, lv_pphase_age)
      t2 <- paste0(group_var_pphase_age, lv_pphase_age, ":", base_term_pphase_age)
      if (t1 %in% names(beta_pphase_age)) idx_pphase_age[t1] <- 1
      if (t2 %in% names(beta_pphase_age)) idx_pphase_age[t2] <- 1
    }

    est <- lincomb_exp_pphase_age(beta_pphase_age, V_pphase_age,
                              idx_pphase_age[names(idx_pphase_age) %in% names(beta_pphase_age)],
                              label_pphase_age = label_pphase_age)

    tibble(level_pphase_age = lv_pphase_age) %>% bind_cols(est)
  })
}

# ===========================================================
# 3) DAILY ASTHMA OUTCOMES (San Juan health region) BY AGE
# ===========================================================

claims_raw_pphase_age <- read_excel(claims_path_pphase_age, sheet = claims_sheet_pphase_age) %>%
  clean_names() %>%
  mutate(from_date = parse_date_any_pphase_age(from_date)) %>%
  ensure_care_setting_pphase_age() %>%
  ensure_asthma_primary_pphase_age() %>%
  ensure_age_group_pphase_age()

claims_pphase_age <- claims_raw_pphase_age %>%
  mutate(
    municipality_desc_pphase_age = toupper(as.character(municipality_desc)),
    year_pphase_age = year(from_date),
    health_region_pphase_age = case_when(
      municipality_desc_pphase_age %in% san_juan_munis_pphase_age ~ "San Juan",
      municipality_desc_pphase_age %in% caguas_munis_pphase_age   ~ "Caguas",
      TRUE                                                ~ NA_character_
    ),
    health_region_pphase_age = factor(health_region_pphase_age, levels = c("San Juan", "Caguas"))
  ) %>%
  filter(year_pphase_age %in% years_all_pphase_age) %>%
  filter(!is.na(health_region_pphase_age))

claims_sj_pphase_age <- claims_pphase_age %>%
  filter(
    health_region_pphase_age == "San Juan",
    asthma_primary_pphase_age == "asthma",
    care_setting_pphase_age %in% clinical_settings_pphase_age,
    !is.na(age_group_pphase_age)
  ) %>%
  mutate(
    acute_flag_pphase_age = as.integer(care_setting_pphase_age %in% c("ED", "Inpatient")),
    age_group_pphase_age = droplevels(age_group_pphase_age)
  )

calendar_pphase_age <- tibble(
  date = seq.Date(
    as.Date(sprintf("%d-01-01", min(years_all_pphase_age))),
    as.Date(sprintf("%d-12-31", max(years_all_pphase_age))),
    by = "day"
  )
) %>%
  mutate(
    year_pphase_age = year(date),
    doy_pphase_age  = yday(date),
    dow_pphase_age  = wday(date, label = TRUE, week_start = 1),
    week_start_pphase_age = floor_date(date, unit = "week", week_start = 1),
    week_of_year_pphase_age = isoweek(date),
    date_num_pphase_age = as.numeric(date)
  ) %>%
  filter(year_pphase_age %in% years_all_pphase_age)

daily_asthma_pphase_age <- claims_sj_pphase_age %>%
  group_by(date = from_date, age_group_pphase_age) %>%
  summarise(
    total_claims_pphase_age = n(),
    acute_claims_pphase_age = sum(acute_flag_pphase_age == 1, na.rm = TRUE),
    nonacute_claims_pphase_age = total_claims_pphase_age - acute_claims_pphase_age,
    .groups = "drop"
  )

age_levels_pphase_age <- levels(claims_sj_pphase_age$age_group_pphase_age)

daily_outcome_pphase_age <- tidyr::expand_grid(
  date = calendar_pphase_age$date,
  age_group_pphase_age = factor(age_levels_pphase_age, levels = age_levels_pphase_age)
) %>%
  left_join(calendar_pphase_age, by = "date") %>%
  left_join(daily_asthma_pphase_age, by = c("date", "age_group_pphase_age")) %>%
  mutate(
    across(c(total_claims_pphase_age, acute_claims_pphase_age, nonacute_claims_pphase_age),
           ~ replace_na(.x, 0L))
  )

# ===========================================================
# 4) DAILY SPORES/POLLEN with weekly QC (San Juan station)
# ===========================================================

spores_raw_pphase_age <- map_dfr(spore_files_sj_pphase_age, read_spores_sj_onefile_pphase_age) %>%
  group_by(date) %>%
  summarise(
    spores_total_raw_pphase_age = safe_mean_pphase_age(spores_total_raw_pphase_age),
    pollen_total_raw_pphase_age = safe_mean_pphase_age(pollen_total_raw_pphase_age),
    .groups = "drop"
  ) %>%
  arrange(date)

spores_daily_pphase_age <- tibble(date = calendar_pphase_age$date) %>%
  left_join(spores_raw_pphase_age, by = "date") %>%
  arrange(date)

spores_daily_pphase_age <- weekly_qc_impute_or_exclude_pphase_age(
  spores_daily_pphase_age,
  value_col_pphase_age = spores_total_raw_pphase_age,
  base_name_pphase_age = "spores_total",
  max_missing_pphase_age = max_missing_days_week_pphase_age
)

spores_daily_pphase_age <- weekly_qc_impute_or_exclude_pphase_age(
  spores_daily_pphase_age,
  value_col_pphase_age = pollen_total_raw_pphase_age,
  base_name_pphase_age = "pollen_total",
  max_missing_pphase_age = max_missing_days_week_pphase_age
)

spores_daily_pphase_age <- spores_daily_pphase_age %>%
  mutate(spores_10k_pphase_age = spores_total_pphase_age / spore_scale_divisor_pphase_age)

# ===========================================================
# 5) WEATHER + MERGE into DAILY dataset + MA7 + pandemic phase
# ===========================================================

weather_daily_pphase_age <- map_dfr(years_all_pphase_age, ~ pull_weather_sju_pphase_age(.x, sju_station_id_pphase_age))

daily_pphase_age <- daily_outcome_pphase_age %>%
  left_join(spores_daily_pphase_age,  by = "date") %>%
  left_join(weather_daily_pphase_age, by = "date") %>%
  mutate(
    pandemic_phase_pphase_age = case_when(
      year_pphase_age < 2020  ~ "pre2020",
      year_pphase_age == 2020 ~ "y2020",
      year_pphase_age >= 2021 ~ "y2021plus"
    ),
    pandemic_phase_pphase_age = factor(pandemic_phase_pphase_age, levels = c("pre2020", "y2020", "y2021plus"))
  ) %>%
  arrange(age_group_pphase_age, date) %>%
  group_by(age_group_pphase_age) %>%
  mutate(
    # MA7 is computed within age group over the full period, so it crosses year boundaries.
    # For a within-year MA7, group_by(age_group_pphase_age, year_pphase_age) instead.
    spores_10k_ma7_pphase_age = zoo::rollapply(spores_10k_pphase_age, width = 7, FUN = mean,
                                          align = "right", fill = NA_real_, partial = FALSE)
  ) %>%
  ungroup()

weather_covars_use_pphase_age <- intersect(weather_covars_pphase_age, names(daily_pphase_age))
if (length(weather_covars_use_pphase_age) == 0) stop("No weather covariates found after weather pull.")

# ===========================================================
# A) BASELINE (2018–2019): weekly NB + weekly logit with spores × year × age
# ===========================================================

weekly_baseline_pphase_age <- daily_pphase_age %>%
  filter(year_pphase_age %in% years_baseline_pphase_age) %>%
  filter(week_of_year_pphase_age <= 52) %>%  # drop week 53
  mutate(week_start_pphase_age = floor_date(date, unit = "week", week_start = 1)) %>%
  group_by(year_pphase_age, week_start_pphase_age, week_of_year_pphase_age, age_group_pphase_age) %>%
  summarise(
    total_claims_wk_pphase_age    = sum(total_claims_pphase_age),
    acute_claims_wk_pphase_age    = sum(acute_claims_pphase_age),
    nonacute_claims_wk_pphase_age = sum(nonacute_claims_pphase_age),

    spores_cum_wk_pphase_age = sum(spores_10k_pphase_age, na.rm = FALSE),

    across(all_of(weather_covars_use_pphase_age), ~ safe_mean_pphase_age(.x), .names = "{.col}_wkmean_pphase_age"),

    spores_week_excluded_pphase_age = any(spores_total_week_excluded_pphase_age, na.rm = TRUE),
    pollen_week_excluded_pphase_age = any(pollen_total_week_excluded_pphase_age, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  filter(!spores_week_excluded_pphase_age) %>%
  filter(!is.na(spores_cum_wk_pphase_age)) %>%
  # Drop weeks with missing weather covariates
  filter(if_all(ends_with("_wkmean_pphase_age"), ~ !is.na(.x))) %>%
  # Drop non-finite spore totals
  filter(is.finite(spores_cum_wk_pphase_age)) %>%
  mutate(
    year_f_pphase_age = factor(year_pphase_age),
    age_group_pphase_age = factor(age_group_pphase_age, levels = age_levels_pphase_age)
  )

# Stop if nothing is left after filtering
if (nrow(weekly_baseline_pphase_age) == 0) {
  stop("No valid weekly baseline data after filtering. Check spore and weather data availability.")
}

# Baseline data summary
cat("\n=== Data Summary for Baseline Models ===\n")
cat("Number of weeks:", nrow(weekly_baseline_pphase_age), "\n")
cat("Age groups:", paste(levels(weekly_baseline_pphase_age$age_group_pphase_age), collapse = ", "), "\n")
cat("Spores range:", range(weekly_baseline_pphase_age$spores_cum_wk_pphase_age, na.rm = TRUE), "\n")

available_years_pphase_age <- sort(unique(weekly_baseline_pphase_age$year_pphase_age))

if (length(available_years_pphase_age) >= 1) {
  weekly_baseline_pphase_age <- weekly_baseline_pphase_age %>%
    mutate(year_f_pphase_age = relevel(factor(year_pphase_age), ref = as.character(available_years_pphase_age[1])))
}

knots_week_pphase_age <- list(week_of_year_pphase_age = c(0.5, 52.5))
wk_weather_cols_pphase_age <- names(weekly_baseline_pphase_age)[stringr::str_ends(names(weekly_baseline_pphase_age), "_wkmean_pphase_age")]

n_weeks_pphase_age <- nrow(weekly_baseline_pphase_age)
k_week_pphase_age <- min(20, max(4, floor(n_weeks_pphase_age / 10)))

if (length(available_years_pphase_age) >= 2) {
  f_wk_nb_pphase_age <- as.formula(
    paste0(
      "total_claims_wk_pphase_age ~ spores_cum_wk_pphase_age * year_f_pphase_age * age_group_pphase_age + ",
      "s(week_of_year_pphase_age, bs='cc', k=", k_week_pphase_age, ")",
      if (length(wk_weather_cols_pphase_age) > 0) paste0(" + ", paste(wk_weather_cols_pphase_age, collapse = " + ")) else ""
    )
  )
} else {
  f_wk_nb_pphase_age <- as.formula(
    paste0(
      "total_claims_wk_pphase_age ~ spores_cum_wk_pphase_age * age_group_pphase_age + ",
      "s(week_of_year_pphase_age, bs='cc', k=", k_week_pphase_age, ")",
      if (length(wk_weather_cols_pphase_age) > 0) paste0(" + ", paste(wk_weather_cols_pphase_age, collapse = " + ")) else ""
    )
  )
}

cat("\n=== Fitting baseline weekly NB model (BY AGE) ===\n")
fit_wk_nb_pphase_age <- gam(
  f_wk_nb_pphase_age,
  data   = weekly_baseline_pphase_age,
  family = nb(),
  method = "REML",
  knots  = knots_week_pphase_age,
  control = gam.control(trace = TRUE, maxit = 100)
)

weekly_baseline_logit_pphase_age <- weekly_baseline_pphase_age %>% 
  filter(total_claims_wk_pphase_age > 0)

# Logit data checks
cat("\n=== Logit Model Data Summary ===\n")
cat("Number of weeks with claims > 0:", nrow(weekly_baseline_logit_pphase_age), "\n")
logit_years_pphase_age <- sort(unique(weekly_baseline_logit_pphase_age$year_pphase_age))
cat("Available years in logit data:", paste(logit_years_pphase_age, collapse = ", "), "\n")

# Check age group representation
age_counts_logit <- weekly_baseline_logit_pphase_age %>%
  count(age_group_pphase_age, year_pphase_age)
cat("\nAge group representation by year:\n")
print(age_counts_logit)

# Check for problematic values
cat("\nChecking for data issues:\n")
cat("Any NA in spores_cum_wk:", any(is.na(weekly_baseline_logit_pphase_age$spores_cum_wk_pphase_age)), "\n")
cat("Any Inf in spores_cum_wk:", any(is.infinite(weekly_baseline_logit_pphase_age$spores_cum_wk_pphase_age)), "\n")
cat("Any NA in weather covars:", any(is.na(weekly_baseline_logit_pphase_age[wk_weather_cols_pphase_age])), "\n")
cat("Spores range:", range(weekly_baseline_logit_pphase_age$spores_cum_wk_pphase_age), "\n")

# Additional filter for complete cases
weekly_baseline_logit_pphase_age <- weekly_baseline_logit_pphase_age %>%
  filter(
    is.finite(spores_cum_wk_pphase_age),
    acute_claims_wk_pphase_age >= 0,
    nonacute_claims_wk_pphase_age >= 0
  )

if (length(wk_weather_cols_pphase_age) > 0) {
  weekly_baseline_logit_pphase_age <- weekly_baseline_logit_pphase_age %>%
    filter(if_all(all_of(wk_weather_cols_pphase_age), ~ is.finite(.x)))
}

cat("\nRows after complete case filtering:", nrow(weekly_baseline_logit_pphase_age), "\n")

# Decide on model complexity based on available data
n_age_groups_logit <- n_distinct(weekly_baseline_logit_pphase_age$age_group_pphase_age)
min_obs_per_group <- weekly_baseline_logit_pphase_age %>%
  count(age_group_pphase_age, year_pphase_age) %>%
  pull(n) %>%
  min()

cat("Number of age groups:", n_age_groups_logit, "\n")
cat("Min observations per age×year:", min_obs_per_group, "\n")

# Build formula based on data availability
if (length(logit_years_pphase_age) >= 2 && n_age_groups_logit >= 2 && min_obs_per_group >= 3) {
  # Full 3-way interaction
  f_wk_logit_pphase_age <- as.formula(
    paste0(
      "cbind(acute_claims_wk_pphase_age, nonacute_claims_wk_pphase_age) ~ spores_cum_wk_pphase_age * year_f_pphase_age * age_group_pphase_age + ",
      "s(week_of_year_pphase_age, bs='cc', k=", k_week_pphase_age, ")",
      if (length(wk_weather_cols_pphase_age) > 0) paste0(" + ", paste(wk_weather_cols_pphase_age, collapse = " + ")) else ""
    )
  )
  cat("\nUsing full 3-way interaction model\n")
} else if (n_age_groups_logit >= 2) {
  # Age only (no year interaction)
  f_wk_logit_pphase_age <- as.formula(
    paste0(
      "cbind(acute_claims_wk_pphase_age, nonacute_claims_wk_pphase_age) ~ spores_cum_wk_pphase_age * age_group_pphase_age + ",
      "s(week_of_year_pphase_age, bs='cc', k=", k_week_pphase_age, ")",
      if (length(wk_weather_cols_pphase_age) > 0) paste0(" + ", paste(wk_weather_cols_pphase_age, collapse = " + ")) else ""
    )
  )
  cat("\nUsing age-stratified model (insufficient data for year interaction)\n")
} else {
  # Minimal model
  f_wk_logit_pphase_age <- as.formula(
    paste0(
      "cbind(acute_claims_wk_pphase_age, nonacute_claims_wk_pphase_age) ~ spores_cum_wk_pphase_age + ",
      "s(week_of_year_pphase_age, bs='cc', k=", k_week_pphase_age, ")",
      if (length(wk_weather_cols_pphase_age) > 0) paste0(" + ", paste(wk_weather_cols_pphase_age, collapse = " + ")) else ""
    )
  )
  cat("\nUsing minimal model (insufficient age group variation)\n")
}


cat("\n=== Fitting baseline weekly logit model (BY AGE) ===\n")

# Skip the logit model if there are fewer than 20 weeks with claims
if (nrow(weekly_baseline_logit_pphase_age) < 20) {
  cat("WARNING: Only", nrow(weekly_baseline_logit_pphase_age), "weeks with claims > 0.\n")
  cat("Skipping logit model due to insufficient data.\n")
  
  # Empty result row so downstream code still runs
  baseline_logit_effects_by_year_age_pphase_age <- tibble(
    age_group_pphase_age = "insufficient_data",
    year_f_pphase_age = "insufficient_data",
    log_est_pphase_age = NA_real_,
    log_se_pphase_age = NA_real_,
    OR_pphase_age = NA_real_,
    low_pphase_age = NA_real_,
    high_pphase_age = NA_real_
  )
  
} else {
  # Try to fit the model with error handling
  fit_wk_logit_pphase_age <- try(
    gam(
      f_wk_logit_pphase_age,
      data   = weekly_baseline_logit_pphase_age,
      family = binomial(link = "logit"),
      method = "REML",
      knots  = knots_week_pphase_age,
      control = gam.control(trace = TRUE, maxit = 100)
    ),
    silent = FALSE
  )
  
  if (inherits(fit_wk_logit_pphase_age, "try-error")) {
    cat("\nLogit model failed to converge. Trying simplified model without spline...\n")
    
    # Simplest possible model
    f_wk_logit_simple <- as.formula(
      "cbind(acute_claims_wk_pphase_age, nonacute_claims_wk_pphase_age) ~ spores_cum_wk_pphase_age + age_group_pphase_age"
    )
    
    fit_wk_logit_pphase_age <- try(
      gam(
        f_wk_logit_simple,
        data   = weekly_baseline_logit_pphase_age,
        family = binomial(link = "logit"),
        method = "REML"
      ),
      silent = FALSE
    )
    
    if (inherits(fit_wk_logit_pphase_age, "try-error")) {
      cat("Even simplified logit model failed. Creating placeholder results.\n")
      baseline_logit_effects_by_year_age_pphase_age <- tibble(
        age_group_pphase_age = "model_failed",
        year_f_pphase_age = "model_failed",
        log_est_pphase_age = NA_real_,
        log_se_pphase_age = NA_real_,
        OR_pphase_age = NA_real_,
        low_pphase_age = NA_real_,
        high_pphase_age = NA_real_
      )
    } else {
      # Extract simple effects
      baseline_logit_effects_by_year_age_pphase_age <- get_group_effects_pphase_age(
        fit_pphase_age = fit_wk_logit_pphase_age,
        base_term_pphase_age = "spores_cum_wk_pphase_age",
        group_var_pphase_age = "age_group_pphase_age",
        label_pphase_age = "OR_pphase_age"
      ) %>%
        rename(age_group_pphase_age = level_pphase_age) %>%
        mutate(year_f_pphase_age = "pooled_simple")
    }
  } else {
    # Model converged: extract effects
    if (length(logit_years_pphase_age) >= 2 && n_age_groups_logit >= 2 && min_obs_per_group >= 3) {
      baseline_logit_effects_by_year_age_pphase_age <- get_group_effects_2way_pphase_age(
        fit_pphase_age = fit_wk_logit_pphase_age,
        base_term_pphase_age = "spores_cum_wk_pphase_age",
        group_var1_pphase_age = "year_f_pphase_age",
        group_var2_pphase_age = "age_group_pphase_age",
        label_pphase_age = "OR_pphase_age"
      )
    } else if (n_age_groups_logit >= 2) {
      baseline_logit_effects_by_year_age_pphase_age <- get_group_effects_pphase_age(
        fit_pphase_age = fit_wk_logit_pphase_age,
        base_term_pphase_age = "spores_cum_wk_pphase_age",
        group_var_pphase_age = "age_group_pphase_age",
        label_pphase_age = "OR_pphase_age"
      ) %>%
        rename(age_group_pphase_age = level_pphase_age) %>%
        mutate(year_f_pphase_age = "pooled")
    } else {
      beta <- coef(fit_wk_logit_pphase_age)
      V <- vcov(fit_wk_logit_pphase_age)
      baseline_logit_effects_by_year_age_pphase_age <- lincomb_exp_pphase_age(
        beta, V, c(spores_cum_wk_pphase_age = 1), "OR_pphase_age"
      ) %>%
        mutate(year_f_pphase_age = "pooled", age_group_pphase_age = "all")
    }
  }
}

baseline_nb_effects_by_year_age_pphase_age <- get_group_effects_2way_pphase_age(
  fit_pphase_age = fit_wk_nb_pphase_age,
  base_term_pphase_age = "spores_cum_wk_pphase_age",
  group_var1_pphase_age = "year_f_pphase_age",
  group_var2_pphase_age = "age_group_pphase_age",
  label_pphase_age = "RR_pphase_age"
)

# ===========================================================
# B) PANDEMIC PHASE: daily NB + daily logit with spores × phase × age
# ===========================================================

daily_phase_pphase_age <- daily_pphase_age %>%
  filter(!is.na(spores_10k_ma7_pphase_age)) %>%
  filter(if_all(all_of(weather_covars_use_pphase_age), ~ !is.na(.x))) %>%
  mutate(pandemic_phase_pphase_age = droplevels(pandemic_phase_pphase_age))

# cyclic spline boundary knots for day-of-year
knots_doy_pphase_age <- list(doy_pphase_age = c(0.5, 366.5))

n_days_pphase_age <- nrow(daily_phase_pphase_age)
k_doy_pphase_age <- min(30, max(6, floor(n_days_pphase_age / 50)))

day_weather_cols_pphase_age <- intersect(weather_covars_use_pphase_age, names(daily_phase_pphase_age))

# Drop incomplete rows before building the formula
daily_phase_pphase_age <- daily_phase_pphase_age %>%
  filter(
    !is.na(spores_10k_ma7_pphase_age),
    is.finite(spores_10k_ma7_pphase_age),
    !is.na(total_claims_pphase_age),
    is.finite(total_claims_pphase_age),
    !is.na(pandemic_phase_pphase_age),
    !is.na(age_group_pphase_age),
    !is.na(doy_pphase_age),
    is.finite(doy_pphase_age),
    !is.na(dow_pphase_age),
    !is.na(date_num_pphase_age),
    is.finite(date_num_pphase_age)
  )

# Filter weather covariates if present
if (length(day_weather_cols_pphase_age) > 0) {
  daily_phase_pphase_age <- daily_phase_pphase_age %>%
    filter(if_all(all_of(day_weather_cols_pphase_age), ~ !is.na(.x) & is.finite(.x)))
}

# Drop any remaining rows with NA or Inf in any column
daily_phase_pphase_age <- daily_phase_pphase_age %>%
  filter(complete.cases(.)) %>%
  filter(if_all(where(is.numeric), is.finite))

# Check the model matrix before fitting
cat("\n=== Pre-flight checks before model fitting ===\n")

# Test model matrix creation
test_formula <- as.formula(
  paste0(
    "total_claims_pphase_age ~ spores_10k_ma7_pphase_age * pandemic_phase_pphase_age * age_group_pphase_age + ",
    "doy_pphase_age + dow_pphase_age + date_num_pphase_age",
    if (length(day_weather_cols_pphase_age) > 0) paste0(" + ", paste(day_weather_cols_pphase_age, collapse = " + ")) else ""
  )
)

test_mm <- try(model.matrix(test_formula, data = daily_phase_pphase_age), silent = TRUE)

if (inherits(test_mm, "try-error")) {
  cat("ERROR: Cannot create model matrix. Check factor levels and data.\n")
  stop("Model matrix creation failed.")
}

# Check for NA/Inf in model matrix
if (any(is.na(test_mm))) {
  cat("WARNING: Model matrix contains NA values\n")
  cat("Columns with NA:", colnames(test_mm)[colSums(is.na(test_mm)) > 0], "\n")
}

if (any(is.infinite(test_mm))) {
  cat("WARNING: Model matrix contains Inf values\n")
  cat("Columns with Inf:", colnames(test_mm)[colSums(is.infinite(test_mm)) > 0], "\n")
}

# Diagnostics AFTER filtering
cat("\n=== Diagnostic checks for daily NB model ===\n")
cat("Rows after complete filtering:", nrow(daily_phase_pphase_age), "\n")

if (nrow(daily_phase_pphase_age) == 0) {
  stop("No valid data remaining after filtering. Check your spore and weather data.")
}

# Check group representation
group_counts <- daily_phase_pphase_age %>%
  count(pandemic_phase_pphase_age, age_group_pphase_age)
cat("\nObservations per phase × age group:\n")
print(group_counts)

if (any(group_counts$n < 10)) {
  warning("Some phase×age combinations have < 10 observations. Model may be unstable.")
}

# Summary statistics for key variables
cat("\nKey variable ranges:\n")
cat("total_claims:", range(daily_phase_pphase_age$total_claims_pphase_age), "\n")
cat("spores_10k_ma7:", range(daily_phase_pphase_age$spores_10k_ma7_pphase_age), "\n")
cat("doy:", range(daily_phase_pphase_age$doy_pphase_age), "\n")
cat("date_num:", range(daily_phase_pphase_age$date_num_pphase_age), "\n")

# Standardize date_num for numerical stability
daily_phase_pphase_age <- daily_phase_pphase_age %>%
  mutate(
    date_num_scaled_pphase_age = scale(date_num_pphase_age)[,1],
    doy_pphase_age = as.numeric(doy_pphase_age)
  )

cat("date_num_scaled range:", range(daily_phase_pphase_age$date_num_scaled_pphase_age), "\n")

# Daily NB uses the standardized date_num
f_day_nb_pphase_age <- as.formula(
  paste0(
    "total_claims_pphase_age ~ spores_10k_ma7_pphase_age * pandemic_phase_pphase_age * age_group_pphase_age + ",
    "s(doy_pphase_age, bs='cc', k=", k_doy_pphase_age, ") + dow_pphase_age + s(date_num_pphase_age, bs='cr', k=10)",
    if (length(day_weather_cols_pphase_age) > 0) paste0(" + ", paste(day_weather_cols_pphase_age, collapse = " + ")) else ""
  )
)

cat("\n=== Fitting pandemic-phase daily NB model (BY AGE) ===\n")
cat("Formula:", deparse(f_day_nb_pphase_age), "\n")

f_day_nb_pphase_age <- as.formula(
  paste0(
    "total_claims_pphase_age ~ spores_10k_ma7_pphase_age * pandemic_phase_pphase_age * age_group_pphase_age + ",
    "s(doy_pphase_age, bs='cc', k=", k_doy_pphase_age, ") + dow_pphase_age + s(date_num_pphase_age, bs='cr', k=10)",
    if (length(day_weather_cols_pphase_age) > 0) paste0(" + ", paste(day_weather_cols_pphase_age, collapse = " + ")) else ""
  )
)

cat("\n=== Fitting pandemic-phase daily NB model (BY AGE) ===\n")
cat("Formula:", deparse(f_day_nb_pphase_age), "\n")

fit_day_nb_pphase_age <- gam(
  f_day_nb_pphase_age,
  data   = daily_phase_pphase_age,
  family = nb(),
  method = "REML",
  knots  = knots_doy_pphase_age,
  control = gam.control(trace = TRUE, maxit = 100)
)

phase_nb_effects_pphase_age <- get_group_effects_2way_pphase_age(
  fit_pphase_age = fit_day_nb_pphase_age,
  base_term_pphase_age = "spores_10k_ma7_pphase_age",
  group_var1_pphase_age = "pandemic_phase_pphase_age",
  group_var2_pphase_age = "age_group_pphase_age",
  label_pphase_age = "RR_pphase_age"
)

# Daily logit model (days with at least one claim)
daily_phase_logit_cc_pphase_age <- daily_phase_pphase_age %>%
  filter(total_claims_pphase_age > 0)

f_day_logit_pphase_age <- as.formula(
  paste0(
    "cbind(acute_claims_pphase_age, nonacute_claims_pphase_age) ~ spores_10k_ma7_pphase_age * pandemic_phase_pphase_age * age_group_pphase_age + ",
    "s(doy_pphase_age, bs='cc', k=", k_doy_pphase_age, ") + dow_pphase_age + s(date_num_pphase_age, bs='cr', k=10)",
    if (length(day_weather_cols_pphase_age) > 0) paste0(" + ", paste(day_weather_cols_pphase_age, collapse = " + ")) else ""
  )
)

cat("\n=== Fitting pandemic-phase daily logit model (BY AGE) ===\n")
fit_day_logit_pphase_age <- gam(
  f_day_logit_pphase_age,
  data   = daily_phase_logit_cc_pphase_age,
  family = binomial(link = "logit"),
  method = "REML",
  knots  = knots_doy_pphase_age,
  control = gam.control(trace = TRUE, maxit = 100)
)

phase_logit_effects_pphase_age <- get_group_effects_2way_pphase_age(
  fit_pphase_age = fit_day_logit_pphase_age,
  base_term_pphase_age = "spores_10k_ma7_pphase_age",
  group_var1_pphase_age = "pandemic_phase_pphase_age",
  group_var2_pphase_age = "age_group_pphase_age",
  label_pphase_age = "OR_pphase_age"
)

write.csv(baseline_nb_effects_by_year_age_pphase_age,
          file.path(out_dir_pphase_age, "baseline_weekly_NB_spores_by_year_age_pphase_age.csv"),
          row.names = FALSE)

write.csv(baseline_logit_effects_by_year_age_pphase_age,
          file.path(out_dir_pphase_age, "baseline_weekly_LOGIT_spores_by_year_age_pphase_age.csv"),
          row.names = FALSE)

write.csv(phase_nb_effects_pphase_age,
          file.path(out_dir_pphase_age, "pooled_daily_NB_spores_ma7_by_pandemic_phase_age_pphase_age.csv"),
          row.names = FALSE)

write.csv(phase_logit_effects_pphase_age,
          file.path(out_dir_pphase_age, "pooled_daily_LOGIT_spores_ma7_by_pandemic_phase_age_pphase_age.csv"),
          row.names = FALSE)

cat("\n=== Summary Statistics by Age Group ===\n")
daily_pphase_age %>%
  group_by(age_group_pphase_age) %>%
  summarise(
    n_days_pphase_age = n(),
    total_claims_pphase_age = sum(total_claims_pphase_age),
    acute_claims_pphase_age = sum(acute_claims_pphase_age),
    acute_prop_pphase_age = if_else(total_claims_pphase_age > 0, acute_claims_pphase_age / total_claims_pphase_age, NA_real_),
    .groups = "drop"
  ) %>%
  print()

cat("\nDone. Age-stratified outputs written to:", out_dir_pphase_age, "\n")
cat("Files created:\n")
cat("  - baseline_weekly_NB_spores_by_year_age_pphase_age.csv\n")
cat("  - baseline_weekly_LOGIT_spores_by_year_age_pphase_age.csv\n")
cat("  - pooled_daily_NB_spores_ma7_by_pandemic_phase_age_pphase_age.csv\n")
cat("  - pooled_daily_LOGIT_spores_ma7_by_pandemic_phase_age_pphase_age.csv\n")
