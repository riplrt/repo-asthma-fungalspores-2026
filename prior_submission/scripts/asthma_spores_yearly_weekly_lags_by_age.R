############################################################
## Spores ↔ Asthma models: weekly and daily lag/MA analyses by age
## (San Juan health region, 2018–2019) - STRATIFIED BY AGE
##
## Models (overall + by season + by age group):
##   - Weekly:  Spearman + partial Spearman (controls weekly mean weather)
##   - Weekly:  Binomial logit (acute fraction) + NB regression (total counts)
##   - Daily :  Binomial logit + NB regression using lags / moving averages
##
## Age groups: 0–17 (children), 18–64 (adults), 65+ (elderly)
##
## Exposure QC (weekly rule; applied separately to spores & pollen):
##   - If a week has <= 3 missing days: impute missing days with that week's mean
##   - If a week has  > 3 missing days OR all values missing: exclude the week
##     (set the entire week to NA so it drops from models)
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(lubridate)
  library(rnoaa)
  library(ppcor)
  library(MASS)
  library(splines)
  library(zoo)
})

# ===========================================================
# 1) SETTINGS
# ===========================================================

claims_path_wlag_age  <- "data/2022_asma_From_date_ver 3.xlsx"
claims_sheet_wlag_age <- 1

# One or more files; will be row-bound
spore_files_sj_wlag_age <- c("data/Datos de esporas San Juan Jan 2018-2022.xlsx")

years_to_run_wlag_age <- 2018:2019

clinical_settings_wlag_age <- c("ED", "Inpatient", "Urgent", "Outpatient")

# Age group definitions (0–17, 18–64, 65+)
age_breaks_wlag_age <- c(0, 18, 65, Inf)
age_labels_wlag_age <- c("0-17", "18-64", "65+")

# Scaling for spores in models (per 10,000 spores)
spore_scale_divisor_wlag_age <- 10000

# Weekly imputation/exclusion rule
max_missing_days_week_wlag_age <- 3

# Minimum thresholds
min_weeks_overall_wlag_age <- 10
min_weeks_season_wlag_age  <- 6
min_days_overall_wlag_age  <- 60
min_days_season_wlag_age   <- 40

# GHCND station ID for SJU area (meteo_pull_monitors expects a GHCND ID)
sju_station_id_wlag_age <- "RQW00011641"
if (toupper(sju_station_id_wlag_age) == "SJU") sju_station_id_wlag_age <- "RQW00011641"

# Weather covariates to control for (from pull_weather_sju_wlag_age output)
weather_covars_wlag_age <- c("tmean_c_wlag_age", "rh_mean_wlag_age", "wind_ms_wlag_age", "prcp_mm_wlag_age", "dewpoint_c_wlag_age")

# Output directory
out_dir_wlag_age <- "statistical_results"
dir.create(out_dir_wlag_age, showWarnings = FALSE, recursive = TRUE)

# ===========================================================
# 2) HELPERS
# ===========================================================

parse_date_any_wlag_age <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30")) # Excel numeric dates
  x <- as.character(x)

  d <- suppressWarnings(ymd(x))
  d <- coalesce(d, suppressWarnings(mdy(x)))
  d <- coalesce(d, suppressWarnings(dmy(x)))
  coalesce(d, suppressWarnings(as.Date(x)))
}

safe_mean_wlag_age <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

assign_season_wlag_age <- function(date) {
  m <- month(date)
  case_when(
    m %in% c(12, 1, 2)  ~ "winter",
    m %in% c(3, 4, 5)   ~ "spring",
    m %in% c(6, 7, 8)   ~ "summer",
    m %in% c(9, 10, 11) ~ "fall",
    TRUE                ~ NA_character_
  )
}

# ---- Claims helpers ----

read_claims_wlag_age <- function(path, sheet) {
  read_excel(path, sheet = sheet) %>%
    clean_names() %>%
    mutate(from_date = parse_date_any_wlag_age(from_date))
}

care_setting_fun_wlag_age <- function(hcg_setting, hcg_mr_line_desc) {
  s <- tolower(as.character(hcg_setting))
  l <- tolower(as.character(hcg_mr_line_desc))

  if (grepl("ed|emergency", s) || grepl("\\bed\\b|\\ber\\b", l)) return("ED")
  if (grepl("inpatient|hospital", s) || grepl("inpatient|hospital", l)) return("Inpatient")
  if (grepl("urgent", s) || grepl("urgent", l)) return("Urgent")
  if (grepl("office|clinic|outpatient", s) || grepl("office|clinic", l)) return("Outpatient")
  if (grepl("dme|supplies", l)) return("DME")
  "Other"
}

ensure_care_setting_wlag_age <- function(df) {
  if ("care_setting_wlag_age" %in% names(df)) return(df)

  if ("care_setting" %in% names(df)) {
    return(df %>% mutate(care_setting_wlag_age = care_setting))
  }

  needed <- c("hcg_setting", "hcg_mr_line_desc")
  if (!all(needed %in% names(df))) {
    stop("No care_setting/care_setting_wlag_age and missing hcg_setting/hcg_mr_line_desc to create it.")
  }

  df %>%
    mutate(care_setting_wlag_age = mapply(care_setting_fun_wlag_age, hcg_setting, hcg_mr_line_desc))
}

ensure_asthma_primary_wlag_age <- function(df) {
  if ("asthma_primary_wlag_age" %in% names(df)) return(df)

  if ("asthma_primary" %in% names(df)) {
    return(df %>% mutate(asthma_primary_wlag_age = asthma_primary))
  }

  if (!"icd_diag_01_primary" %in% names(df)) {
    stop("Missing asthma_primary and icd_diag_01_primary; create asthma_primary before running.")
  }

  df %>%
    mutate(
      icd_diag_01_primary = str_to_upper(str_trim(as.character(icd_diag_01_primary))),
      asthma_primary_wlag_age = if_else(str_detect(icd_diag_01_primary, "^J45|^J46"), "asthma", "other")
    )
}

# ---- Age processing helper ----
ensure_age_group_wlag_age <- function(df) {
  if ("age_group_wlag_age" %in% names(df)) return(df)

  if (!"age" %in% names(df)) {
    stop("Missing 'age' column; cannot create age groups.")
  }

  df %>%
    mutate(
      age_raw_wlag_age = str_trim(as.character(age)),
      age_num_wlag_age = case_when(
        str_detect(age_raw_wlag_age, "\\+$") ~ suppressWarnings(as.numeric(str_remove(age_raw_wlag_age, "\\+"))),
        TRUE ~ suppressWarnings(as.numeric(age_raw_wlag_age))
      ),
      age_group_wlag_age = cut(
        age_num_wlag_age,
        breaks = age_breaks_wlag_age,
        labels = age_labels_wlag_age,
        right = FALSE,
        include.lowest = TRUE
      )
    ) %>%
    dplyr::select(-age_raw_wlag_age, -age_num_wlag_age)
}

# ---- Spores/pollen helpers ----

read_spores_sj_onefile_wlag_age <- function(path) {
  x <- read_excel(path) %>% clean_names()

  date_col <- if ("date" %in% names(x)) "date" else names(x)[1]

  # spores
  spore_candidates <- c("spores_total", "total_spores", "esporas_de_hongos")
  spore_col <- intersect(spore_candidates, names(x))[1]
  if (is.na(spore_col) || length(spore_col) == 0) {
    spore_col <- grep("esporas|spore", names(x), ignore.case = TRUE, value = TRUE)[1]
  }
  if (is.na(spore_col) || length(spore_col) == 0) stop("Could not identify a spores column in: ", path)

  # pollen (optional)
  pollen_candidates <- c("pollen_total", "polen_san_juan", "polen")
  pollen_col <- intersect(pollen_candidates, names(x))[1]
  if (is.na(pollen_col) || length(pollen_col) == 0) {
    pollen_col <- grep("polen|pollen", names(x), ignore.case = TRUE, value = TRUE)[1]
  }

  x %>%
    transmute(
      date = parse_date_any_wlag_age(.data[[date_col]]),
      spores_total_raw = suppressWarnings(as.numeric(.data[[spore_col]])),
      pollen_total_raw = if (!is.na(pollen_col) && length(pollen_col) > 0)
        suppressWarnings(as.numeric(.data[[pollen_col]])) else NA_real_
    ) %>%
    filter(!is.na(date))
}

weekly_qc_impute_or_exclude_wlag_age <- function(df, value_col, out_col, max_missing = 3) {
  value_sym <- rlang::ensym(value_col)

  df %>%
    mutate(week_start_tmp = floor_date(date, unit = "week", week_start = 1)) %>%
    group_by(week_start_tmp) %>%
    mutate(
      n_missing_week = sum(is.na(!!value_sym)),
      week_mean      = safe_mean_wlag_age(!!value_sym),
      week_ok        = (n_missing_week <= max_missing) & !is.na(week_mean),
      week_excluded  = !week_ok,
      imputed_day    = week_ok & is.na(!!value_sym),
      value_clean    = case_when(
        week_ok & is.na(!!value_sym) ~ week_mean,
        week_ok ~ as.numeric(!!value_sym),
        TRUE ~ NA_real_
      )
    ) %>%
    ungroup() %>%
    mutate(
      !!out_col := value_clean,
      !!paste0(out_col, "_imputed_wlag_age") := imputed_day,
      !!paste0(out_col, "_week_excluded_wlag_age") := week_excluded,
      !!paste0(out_col, "_n_missing_week_wlag_age") := n_missing_week,
      !!paste0(out_col, "_week_mean_wlag_age") := week_mean
    ) %>%
    dplyr::select(-week_start_tmp, -n_missing_week, -week_mean, -week_ok, -week_excluded, -imputed_day, -value_clean)
}

# ---- Weather helper (GHCND via rnoaa) ----
pull_weather_sju_wlag_age <- function(year_wlag_age, station_id_wlag_age) {
  vars <- c("PRCP", "TMAX", "TMIN", "AWND", "ADPT", "RHAV")

  w <- meteo_pull_monitors(
    monitors   = station_id_wlag_age,
    var        = vars,
    date_min   = sprintf("%d-01-01", year_wlag_age),
    date_max   = sprintf("%d-12-31", year_wlag_age),
    keep_flags = FALSE
  ) %>%
    as_tibble() %>%
    mutate(date = as.Date(date))

  for (nm in c("prcp","tmax","tmin","awnd","adpt","rhav")) {
    if (!nm %in% names(w)) w[[nm]] <- NA_real_
  }

  w %>%
    transmute(
      date = date,
      prcp_mm_wlag_age    = prcp / 10,                        # tenths mm -> mm
      tmax_c_wlag_age     = tmax / 10,                        # tenths °C -> °C
      tmin_c_wlag_age     = tmin / 10,
      tmean_c_wlag_age    = (tmax_c_wlag_age + tmin_c_wlag_age) / 2,
      wind_ms_wlag_age    = awnd / 10,                        # tenths m/s -> m/s
      rh_mean_wlag_age    = rhav,                             # %
      dewpoint_c_wlag_age = adpt / 10                         # tenths °C -> °C
    )
}

# ---- Extract exponentiated effect (OR/RR) + Wald CI ----
extract_effect_wlag_age <- function(model, term, exp_scale = TRUE) {
  sm <- summary(model)
  coefs <- sm$coefficients

  if (!term %in% rownames(coefs)) {
    return(tibble(estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
  }

  est <- coefs[term, "Estimate"]
  se  <- coefs[term, "Std. Error"]
  pcol <- intersect(colnames(coefs), c("Pr(>|z|)", "Pr(>|t|)"))
  pval <- if (length(pcol) == 1) coefs[term, pcol] else NA_real_

  ci <- est + c(-1, 1) * 1.96 * se

  if (exp_scale) {
    tibble(estimate = exp(est), ci_low = exp(ci[1]), ci_high = exp(ci[2]), p = pval)
  } else {
    tibble(estimate = est, ci_low = ci[1], ci_high = ci[2], p = pval)
  }
}

# ===========================================================
# 3) BUILD DAILY ASTHMA OUTCOMES (San Juan health region) - BY AGE GROUP
# ===========================================================

claims_df_wlag_age <- read_claims_wlag_age(claims_path_wlag_age, claims_sheet_wlag_age) %>%
  ensure_care_setting_wlag_age() %>%
  ensure_asthma_primary_wlag_age() %>%
  ensure_age_group_wlag_age() %>%
  mutate(
    municipality_desc = toupper(as.character(municipality_desc)),
    health_region_wlag_age = case_when(
      municipality_desc %in% c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA") ~ "San Juan",
      municipality_desc %in% c("CAGUAS", "GURABO", "JUNCOS", "SAN LORENZO") ~ "Caguas",
      TRUE ~ "Other"
    ),
    health_region_wlag_age = factor(health_region_wlag_age, levels = c("San Juan", "Caguas")),
    year_wlag_age = year(from_date)
  ) %>%
  filter(year_wlag_age %in% years_to_run_wlag_age)

claims_sj_region_wlag_age <- claims_df_wlag_age %>%
  filter(
    health_region_wlag_age == "San Juan",
    asthma_primary_wlag_age == "asthma",
    care_setting_wlag_age %in% clinical_settings_wlag_age,
    !is.na(age_group_wlag_age)
  ) %>%
  mutate(acute_flag_wlag_age = as.integer(care_setting_wlag_age %in% c("ED", "Inpatient")))

calendar_all_wlag_age <- tibble(
  date = seq.Date(
    as.Date(sprintf("%d-01-01", min(years_to_run_wlag_age))),
    as.Date(sprintf("%d-12-31", max(years_to_run_wlag_age))),
    by = "day"
  )
) %>%
  mutate(
    year_wlag_age = year(date),
    doy_wlag_age  = yday(date),
    dow_wlag_age  = wday(date, label = TRUE, week_start = 1),
    week_start_wlag_age = floor_date(date, unit = "week", week_start = 1),
    season_wlag_age = assign_season_wlag_age(date)
  ) %>%
  filter(year_wlag_age %in% years_to_run_wlag_age)

# Summarize by date AND age group
daily_asthma_sj_wlag_age <- claims_sj_region_wlag_age %>%
  group_by(date = from_date, age_group_wlag_age) %>%
  summarise(
    total_claims_wlag_age    = n(),
    acute_claims_wlag_age    = sum(acute_flag_wlag_age, na.rm = TRUE),
    nonacute_claims_wlag_age = total_claims_wlag_age - acute_claims_wlag_age,
    ed_claims_wlag_age       = sum(care_setting_wlag_age == "ED", na.rm = TRUE),
    inpt_claims_wlag_age     = sum(care_setting_wlag_age == "Inpatient", na.rm = TRUE),
    outpt_claims_wlag_age    = sum(care_setting_wlag_age == "Outpatient", na.rm = TRUE),
    urgent_claims_wlag_age   = sum(care_setting_wlag_age == "Urgent", na.rm = TRUE),
    .groups = "drop"
  )

# Create complete calendar grid for each age group
daily_sj_wlag_age <- calendar_all_wlag_age %>%
  crossing(age_group_wlag_age = factor(age_labels_wlag_age, levels = age_labels_wlag_age)) %>%
  left_join(daily_asthma_sj_wlag_age, by = c("date", "age_group_wlag_age")) %>%
  mutate(
    across(
      c(total_claims_wlag_age, acute_claims_wlag_age, nonacute_claims_wlag_age,
        ed_claims_wlag_age, inpt_claims_wlag_age, outpt_claims_wlag_age, urgent_claims_wlag_age),
      ~ replace_na(.x, 0L)
    )
  )

# ===========================================================
# 4) SPORES/POLLEN (San Juan station) + WEATHER
# ===========================================================

# Read and combine spores/pollen
spores_sj_raw_wlag_age <- map_dfr(spore_files_sj_wlag_age, read_spores_sj_onefile_wlag_age) %>%
  group_by(date) %>%
  summarise(
    spores_total_raw = safe_mean_wlag_age(spores_total_raw),
    pollen_total_raw = safe_mean_wlag_age(pollen_total_raw),
    .groups = "drop"
  ) %>%
  arrange(date)

# Expand to full daily calendar so missing dates count as missing
calendar_exp_wlag_age <- tibble(date = calendar_all_wlag_age$date)

spores_sj_daily_wlag_age <- calendar_exp_wlag_age %>%
  left_join(spores_sj_raw_wlag_age, by = "date") %>%
  arrange(date)

# Apply weekly QC separately (creates spores_total_wlag_age and pollen_total_wlag_age)
spores_sj_daily_wlag_age <- weekly_qc_impute_or_exclude_wlag_age(
  spores_sj_daily_wlag_age,
  value_col = spores_total_raw,
  out_col   = "spores_total_wlag_age",
  max_missing = max_missing_days_week_wlag_age
)

spores_sj_daily_wlag_age <- weekly_qc_impute_or_exclude_wlag_age(
  spores_sj_daily_wlag_age,
  value_col = pollen_total_raw,
  out_col   = "pollen_total_wlag_age",
  max_missing = max_missing_days_week_wlag_age
)

# Pull weather for each year and bind
weather_sj_wlag_age <- map_dfr(years_to_run_wlag_age, ~ pull_weather_sju_wlag_age(.x, sju_station_id_wlag_age))

# Merge exposures/weather with age-stratified outcomes
daily_merged_wlag_age <- daily_sj_wlag_age %>%
  left_join(
    spores_sj_daily_wlag_age %>%
      dplyr::select(date,
             spores_total_wlag_age, pollen_total_wlag_age,
             spores_total_wlag_age_imputed_wlag_age, spores_total_wlag_age_week_excluded_wlag_age, spores_total_wlag_age_n_missing_week_wlag_age,
             pollen_total_wlag_age_imputed_wlag_age, pollen_total_wlag_age_week_excluded_wlag_age, pollen_total_wlag_age_n_missing_week_wlag_age),
    by = "date"
  ) %>%
  left_join(weather_sj_wlag_age, by = "date") %>%
  arrange(age_group_wlag_age, date) %>%
  group_by(age_group_wlag_age, year_wlag_age) %>%
  mutate(
    spores_10k_wlag_age = spores_total_wlag_age / spore_scale_divisor_wlag_age,

    # lags (within year + age group ordering)
    spores_10k_lag0_wlag_age = spores_10k_wlag_age,
    spores_10k_lag1_wlag_age = lag(spores_10k_wlag_age, 1),
    spores_10k_lag2_wlag_age = lag(spores_10k_wlag_age, 2),
    spores_10k_lag3_wlag_age = lag(spores_10k_wlag_age, 3),

    # moving averages
    spores_10k_ma3_wlag_age = zoo::rollapply(spores_10k_wlag_age, width = 3, FUN = mean,
                                         align = "right", fill = NA_real_, partial = FALSE),
    spores_10k_ma7_wlag_age = zoo::rollapply(spores_10k_wlag_age, width = 7, FUN = mean,
                                         align = "right", fill = NA_real_, partial = FALSE)
  ) %>%
  ungroup()

weather_covars_use_wlag_age <- intersect(weather_covars_wlag_age, names(daily_merged_wlag_age))
if (length(weather_covars_use_wlag_age) == 0) stop("No requested weather covariates found after weather pull.")

# ===========================================================
# 5) DIAGNOSTICS
# ===========================================================

cat("\n=== Exposure QC diagnostics (daily) by age group ===\n")
daily_merged_wlag_age %>%
  group_by(year_wlag_age, age_group_wlag_age) %>%
  summarise(
    n_days = n(),
    n_spores_nonmiss = sum(!is.na(spores_10k_wlag_age)),
    n_weather_complete = sum(if_all(all_of(weather_covars_use_wlag_age), ~ !is.na(.x))),
    n_weeks_spores_excluded = n_distinct(week_start_wlag_age[spores_total_wlag_age_week_excluded_wlag_age %in% TRUE]),
    .groups = "drop"
  ) %>%
  print()

cat("\n=== Sample sizes by age group ===\n")
claims_sj_region_wlag_age %>%
  count(year_wlag_age, age_group_wlag_age) %>%
  print()

# ===========================================================
# 6) MODELS (overall + by season + by age group)
# ===========================================================

run_models_one_year_season_age_wlag_age <- function(
  year_value,
  season_value,
  age_value,
  dat,
  weather_covars_use,
  min_weeks,
  min_days
) {

  # Subset year + age + season
  dat0 <- dat %>% filter(year_wlag_age == year_value, age_group_wlag_age == age_value)
  if (season_value != "All seasons") dat0 <- dat0 %>% filter(season_wlag_age == season_value)

  # -------------------------------------------------------
  # A) WEEKLY data
  # -------------------------------------------------------
  dat_weekly <- dat0 %>%
    filter(!is.na(spores_10k_wlag_age)) %>%
    filter(if_all(all_of(weather_covars_use), ~ !is.na(.x))) %>%
    group_by(week_start_wlag_age) %>%
    summarise(
      total_claims_wlag_age    = sum(total_claims_wlag_age),
      acute_claims_wlag_age    = sum(acute_claims_wlag_age),
      nonacute_claims_wlag_age = sum(nonacute_claims_wlag_age),
      spores_mean_wlag_age     = mean(spores_10k_wlag_age, na.rm = TRUE),
      across(all_of(weather_covars_use), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
      .groups = "drop"
    )

  # A1) Correlations
  cor_row <- tibble(
    year_wlag_age = year_value,
    season_wlag_age = season_value,
    age_group_wlag_age = age_value,
    n_weeks_wlag_age = nrow(dat_weekly),
    spearman_rho_wlag_age = NA_real_,
    spearman_p_wlag_age   = NA_real_,
    partial_rho_wlag_age  = NA_real_,
    partial_p_wlag_age    = NA_real_
  )

  if (nrow(dat_weekly) >= min_weeks) {
    ct <- suppressWarnings(cor.test(dat_weekly$total_claims_wlag_age, dat_weekly$spores_mean_wlag_age, method = "spearman"))
    cor_row$spearman_rho_wlag_age <- unname(ct$estimate)
    cor_row$spearman_p_wlag_age   <- ct$p.value

    zmat <- dat_weekly %>% dplyr::select(starts_with("mean_")) %>% as.data.frame()
    ok <- complete.cases(dat_weekly$total_claims_wlag_age, dat_weekly$spores_mean_wlag_age, zmat)

    if (sum(ok) >= min_weeks && ncol(zmat) > 0) {
      pt <- tryCatch(
        suppressWarnings(ppcor::pcor.test(
          x = dat_weekly$spores_mean_wlag_age[ok],
          y = dat_weekly$total_claims_wlag_age[ok],
          z = zmat[ok, , drop = FALSE],
          method = "spearman"
        )),
        error = function(e) NULL
      )
      if (!is.null(pt)) {
        cor_row$partial_rho_wlag_age <- as.numeric(pt$estimate)
        cor_row$partial_p_wlag_age   <- as.numeric(pt$p.value)
      }
    }
  }

  # A2) Weekly logit: acute fraction
  weekly_logit_row <- tibble(
    year_wlag_age = year_value,
    season_wlag_age = season_value,
    age_group_wlag_age = age_value,
    exposure_var_wlag_age = "spores_mean_wlag_age",
    n_weeks_wlag_age = nrow(dat_weekly),
    spores_or_wlag_age = NA_real_,
    spores_or_low_wlag_age = NA_real_,
    spores_or_high_wlag_age = NA_real_,
    spores_p_wlag_age = NA_real_
  )

  if (nrow(dat_weekly) >= min_weeks) {
    dat_week_logit <- dat_weekly %>% filter(total_claims_wlag_age > 0)
    if (nrow(dat_week_logit) >= min_weeks) {
      weather_week_cols <- names(dat_week_logit)[str_starts(names(dat_week_logit), "mean_")]

      f_week_logit <- as.formula(
        paste0(
          "cbind(acute_claims_wlag_age, nonacute_claims_wlag_age) ~ spores_mean_wlag_age",
          if (length(weather_week_cols) > 0) paste0(" + ", paste(weather_week_cols, collapse = " + ")) else ""
        )
      )

      fit_week_logit <- tryCatch(glm(f_week_logit, data = dat_week_logit, family = binomial()),
                                 error = function(e) NULL)
      if (!is.null(fit_week_logit)) {
        eff <- extract_effect_wlag_age(fit_week_logit, "spores_mean_wlag_age", exp_scale = TRUE)
        weekly_logit_row$spores_or_wlag_age      <- eff$estimate
        weekly_logit_row$spores_or_low_wlag_age  <- eff$ci_low
        weekly_logit_row$spores_or_high_wlag_age <- eff$ci_high
        weekly_logit_row$spores_p_wlag_age       <- eff$p
      }
    }
  }

  # A3) Weekly NB: total counts
  weekly_reg_row <- tibble(
    year_wlag_age = year_value,
    season_wlag_age = season_value,
    age_group_wlag_age = age_value,
    exposure_var_wlag_age = "spores_mean_wlag_age",
    n_weeks_wlag_age = nrow(dat_weekly),
    spores_rr_wlag_age = NA_real_,
    spores_rr_low_wlag_age = NA_real_,
    spores_rr_high_wlag_age = NA_real_,
    spores_p_wlag_age = NA_real_
  )

  if (nrow(dat_weekly) >= min_weeks) {
    weather_week_cols <- names(dat_weekly)[str_starts(names(dat_weekly), "mean_")]

    f_week_nb <- as.formula(
      paste0(
        "total_claims_wlag_age ~ spores_mean_wlag_age",
        if (length(weather_week_cols) > 0) paste0(" + ", paste(weather_week_cols, collapse = " + ")) else ""
      )
    )

    fit_week_nb <- tryCatch(MASS::glm.nb(f_week_nb, data = dat_weekly),
                            error = function(e) NULL)
    if (!is.null(fit_week_nb)) {
      eff <- extract_effect_wlag_age(fit_week_nb, "spores_mean_wlag_age", exp_scale = TRUE)
      weekly_reg_row$spores_rr_wlag_age      <- eff$estimate
      weekly_reg_row$spores_rr_low_wlag_age  <- eff$ci_low
      weekly_reg_row$spores_rr_high_wlag_age <- eff$ci_high
      weekly_reg_row$spores_p_wlag_age       <- eff$p
    }
  }

  # -------------------------------------------------------
  # B) DAILY models (lags / moving averages)
  # -------------------------------------------------------
  exposure_vars <- c(
    "spores_10k_lag0_wlag_age",
    "spores_10k_lag1_wlag_age",
    "spores_10k_lag2_wlag_age",
    "spores_10k_lag3_wlag_age",
    "spores_10k_ma3_wlag_age",
    "spores_10k_ma7_wlag_age"
  )

  # Daily logit: claim-days only
  dat_daily_logit <- dat0 %>%
    filter(total_claims_wlag_age > 0) %>%
    filter(if_all(all_of(weather_covars_use), ~ !is.na(.x)))

  logit_daily_tbl <- map_dfr(exposure_vars, function(ex_var) {
    if (!ex_var %in% names(dat_daily_logit)) return(NULL)

    d <- dat_daily_logit %>% filter(!is.na(.data[[ex_var]]))

    row <- tibble(
      year_wlag_age = year_value,
      season_wlag_age = season_value,
      age_group_wlag_age = age_value,
      exposure_var_wlag_age = ex_var,
      n_days_wlag_age = nrow(d),
      spores_or_wlag_age = NA_real_,
      spores_or_low_wlag_age = NA_real_,
      spores_or_high_wlag_age = NA_real_,
      spores_p_wlag_age = NA_real_
    )

    if (nrow(d) < min_days) return(row)

    f <- as.formula(
      paste0(
        "cbind(acute_claims_wlag_age, nonacute_claims_wlag_age) ~ ",
        ex_var, " + ",
        paste(weather_covars_use, collapse = " + "),
        " + ns(doy_wlag_age, df = 6) + factor(dow_wlag_age)"
      )
    )

    fit <- tryCatch(glm(f, data = d, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(row)

    eff <- extract_effect_wlag_age(fit, ex_var, exp_scale = TRUE)
    row$spores_or_wlag_age      <- eff$estimate
    row$spores_or_low_wlag_age  <- eff$ci_low
    row$spores_or_high_wlag_age <- eff$ci_high
    row$spores_p_wlag_age       <- eff$p
    row
  })

  # Daily NB: total daily counts
  dat_daily_nb <- dat0 %>%
    filter(if_all(all_of(weather_covars_use), ~ !is.na(.x)))

  reg_daily_tbl <- map_dfr(exposure_vars, function(ex_var) {
    if (!ex_var %in% names(dat_daily_nb)) return(NULL)

    d <- dat_daily_nb %>% filter(!is.na(.data[[ex_var]]))

    row <- tibble(
      year_wlag_age = year_value,
      season_wlag_age = season_value,
      age_group_wlag_age = age_value,
      exposure_var_wlag_age = ex_var,
      n_days_wlag_age = nrow(d),
      spores_rr_wlag_age = NA_real_,
      spores_rr_low_wlag_age = NA_real_,
      spores_rr_high_wlag_age = NA_real_,
      spores_p_wlag_age = NA_real_
    )

    if (nrow(d) < min_days) return(row)

    f <- as.formula(
      paste0(
        "total_claims_wlag_age ~ ",
        ex_var, " + ",
        paste(weather_covars_use, collapse = " + "),
        " + ns(doy_wlag_age, df = 6) + factor(dow_wlag_age)"
      )
    )

    fit <- tryCatch(MASS::glm.nb(f, data = d), error = function(e) NULL)
    if (is.null(fit)) return(row)

    eff <- extract_effect_wlag_age(fit, ex_var, exp_scale = TRUE)
    row$spores_rr_wlag_age      <- eff$estimate
    row$spores_rr_low_wlag_age  <- eff$ci_low
    row$spores_rr_high_wlag_age <- eff$ci_high
    row$spores_p_wlag_age       <- eff$p
    row
  })

  list(
    cor_wlag_age          = cor_row,
    logit_weekly_wlag_age = weekly_logit_row,
    reg_weekly_wlag_age   = weekly_reg_row,
    logit_daily_wlag_age  = logit_daily_tbl,
    reg_daily_wlag_age    = reg_daily_tbl
  )
}

season_grid_wlag_age <- tibble(
  season_wlag_age = c("All seasons", "winter", "spring", "summer", "fall"),
  min_weeks = c(min_weeks_overall_wlag_age, rep(min_weeks_season_wlag_age, 4)),
  min_days  = c(min_days_overall_wlag_age,  rep(min_days_season_wlag_age, 4))
)

# Add age groups to the analysis grid
run_grid_wlag_age <- crossing(
  year_wlag_age = years_to_run_wlag_age,
  age_group_wlag_age = age_labels_wlag_age,
  season_grid_wlag_age
)

cat(sprintf("\n=== Running %d age-stratified analyses ===\n", nrow(run_grid_wlag_age)))

results_list_wlag_age <- pmap(
  run_grid_wlag_age,
  function(year_wlag_age, age_group_wlag_age, season_wlag_age, min_weeks, min_days) {
    run_models_one_year_season_age_wlag_age(
      year_value = year_wlag_age,
      season_value = season_wlag_age,
      age_value = age_group_wlag_age,
      dat = daily_merged_wlag_age,
      weather_covars_use = weather_covars_use_wlag_age,
      min_weeks = min_weeks,
      min_days  = min_days
    )
  }
)

cor_table_wlag_age          <- bind_rows(map(results_list_wlag_age, "cor_wlag_age"))
logit_weekly_table_wlag_age <- bind_rows(map(results_list_wlag_age, "logit_weekly_wlag_age"))
reg_weekly_table_wlag_age   <- bind_rows(map(results_list_wlag_age, "reg_weekly_wlag_age"))
logit_daily_table_wlag_age  <- bind_rows(map(results_list_wlag_age, "logit_daily_wlag_age"))
reg_daily_table_wlag_age    <- bind_rows(map(results_list_wlag_age, "reg_daily_wlag_age"))

write.csv(cor_table_wlag_age,
          file.path(out_dir_wlag_age, "SJ_region_spores_asthma_correlations_by_year_season_age_wlag_age.csv"),
          row.names = FALSE)

write.csv(logit_daily_table_wlag_age,
          file.path(out_dir_wlag_age, "SJ_region_spores_asthma_logit_DAILY_lags_by_year_season_age_wlag_age.csv"),
          row.names = FALSE)

write.csv(reg_daily_table_wlag_age,
          file.path(out_dir_wlag_age, "SJ_region_spores_asthma_regression_DAILY_lags_by_year_season_age_wlag_age.csv"),
          row.names = FALSE)

write.csv(logit_weekly_table_wlag_age,
          file.path(out_dir_wlag_age, "SJ_region_spores_asthma_logit_WEEKLY_by_year_season_age_wlag_age.csv"),
          row.names = FALSE)

write.csv(reg_weekly_table_wlag_age,
          file.path(out_dir_wlag_age, "SJ_region_spores_asthma_regression_WEEKLY_by_year_season_age_wlag_age.csv"),
          row.names = FALSE)

print(cor_table_wlag_age)
print(logit_weekly_table_wlag_age)
print(reg_weekly_table_wlag_age)
print(logit_daily_table_wlag_age)
print(reg_daily_table_wlag_age)

cat("\nDone. Age-stratified outputs written to:", out_dir_wlag_age, "\n")
cat("Output files include age_group_wlag_age column for stratification.\n")
