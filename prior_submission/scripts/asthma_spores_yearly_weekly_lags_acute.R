############################################################
## Spores ↔ Asthma models: weekly and daily lag/MA analyses (acute)
## (San Juan health region, 2018–2019)
##
## Models (overall + by season):
##   - Weekly: Spearman + partial Spearman (controls weekly mean weather)
##   - Weekly: Binomial logit (acute fraction) + NB regression (total counts)
##   - Daily : Binomial logit + NB regression using lags / moving averages
##
## Exposure QC (weekly rule; applied separately to spores & pollen):
##   - If a week has <= 3 missing days: impute missing days with that week’s mean
##   - If a week has  > 3 missing days OR all values missing: exclude the week
##     (set the entire week to NA so it drops from models)
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(rnoaa)
  library(ppcor)
  library(MASS)
  library(splines)
  library(zoo)
})

# 1) SETTINGS ------

claims_path_wlag_acute  <- "data/2022_asma_From_date_ver 3.xlsx"
claims_sheet_wlag_acute <- 1

# One or more files; will be row-bound
spore_files_sj_wlag_acute <- c("data/Datos de esporas San Juan Jan 2018-2022.xlsx")

years_to_run_wlag_acute <- 2018:2019

clinical_settings_wlag_acute <- c("ED", "Inpatient", "Urgent", "Outpatient")

# Scaling for spores in models (per 10,000 spores)
spore_scale_divisor_wlag_acute <- 10000

# Weekly imputation/exclusion rule
max_missing_days_week_wlag_acute <- 3

# Minimum thresholds
min_weeks_overall_wlag_acute <- 10
min_weeks_season_wlag_acute  <- 6
min_days_overall_wlag_acute  <- 60
min_days_season_wlag_acute   <- 40

# GHCND station ID for SJU area (edit if you use a different station)
# NOTE: meteo_pull_monitors expects a GHCND ID (not "SJU")
sju_station_id_wlag_acute <- "RQW00011641"
if (toupper(sju_station_id_wlag_acute) == "SJU") sju_station_id_wlag_acute <- "RQW00011641"

# Weather covariates to control for (from pull_weather_sju_wlag_acute output)
weather_covars_wlag_acute <- c("tmean_c_wlag_acute", "rh_mean_wlag_acute", "wind_ms_wlag_acute", "prcp_mm_wlag_acute", "dewpoint_c_wlag_acute")

# Output directory
out_dir_wlag_acute <- "statistical_results"
dir.create(out_dir_wlag_acute, showWarnings = FALSE, recursive = TRUE)

# 2) HELPERS ------

parse_date_any_wlag_acute <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30")) # Excel numeric dates
  x <- as.character(x)

  d <- suppressWarnings(ymd(x))
  d <- coalesce(d, suppressWarnings(mdy(x)))
  d <- coalesce(d, suppressWarnings(dmy(x)))
  coalesce(d, suppressWarnings(as.Date(x)))
}

safe_mean_wlag_acute <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

assign_season_wlag_acute <- function(date) {
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

read_claims_wlag_acute <- function(path, sheet) {
  read_excel(path, sheet = sheet) %>%
    clean_names() %>%
    mutate(from_date = parse_date_any_wlag_acute(from_date))
}

care_setting_fun_wlag_acute <- function(hcg_setting, hcg_mr_line_desc) {
  s <- tolower(as.character(hcg_setting))
  l <- tolower(as.character(hcg_mr_line_desc))

  if (grepl("ed|emergency", s) || grepl("\\bed\\b|\\ber\\b", l)) return("ED")
  if (grepl("inpatient|hospital", s) || grepl("inpatient|hospital", l)) return("Inpatient")
  if (grepl("urgent", s) || grepl("urgent", l)) return("Urgent")
  if (grepl("office|clinic|outpatient", s) || grepl("office|clinic", l)) return("Outpatient")
  if (grepl("dme|supplies", l)) return("DME")
  "Other"
}

ensure_care_setting_wlag_acute <- function(df) {
  if ("care_setting_wlag_acute" %in% names(df)) return(df)

  if ("care_setting" %in% names(df)) {
    return(df %>% mutate(care_setting_wlag_acute = care_setting))
  }

  needed <- c("hcg_setting", "hcg_mr_line_desc")
  if (!all(needed %in% names(df))) {
    stop("No care_setting/care_setting_wlag_acute and missing hcg_setting/hcg_mr_line_desc to create it.")
  }

  df %>%
    mutate(care_setting_wlag_acute = mapply(care_setting_fun_wlag_acute, hcg_setting, hcg_mr_line_desc))
}

ensure_asthma_primary_wlag_acute <- function(df) {
  if ("asthma_primary_wlag_acute" %in% names(df)) return(df)

  if ("asthma_primary" %in% names(df)) {
    return(df %>% mutate(asthma_primary_wlag_acute = asthma_primary))
  }

  if (!"icd_diag_01_primary" %in% names(df)) {
    stop("Missing asthma_primary and icd_diag_01_primary; create asthma_primary before running.")
  }

  df %>%
    mutate(
      icd_diag_01_primary = str_to_upper(str_trim(as.character(icd_diag_01_primary))),
      asthma_primary_wlag_acute = if_else(str_detect(icd_diag_01_primary, "^J45|^J46"), "asthma", "other")
    )
}

# ---- Spores/pollen helpers ----

read_spores_sj_onefile_wlag_acute <- function(path) {
  x <- read_excel(path) %>% clean_names()

  date_col <- if ("date" %in% names(x)) "date" else names(x)[1]

  # spores
  spore_candidates <- c("spores_total", "total_spores", "esporas_de_hongos")
  spore_col <- spore_candidates[spore_candidates %in% names(x)][1]
  if (is.na(spore_col)) spore_col <- grep("esporas|spore", names(x), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(spore_col)) stop("Could not identify a spores column in: ", path)

  # pollen (optional)
  pollen_candidates <- c("pollen_total", "polen_san_juan", "polen")
  pollen_col <- pollen_candidates[pollen_candidates %in% names(x)][1]
  if (is.na(pollen_col)) pollen_col <- grep("polen|pollen", names(x), ignore.case = TRUE, value = TRUE)[1]

  x %>%
    transmute(
      date = parse_date_any_wlag_acute(.data[[date_col]]),
      spores_total_raw = suppressWarnings(as.numeric(.data[[spore_col]])),
      pollen_total_raw = if (!is.na(pollen_col)) suppressWarnings(as.numeric(.data[[pollen_col]])) else NA_real_
    ) %>%
    filter(!is.na(date))
}

weekly_qc_impute_or_exclude_wlag_acute <- function(df, value_col, out_col, max_missing = 3) {
  value_sym <- rlang::ensym(value_col)

  df %>%
    mutate(week_start_tmp = floor_date(date, unit = "week", week_start = 1)) %>%
    group_by(week_start_tmp) %>%
    mutate(
      n_missing_week = sum(is.na(!!value_sym)),
      week_mean      = safe_mean_wlag_acute(!!value_sym),
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
      !!paste0(out_col, "_imputed_wlag_acute") := imputed_day,
      !!paste0(out_col, "_week_excluded_wlag_acute") := week_excluded,
      !!paste0(out_col, "_n_missing_week_wlag_acute") := n_missing_week,
      !!paste0(out_col, "_week_mean_wlag_acute") := week_mean
    ) %>%
    dplyr::select(-week_start_tmp, -n_missing_week, -week_mean, -week_ok, -week_excluded, -imputed_day, -value_clean)
}

# ---- Weather helper (GHCND via rnoaa) ----
pull_weather_sju_wlag_acute <- function(year_wlag_acute, station_id_wlag_acute) {
  vars <- c("PRCP", "TMAX", "TMIN", "AWND", "ADPT", "RHAV")

  w <- meteo_pull_monitors(
    monitors   = station_id_wlag_acute,
    var        = vars,
    date_min   = sprintf("%d-01-01", year_wlag_acute),
    date_max   = sprintf("%d-12-31", year_wlag_acute),
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
      prcp_mm_wlag_acute    = prcp / 10,                        # tenths mm -> mm
      tmax_c_wlag_acute     = tmax / 10,                        # tenths °C -> °C
      tmin_c_wlag_acute     = tmin / 10,
      tmean_c_wlag_acute    = (tmax_c_wlag_acute + tmin_c_wlag_acute) / 2,
      wind_ms_wlag_acute    = awnd / 10,                        # tenths m/s -> m/s
      rh_mean_wlag_acute    = rhav,                             # %
      dewpoint_c_wlag_acute = adpt / 10                         # tenths °C -> °C
    )
}

# ---- Extract exponentiated effect (OR/RR) + Wald CI ----
extract_effect_wlag_acute <- function(model, term, exp_scale = TRUE) {
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

# 3) BUILD DAILY ASTHMA OUTCOMES (San Juan health region) ------

claims_df_wlag_acute <- read_claims_wlag_acute(claims_path_wlag_acute, claims_sheet_wlag_acute) %>%
  ensure_care_setting_wlag_acute() %>%
  ensure_asthma_primary_wlag_acute() %>%
  mutate(
    municipality_desc = toupper(as.character(municipality_desc)),
    health_region_wlag_acute = case_when(
      municipality_desc %in% c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA") ~ "San Juan",
      municipality_desc %in% c("CAGUAS", "GURABO", "JUNCOS", "SAN LORENZO") ~ "Caguas",
      TRUE ~ "Other"
    ),
    health_region_wlag_acute = factor(health_region_wlag_acute, levels = c("San Juan", "Caguas")),
    year_wlag_acute = year(from_date)
  ) %>%
  filter(year_wlag_acute %in% years_to_run_wlag_acute)

claims_sj_region_wlag_acute <- claims_df_wlag_acute %>%
  filter(
    health_region_wlag_acute == "San Juan",
    asthma_primary_wlag_acute == "asthma",
    care_setting_wlag_acute %in% clinical_settings_wlag_acute
  ) %>%
  mutate(acute_flag_wlag_acute = as.integer(care_setting_wlag_acute %in% c("ED", "Inpatient")))

calendar_all_wlag_acute <- tibble(
  date = seq.Date(
    as.Date(sprintf("%d-01-01", min(years_to_run_wlag_acute))),
    as.Date(sprintf("%d-12-31", max(years_to_run_wlag_acute))),
    by = "day"
  )
) %>%
  mutate(
    year_wlag_acute = year(date),
    doy_wlag_acute  = yday(date),
    dow_wlag_acute  = wday(date, label = TRUE, week_start = 1),
    week_start_wlag_acute = floor_date(date, unit = "week", week_start = 1),
    season_wlag_acute = assign_season_wlag_acute(date)
  ) %>%
  filter(year_wlag_acute %in% years_to_run_wlag_acute)

daily_asthma_sj_wlag_acute <- claims_sj_region_wlag_acute %>%
  group_by(date = from_date) %>%
  summarise(
    total_claims_wlag_acute    = n(),
    acute_claims_wlag_acute    = sum(acute_flag_wlag_acute, na.rm = TRUE),
    nonacute_claims_wlag_acute = total_claims_wlag_acute - acute_claims_wlag_acute,
    ed_claims_wlag_acute       = sum(care_setting_wlag_acute == "ED", na.rm = TRUE),
    inpt_claims_wlag_acute     = sum(care_setting_wlag_acute == "Inpatient", na.rm = TRUE),
    outpt_claims_wlag_acute    = sum(care_setting_wlag_acute == "Outpatient", na.rm = TRUE),
    urgent_claims_wlag_acute   = sum(care_setting_wlag_acute == "Urgent", na.rm = TRUE),
    .groups = "drop"
  )

daily_sj_wlag_acute <- calendar_all_wlag_acute %>%
  left_join(daily_asthma_sj_wlag_acute, by = "date") %>%
  mutate(
    across(
      c(total_claims_wlag_acute, acute_claims_wlag_acute, nonacute_claims_wlag_acute,
        ed_claims_wlag_acute, inpt_claims_wlag_acute, outpt_claims_wlag_acute, urgent_claims_wlag_acute),
      ~ replace_na(.x, 0L)
    )
  )

# 4) SPORES/POLLEN (San Juan station) + WEATHER -----

spores_sj_raw_wlag_acute <- map_dfr(spore_files_sj_wlag_acute, read_spores_sj_onefile_wlag_acute) %>%
  group_by(date) %>%
  summarise(
    spores_total_raw = safe_mean_wlag_acute(spores_total_raw),
    pollen_total_raw = safe_mean_wlag_acute(pollen_total_raw),
    .groups = "drop"
  ) %>%
  arrange(date)

calendar_exp_wlag_acute <- tibble(date = calendar_all_wlag_acute$date)

spores_sj_daily_wlag_acute <- calendar_exp_wlag_acute %>%
  left_join(spores_sj_raw_wlag_acute, by = "date") %>%
  arrange(date)

# Apply weekly QC separately (creates spores_total_wlag_acute and pollen_total_wlag_acute)
spores_sj_daily_wlag_acute <- weekly_qc_impute_or_exclude_wlag_acute(
  spores_sj_daily_wlag_acute,
  value_col = spores_total_raw,
  out_col   = "spores_total_wlag_acute",
  max_missing = max_missing_days_week_wlag_acute
)

spores_sj_daily_wlag_acute <- weekly_qc_impute_or_exclude_wlag_acute(
  spores_sj_daily_wlag_acute,
  value_col = pollen_total_raw,
  out_col   = "pollen_total_wlag_acute",
  max_missing = max_missing_days_week_wlag_acute
)

weather_sj_wlag_acute <- map_dfr(years_to_run_wlag_acute, ~ pull_weather_sju_wlag_acute(.x, sju_station_id_wlag_acute))

daily_merged_wlag_acute <- daily_sj_wlag_acute %>%
  left_join(
    spores_sj_daily_wlag_acute %>%
      dplyr::select(date,
             spores_total_wlag_acute, pollen_total_wlag_acute,
             spores_total_wlag_acute_imputed_wlag_acute, spores_total_wlag_acute_week_excluded_wlag_acute, spores_total_wlag_acute_n_missing_week_wlag_acute,
             pollen_total_wlag_acute_imputed_wlag_acute, pollen_total_wlag_acute_week_excluded_wlag_acute, pollen_total_wlag_acute_n_missing_week_wlag_acute),
    by = "date"
  ) %>%
  left_join(weather_sj_wlag_acute, by = "date") %>%
  arrange(date) %>%
  group_by(year_wlag_acute) %>%
  mutate(
    spores_10k_wlag_acute = spores_total_wlag_acute / spore_scale_divisor_wlag_acute,

    spores_10k_lag0_wlag_acute = spores_10k_wlag_acute,
    spores_10k_lag1_wlag_acute = lag(spores_10k_wlag_acute, 1),
    spores_10k_lag2_wlag_acute = lag(spores_10k_wlag_acute, 2),
    spores_10k_lag3_wlag_acute = lag(spores_10k_wlag_acute, 3),

    spores_10k_ma3_wlag_acute = zoo::rollapply(spores_10k_wlag_acute, width = 3, FUN = mean,
                                         align = "right", fill = NA_real_, partial = FALSE),
    spores_10k_ma7_wlag_acute = zoo::rollapply(spores_10k_wlag_acute, width = 7, FUN = mean,
                                         align = "right", fill = NA_real_, partial = FALSE)
  ) %>%
  ungroup()

weather_covars_use_wlag_acute <- intersect(weather_covars_wlag_acute, names(daily_merged_wlag_acute))
if (length(weather_covars_use_wlag_acute) == 0) stop("No requested weather covariates found after weather pull.")

# 5) DIAGNOSTICS -----

cat("\n=== Exposure QC diagnostics (daily) ===\n")
daily_merged_wlag_acute %>%
  group_by(year_wlag_acute) %>%
  summarise(
    n_days = n(),
    n_spores_nonmiss = sum(!is.na(spores_10k_wlag_acute)),
    n_weather_complete = sum(if_all(all_of(weather_covars_use_wlag_acute), ~ !is.na(.x))),
    n_weeks_spores_excluded = n_distinct(week_start_wlag_acute[spores_total_wlag_acute_week_excluded_wlag_acute %in% TRUE]),
    .groups = "drop"
  ) %>%
  print()

# 6) MODELS (overall + by season) -----

run_models_one_year_season_wlag_acute <- function(
  year_value,
  season_value,
  dat,
  weather_covars_use,
  min_weeks,
  min_days
) {

  dat0 <- dat %>% filter(year_wlag_acute == year_value)
  if (season_value != "All seasons") dat0 <- dat0 %>% filter(season_wlag_acute == season_value)

  # -------------------------
  # A) WEEKLY summary
  # -------------------------
  dat_weekly <- dat0 %>%
    filter(!is.na(spores_10k_wlag_acute)) %>%
    filter(if_all(all_of(weather_covars_use), ~ !is.na(.x))) %>%
    group_by(week_start_wlag_acute) %>%
    summarise(
      total_claims_wlag_acute    = sum(total_claims_wlag_acute),
      acute_claims_wlag_acute    = sum(acute_claims_wlag_acute),
      nonacute_claims_wlag_acute = sum(nonacute_claims_wlag_acute),
      spores_mean_wlag_acute     = mean(spores_10k_wlag_acute, na.rm = TRUE),
      across(all_of(weather_covars_use), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
      .groups = "drop"
    )

  # A1) correlations
  cor_row <- tibble(
    year_wlag_acute = year_value,
    season_wlag_acute = season_value,
    n_weeks_wlag_acute = nrow(dat_weekly),
    spearman_rho_wlag_acute = NA_real_,
    spearman_p_wlag_acute   = NA_real_,
    partial_rho_wlag_acute  = NA_real_,
    partial_p_wlag_acute    = NA_real_
  )

  if (nrow(dat_weekly) >= min_weeks) {
    ct <- suppressWarnings(cor.test(dat_weekly$total_claims_wlag_acute, dat_weekly$spores_mean_wlag_acute, method = "spearman"))
    cor_row$spearman_rho_wlag_acute <- unname(ct$estimate)
    cor_row$spearman_p_wlag_acute   <- ct$p.value

    zmat <- dat_weekly %>% dplyr::select(starts_with("mean_")) %>% as.data.frame()
    ok <- complete.cases(dat_weekly$total_claims_wlag_acute, dat_weekly$spores_mean_wlag_acute, zmat)

    if (sum(ok) >= min_weeks && ncol(zmat) > 0) {
      pt <- tryCatch(
        suppressWarnings(ppcor::pcor.test(
          x = dat_weekly$spores_mean_wlag_acute[ok],
          y = dat_weekly$total_claims_wlag_acute[ok],
          z = zmat[ok, , drop = FALSE],
          method = "spearman"
        )),
        error = function(e) NULL
      )
      if (!is.null(pt)) {
        cor_row$partial_rho_wlag_acute <- pt$estimate
        cor_row$partial_p_wlag_acute   <- pt$p.value
      }
    }
  }

  # A2) weekly logit (acute fraction)
  weekly_logit_row <- tibble(
    year_wlag_acute = year_value,
    season_wlag_acute = season_value,
    exposure_var_wlag_acute = "spores_mean_wlag_acute",
    n_weeks_wlag_acute = nrow(dat_weekly),
    spores_or_wlag_acute = NA_real_,
    spores_or_low_wlag_acute = NA_real_,
    spores_or_high_wlag_acute = NA_real_,
    spores_p_wlag_acute = NA_real_
  )

  if (nrow(dat_weekly) >= min_weeks) {
    dat_week_logit <- dat_weekly %>% filter(total_claims_wlag_acute > 0)
    if (nrow(dat_week_logit) >= min_weeks) {
      weather_week_cols <- names(dat_week_logit)[str_starts(names(dat_week_logit), "mean_")]

      f_week_logit <- as.formula(
        paste0(
          "cbind(acute_claims_wlag_acute, nonacute_claims_wlag_acute) ~ spores_mean_wlag_acute",
          if (length(weather_week_cols) > 0) paste0(" + ", paste(weather_week_cols, collapse = " + ")) else ""
        )
      )

      fit_week_logit <- tryCatch(glm(f_week_logit, data = dat_week_logit, family = binomial()), error = function(e) NULL)
      if (!is.null(fit_week_logit)) {
        eff <- extract_effect_wlag_acute(fit_week_logit, "spores_mean_wlag_acute", exp_scale = TRUE)
        weekly_logit_row$spores_or_wlag_acute      <- eff$estimate
        weekly_logit_row$spores_or_low_wlag_acute  <- eff$ci_low
        weekly_logit_row$spores_or_high_wlag_acute <- eff$ci_high
        weekly_logit_row$spores_p_wlag_acute       <- eff$p
      }
    }
  }

  # A3) weekly NB (total counts)
  weekly_reg_row <- tibble(
    year_wlag_acute = year_value,
    season_wlag_acute = season_value,
    exposure_var_wlag_acute = "spores_mean_wlag_acute",
    n_weeks_wlag_acute = nrow(dat_weekly),
    spores_rr_wlag_acute = NA_real_,
    spores_rr_low_wlag_acute = NA_real_,
    spores_rr_high_wlag_acute = NA_real_,
    spores_p_wlag_acute = NA_real_
  )

  if (nrow(dat_weekly) >= min_weeks) {
    weather_week_cols <- names(dat_weekly)[str_starts(names(dat_weekly), "mean_")]

    f_week_nb <- as.formula(
      paste0(
        "total_claims_wlag_acute ~ spores_mean_wlag_acute",
        if (length(weather_week_cols) > 0) paste0(" + ", paste(weather_week_cols, collapse = " + ")) else ""
      )
    )

    fit_week_nb <- tryCatch(MASS::glm.nb(f_week_nb, data = dat_weekly), error = function(e) NULL)
    if (!is.null(fit_week_nb)) {
      eff <- extract_effect_wlag_acute(fit_week_nb, "spores_mean_wlag_acute", exp_scale = TRUE)
      weekly_reg_row$spores_rr_wlag_acute      <- eff$estimate
      weekly_reg_row$spores_rr_low_wlag_acute  <- eff$ci_low
      weekly_reg_row$spores_rr_high_wlag_acute <- eff$ci_high
      weekly_reg_row$spores_p_wlag_acute       <- eff$p
    }
  }

  # B) DAILY models (lags/MA) ------
  exposure_vars <- c(
    "spores_10k_lag0_wlag_acute",
    "spores_10k_lag1_wlag_acute",
    "spores_10k_lag2_wlag_acute",
    "spores_10k_lag3_wlag_acute",
    "spores_10k_ma3_wlag_acute",
    "spores_10k_ma7_wlag_acute"
  )

  # Daily logit: claim-days only
  dat_daily_logit <- dat0 %>%
    filter(total_claims_wlag_acute > 0) %>%
    filter(if_all(all_of(weather_covars_use), ~ !is.na(.x)))

  logit_daily_tbl <- map_dfr(exposure_vars, function(ex_var) {
    if (!ex_var %in% names(dat_daily_logit)) return(NULL)

    d <- dat_daily_logit %>% filter(!is.na(.data[[ex_var]]))

    row <- tibble(
      year_wlag_acute = year_value,
      season_wlag_acute = season_value,
      exposure_var_wlag_acute = ex_var,
      n_days_wlag_acute = nrow(d),
      spores_or_wlag_acute = NA_real_,
      spores_or_low_wlag_acute = NA_real_,
      spores_or_high_wlag_acute = NA_real_,
      spores_p_wlag_acute = NA_real_
    )

    if (nrow(d) < min_days) return(row)

    f <- as.formula(
      paste0(
        "cbind(acute_claims_wlag_acute, nonacute_claims_wlag_acute) ~ ",
        ex_var, " + ",
        paste(weather_covars_use, collapse = " + "),
        " + ns(doy_wlag_acute, df = 6) + factor(dow_wlag_acute)"
      )
    )

    fit <- tryCatch(glm(f, data = d, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(row)

    eff <- extract_effect_wlag_acute(fit, ex_var, exp_scale = TRUE)
    row$spores_or_wlag_acute      <- eff$estimate
    row$spores_or_low_wlag_acute  <- eff$ci_low
    row$spores_or_high_wlag_acute <- eff$ci_high
    row$spores_p_wlag_acute       <- eff$p
    row
  })

  # Daily NB: total counts
  dat_daily_nb <- dat0 %>%
    filter(if_all(all_of(weather_covars_use), ~ !is.na(.x)))

  reg_daily_tbl <- map_dfr(exposure_vars, function(ex_var) {
    if (!ex_var %in% names(dat_daily_nb)) return(NULL)

    d <- dat_daily_nb %>% filter(!is.na(.data[[ex_var]]))

    row <- tibble(
      year_wlag_acute = year_value,
      season_wlag_acute = season_value,
      exposure_var_wlag_acute = ex_var,
      n_days_wlag_acute = nrow(d),
      spores_rr_wlag_acute = NA_real_,
      spores_rr_low_wlag_acute = NA_real_,
      spores_rr_high_wlag_acute = NA_real_,
      spores_p_wlag_acute = NA_real_
    )

    if (nrow(d) < min_days) return(row)

    f <- as.formula(
      paste0(
        "total_claims_wlag_acute ~ ",
        ex_var, " + ",
        paste(weather_covars_use, collapse = " + "),
        " + ns(doy_wlag_acute, df = 6) + factor(dow_wlag_acute)"
      )
    )

    fit <- tryCatch(MASS::glm.nb(f, data = d), error = function(e) NULL)
    if (is.null(fit)) return(row)

    eff <- extract_effect_wlag_acute(fit, ex_var, exp_scale = TRUE)
    row$spores_rr_wlag_acute      <- eff$estimate
    row$spores_rr_low_wlag_acute  <- eff$ci_low
    row$spores_rr_high_wlag_acute <- eff$ci_high
    row$spores_p_wlag_acute       <- eff$p
    row
  })

  list(
    cor_wlag_acute          = cor_row,
    logit_weekly_wlag_acute = weekly_logit_row,
    reg_weekly_wlag_acute   = weekly_reg_row,
    logit_daily_wlag_acute  = logit_daily_tbl,
    reg_daily_wlag_acute    = reg_daily_tbl
  )
}

season_grid_wlag_acute <- tibble(
  season_wlag_acute = c("All seasons", "winter", "spring", "summer", "fall"),
  min_weeks = c(min_weeks_overall_wlag_acute, rep(min_weeks_season_wlag_acute, 4)),
  min_days  = c(min_days_overall_wlag_acute,  rep(min_days_season_wlag_acute, 4))
)

run_grid_wlag_acute <- crossing(
  year_wlag_acute = years_to_run_wlag_acute,
  season_grid_wlag_acute
)

results_list_wlag_acute <- pmap(
  run_grid_wlag_acute,
  function(year_wlag_acute, season_wlag_acute, min_weeks, min_days) {
    run_models_one_year_season_wlag_acute(
      year_value = year_wlag_acute,
      season_value = season_wlag_acute,
      dat = daily_merged_wlag_acute,
      weather_covars_use = weather_covars_use_wlag_acute,
      min_weeks = min_weeks,
      min_days  = min_days
    )
  }
)

cor_table_wlag_acute          <- bind_rows(map(results_list_wlag_acute, "cor_wlag_acute"))
logit_weekly_table_wlag_acute <- bind_rows(map(results_list_wlag_acute, "logit_weekly_wlag_acute"))
reg_weekly_table_wlag_acute   <- bind_rows(map(results_list_wlag_acute, "reg_weekly_wlag_acute"))
logit_daily_table_wlag_acute  <- bind_rows(map(results_list_wlag_acute, "logit_daily_wlag_acute"))
reg_daily_table_wlag_acute    <- bind_rows(map(results_list_wlag_acute, "reg_daily_wlag_acute"))

write.csv(cor_table_wlag_acute,
          file.path(out_dir_wlag_acute, "SJ_region_spores_asthma_correlations_by_year_season_wlag_acute.csv"),
          row.names = FALSE)

write.csv(logit_daily_table_wlag_acute,
          file.path(out_dir_wlag_acute, "SJ_region_spores_asthma_logit_DAILY_lags_by_year_season_wlag_acute.csv"),
          row.names = FALSE)

write.csv(reg_daily_table_wlag_acute,
          file.path(out_dir_wlag_acute, "SJ_region_spores_asthma_regression_DAILY_lags_by_year_season_wlag_acute.csv"),
          row.names = FALSE)

write.csv(logit_weekly_table_wlag_acute,
          file.path(out_dir_wlag_acute, "SJ_region_spores_asthma_logit_WEEKLY_by_year_season_wlag_acute.csv"),
          row.names = FALSE)

write.csv(reg_weekly_table_wlag_acute,
          file.path(out_dir_wlag_acute, "SJ_region_spores_asthma_regression_WEEKLY_by_year_season_wlag_acute.csv"),
          row.names = FALSE)

print(cor_table_wlag_acute)
print(logit_weekly_table_wlag_acute)
print(reg_weekly_table_wlag_acute)
print(logit_daily_table_wlag_acute)
print(reg_daily_table_wlag_acute)

cat("\nDone. Outputs written to:", out_dir_wlag_acute, "\n")
