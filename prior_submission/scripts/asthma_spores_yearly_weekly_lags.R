############################################################
## Spores ↔ Asthma models: weekly and daily lag/MA analyses
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

claims_path_wlag  <- "data/2022_asma_From_date_ver 3.xlsx"
claims_sheet_wlag <- 1

# One or more files; will be row-bound
spore_files_sj_wlag <- c("data/Datos de esporas San Juan Jan 2018-2022.xlsx")

years_to_run_wlag <- 2018:2019

clinical_settings_wlag <- c("ED", "Inpatient", "Urgent", "Outpatient")

# Scaling for spores in models (per 10,000 spores)
spore_scale_divisor_wlag <- 10000

# Weekly imputation/exclusion rule
max_missing_days_week_wlag <- 3

# Minimum thresholds
min_weeks_overall_wlag <- 10
min_weeks_season_wlag  <- 6
min_days_overall_wlag  <- 60
min_days_season_wlag   <- 40

# GHCND station ID for SJU area (edit if you use a different station)
# NOTE: meteo_pull_monitors expects a GHCND ID (not "SJU")
sju_station_id_wlag <- "RQW00011641"
if (toupper(sju_station_id_wlag) == "SJU") sju_station_id_wlag <- "RQW00011641"

# Weather covariates to control for (from pull_weather_sju_wlag output)
weather_covars_wlag <- c("tmean_c_wlag", "rh_mean_wlag", "wind_ms_wlag", "prcp_mm_wlag", "dewpoint_c_wlag")

# Output directory
out_dir_wlag <- "statistical_results"
dir.create(out_dir_wlag, showWarnings = FALSE, recursive = TRUE)

# 2) HELPERS ------

parse_date_any_wlag <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30")) # Excel numeric dates
  x <- as.character(x)

  d <- suppressWarnings(ymd(x))
  d <- coalesce(d, suppressWarnings(mdy(x)))
  d <- coalesce(d, suppressWarnings(dmy(x)))
  coalesce(d, suppressWarnings(as.Date(x)))
}

safe_mean_wlag <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

assign_season_wlag <- function(date) {
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

read_claims_wlag <- function(path, sheet) {
  read_excel(path, sheet = sheet) %>%
    clean_names() %>%
    mutate(from_date = parse_date_any_wlag(from_date))
}

care_setting_fun_wlag <- function(hcg_setting, hcg_mr_line_desc) {
  s <- tolower(as.character(hcg_setting))
  l <- tolower(as.character(hcg_mr_line_desc))

  if (grepl("ed|emergency", s) || grepl("\\bed\\b|\\ber\\b", l)) return("ED")
  if (grepl("inpatient|hospital", s) || grepl("inpatient|hospital", l)) return("Inpatient")
  if (grepl("urgent", s) || grepl("urgent", l)) return("Urgent")
  if (grepl("office|clinic|outpatient", s) || grepl("office|clinic", l)) return("Outpatient")
  if (grepl("dme|supplies", l)) return("DME")
  "Other"
}

ensure_care_setting_wlag <- function(df) {
  if ("care_setting_wlag" %in% names(df)) return(df)

  if ("care_setting" %in% names(df)) {
    return(df %>% mutate(care_setting_wlag = care_setting))
  }

  needed <- c("hcg_setting", "hcg_mr_line_desc")
  if (!all(needed %in% names(df))) {
    stop("No care_setting/care_setting_wlag and missing hcg_setting/hcg_mr_line_desc to create it.")
  }

  df %>%
    mutate(care_setting_wlag = mapply(care_setting_fun_wlag, hcg_setting, hcg_mr_line_desc))
}

ensure_asthma_primary_wlag <- function(df) {
  if ("asthma_primary_wlag" %in% names(df)) return(df)

  if ("asthma_primary" %in% names(df)) {
    return(df %>% mutate(asthma_primary_wlag = asthma_primary))
  }

  if (!"icd_diag_01_primary" %in% names(df)) {
    stop("Missing asthma_primary and icd_diag_01_primary; create asthma_primary before running.")
  }

  df %>%
    mutate(
      icd_diag_01_primary = str_to_upper(str_trim(as.character(icd_diag_01_primary))),
      asthma_primary_wlag = if_else(str_detect(icd_diag_01_primary, "^J45|^J46"), "asthma", "other")
    )
}

# ---- Spores/pollen helpers ----

read_spores_sj_onefile_wlag <- function(path) {
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
      date = parse_date_any_wlag(.data[[date_col]]),
      spores_total_raw = suppressWarnings(as.numeric(.data[[spore_col]])),
      pollen_total_raw = if (!is.na(pollen_col)) suppressWarnings(as.numeric(.data[[pollen_col]])) else NA_real_
    ) %>%
    filter(!is.na(date))
}

weekly_qc_impute_or_exclude_wlag <- function(df, value_col, out_col, max_missing = 3) {
  value_sym <- rlang::ensym(value_col)

  df %>%
    mutate(week_start_tmp = floor_date(date, unit = "week", week_start = 1)) %>%
    group_by(week_start_tmp) %>%
    mutate(
      n_missing_week = sum(is.na(!!value_sym)),
      week_mean      = safe_mean_wlag(!!value_sym),
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
      !!paste0(out_col, "_imputed_wlag") := imputed_day,
      !!paste0(out_col, "_week_excluded_wlag") := week_excluded,
      !!paste0(out_col, "_n_missing_week_wlag") := n_missing_week,
      !!paste0(out_col, "_week_mean_wlag") := week_mean
    ) %>%
    dplyr::select(-week_start_tmp, -n_missing_week, -week_mean, -week_ok, -week_excluded, -imputed_day, -value_clean)
}

# ---- Weather helper (GHCND via rnoaa) ----
pull_weather_sju_wlag <- function(year_wlag, station_id_wlag) {
  vars <- c("PRCP", "TMAX", "TMIN", "AWND", "ADPT", "RHAV")

  w <- meteo_pull_monitors(
    monitors   = station_id_wlag,
    var        = vars,
    date_min   = sprintf("%d-01-01", year_wlag),
    date_max   = sprintf("%d-12-31", year_wlag),
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
      prcp_mm_wlag    = prcp / 10,                        # tenths mm -> mm
      tmax_c_wlag     = tmax / 10,                        # tenths °C -> °C
      tmin_c_wlag     = tmin / 10,
      tmean_c_wlag    = (tmax_c_wlag + tmin_c_wlag) / 2,
      wind_ms_wlag    = awnd / 10,                        # tenths m/s -> m/s
      rh_mean_wlag    = rhav,                             # %
      dewpoint_c_wlag = adpt / 10                         # tenths °C -> °C
    )
}

# ---- Extract exponentiated effect (OR/RR) + Wald CI ----
extract_effect_wlag <- function(model, term, exp_scale = TRUE) {
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

claims_df_wlag <- read_claims_wlag(claims_path_wlag, claims_sheet_wlag) %>%
  ensure_care_setting_wlag() %>%
  ensure_asthma_primary_wlag() %>%
  mutate(
    municipality_desc = toupper(as.character(municipality_desc)),
    health_region_wlag = case_when(
      municipality_desc %in% c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA") ~ "San Juan",
      municipality_desc %in% c("CAGUAS", "GURABO", "JUNCOS", "SAN LORENZO") ~ "Caguas",
      TRUE ~ "Other"
    ),
    health_region_wlag = factor(health_region_wlag, levels = c("San Juan", "Caguas")),
    year_wlag = year(from_date)
  ) %>%
  filter(year_wlag %in% years_to_run_wlag)

claims_sj_region_wlag <- claims_df_wlag %>%
  filter(
    health_region_wlag == "San Juan",
    asthma_primary_wlag == "asthma",
    care_setting_wlag %in% clinical_settings_wlag
  ) %>%
  mutate(acute_flag_wlag = as.integer(care_setting_wlag %in% c("ED", "Inpatient")))

calendar_all_wlag <- tibble(
  date = seq.Date(
    as.Date(sprintf("%d-01-01", min(years_to_run_wlag))),
    as.Date(sprintf("%d-12-31", max(years_to_run_wlag))),
    by = "day"
  )
) %>%
  mutate(
    year_wlag = year(date),
    doy_wlag  = yday(date),
    dow_wlag  = wday(date, label = TRUE, week_start = 1),
    week_start_wlag = floor_date(date, unit = "week", week_start = 1),
    season_wlag = assign_season_wlag(date)
  ) %>%
  filter(year_wlag %in% years_to_run_wlag)

daily_asthma_sj_wlag <- claims_sj_region_wlag %>%
  group_by(date = from_date) %>%
  summarise(
    total_claims_wlag    = n(),
    acute_claims_wlag    = sum(acute_flag_wlag, na.rm = TRUE),
    nonacute_claims_wlag = total_claims_wlag - acute_claims_wlag,
    ed_claims_wlag       = sum(care_setting_wlag == "ED", na.rm = TRUE),
    inpt_claims_wlag     = sum(care_setting_wlag == "Inpatient", na.rm = TRUE),
    outpt_claims_wlag    = sum(care_setting_wlag == "Outpatient", na.rm = TRUE),
    urgent_claims_wlag   = sum(care_setting_wlag == "Urgent", na.rm = TRUE),
    .groups = "drop"
  )

daily_sj_wlag <- calendar_all_wlag %>%
  left_join(daily_asthma_sj_wlag, by = "date") %>%
  mutate(
    across(
      c(total_claims_wlag, acute_claims_wlag, nonacute_claims_wlag,
        ed_claims_wlag, inpt_claims_wlag, outpt_claims_wlag, urgent_claims_wlag),
      ~ replace_na(.x, 0L)
    )
  )

# 4) SPORES/POLLEN (San Juan station) + WEATHER -----

spores_sj_raw_wlag <- map_dfr(spore_files_sj_wlag, read_spores_sj_onefile_wlag) %>%
  group_by(date) %>%
  summarise(
    spores_total_raw = safe_mean_wlag(spores_total_raw),
    pollen_total_raw = safe_mean_wlag(pollen_total_raw),
    .groups = "drop"
  ) %>%
  arrange(date)

calendar_exp_wlag <- tibble(date = calendar_all_wlag$date)

spores_sj_daily_wlag <- calendar_exp_wlag %>%
  left_join(spores_sj_raw_wlag, by = "date") %>%
  arrange(date)

# Apply weekly QC separately (creates spores_total_wlag and pollen_total_wlag)
spores_sj_daily_wlag <- weekly_qc_impute_or_exclude_wlag(
  spores_sj_daily_wlag,
  value_col = spores_total_raw,
  out_col   = "spores_total_wlag",
  max_missing = max_missing_days_week_wlag
)

spores_sj_daily_wlag <- weekly_qc_impute_or_exclude_wlag(
  spores_sj_daily_wlag,
  value_col = pollen_total_raw,
  out_col   = "pollen_total_wlag",
  max_missing = max_missing_days_week_wlag
)

weather_sj_wlag <- map_dfr(years_to_run_wlag, ~ pull_weather_sju_wlag(.x, sju_station_id_wlag))

daily_merged_wlag <- daily_sj_wlag %>%
  left_join(
    spores_sj_daily_wlag %>%
      dplyr::select(date,
             spores_total_wlag, pollen_total_wlag,
             spores_total_wlag_imputed_wlag, spores_total_wlag_week_excluded_wlag, spores_total_wlag_n_missing_week_wlag,
             pollen_total_wlag_imputed_wlag, pollen_total_wlag_week_excluded_wlag, pollen_total_wlag_n_missing_week_wlag),
    by = "date"
  ) %>%
  left_join(weather_sj_wlag, by = "date") %>%
  arrange(date) %>%
  group_by(year_wlag) %>%
  mutate(
    spores_10k_wlag = spores_total_wlag / spore_scale_divisor_wlag,

    spores_10k_lag0_wlag = spores_10k_wlag,
    spores_10k_lag1_wlag = lag(spores_10k_wlag, 1),
    spores_10k_lag2_wlag = lag(spores_10k_wlag, 2),
    spores_10k_lag3_wlag = lag(spores_10k_wlag, 3),

    spores_10k_ma3_wlag = zoo::rollapply(spores_10k_wlag, width = 3, FUN = mean,
                                         align = "right", fill = NA_real_, partial = FALSE),
    spores_10k_ma7_wlag = zoo::rollapply(spores_10k_wlag, width = 7, FUN = mean,
                                         align = "right", fill = NA_real_, partial = FALSE)
  ) %>%
  ungroup()

weather_covars_use_wlag <- intersect(weather_covars_wlag, names(daily_merged_wlag))
if (length(weather_covars_use_wlag) == 0) stop("No requested weather covariates found after weather pull.")

# 5) DIAGNOSTICS -----

cat("\n=== Exposure QC diagnostics (daily) ===\n")
daily_merged_wlag %>%
  group_by(year_wlag) %>%
  summarise(
    n_days = n(),
    n_spores_nonmiss = sum(!is.na(spores_10k_wlag)),
    n_weather_complete = sum(if_all(all_of(weather_covars_use_wlag), ~ !is.na(.x))),
    n_weeks_spores_excluded = n_distinct(week_start_wlag[spores_total_wlag_week_excluded_wlag %in% TRUE]),
    .groups = "drop"
  ) %>%
  print()

# 6) MODELS (overall + by season) -----

run_models_one_year_season_wlag <- function(
  year_value,
  season_value,
  dat,
  weather_covars_use,
  min_weeks,
  min_days
) {

  dat0 <- dat %>% filter(year_wlag == year_value)
  if (season_value != "All seasons") dat0 <- dat0 %>% filter(season_wlag == season_value)

  # -------------------------
  # A) WEEKLY summary
  # -------------------------
  dat_weekly <- dat0 %>%
    filter(!is.na(spores_10k_wlag)) %>%
    filter(if_all(all_of(weather_covars_use), ~ !is.na(.x))) %>%
    group_by(week_start_wlag) %>%
    summarise(
      total_claims_wlag    = sum(total_claims_wlag),
      acute_claims_wlag    = sum(acute_claims_wlag),
      nonacute_claims_wlag = sum(nonacute_claims_wlag),
      spores_mean_wlag     = mean(spores_10k_wlag, na.rm = TRUE),
      across(all_of(weather_covars_use), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
      .groups = "drop"
    )

  # A1) correlations
  cor_row <- tibble(
    year_wlag = year_value,
    season_wlag = season_value,
    n_weeks_wlag = nrow(dat_weekly),
    spearman_rho_wlag = NA_real_,
    spearman_p_wlag   = NA_real_,
    partial_rho_wlag  = NA_real_,
    partial_p_wlag    = NA_real_
  )

  if (nrow(dat_weekly) >= min_weeks) {
    ct <- suppressWarnings(cor.test(dat_weekly$total_claims_wlag, dat_weekly$spores_mean_wlag, method = "spearman"))
    cor_row$spearman_rho_wlag <- unname(ct$estimate)
    cor_row$spearman_p_wlag   <- ct$p.value

    zmat <- dat_weekly %>% dplyr::select(starts_with("mean_")) %>% as.data.frame()
    ok <- complete.cases(dat_weekly$total_claims_wlag, dat_weekly$spores_mean_wlag, zmat)

    if (sum(ok) >= min_weeks && ncol(zmat) > 0) {
      pt <- tryCatch(
        suppressWarnings(ppcor::pcor.test(
          x = dat_weekly$spores_mean_wlag[ok],
          y = dat_weekly$total_claims_wlag[ok],
          z = zmat[ok, , drop = FALSE],
          method = "spearman"
        )),
        error = function(e) NULL
      )
      if (!is.null(pt)) {
        cor_row$partial_rho_wlag <- pt$estimate
        cor_row$partial_p_wlag   <- pt$p.value
      }
    }
  }

  # A2) weekly logit (acute fraction)
  weekly_logit_row <- tibble(
    year_wlag = year_value,
    season_wlag = season_value,
    exposure_var_wlag = "spores_mean_wlag",
    n_weeks_wlag = nrow(dat_weekly),
    spores_or_wlag = NA_real_,
    spores_or_low_wlag = NA_real_,
    spores_or_high_wlag = NA_real_,
    spores_p_wlag = NA_real_
  )

  if (nrow(dat_weekly) >= min_weeks) {
    dat_week_logit <- dat_weekly %>% filter(total_claims_wlag > 0)
    if (nrow(dat_week_logit) >= min_weeks) {
      weather_week_cols <- names(dat_week_logit)[str_starts(names(dat_week_logit), "mean_")]

      f_week_logit <- as.formula(
        paste0(
          "cbind(acute_claims_wlag, nonacute_claims_wlag) ~ spores_mean_wlag",
          if (length(weather_week_cols) > 0) paste0(" + ", paste(weather_week_cols, collapse = " + ")) else ""
        )
      )

      fit_week_logit <- tryCatch(glm(f_week_logit, data = dat_week_logit, family = binomial()), error = function(e) NULL)
      if (!is.null(fit_week_logit)) {
        eff <- extract_effect_wlag(fit_week_logit, "spores_mean_wlag", exp_scale = TRUE)
        weekly_logit_row$spores_or_wlag      <- eff$estimate
        weekly_logit_row$spores_or_low_wlag  <- eff$ci_low
        weekly_logit_row$spores_or_high_wlag <- eff$ci_high
        weekly_logit_row$spores_p_wlag       <- eff$p
      }
    }
  }

  # A3) weekly NB (total counts)
  weekly_reg_row <- tibble(
    year_wlag = year_value,
    season_wlag = season_value,
    exposure_var_wlag = "spores_mean_wlag",
    n_weeks_wlag = nrow(dat_weekly),
    spores_rr_wlag = NA_real_,
    spores_rr_low_wlag = NA_real_,
    spores_rr_high_wlag = NA_real_,
    spores_p_wlag = NA_real_
  )

  if (nrow(dat_weekly) >= min_weeks) {
    weather_week_cols <- names(dat_weekly)[str_starts(names(dat_weekly), "mean_")]

    f_week_nb <- as.formula(
      paste0(
        "total_claims_wlag ~ spores_mean_wlag",
        if (length(weather_week_cols) > 0) paste0(" + ", paste(weather_week_cols, collapse = " + ")) else ""
      )
    )

    fit_week_nb <- tryCatch(MASS::glm.nb(f_week_nb, data = dat_weekly), error = function(e) NULL)
    if (!is.null(fit_week_nb)) {
      eff <- extract_effect_wlag(fit_week_nb, "spores_mean_wlag", exp_scale = TRUE)
      weekly_reg_row$spores_rr_wlag      <- eff$estimate
      weekly_reg_row$spores_rr_low_wlag  <- eff$ci_low
      weekly_reg_row$spores_rr_high_wlag <- eff$ci_high
      weekly_reg_row$spores_p_wlag       <- eff$p
    }
  }

  # B) DAILY models (lags/MA) ------
  exposure_vars <- c(
    "spores_10k_lag0_wlag",
    "spores_10k_lag1_wlag",
    "spores_10k_lag2_wlag",
    "spores_10k_lag3_wlag",
    "spores_10k_ma3_wlag",
    "spores_10k_ma7_wlag"
  )

  # Daily logit: claim-days only
  dat_daily_logit <- dat0 %>%
    filter(total_claims_wlag > 0) %>%
    filter(if_all(all_of(weather_covars_use), ~ !is.na(.x)))

  logit_daily_tbl <- map_dfr(exposure_vars, function(ex_var) {
    if (!ex_var %in% names(dat_daily_logit)) return(NULL)

    d <- dat_daily_logit %>% filter(!is.na(.data[[ex_var]]))

    row <- tibble(
      year_wlag = year_value,
      season_wlag = season_value,
      exposure_var_wlag = ex_var,
      n_days_wlag = nrow(d),
      spores_or_wlag = NA_real_,
      spores_or_low_wlag = NA_real_,
      spores_or_high_wlag = NA_real_,
      spores_p_wlag = NA_real_
    )

    if (nrow(d) < min_days) return(row)

    f <- as.formula(
      paste0(
        "cbind(acute_claims_wlag, nonacute_claims_wlag) ~ ",
        ex_var, " + ",
        paste(weather_covars_use, collapse = " + "),
        " + ns(doy_wlag, df = 6) + factor(dow_wlag)"
      )
    )

    fit <- tryCatch(glm(f, data = d, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(row)

    eff <- extract_effect_wlag(fit, ex_var, exp_scale = TRUE)
    row$spores_or_wlag      <- eff$estimate
    row$spores_or_low_wlag  <- eff$ci_low
    row$spores_or_high_wlag <- eff$ci_high
    row$spores_p_wlag       <- eff$p
    row
  })

  # Daily NB: total counts
  dat_daily_nb <- dat0 %>%
    filter(if_all(all_of(weather_covars_use), ~ !is.na(.x)))

  reg_daily_tbl <- map_dfr(exposure_vars, function(ex_var) {
    if (!ex_var %in% names(dat_daily_nb)) return(NULL)

    d <- dat_daily_nb %>% filter(!is.na(.data[[ex_var]]))

    row <- tibble(
      year_wlag = year_value,
      season_wlag = season_value,
      exposure_var_wlag = ex_var,
      n_days_wlag = nrow(d),
      spores_rr_wlag = NA_real_,
      spores_rr_low_wlag = NA_real_,
      spores_rr_high_wlag = NA_real_,
      spores_p_wlag = NA_real_
    )

    if (nrow(d) < min_days) return(row)

    f <- as.formula(
      paste0(
        "total_claims_wlag ~ ",
        ex_var, " + ",
        paste(weather_covars_use, collapse = " + "),
        " + ns(doy_wlag, df = 6) + factor(dow_wlag)"
      )
    )

    fit <- tryCatch(MASS::glm.nb(f, data = d), error = function(e) NULL)
    if (is.null(fit)) return(row)

    eff <- extract_effect_wlag(fit, ex_var, exp_scale = TRUE)
    row$spores_rr_wlag      <- eff$estimate
    row$spores_rr_low_wlag  <- eff$ci_low
    row$spores_rr_high_wlag <- eff$ci_high
    row$spores_p_wlag       <- eff$p
    row
  })

  list(
    cor_wlag          = cor_row,
    logit_weekly_wlag = weekly_logit_row,
    reg_weekly_wlag   = weekly_reg_row,
    logit_daily_wlag  = logit_daily_tbl,
    reg_daily_wlag    = reg_daily_tbl
  )
}

season_grid_wlag <- tibble(
  season_wlag = c("All seasons", "winter", "spring", "summer", "fall"),
  min_weeks = c(min_weeks_overall_wlag, rep(min_weeks_season_wlag, 4)),
  min_days  = c(min_days_overall_wlag,  rep(min_days_season_wlag, 4))
)

run_grid_wlag <- crossing(
  year_wlag = years_to_run_wlag,
  season_grid_wlag
)

results_list_wlag <- pmap(
  run_grid_wlag,
  function(year_wlag, season_wlag, min_weeks, min_days) {
    run_models_one_year_season_wlag(
      year_value = year_wlag,
      season_value = season_wlag,
      dat = daily_merged_wlag,
      weather_covars_use = weather_covars_use_wlag,
      min_weeks = min_weeks,
      min_days  = min_days
    )
  }
)

cor_table_wlag          <- bind_rows(map(results_list_wlag, "cor_wlag"))
logit_weekly_table_wlag <- bind_rows(map(results_list_wlag, "logit_weekly_wlag"))
reg_weekly_table_wlag   <- bind_rows(map(results_list_wlag, "reg_weekly_wlag"))
logit_daily_table_wlag  <- bind_rows(map(results_list_wlag, "logit_daily_wlag"))
reg_daily_table_wlag    <- bind_rows(map(results_list_wlag, "reg_daily_wlag"))

write.csv(cor_table_wlag,
          file.path(out_dir_wlag, "SJ_region_spores_asthma_correlations_by_year_season_wlag.csv"),
          row.names = FALSE)

write.csv(logit_daily_table_wlag,
          file.path(out_dir_wlag, "SJ_region_spores_asthma_logit_DAILY_lags_by_year_season_wlag.csv"),
          row.names = FALSE)

write.csv(reg_daily_table_wlag,
          file.path(out_dir_wlag, "SJ_region_spores_asthma_regression_DAILY_lags_by_year_season_wlag.csv"),
          row.names = FALSE)

write.csv(logit_weekly_table_wlag,
          file.path(out_dir_wlag, "SJ_region_spores_asthma_logit_WEEKLY_by_year_season_wlag.csv"),
          row.names = FALSE)

write.csv(reg_weekly_table_wlag,
          file.path(out_dir_wlag, "SJ_region_spores_asthma_regression_WEEKLY_by_year_season_wlag.csv"),
          row.names = FALSE)

print(cor_table_wlag)
print(logit_weekly_table_wlag)
print(reg_weekly_table_wlag)
print(logit_daily_table_wlag)
print(reg_daily_table_wlag)

cat("\nDone. Outputs written to:", out_dir_wlag, "\n")
