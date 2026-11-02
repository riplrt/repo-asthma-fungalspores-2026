############################################################
## San Juan Health Region: yearly spores ↔ asthma models (v2_a_adj + age)
## (correlations, logit, regression) controlled by weather AND age
##
## WEEKLY EXPOSURE QC (spores + pollen, separately):
##   - For weeks with <= 3 missing days: impute missing daily values
##     using that week's mean of observed values
##   - For weeks with > 3 missing days (or all missing): exclude that week
##     by setting exposure values to NA for all days in that week
##
## AGE ADJUSTMENT:
##   - Creates age_group_v2_a_adj from age (or uses existing age_group)
##   - Age-stratified weekly correlations (Spearman + partial Spearman)
##   - Daily logit (acute fraction) and NB regression include age_group_v2_a_adj
##
## Outputs:
##   1) statistical_results/SJ_region_spores_asthma_correlations_by_year_age_v2_a_adj.csv
##   2) statistical_results/SJ_region_spores_asthma_logit_by_year_age_v2_a_adj.csv
##   3) statistical_results/SJ_region_spores_asthma_regression_by_year_age_v2_a_adj.csv
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(lubridate)
  library(zoo)
  library(rnoaa)
  library(ppcor)
  library(MASS)
  library(splines)
  library(broom)
})

# ===========================================================
# 1) SETTINGS (EDIT)
# ===========================================================

claims_path_v2_a_adj  <- "data/2022_asma_From_date_ver 3.xlsx"
claims_sheet_v2_a_adj <- 1

# San Juan station spores/pollen files (can be multiple; row-bound)
spore_files_sj_v2_a_adj <- c("data/Datos de esporas San Juan Jan 2018-2022.xlsx")

years_to_run_v2_a_adj <- 2018:2021

# NOAA token (recommended: set env var NOAATOKEN)
# Sys.setenv(NOAATOKEN="YOUR_TOKEN"); options(noaakey = Sys.getenv("NOAATOKEN"))
options(noaakey = Sys.getenv("NOAATOKEN"))

# GHCND station ID for SJU area
sju_station_id_v2_a_adj <- "RQW00011641"

clinical_settings_v2_a_adj <- c("ED", "Inpatient", "Outpatient", "Urgent")

# Weather covariates to adjust for (must match pull_weather_sju_v2_a_adj outputs)
weather_covars_v2_a_adj <- c("tmean_c_v2_a_adj", "rh_mean_v2_a_adj", "wind_ms_v2_a_adj", "prcp_mm_v2_a_adj", "dewpoint_c_v2_a_adj")

# Exposure scaling
spore_scale_divisor_v2_a_adj <- 10000

# Weekly QC rule
max_missing_days_week_v2_a_adj <- 3

# Minimum data thresholds (tune as needed)
min_weeks_for_cor_v2_a_adj <- 10
min_rows_for_models_v2_a_adj <- 60

# Output directory
out_dir_v2_a_adj <- "statistical_results"
dir.create(out_dir_v2_a_adj, showWarnings = FALSE, recursive = TRUE)

# ===========================================================
# 2) HELPERS
# ===========================================================

parse_date_any_v2_a_adj <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))  # Excel dates
  x <- as.character(x)
  d <- suppressWarnings(lubridate::ymd(x))
  d <- dplyr::coalesce(d, suppressWarnings(lubridate::mdy(x)))
  d <- dplyr::coalesce(d, suppressWarnings(lubridate::dmy(x)))
  dplyr::coalesce(d, suppressWarnings(as.Date(x)))
}

care_setting_fun_v2_a_adj <- function(hcg_setting, hcg_mr_line_desc) {
  s <- tolower(as.character(hcg_setting))
  l <- tolower(as.character(hcg_mr_line_desc))

  if (grepl("ed|emergency", s) || grepl("\\bed\\b|\\ber\\b", l)) return("ED")
  if (grepl("inpatient|hospital", s) || grepl("inpatient|hospital", l)) return("Inpatient")
  if (grepl("urgent", s) || grepl("urgent", l)) return("Urgent")
  if (grepl("office|clinic|outpatient", s) || grepl("office|clinic", l)) return("Outpatient")
  if (grepl("dme|supplies", l)) return("DME")
  "Other"
}

ensure_care_setting_v2_a_adj <- function(dat) {
  if ("care_setting_v2_a_adj" %in% names(dat)) return(dat)

  if ("care_setting" %in% names(dat)) {
    return(dat %>% mutate(care_setting_v2_a_adj = care_setting))
  }

  needed <- c("hcg_setting", "hcg_mr_line_desc")
  if (!all(needed %in% names(dat))) {
    stop("No care_setting/care_setting_v2_a_adj and missing hcg_setting/hcg_mr_line_desc to create it.")
  }

  dat %>%
    mutate(care_setting_v2_a_adj = mapply(care_setting_fun_v2_a_adj, hcg_setting, hcg_mr_line_desc))
}

ensure_asthma_primary_v2_a_adj <- function(dat) {
  if ("asthma_primary_v2_a_adj" %in% names(dat)) return(dat)

  if ("asthma_primary" %in% names(dat)) {
    return(dat %>% mutate(asthma_primary_v2_a_adj = asthma_primary))
  }

  if (!"icd_diag_01_primary" %in% names(dat)) {
    stop("Missing asthma_primary/asthma_primary_v2_a_adj and icd_diag_01_primary; create asthma_primary before modeling.")
  }

  dat %>%
    mutate(
      icd_diag_01_primary = stringr::str_to_upper(stringr::str_trim(as.character(icd_diag_01_primary))),
      asthma_primary_v2_a_adj = if_else(stringr::str_detect(icd_diag_01_primary, "^J45|^J46"),
                                 "asthma", "non-asthma")
    )
}

create_age_groups_v2_a_adj <- function(dat) {
  if ("age_group_v2_a_adj" %in% names(dat)) return(dat)

  if ("age_group" %in% names(dat)) {
    return(dat %>% mutate(age_group_v2_a_adj = as.factor(age_group)))
  }

  if (!"age" %in% names(dat)) {
    warning("No age column found; cannot create age_group_v2_a_adj. Setting to NA.")
    return(dat %>% mutate(age_group_v2_a_adj = NA_character_))
  }

  dat %>%
    mutate(
      age_raw_v2_a_adj = stringr::str_trim(as.character(age)),
      age_num_v2_a_adj = dplyr::case_when(
        stringr::str_detect(age_raw_v2_a_adj, "\\+$") ~ suppressWarnings(as.numeric(stringr::str_remove(age_raw_v2_a_adj, "\\+"))),
        TRUE ~ suppressWarnings(as.numeric(age_raw_v2_a_adj))
      ),
      age_group_v2_a_adj = cut(
        age_num_v2_a_adj,
        breaks = c(0, 4, 11, 17, 39, 64, 74, 84, Inf),
        labels = c("0-4", "5-11", "12-17", "18-39", "40-64", "65-74", "75-84", "85+"),
        right  = TRUE,
        include.lowest = TRUE
      )
    ) %>%
    dplyr::select(-age_raw_v2_a_adj, -age_num_v2_a_adj)
}

read_claims_v2_a_adj <- function(path, sheet = 1) {
  dat <- read_excel(path, sheet = sheet) %>%
    clean_names()

  if (!"from_date" %in% names(dat)) stop("Missing from_date (after clean_names()).")

  dat %>%
    mutate(
      from_date = parse_date_any_v2_a_adj(from_date),
      year_v2_a_adj   = lubridate::year(from_date)
    )
}

read_spores_sj_onefile_v2_a_adj <- function(file) {
  x <- read_excel(file) %>% clean_names()

  date_col <- names(x)[1]
  spore_col <- names(x)[stringr::str_detect(names(x), "esporas|spores")][1]
  pollen_col <- names(x)[stringr::str_detect(names(x), "polen|pollen")][1]

  if (is.na(spore_col)) stop("Could not find spores column in: ", file)

  x %>%
    transmute(
      date         = parse_date_any_v2_a_adj(.data[[date_col]]),
      spores_total = suppressWarnings(as.numeric(.data[[spore_col]])),
      pollen_total = if (!is.na(pollen_col)) suppressWarnings(as.numeric(.data[[pollen_col]])) else NA_real_
    ) %>%
    filter(!is.na(date))
}

weekly_mean_impute_or_exclude_v2_a_adj <- function(df, value_col, max_missing_days_week = 3) {
  value_sym  <- rlang::ensym(value_col)
  value_name <- rlang::as_string(value_sym)

  out_name      <- paste0(value_name, "_v2_a_adj")
  mean_name     <- paste0(value_name, "_week_mean_v2_a_adj")
  missing_name  <- paste0(value_name, "_n_missing_week_v2_a_adj")
  imputed_name  <- paste0(value_name, "_imputed_v2_a_adj")
  excluded_name <- paste0(value_name, "_week_excluded_v2_a_adj")

  df %>%
    mutate(week_start_v2_a_adj = lubridate::floor_date(date_v2_a_adj, unit = "week", week_start = 1)) %>%
    group_by(week_start_v2_a_adj) %>%
    mutate(
      n_missing_tmp_v2_a_adj = sum(is.na(!!value_sym)),
      week_mean_tmp_v2_a_adj = mean(!!value_sym, na.rm = TRUE),
      week_mean_tmp_v2_a_adj = if_else(is.nan(week_mean_tmp_v2_a_adj), NA_real_, week_mean_tmp_v2_a_adj),
      week_ok_tmp_v2_a_adj   = (n_missing_tmp_v2_a_adj <= max_missing_days_week) & !is.na(week_mean_tmp_v2_a_adj),

      !!missing_name  := n_missing_tmp_v2_a_adj,
      !!mean_name     := week_mean_tmp_v2_a_adj,
      !!excluded_name := !week_ok_tmp_v2_a_adj,
      !!imputed_name  := week_ok_tmp_v2_a_adj & is.na(!!value_sym),

      !!out_name := case_when(
        week_ok_tmp_v2_a_adj & is.na(!!value_sym) ~ week_mean_tmp_v2_a_adj,
        week_ok_tmp_v2_a_adj ~ as.numeric(!!value_sym),
        TRUE ~ NA_real_
      )
    ) %>%
    ungroup() %>%
    dplyr::select(-n_missing_tmp_v2_a_adj, -week_mean_tmp_v2_a_adj, -week_ok_tmp_v2_a_adj)
}

pull_weather_sju_v2_a_adj <- function(y, station_id) {
  vars <- c("PRCP", "TMAX", "TMIN", "AWND", "ADPT", "RHAV")

  w <- meteo_pull_monitors(
    monitors   = station_id,
    var        = vars,
    date_min   = sprintf("%d-01-01", y),
    date_max   = sprintf("%d-12-31", y),
    keep_flags = FALSE
  )

  for (nm in c("prcp","tmax","tmin","awnd","adpt","rhav")) {
    if (!nm %in% names(w)) w[[nm]] <- NA_real_
  }

  w %>%
    mutate(
      date = as.Date(date),
      prcp_mm_v2_a_adj = prcp / 10,
      tmax_c_v2_a_adj  = tmax / 10,
      tmin_c_v2_a_adj  = tmin / 10,
      tmean_c_v2_a_adj = (tmax_c_v2_a_adj + tmin_c_v2_a_adj) / 2,
      wind_ms_v2_a_adj = awnd,
      dewpoint_c_v2_a_adj = adpt / 10,
      rh_mean_v2_a_adj = rhav
    ) %>%
    dplyr::select(date, prcp_mm_v2_a_adj, tmax_c_v2_a_adj, tmin_c_v2_a_adj, tmean_c_v2_a_adj, wind_ms_v2_a_adj, dewpoint_c_v2_a_adj, rh_mean_v2_a_adj)
}

extract_term_exp_v2_a_adj <- function(model, term) {
  sm <- summary(model)
  coefs <- sm$coefficients
  if (!term %in% rownames(coefs)) {
    return(tibble(effect = NA_real_, low = NA_real_, high = NA_real_, p = NA_real_))
  }
  est <- coefs[term, "Estimate"]
  se  <- coefs[term, "Std. Error"]
  pcol <- intersect(colnames(coefs), c("Pr(>|z|)", "Pr(>|t|)"))
  pval <- if (length(pcol) == 1) coefs[term, pcol] else NA_real_
  ci <- est + c(-1, 1) * 1.96 * se
  tibble(effect = exp(est), low = exp(ci[1]), high = exp(ci[2]), p = pval)
}

safe_spearman_v2_a_adj <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(tibble(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  tibble(rho = unname(ct$estimate), p = ct$p.value)
}

safe_partial_spearman_v2_a_adj <- function(x, y, zmat) {
  ok <- complete.cases(x, y, zmat)
  if (sum(ok) < 5 || ncol(zmat) < 1) return(tibble(partial_rho = NA_real_, partial_p = NA_real_))
  pt <- tryCatch(
    suppressWarnings(ppcor::pcor.test(x[ok], y[ok], z = zmat[ok, , drop = FALSE], method = "spearman")),
    error = function(e) NULL
  )
  if (is.null(pt)) return(tibble(partial_rho = NA_real_, partial_p = NA_real_))
  tibble(partial_rho = as.numeric(pt$estimate), partial_p = as.numeric(pt$p.value))
}

# ===========================================================
# 3) READ + PREP CLAIMS (San Juan health region)
# ===========================================================

claims_raw_v2_a_adj <- read_claims_v2_a_adj(claims_path_v2_a_adj, claims_sheet_v2_a_adj) %>%
  ensure_care_setting_v2_a_adj() %>%
  ensure_asthma_primary_v2_a_adj() %>%
  create_age_groups_v2_a_adj()

claims_df_v2_a_adj <- claims_raw_v2_a_adj %>%
  mutate(
    municipality_desc = as.character(municipality_desc),
    health_region_v2_a_adj = if_else(municipality_desc == "Caguas", "Caguas", "San Juan"),
    health_region_v2_a_adj = factor(health_region_v2_a_adj, levels = c("San Juan", "Caguas"))
  )

claims_sj_region_v2_a_adj <- claims_df_v2_a_adj %>%
  filter(
    year_v2_a_adj %in% years_to_run_v2_a_adj,
    health_region_v2_a_adj == "San Juan",
    asthma_primary_v2_a_adj == "asthma",
    care_setting_v2_a_adj %in% clinical_settings_v2_a_adj,
    !is.na(age_group_v2_a_adj)
  ) %>%
  mutate(
    acute_flag_v2_a_adj = as.integer(care_setting_v2_a_adj %in% c("ED", "Inpatient")),
    date_v2_a_adj = as.Date(from_date),
    age_group_v2_a_adj = droplevels(age_group_v2_a_adj)
  )

daily_counts_age_v2_a_adj <- claims_sj_region_v2_a_adj %>%
  group_by(date_v2_a_adj, age_group_v2_a_adj) %>%
  summarise(
    total_claims_v2_a_adj    = n(),
    acute_claims_v2_a_adj    = sum(acute_flag_v2_a_adj == 1, na.rm = TRUE),
    nonacute_claims_v2_a_adj = total_claims_v2_a_adj - acute_claims_v2_a_adj,
    .groups = "drop"
  )

calendar_v2_a_adj <- tibble(
  date_v2_a_adj = seq.Date(as.Date(sprintf("%d-01-01", min(years_to_run_v2_a_adj))),
                     as.Date(sprintf("%d-12-31", max(years_to_run_v2_a_adj))),
                     by = "day")
) %>%
  mutate(
    year_v2_a_adj = lubridate::year(date_v2_a_adj),
    doy_v2_a_adj  = lubridate::yday(date_v2_a_adj),
    dow_v2_a_adj  = lubridate::wday(date_v2_a_adj, label = TRUE, week_start = 1),
    week_start_v2_a_adj = lubridate::floor_date(date_v2_a_adj, unit = "week", week_start = 1),
    week_of_year_v2_a_adj = lubridate::isoweek(date_v2_a_adj)
  ) %>%
  filter(year_v2_a_adj %in% years_to_run_v2_a_adj)

age_levels_v2_a_adj <- levels(claims_sj_region_v2_a_adj$age_group_v2_a_adj)

daily_outcomes_age_v2_a_adj <- tidyr::expand_grid(
  date_v2_a_adj = calendar_v2_a_adj$date_v2_a_adj,
  age_group_v2_a_adj = factor(age_levels_v2_a_adj, levels = age_levels_v2_a_adj)
) %>%
  left_join(calendar_v2_a_adj, by = "date_v2_a_adj") %>%
  left_join(daily_counts_age_v2_a_adj, by = c("date_v2_a_adj", "age_group_v2_a_adj")) %>%
  mutate(
    across(c(total_claims_v2_a_adj, acute_claims_v2_a_adj, nonacute_claims_v2_a_adj),
           ~ tidyr::replace_na(.x, 0L))
  )

# ===========================================================
# 4) READ + PREP SPORES/POLLEN (weekly QC)
# ===========================================================

spores_sj_raw_v2_a_adj <- purrr::map_dfr(spore_files_sj_v2_a_adj, read_spores_sj_onefile_v2_a_adj) %>%
  group_by(date) %>%
  summarise(
    spores_total = mean(spores_total, na.rm = TRUE),
    pollen_total = mean(pollen_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(date)

spores_sj_daily_v2_a_adj <- tibble(date_v2_a_adj = calendar_v2_a_adj$date_v2_a_adj) %>%
  left_join(spores_sj_raw_v2_a_adj, by = c("date_v2_a_adj" = "date"))

spores_sj_daily_v2_a_adj <- weekly_mean_impute_or_exclude_v2_a_adj(spores_sj_daily_v2_a_adj, spores_total, max_missing_days_week_v2_a_adj)
spores_sj_daily_v2_a_adj <- weekly_mean_impute_or_exclude_v2_a_adj(spores_sj_daily_v2_a_adj, pollen_total, max_missing_days_week_v2_a_adj)

spores_sj_daily_v2_a_adj <- spores_sj_daily_v2_a_adj %>%
  mutate(spores_10k_v2_a_adj = spores_total_v2_a_adj / spore_scale_divisor_v2_a_adj)

# ===========================================================
# 5) WEATHER + MERGE DAILY MODEL DATA
# ===========================================================

weather_sj_v2_a_adj <- purrr::map_dfr(years_to_run_v2_a_adj, ~ pull_weather_sju_v2_a_adj(.x, sju_station_id_v2_a_adj))

daily_model_df_v2_a_adj <- daily_outcomes_age_v2_a_adj %>%
  left_join(
    spores_sj_daily_v2_a_adj %>%
      dplyr::select(date_v2_a_adj, spores_10k_v2_a_adj,
             spores_total_v2_a_adj, spores_total_week_excluded_v2_a_adj,
             pollen_total_v2_a_adj, pollen_total_week_excluded_v2_a_adj),
    by = "date_v2_a_adj"
  ) %>%
  left_join(weather_sj_v2_a_adj, by = c("date_v2_a_adj" = "date"))

weather_covars_use_v2_a_adj <- intersect(weather_covars_v2_a_adj, names(daily_model_df_v2_a_adj))
if (length(weather_covars_use_v2_a_adj) == 0) stop("None of the requested weather covariates were found after weather pull.")

# ===========================================================
# 6) CORRELATIONS (weekly, age-stratified)
# ===========================================================

calc_correlations_year_age_v2_a_adj <- function(y, dat) {
  dat_y <- dat %>%
    filter(year_v2_a_adj == y) %>%
    filter(!is.na(spores_10k_v2_a_adj)) %>%
    filter(if_all(all_of(weather_covars_use_v2_a_adj), ~ !is.na(.x))) %>%
    filter(week_of_year_v2_a_adj <= 52)

  weekly_y <- dat_y %>%
    group_by(age_group_v2_a_adj, week_start_v2_a_adj) %>%
    summarise(
      total_claims_wk_v2_a_adj = sum(total_claims_v2_a_adj),
      spores_mean_wk_v2_a_adj  = mean(spores_10k_v2_a_adj, na.rm = TRUE),
      across(all_of(weather_covars_use_v2_a_adj), ~ mean(.x, na.rm = TRUE), .names = "{.col}_wkmean_v2_a_adj"),
      .groups = "drop"
    )

  weather_wk_cols <- names(weekly_y)[stringr::str_ends(names(weekly_y), "_wkmean_v2_a_adj")]

  weekly_y %>%
    group_by(age_group_v2_a_adj) %>%
    group_modify(~{
      df_g <- .x
      n_weeks <- nrow(df_g)

      spe <- if (n_weeks >= min_weeks_for_cor_v2_a_adj) {
        safe_spearman_v2_a_adj(df_g$total_claims_wk_v2_a_adj, df_g$spores_mean_wk_v2_a_adj)
      } else tibble(rho = NA_real_, p = NA_real_)

      zmat <- if (length(weather_wk_cols) > 0) as.matrix(df_g[, weather_wk_cols, drop = FALSE]) else matrix(nrow=nrow(df_g), ncol=0)

      part <- if (n_weeks >= min_weeks_for_cor_v2_a_adj && ncol(zmat) > 0) {
        safe_partial_spearman_v2_a_adj(df_g$spores_mean_wk_v2_a_adj, df_g$total_claims_wk_v2_a_adj, zmat)
      } else tibble(partial_rho = NA_real_, partial_p = NA_real_)

      tibble(
        year_v2_a_adj = y,
        n_weeks_v2_a_adj = n_weeks,
        spearman_rho_v2_a_adj = spe$rho,
        spearman_p_v2_a_adj   = spe$p,
        partial_rho_v2_a_adj  = part$partial_rho,
        partial_p_v2_a_adj    = part$partial_p
      )
    }) %>%
    ungroup()
}

cor_table_v2_a_adj <- purrr::map_dfr(years_to_run_v2_a_adj, ~ calc_correlations_year_age_v2_a_adj(.x, daily_model_df_v2_a_adj))

write.csv(cor_table_v2_a_adj,
          file.path(out_dir_v2_a_adj, "SJ_region_spores_asthma_correlations_by_year_age_v2_a_adj.csv"),
          row.names = FALSE)

# ===========================================================
# 7) YEAR-SPECIFIC MODELS (logit + NB), age-adjusted
# ===========================================================

fit_year_models_v2_a_adj <- function(y, dat) {
  dat_y <- dat %>%
    filter(year_v2_a_adj == y) %>%
    filter(!is.na(spores_10k_v2_a_adj)) %>%
    filter(if_all(all_of(weather_covars_use_v2_a_adj), ~ !is.na(.x)))

  dat_logit <- dat_y %>% filter(total_claims_v2_a_adj > 0)

  logit_row <- tibble(
    year_v2_a_adj = y,
    n_rows_v2_a_adj = nrow(dat_logit),
    spores_or_v2_a_adj = NA_real_,
    spores_or_low_v2_a_adj = NA_real_,
    spores_or_high_v2_a_adj = NA_real_,
    spores_p_v2_a_adj = NA_real_
  )

  if (nrow(dat_logit) >= min_rows_for_models_v2_a_adj) {
    f_logit <- as.formula(
      paste0(
        "cbind(acute_claims_v2_a_adj, nonacute_claims_v2_a_adj) ~ spores_10k_v2_a_adj + age_group_v2_a_adj + ",
        paste(weather_covars_use_v2_a_adj, collapse = " + "),
        " + ns(doy_v2_a_adj, df = 6) + factor(dow_v2_a_adj)"
      )
    )

    fit_logit <- tryCatch(glm(f_logit, data = dat_logit, family = binomial(link = "logit")),
                          error = function(e) NULL)

    if (!is.null(fit_logit)) {
      eff <- extract_term_exp_v2_a_adj(fit_logit, "spores_10k_v2_a_adj")
      logit_row$spores_or_v2_a_adj      <- eff$effect
      logit_row$spores_or_low_v2_a_adj  <- eff$low
      logit_row$spores_or_high_v2_a_adj <- eff$high
      logit_row$spores_p_v2_a_adj       <- eff$p
    }
  }

  reg_row <- tibble(
    year_v2_a_adj = y,
    n_rows_v2_a_adj = nrow(dat_y),
    spores_rr_v2_a_adj = NA_real_,
    spores_rr_low_v2_a_adj = NA_real_,
    spores_rr_high_v2_a_adj = NA_real_,
    spores_p_v2_a_adj = NA_real_
  )

  if (nrow(dat_y) >= min_rows_for_models_v2_a_adj) {
    f_nb <- as.formula(
      paste0(
        "total_claims_v2_a_adj ~ spores_10k_v2_a_adj + age_group_v2_a_adj + ",
        paste(weather_covars_use_v2_a_adj, collapse = " + "),
        " + ns(doy_v2_a_adj, df = 6) + factor(dow_v2_a_adj)"
      )
    )

    fit_nb <- tryCatch(MASS::glm.nb(f_nb, data = dat_y),
                       error = function(e) NULL)

    if (!is.null(fit_nb)) {
      eff <- extract_term_exp_v2_a_adj(fit_nb, "spores_10k_v2_a_adj")
      reg_row$spores_rr_v2_a_adj      <- eff$effect
      reg_row$spores_rr_low_v2_a_adj  <- eff$low
      reg_row$spores_rr_high_v2_a_adj <- eff$high
      reg_row$spores_p_v2_a_adj       <- eff$p
    }
  }

  list(logit = logit_row, reg = reg_row)
}

model_list_v2_a_adj <- purrr::map(years_to_run_v2_a_adj, ~ fit_year_models_v2_a_adj(.x, daily_model_df_v2_a_adj))

logit_table_v2_a_adj <- bind_rows(purrr::map(model_list_v2_a_adj, "logit"))
reg_table_v2_a_adj   <- bind_rows(purrr::map(model_list_v2_a_adj, "reg"))

write.csv(logit_table_v2_a_adj,
          file.path(out_dir_v2_a_adj, "SJ_region_spores_asthma_logit_by_year_age_v2_a_adj.csv"),
          row.names = FALSE)

write.csv(reg_table_v2_a_adj,
          file.path(out_dir_v2_a_adj, "SJ_region_spores_asthma_regression_by_year_age_v2_a_adj.csv"),
          row.names = FALSE)

cat("\nDone. Outputs written to: ", out_dir_v2_a_adj, "\n", sep = "")
