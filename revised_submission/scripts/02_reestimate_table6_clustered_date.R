#!/usr/bin/env Rscript
############################################################
## Re-estimate Table 6 Family with Date-Clustered SEs
## Focus:  Pandemic phase x age, encounter volume and acute fraction
##   - Daily NB GAM + binomial GAM
##   - Key term: spores_10k_ma7 x pandemic_phase x age_group
##   - Controls: cyclic day-of-year spline + day-of-week + smooth time trend + weather
##   - BH adjustment across the 48-test family
## Note: uses the original corrected spore file and rnoaa weather; see
##   01_spore_age_period_interaction_tests.R for the San Juan series rerun
##
## Outputs (CSV) written to: R2_multiplicity_clusteredSE/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(lubridate)
  library(rnoaa)
  library(mgcv)
  library(zoo)
  library(sandwich)
  library(lmtest)
})

claims_path <- "data/2022_asma_From_date_ver 3.xlsx"
claims_sheet <- 1
spore_files_sj <- c("data/Datos esporas y polen 24-08-19 a 31-12-22 Corregidas 12-09-25.xlsx")
out_dir <- "Journal_Aerobiology_Submission_response/R2_multiplicity_clusteredSE/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

options(noaakey = Sys.getenv("NOAATOKEN"))
sju_station_id <- "RQW00011641"
years_all <- 2018:2022
clinical_settings <- c("ED", "Inpatient", "Outpatient", "Urgent")
max_missing_days_week <- 3
spore_scale_divisor <- 10000
weather_covars <- c("tmean_c", "rh_mean", "wind_ms", "prcp_mm", "dewpoint_c")

age_breaks <- c(0, 4, 11, 17, 39, 64, 74, 84, Inf)
age_labels <- c("0-4", "5-11", "12-17", "18-39", "40-64", "65-74", "75-84", "85+")

san_juan_munis <- c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA")
caguas_munis <- c("CAGUAS", "GURABO", "JUNCOS", "SAN LORENZO")

parse_date_any <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x <- as.character(x)
  d <- suppressWarnings(ymd(x))
  d <- dplyr::coalesce(d, suppressWarnings(mdy(x)))
  d <- dplyr::coalesce(d, suppressWarnings(dmy(x)))
  dplyr::coalesce(d, suppressWarnings(as.Date(x)))
}

safe_mean <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

care_setting_fun <- function(hcg_setting, hcg_mr_line_desc) {
  s <- tolower(as.character(hcg_setting))
  l <- tolower(as.character(hcg_mr_line_desc))
  if (grepl("ed|emergency", s) || grepl("\\bed\\b|\\ber\\b", l)) return("ED")
  if (grepl("inpatient|hospital", s) || grepl("inpatient|hospital", l)) return("Inpatient")
  if (grepl("urgent", s) || grepl("urgent", l)) return("Urgent")
  if (grepl("office|clinic|outpatient", s) || grepl("office|clinic", l)) return("Outpatient")
  if (grepl("dme|supplies", l)) return("DME")
  "Other"
}

ensure_care_setting <- function(df) {
  if ("care_setting" %in% names(df)) return(df)
  df |> mutate(care_setting = mapply(care_setting_fun, hcg_setting, hcg_mr_line_desc))
}

ensure_asthma_primary <- function(df) {
  if ("asthma_primary" %in% names(df)) return(df)
  df |>
    mutate(
      icd_diag_01_primary = stringr::str_to_upper(stringr::str_trim(as.character(icd_diag_01_primary))),
      asthma_primary = if_else(stringr::str_detect(icd_diag_01_primary, "^J45|^J46"), "asthma", "other")
    )
}

ensure_age_group <- function(df) {
  if ("age_group" %in% names(df)) return(df |> mutate(age_group = as.factor(age_group)))
  df |>
    mutate(
      age_raw = stringr::str_trim(as.character(age)),
      age_num = dplyr::case_when(
        stringr::str_detect(age_raw, "\\+$") ~ suppressWarnings(as.numeric(stringr::str_remove(age_raw, "\\+"))),
        TRUE ~ suppressWarnings(as.numeric(age_raw))
      ),
      age_group = cut(age_num, breaks = age_breaks, labels = age_labels, right = TRUE, include.lowest = TRUE)
    ) |>
    select(-age_raw, -age_num)
}

weekly_qc_impute_or_exclude <- function(df, value_col, base_name, max_missing = 3) {
  value_sym <- rlang::ensym(value_col)
  out_col <- paste0(base_name)
  week_excl_col <- paste0(base_name, "_week_excluded")

  df |>
    mutate(week_start = floor_date(date, unit = "week", week_start = 1)) |>
    group_by(week_start) |>
    mutate(
      n_missing_week = sum(is.na(!!value_sym)),
      week_mean = safe_mean(!!value_sym),
      week_ok = (n_missing_week <= max_missing) & !is.na(week_mean),
      !!week_excl_col := !week_ok,
      !!out_col := case_when(
        week_ok & is.na(!!value_sym) ~ week_mean,
        week_ok ~ as.numeric(!!value_sym),
        TRUE ~ NA_real_
      )
    ) |>
    ungroup() |>
    select(-n_missing_week, -week_mean, -week_ok)
}

needs_div10 <- function(x, thr) {
  x <- suppressWarnings(as.numeric(x))
  q <- suppressWarnings(quantile(abs(x), probs = 0.99, na.rm = TRUE, names = FALSE))
  if (is.na(q)) return(FALSE)
  q > thr
}

pull_weather_sju <- function(year_value, station_id) {
  vars <- c("PRCP", "TMAX", "TMIN", "AWND", "ADPT", "RHAV")
  w <- meteo_pull_monitors(
    monitors = station_id,
    var = vars,
    date_min = sprintf("%d-01-01", year_value),
    date_max = sprintf("%d-12-31", year_value),
    keep_flags = FALSE
  ) |>
    as_tibble() |>
    mutate(date = as.Date(date))

  for (nm in c("prcp", "tmax", "tmin", "awnd", "adpt", "rhav")) {
    if (!nm %in% names(w)) w[[nm]] <- NA_real_
  }

  prcp_mm <- if (needs_div10(w$prcp, 300)) w$prcp / 10 else w$prcp
  tmax_c <- if (needs_div10(w$tmax, 80)) w$tmax / 10 else w$tmax
  tmin_c <- if (needs_div10(w$tmin, 80)) w$tmin / 10 else w$tmin
  dewpoint_c <- if (needs_div10(w$adpt, 80)) w$adpt / 10 else w$adpt
  wind_ms <- if (needs_div10(w$awnd, 30)) w$awnd / 10 else w$awnd

  tibble(
    date = w$date,
    prcp_mm = prcp_mm,
    tmax_c = tmax_c,
    tmin_c = tmin_c,
    tmean_c = (tmax_c + tmin_c) / 2,
    dewpoint_c = dewpoint_c,
    rh_mean = w$rhav,
    wind_ms = wind_ms
  )
}

lincomb_stats <- function(beta, V, idx) {
  b <- sum(beta[names(idx)] * idx, na.rm = TRUE)
  v <- as.numeric(t(idx) %*% V[names(idx), names(idx), drop = FALSE] %*% idx)
  se <- sqrt(v)
  c(est = b, se = se)
}

get_group_effects_clustered <- function(fit, base_term, group_var1, group_var2, vcov_mat, label = "effect") {
  beta <- coef(fit)
  lvls1 <- levels(model.frame(fit)[[group_var1]])
  lvls2 <- levels(model.frame(fit)[[group_var2]])
  ref1 <- lvls1[1]
  ref2 <- lvls2[1]

  out <- map_dfr(lvls1, function(lv1) {
    map_dfr(lvls2, function(lv2) {
      idx <- rep(0, length(beta))
      names(idx) <- names(beta)
      idx[base_term] <- 1

      if (lv1 != ref1) {
        t1 <- paste0(base_term, ":", group_var1, lv1)
        t2 <- paste0(group_var1, lv1, ":", base_term)
        if (t1 %in% names(beta)) idx[t1] <- 1
        if (t2 %in% names(beta)) idx[t2] <- 1
      }

      if (lv2 != ref2) {
        t1 <- paste0(base_term, ":", group_var2, lv2)
        t2 <- paste0(group_var2, lv2, ":", base_term)
        if (t1 %in% names(beta)) idx[t1] <- 1
        if (t2 %in% names(beta)) idx[t2] <- 1
      }

      if (lv1 != ref1 && lv2 != ref2) {
        patterns <- c(
          paste0(base_term, ":", group_var1, lv1, ":", group_var2, lv2),
          paste0(base_term, ":", group_var2, lv2, ":", group_var1, lv1),
          paste0(group_var1, lv1, ":", base_term, ":", group_var2, lv2),
          paste0(group_var2, lv2, ":", base_term, ":", group_var1, lv1)
        )
        for (pat in patterns) if (pat %in% names(beta)) idx[pat] <- 1
      }

      idx_use <- idx[names(idx) %in% names(beta)]
      st <- lincomb_stats(beta, vcov_mat, idx_use)
      z <- st["est"] / st["se"]
      p <- 2 * pnorm(-abs(z))
      ci <- st["est"] + c(-1, 1) * 1.96 * st["se"]

      tibble(
        !!group_var1 := lv1,
        !!group_var2 := lv2,
        log_est = unname(st["est"]),
        log_se = unname(st["se"]),
        !!label := exp(unname(st["est"])),
        low = exp(ci[1]),
        high = exp(ci[2]),
        p_value = p
      )
    })
  })

  out
}

claims <- read_excel(claims_path, sheet = claims_sheet) |>
  clean_names() |>
  mutate(from_date = parse_date_any(from_date)) |>
  ensure_care_setting() |>
  ensure_asthma_primary() |>
  ensure_age_group() |>
  mutate(
    municipality_desc = toupper(as.character(municipality_desc)),
    year = year(from_date),
    health_region = case_when(
      municipality_desc %in% san_juan_munis ~ "San Juan",
      municipality_desc %in% caguas_munis ~ "Caguas",
      TRUE ~ NA_character_
    )
  ) |>
  filter(year %in% years_all, !is.na(health_region))

claims_sj <- claims |>
  filter(
    health_region == "San Juan",
    asthma_primary == "asthma",
    care_setting %in% clinical_settings,
    !is.na(age_group)
  ) |>
  mutate(acute_flag = as.integer(care_setting %in% c("ED", "Inpatient")))

calendar <- tibble(
  date = seq.Date(as.Date(sprintf("%d-01-01", min(years_all))), as.Date(sprintf("%d-12-31", max(years_all))), by = "day")
) |>
  mutate(
    year = year(date),
    doy = yday(date),
    dow = wday(date, label = TRUE, week_start = 1),
    date_num = as.numeric(date)
  ) |>
  filter(year %in% years_all)

daily_asthma <- claims_sj |>
  group_by(date = from_date, age_group) |>
  summarise(
    total_claims = n(),
    acute_claims = sum(acute_flag == 1, na.rm = TRUE),
    nonacute_claims = total_claims - acute_claims,
    .groups = "drop"
  )

age_levels <- levels(droplevels(claims_sj$age_group))

daily_outcome <- tidyr::expand_grid(
  date = calendar$date,
  age_group = factor(age_levels, levels = age_levels)
) |>
  left_join(calendar, by = "date") |>
  left_join(daily_asthma, by = c("date", "age_group")) |>
  mutate(across(c(total_claims, acute_claims, nonacute_claims), ~ replace_na(.x, 0L)))

read_spores_sj_onefile <- function(path) {
  x <- read_excel(path) |> clean_names()
  date_col <- if ("date" %in% names(x)) "date" else names(x)[1]
  spore_col <- grep("esporas|spore", names(x), ignore.case = TRUE, value = TRUE)[1]
  pollen_col <- grep("polen|pollen", names(x), ignore.case = TRUE, value = TRUE)[1]
  x |>
    transmute(
      date = parse_date_any(.data[[date_col]]),
      spores_total_raw = suppressWarnings(as.numeric(.data[[spore_col]])),
      pollen_total_raw = if (!is.na(pollen_col)) suppressWarnings(as.numeric(.data[[pollen_col]])) else NA_real_
    ) |>
    filter(!is.na(date))
}

spores_raw <- map_dfr(spore_files_sj, read_spores_sj_onefile) |>
  group_by(date) |>
  summarise(
    spores_total_raw = safe_mean(spores_total_raw),
    pollen_total_raw = safe_mean(pollen_total_raw),
    .groups = "drop"
  ) |>
  arrange(date)

spores_daily <- tibble(date = calendar$date) |>
  left_join(spores_raw, by = "date") |>
  arrange(date)

spores_daily <- weekly_qc_impute_or_exclude(spores_daily, spores_total_raw, "spores_total", max_missing_days_week)
spores_daily <- weekly_qc_impute_or_exclude(spores_daily, pollen_total_raw, "pollen_total", max_missing_days_week)
spores_daily <- spores_daily |> mutate(spores_10k = spores_total / spore_scale_divisor)

weather_daily <- map_dfr(years_all, ~ pull_weather_sju(.x, sju_station_id))

daily <- daily_outcome |>
  left_join(spores_daily, by = "date") |>
  left_join(weather_daily, by = "date") |>
  mutate(
    pandemic_phase = case_when(
      year < 2020 ~ "pre2020",
      year == 2020 ~ "y2020",
      TRUE ~ "y2021plus"
    ),
    pandemic_phase = factor(pandemic_phase, levels = c("pre2020", "y2020", "y2021plus"))
  ) |>
  arrange(age_group, date) |>
  group_by(age_group) |>
  mutate(
    spores_10k_ma7 = zoo::rollapply(spores_10k, width = 7, FUN = mean, align = "right", fill = NA_real_, partial = FALSE)
  ) |>
  ungroup()

weather_covars_use <- intersect(weather_covars, names(daily))
k_doy <- min(25, max(6, floor(nrow(daily) / 400)))
knots_doy <- list(doy = c(0.5, 366.5))

daily_mod <- daily |>
  filter(!is.na(spores_10k_ma7)) |>
  filter(if_all(all_of(weather_covars_use), ~ is.finite(.x)))

f_nb <- as.formula(
  paste0(
    "total_claims ~ spores_10k_ma7 * pandemic_phase * age_group + ",
    "s(doy, bs='cc', k=", k_doy, ") + dow + s(date_num, bs='cr', k=10)",
    if (length(weather_covars_use) > 0) paste0(" + ", paste(weather_covars_use, collapse = " + ")) else ""
  )
)

fit_nb <- gam(
  f_nb,
  data = daily_mod,
  family = nb(),
  method = "REML",
  knots = knots_doy
)

daily_logit <- daily_mod |> filter(total_claims > 0)

f_logit <- as.formula(
  paste0(
    "cbind(acute_claims, nonacute_claims) ~ spores_10k_ma7 * pandemic_phase * age_group + ",
    "s(doy, bs='cc', k=", k_doy, ") + dow + s(date_num, bs='cr', k=10)",
    if (length(weather_covars_use) > 0) paste0(" + ", paste(weather_covars_use, collapse = " + ")) else ""
  )
)

fit_logit <- gam(
  f_logit,
  data = daily_logit,
  family = binomial(link = "logit"),
  method = "REML",
  knots = knots_doy
)

V_nb_cl <- sandwich::vcovCL(fit_nb, cluster = daily_mod$date)
V_logit_cl <- sandwich::vcovCL(fit_logit, cluster = daily_logit$date)

t6a_clustered <- get_group_effects_clustered(
  fit = fit_nb,
  base_term = "spores_10k_ma7",
  group_var1 = "pandemic_phase",
  group_var2 = "age_group",
  vcov_mat = V_nb_cl,
  label = "RR_clustered"
)

t6b_clustered <- get_group_effects_clustered(
  fit = fit_logit,
  base_term = "spores_10k_ma7",
  group_var1 = "pandemic_phase",
  group_var2 = "age_group",
  vcov_mat = V_logit_cl,
  label = "OR_clustered"
)

joint_p <- c(t6a_clustered$p_value, t6b_clustered$p_value)
joint_bh <- p.adjust(joint_p, method = "BH")

t6a_clustered <- t6a_clustered |>
  mutate(
    p_bh_48test_family = joint_bh[seq_len(n())],
    significant_nominal = p_value < 0.05,
    significant_bh = p_bh_48test_family < 0.05
  )

t6b_clustered <- t6b_clustered |>
  mutate(
    p_bh_48test_family = joint_bh[nrow(t6a_clustered) + seq_len(n())],
    significant_nominal = p_value < 0.05,
    significant_bh = p_bh_48test_family < 0.05
  )

write_csv(t6a_clustered, file.path(out_dir, "Table_S5_pandemic_phase_age_volume_clustered_by_date_draft.csv"))
write_csv(t6b_clustered, file.path(out_dir, "Table_S6_pandemic_phase_age_acute_fraction_clustered_by_date_draft.csv"))

interaction_tests <- tibble(
  model = c("NB_volume", "LOGIT_acute_fraction"),
  interaction = c("spores_10k_ma7:age_group", "spores_10k_ma7:age_group"),
  p_value = NA_real_,
  note = "Placeholder: add explicit nested-model interaction tests in a dedicated script"
)

write_csv(interaction_tests, file.path(out_dir, "interaction_tests_placeholder_draft.csv"))

message("Clustered Table 6 drafts written to: ", out_dir)
