#!/usr/bin/env Rscript
############################################################
## Inverse-Association Decomposition (Age 65+, Pre-Pandemic)
## Focus:  Weekly ED/inpatient vs other asthma encounters, 2018-2019
## Tests whether the inverse spore-fraction association reflects numerator or
##   denominator changes by modeling the component counts separately
##   - Weekly NB GAMs (total, acute, other) + binomial GAM (acute fraction)
##   - Key term: spores_10k x year
##   - Controls: cyclic week-of-year spline + weekly mean weather
## Also (2018): spores x 65+ age-subgroup interaction; J44 (COPD) sensitivity
##
## Outputs (CSV) written to: R2_inverse_associations/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(mgcv)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_inverse_associations/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

san_juan_munis <- c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA")
clinical_settings <- c("ED", "Inpatient", "Urgent", "Outpatient")

claims <- read_csv("data/2022_asma_From_date_ver 3.csv", show_col_types = FALSE) |>
  mutate(
    from_date = suppressWarnings(mdy(FROM_DATE)),
    year = year(from_date),
    muni = toupper(MUNICIPALITY_DESC),
    icd1 = toupper(trimws(ICD_DIAG_01_PRIMARY)),
    care_setting = case_when(
      grepl("ed|emergency", tolower(HCG_SETTING)) | grepl("\\bed\\b|\\ber\\b", tolower(HCG_MR_LINE_DESC)) ~ "ED",
      grepl("inpatient|hospital", tolower(HCG_SETTING)) | grepl("inpatient|hospital", tolower(HCG_MR_LINE_DESC)) ~ "Inpatient",
      grepl("urgent", tolower(HCG_SETTING)) | grepl("urgent", tolower(HCG_MR_LINE_DESC)) ~ "Urgent",
      grepl("office|clinic|outpatient", tolower(HCG_SETTING)) | grepl("office|clinic", tolower(HCG_MR_LINE_DESC)) ~ "Outpatient",
      grepl("dme|supplies", tolower(HCG_MR_LINE_DESC)) ~ "DME",
      TRUE ~ "Other"
    ),
    age_num = suppressWarnings(as.numeric(gsub("\\+", "", as.character(Age)))),
    age65grp = cut(age_num, breaks = c(65, 74, 84, Inf), labels = c("65-74", "75-84", "85+"), right = TRUE, include.lowest = TRUE)
  ) |>
  filter(!is.na(from_date), year %in% c(2018, 2019), muni %in% san_juan_munis)

sp <- read_excel("data/Datos de esporas San Juan Jan 2018-2022.xlsx") |>
  transmute(date = as.Date(date), spores_10k = as.numeric(spores_total) / 10000)

wx <- read_csv("data/sju_weather_2018_2019.csv", show_col_types = FALSE) |>
  transmute(
    date = as.Date(date),
    tmean_c = as.numeric((tmax_c + tmin_c) / 2),
    rh_mean = as.numeric(rh_mean),
    wind_ms = as.numeric(wind_ms),
    prcp_mm = as.numeric(prcp_mm),
    dewpoint_c = as.numeric(dewpoint_c)
  )

asthma65 <- claims |>
  filter(grepl("^J45|^J46", icd1), care_setting %in% clinical_settings, age_num >= 65) |>
  mutate(acute = as.integer(care_setting %in% c("ED", "Inpatient")))

weekly <- asthma65 |>
  group_by(date = from_date) |>
  summarise(total = n(), acute = sum(acute), other = total - acute, .groups = "drop") |>
  right_join(tibble(date = seq.Date(as.Date("2018-01-01"), as.Date("2019-12-31"), by = "day")), by = "date") |>
  mutate(across(c(total, acute, other), ~ tidyr::replace_na(.x, 0L))) |>
  mutate(year = year(date), week = isoweek(date)) |>
  left_join(sp, by = "date") |>
  left_join(wx, by = "date") |>
  group_by(year, week) |>
  summarise(
    total = sum(total),
    acute = sum(acute),
    other = sum(other),
    spores_10k = mean(spores_10k, na.rm = TRUE),
    tmean_c = mean(tmean_c, na.rm = TRUE),
    rh_mean = mean(rh_mean, na.rm = TRUE),
    wind_ms = mean(wind_ms, na.rm = TRUE),
    prcp_mm = mean(prcp_mm, na.rm = TRUE),
    dewpoint_c = mean(dewpoint_c, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(week <= 52) |>
  mutate(year = factor(year)) |>
  filter(if_all(c(spores_10k, tmean_c, rh_mean, wind_ms, prcp_mm, dewpoint_c), is.finite))

k <- min(12, max(4, floor(nrow(weekly) / 12)))
kn <- list(week = c(0.5, 52.5))

fit_total <- gam(total ~ spores_10k * year + s(week, bs = "cc", k = k) + tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
                 family = nb(), data = weekly, method = "REML", knots = kn)
fit_acute <- gam(acute ~ spores_10k * year + s(week, bs = "cc", k = k) + tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
                 family = nb(), data = weekly, method = "REML", knots = kn)
fit_other <- gam(other ~ spores_10k * year + s(week, bs = "cc", k = k) + tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
                 family = nb(), data = weekly, method = "REML", knots = kn)
fit_frac <- gam(cbind(acute, other) ~ spores_10k * year + s(week, bs = "cc", k = k) + tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
                family = binomial(), data = weekly, method = "REML", knots = kn)

extract_year <- function(fit, year_level = "2018") {
  b <- coef(fit)
  V <- vcov(fit)
  idx <- rep(0, length(b)); names(idx) <- names(b)
  idx["spores_10k"] <- 1
  if (year_level != levels(weekly$year)[1]) {
    term <- paste0("spores_10k:year", year_level)
    term2 <- paste0("year", year_level, ":spores_10k")
    if (term %in% names(b)) idx[term] <- 1
    if (term2 %in% names(b)) idx[term2] <- 1
  }
  u <- names(idx)[idx != 0]
  est <- sum(b[u] * idx[u])
  se <- sqrt(as.numeric(t(idx[u]) %*% V[u, u, drop = FALSE] %*% idx[u]))
  tibble(
    estimate = exp(est),
    ci_low = exp(est - 1.96 * se),
    ci_high = exp(est + 1.96 * se),
    p_value = 2 * pnorm(-abs(est / se)),
    year = year_level
  )
}

tbl_sep <- bind_rows(
  extract_year(fit_total, "2018") |> mutate(model = "Total >=65 (NB)", outcome = "total_encounters"),
  extract_year(fit_acute, "2018") |> mutate(model = "ED/Inpatient >=65 (NB)", outcome = "acute_encounters"),
  extract_year(fit_other, "2018") |> mutate(model = "Other encounters >=65 (NB)", outcome = "other_encounters"),
  extract_year(fit_frac, "2018") |> mutate(model = "Acute fraction >=65 (logit)", outcome = "acute_fraction")
)

write_csv(tbl_sep, file.path(out_dir, "inverse_assoc_separate_counts_age65_2018_draft.csv"))

weekly_mix <- claims |>
  filter(!is.na(from_date), year(from_date) == 2018, muni %in% san_juan_munis,
         grepl("^J45|^J46", icd1), care_setting %in% clinical_settings,
         age_num >= 65, !is.na(age65grp)) |>
  mutate(acute = as.integer(care_setting %in% c("ED", "Inpatient"))) |>
  group_by(date = from_date, age65grp) |>
  summarise(acute = sum(acute), other = n() - acute, .groups = "drop") |>
  tidyr::complete(date = seq.Date(as.Date("2018-01-01"), as.Date("2018-12-31"), by = "day"),
                  age65grp, fill = list(acute = 0L, other = 0L)) |>
  mutate(week = isoweek(date)) |>
  left_join(sp, by = "date") |>
  left_join(wx, by = "date") |>
  group_by(age65grp, week) |>
  summarise(
    acute = sum(acute),
    other = sum(other),
    spores_10k = mean(spores_10k, na.rm = TRUE),
    tmean_c = mean(tmean_c, na.rm = TRUE),
    rh_mean = mean(rh_mean, na.rm = TRUE),
    wind_ms = mean(wind_ms, na.rm = TRUE),
    prcp_mm = mean(prcp_mm, na.rm = TRUE),
    dewpoint_c = mean(dewpoint_c, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(week <= 52) |>
  filter(if_all(c(spores_10k, tmean_c, rh_mean, wind_ms, prcp_mm, dewpoint_c), is.finite))

fit_mix <- gam(cbind(acute, other) ~ spores_10k * age65grp + s(week, bs = "cc", k = 8) +
                 tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
               family = binomial(), data = weekly_mix, method = "REML",
               knots = list(week = c(0.5, 52.5)))

mix_coefs <- as.data.frame(summary(fit_mix)$p.table) |>
  tibble::rownames_to_column("term") |>
  filter(grepl("spores_10k", term))

write_csv(mix_coefs, file.path(out_dir, "inverse_assoc_age65_mix_terms_2018_draft.csv"))

weekly_j44 <- claims |>
  filter(!is.na(from_date), year(from_date) == 2018, muni %in% san_juan_munis,
         grepl("^J45|^J46|^J44", icd1), care_setting %in% clinical_settings,
         age_num >= 65) |>
  mutate(acute = as.integer(care_setting %in% c("ED", "Inpatient"))) |>
  group_by(date = from_date) |>
  summarise(acute = sum(acute), other = n() - acute, .groups = "drop") |>
  right_join(tibble(date = seq.Date(as.Date("2018-01-01"), as.Date("2018-12-31"), by = "day")), by = "date") |>
  mutate(across(c(acute, other), ~ tidyr::replace_na(.x, 0L))) |>
  mutate(week = isoweek(date)) |>
  left_join(sp, by = "date") |>
  left_join(wx, by = "date") |>
  group_by(week) |>
  summarise(
    acute = sum(acute),
    other = sum(other),
    spores_10k = mean(spores_10k, na.rm = TRUE),
    tmean_c = mean(tmean_c, na.rm = TRUE),
    rh_mean = mean(rh_mean, na.rm = TRUE),
    wind_ms = mean(wind_ms, na.rm = TRUE),
    prcp_mm = mean(prcp_mm, na.rm = TRUE),
    dewpoint_c = mean(dewpoint_c, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(week <= 52) |>
  filter(if_all(c(spores_10k, tmean_c, rh_mean, wind_ms, prcp_mm, dewpoint_c), is.finite))

fit_j44 <- gam(cbind(acute, other) ~ spores_10k + s(week, bs = "cc", k = 8) +
                 tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c,
               family = binomial(), data = weekly_j44, method = "REML",
               knots = list(week = c(0.5, 52.5)))

j44_coef <- as.data.frame(summary(fit_j44)$p.table) |>
  tibble::rownames_to_column("term") |>
  filter(term == "spores_10k")

write_csv(j44_coef, file.path(out_dir, "inverse_assoc_j44_sensitivity_age65_2018_draft.csv"))

message("Inverse-association draft tables written to: ", out_dir)
