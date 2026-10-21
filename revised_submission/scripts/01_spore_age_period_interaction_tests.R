#!/usr/bin/env Rscript
############################################################
## Spore x Age and Spore x Period Interaction Tests
## Focus:  Phase x age MA7 models behind Tables 6a/6b
## Formal heterogeneity tests with standard errors clustered by date
##   - Daily NB (encounter volume) + binomial logit (ED/inpatient fraction)
##   - One spore slope per phase x age cell; joint Wald tests for
##     spore x age (7 df), spore x period (2 df), three-way (14 df)
##
## Departures from asthma_spores_pandemic_phase_age_adjusted.R, and why:
##   - Exposure: San Juan station series ("Datos de esporas San Juan Jan 2018-2022.xlsx").
##     The original reader selected the first "esporas" column of the corrected file,
##     which is esporas_caguas (Caguas station, from 2019-08-24 only).
##   - Weather: local data/sju_weather_2018_2021.csv (rnoaa unavailable), so the
##     analysis period is 2018-2021 and the "y2021plus" phase contains 2021 only.
##   - Seasonality/trend: unpenalized terms (3 Fourier pairs for day of year,
##     ns(date, 8) for trend) instead of penalized mgcv smooths, so that the
##     date-clustered sandwich variance is well defined.
##
## Outputs (CSV) written to: R2_interaction_tests/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(zoo)
  library(MASS)
  library(splines)
  library(sandwich)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_interaction_tests/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

san_juan_munis <- c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA")
clinical_settings <- c("ED", "Inpatient", "Outpatient", "Urgent")
age_breaks <- c(0, 4, 11, 17, 39, 64, 74, 84, Inf)
age_labels <- c("0-4", "5-11", "12-17", "18-39", "40-64", "65-74", "75-84", "85+")
phase_levels <- c("pre2020", "y2020", "y2021plus")
years_all <- 2018:2021

care_setting_fun <- function(hcg_setting, hcg_mr_line_desc) {
  s <- tolower(as.character(hcg_setting))
  l <- tolower(as.character(hcg_mr_line_desc))
  if (grepl("ed|emergency", s) || grepl("\\bed\\b|\\ber\\b", l)) return("ED")
  if (grepl("inpatient|hospital", s) || grepl("inpatient|hospital", l)) return("Inpatient")
  if (grepl("urgent", s) || grepl("urgent", l)) return("Urgent")
  if (grepl("office|clinic|outpatient", s) || grepl("office|clinic", l)) return("Outpatient")
  "Other"
}

# ---- Outcomes: day x age-group panel -------------------------------------
claims <- read_csv("data/2022_asma_From_date_ver 3.csv", show_col_types = FALSE) |>
  mutate(
    date = suppressWarnings(mdy(FROM_DATE)),
    municipality = toupper(trimws(MUNICIPALITY_DESC)),
    icd1 = toupper(trimws(ICD_DIAG_01_PRIMARY)),
    care_setting = mapply(care_setting_fun, HCG_SETTING, HCG_MR_LINE_DESC),
    age_num = suppressWarnings(as.numeric(gsub("\\+", "", as.character(Age)))),
    age_group = cut(age_num, breaks = age_breaks, labels = age_labels,
                    right = TRUE, include.lowest = TRUE)
  ) |>
  filter(!is.na(date), year(date) %in% years_all, municipality %in% san_juan_munis,
         grepl("^J45|^J46", icd1), care_setting %in% clinical_settings, !is.na(age_group))

panel <- claims |>
  group_by(date, age_group) |>
  summarise(
    total = n(),
    acute = sum(care_setting %in% c("ED", "Inpatient")),
    .groups = "drop"
  ) |>
  complete(
    date = seq.Date(as.Date("2018-01-01"), as.Date("2021-12-31"), by = "day"),
    age_group = factor(age_labels, levels = age_labels),
    fill = list(total = 0L, acute = 0L)
  ) |>
  mutate(other = total - acute)

# ---- Exposure: weekly QC rule, then right-aligned MA7 ---------------------
spores <- read_excel("data/Datos de esporas San Juan Jan 2018-2022.xlsx") |>
  transmute(date = as.Date(date), spores_raw = as.numeric(spores_total)) |>
  right_join(tibble(date = seq.Date(as.Date("2018-01-01"), as.Date("2021-12-31"), by = "day")),
             by = "date") |>
  arrange(date) |>
  mutate(week_start = floor_date(date, "week", week_start = 1)) |>
  group_by(week_start) |>
  mutate(
    n_miss = sum(is.na(spores_raw)),
    week_mean = if (all(is.na(spores_raw))) NA_real_ else mean(spores_raw, na.rm = TRUE),
    spores_qc = case_when(
      n_miss > 3 | is.na(week_mean) ~ NA_real_,
      is.na(spores_raw) ~ week_mean,
      TRUE ~ spores_raw
    )
  ) |>
  ungroup() |>
  mutate(spores_10k_ma7 = rollapply(spores_qc / 10000, 7, mean, align = "right",
                                    fill = NA_real_, partial = FALSE)) |>
  dplyr::select(date, spores_10k_ma7)

weather <- read_csv("data/sju_weather_2018_2021.csv", show_col_types = FALSE) |>
  transmute(
    date = as.Date(date),
    tmean_c = (tmax_c + tmin_c) / 2,
    rh_mean, wind_ms, prcp_mm, dewpoint_c
  )

dat <- panel |>
  left_join(spores, by = "date") |>
  left_join(weather, by = "date") |>
  mutate(
    phase = factor(case_when(year(date) < 2020 ~ "pre2020",
                             year(date) == 2020 ~ "y2020",
                             TRUE ~ "y2021plus"), levels = phase_levels),
    pa = interaction(phase, age_group, sep = "|", lex.order = TRUE),
    doy = yday(date),
    dow = wday(date, label = TRUE, week_start = 1),
    date_num = as.numeric(date),
    date_id = as.integer(factor(date))
  ) |>
  filter(!is.na(spores_10k_ma7),
         if_all(c(tmean_c, rh_mean, wind_ms, prcp_mm, dewpoint_c), is.finite))

dat_logit <- dat |> filter(total > 0)

# ---- Models: one spore slope per phase x age cell -------------------------
rhs <- paste(
  "phase * age_group + spores_10k_ma7:pa",
  "sin(2*pi*doy/365.25) + cos(2*pi*doy/365.25)",
  "sin(4*pi*doy/365.25) + cos(4*pi*doy/365.25)",
  "sin(6*pi*doy/365.25) + cos(6*pi*doy/365.25)",
  "ns(date_num, 8) + dow + tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c",
  sep = " + "
)

# Preferred specification: seasonality and trend allowed to differ by age group.
# A shared seasonal curve forces age-specific seasonal patterns (e.g. school-year
# effects in children) into the age-specific spore slopes.
rhs_age <- paste(
  "phase * age_group + spores_10k_ma7:pa",
  paste0("age_group:(sin(2*pi*doy/365.25) + cos(2*pi*doy/365.25) + ",
         "sin(4*pi*doy/365.25) + cos(4*pi*doy/365.25) + ",
         "sin(6*pi*doy/365.25) + cos(6*pi*doy/365.25))"),
  "age_group:ns(date_num, 8) + dow + tmean_c + rh_mean + wind_ms + prcp_mm + dewpoint_c",
  sep = " + "
)

fit_nb <- glm.nb(as.formula(paste("total ~", rhs)), data = dat)
fit_lg <- glm(as.formula(paste("cbind(acute, other) ~", rhs)), family = binomial(), data = dat_logit)
fit_nb_a <- glm.nb(as.formula(paste("total ~", rhs_age)), data = dat)
fit_lg_a <- glm(as.formula(paste("cbind(acute, other) ~", rhs_age)), family = binomial(), data = dat_logit)

slope_names <- function(fit) grep("^spores_10k_ma7:pa", names(coef(fit)), value = TRUE)

cell_table <- function(fit, V, outcome) {
  nm <- slope_names(fit)
  b <- coef(fit)[nm]
  se <- sqrt(diag(V)[nm])
  cells <- sub("^spores_10k_ma7:pa", "", nm)
  tibble(
    outcome = outcome,
    phase = sub("\\|.*$", "", cells),
    age_group = sub("^.*\\|", "", cells),
    estimate = exp(b),
    ci_low = exp(b - 1.96 * se),
    ci_high = exp(b + 1.96 * se),
    p_value = 2 * pnorm(-abs(b / se))
  )
}

# Joint Wald test that a set of linear contrasts of the slopes equals zero
wald <- function(fit, V, C) {
  nm <- slope_names(fit)
  b <- coef(fit)[nm]
  Vs <- V[nm, nm]
  est <- C %*% b
  stat <- as.numeric(t(est) %*% solve(C %*% Vs %*% t(C)) %*% est)
  c(chisq = stat, df = nrow(C), p = pchisq(stat, nrow(C), lower.tail = FALSE))
}

cells <- expand_grid(phase = phase_levels, age_group = age_labels) |>
  mutate(cell = paste(phase, age_group, sep = "|"))

# Contrast rows comparing each level to the first level within groups
diff_contrasts <- function(group_var) {
  rows <- list()
  for (g in unique(cells[[group_var]])) {
    idx <- which(cells[[group_var]] == g)
    ref <- idx[1]
    for (j in idx[-1]) {
      r <- numeric(nrow(cells)); r[j] <- 1; r[ref] <- -1
      rows[[length(rows) + 1]] <- r
    }
  }
  do.call(rbind, rows)
}

# Spore x age within one phase (7 df) and spore x period within one age group (2 df)
within_phase <- function(ph) {
  idx <- which(cells$phase == ph)
  t(sapply(idx[-1], function(j) { r <- numeric(nrow(cells)); r[j] <- 1; r[idx[1]] <- -1; r }))
}
within_age <- function(ag) {
  idx <- which(cells$age_group == ag)
  t(sapply(idx[-1], function(j) { r <- numeric(nrow(cells)); r[j] <- 1; r[idx[1]] <- -1; r }))
}

# Three-way: does the age pattern of slopes differ between phases? (14 df)
three_way <- function() {
  rows <- list()
  for (ph in phase_levels[-1]) for (ag in age_labels[-1]) {
    r <- numeric(nrow(cells))
    r[cells$cell == paste(ph, ag, sep = "|")] <- 1
    r[cells$cell == paste(ph, age_labels[1], sep = "|")] <- -1
    r[cells$cell == paste(phase_levels[1], ag, sep = "|")] <- -1
    r[cells$cell == paste(phase_levels[1], age_labels[1], sep = "|")] <- 1
    rows[[length(rows) + 1]] <- r
  }
  do.call(rbind, rows)
}

tests <- list(
  list(test = "spore x age (all phases)", C = diff_contrasts("phase")),
  list(test = "spore x period (all ages)", C = diff_contrasts("age_group")),
  list(test = "spore x age x period", C = three_way())
)
tests <- c(
  tests,
  lapply(phase_levels, function(ph) list(test = paste0("spore x age within ", ph), C = within_phase(ph))),
  lapply(age_labels, function(ag) list(test = paste0("spore x period within ", ag), C = within_age(ag)))
)

stopifnot(nobs(fit_nb) == nrow(dat), nrow(model.frame(fit_lg)) == nrow(dat_logit),
          nobs(fit_nb_a) == nrow(dat), nrow(model.frame(fit_lg_a)) == nrow(dat_logit))

run_tests <- function(fit_vol, fit_frac, spec) {
  fl <- list(
    list(outcome = "Encounter volume (NB)", fit = fit_vol, cluster = dat$date_id),
    list(outcome = "ED/inpatient encounter fraction (logit)", fit = fit_frac, cluster = dat_logit$date_id)
  )
  # Coefficient order in the model must match the cells table
  for (f in fl) stopifnot(identical(sub("^spores_10k_ma7:pa", "", slope_names(f$fit)), cells$cell))
  for (i in seq_along(fl)) {
    fl[[i]]$V_model <- vcov(fl[[i]]$fit)
    fl[[i]]$V_cl <- vcovCL(fl[[i]]$fit, cluster = fl[[i]]$cluster, type = "HC0")
  }

  tt <- bind_rows(lapply(fl, function(f) bind_rows(lapply(tests, function(t) {
    cl <- wald(f$fit, f$V_cl, t$C)
    md <- wald(f$fit, f$V_model, t$C)
    tibble(spec = spec, outcome = f$outcome, test = t$test, df = cl[["df"]],
           chisq_clustered = cl[["chisq"]], p_clustered = cl[["p"]],
           chisq_model = md[["chisq"]], p_model = md[["p"]])
  }))))

  # BH applied across all 48 phase x age x outcome slopes
  ct <- bind_rows(lapply(fl, function(f) {
    cl <- cell_table(f$fit, f$V_cl, f$outcome)
    md <- cell_table(f$fit, f$V_model, f$outcome) |>
      dplyr::select(outcome, phase, age_group, ci_low_model = ci_low, ci_high_model = ci_high, p_model = p_value)
    left_join(cl, md, by = c("outcome", "phase", "age_group"))
  })) |>
    rename(p_clustered = p_value, ci_low_clustered = ci_low, ci_high_clustered = ci_high) |>
    mutate(spec = spec,
           p_bh_48_clustered = p.adjust(p_clustered, method = "BH"),
           p_bh_48_model = p.adjust(p_model, method = "BH"))

  list(tests = tt, cells = ct)
}

r_common <- run_tests(fit_nb, fit_lg, "common seasonality/trend")
r_age <- run_tests(fit_nb_a, fit_lg_a, "age-specific seasonality/trend")

test_tbl <- bind_rows(r_age$tests, r_common$tests)
cell_tbl <- bind_rows(r_age$cells, r_common$cells)

fit_tbl <- tibble(
  spec = rep(c("common seasonality/trend", "age-specific seasonality/trend"), each = 2),
  outcome = rep(c("Encounter volume (NB)", "ED/inpatient encounter fraction (logit)"), 2),
  aic = c(AIC(fit_nb), AIC(fit_lg), AIC(fit_nb_a), AIC(fit_lg_a))
)

sample_tbl <- tibble(
  model = c("Encounter volume (NB)", "ED/inpatient encounter fraction (logit)"),
  n_rows = c(nrow(dat), nrow(dat_logit)),
  n_dates = c(n_distinct(dat$date), n_distinct(dat_logit$date)),
  first_date = c(min(dat$date), min(dat_logit$date)),
  last_date = c(max(dat$date), max(dat_logit$date))
)

days_by_phase <- dat |> distinct(date, phase) |> count(phase, name = "n_days_in_model")

write_csv(test_tbl, file.path(out_dir, "interaction_tests_spore_age_period_draft.csv"))
write_csv(cell_tbl, file.path(out_dir, "phase_age_slopes_clustered_by_date_draft.csv"))
write_csv(fit_tbl, file.path(out_dir, "interaction_model_aic_draft.csv"))
write_csv(sample_tbl, file.path(out_dir, "interaction_model_sample_draft.csv"))
write_csv(days_by_phase, file.path(out_dir, "interaction_model_days_by_phase_draft.csv"))

message("Interaction test drafts written to: ", out_dir)
