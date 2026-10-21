#!/usr/bin/env Rscript
############################################################
## Draft Supplementary Tables S3-S6 (Multiplicity Inventory)
## Focus:  Test-family counts and full-strata tables from saved model output
## Builds S3 (multiplicity inventory), S4 (pre-pandemic daily severity MA7),
##   S5/S6 (pandemic phase x age: encounter volume / acute fraction)
## Reads existing CSVs in statistical_results/; no models are refit
##
## Outputs (CSV) written to: R2_multiplicity_clusteredSE/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_multiplicity_clusteredSE/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wald_p <- function(est, se) 2 * pnorm(-abs(est / se))

t2_rr <- read_csv("statistical_results/SJ_region_spores_asthma_regression_by_year_age_v2_a_adj.csv", show_col_types = FALSE)
t2_or <- read_csv("statistical_results/SJ_region_spores_asthma_logit_by_year_age_v2_a_adj.csv", show_col_types = FALSE)
t3 <- read_csv("statistical_results/SJ_region_spores_asthma_correlations_by_year_age_v2_a_adj.csv", show_col_types = FALSE)
t4_rr <- read_csv("statistical_results/SJ_region_spores_asthma_regression_WEEKLY_by_year_season_age_wlag_age.csv", show_col_types = FALSE)
t4_or <- read_csv("statistical_results/SJ_region_spores_asthma_logit_WEEKLY_by_year_season_age_wlag_age.csv", show_col_types = FALSE)
t5 <- read_csv("statistical_results/SJ_region_spores_asthma_logit_DAILY_lags_by_year_season_age_wlag_age.csv", show_col_types = FALSE)
t6_rr <- read_csv("statistical_results/pooled_daily_NB_spores_ma7_by_pandemic_phase_age_pphase_age.csv", show_col_types = FALSE)
t6_or <- read_csv("statistical_results/pooled_daily_LOGIT_spores_ma7_by_pandemic_phase_age_pphase_age.csv", show_col_types = FALSE)
t7 <- read_csv("statistical_results/case_crossover_ed_inpatient_results_cases.csv", show_col_types = FALSE)

count_t4 <- nrow(t4_rr |> filter(season_wlag_age == "All seasons")) +
  nrow(t4_or |> filter(season_wlag_age == "All seasons"))

count_t5_ma7 <- nrow(t5 |> filter(exposure_var_wlag_age == "spores_10k_ma7_wlag_age", season_wlag_age == "All seasons"))

count_t5_ma3_ma7 <- nrow(
  t5 |>
    filter(exposure_var_wlag_age %in% c("spores_10k_ma3_wlag_age", "spores_10k_ma7_wlag_age"),
           season_wlag_age == "All seasons")
)

s3 <- tibble(
  family_id = c("F1", "F2", "F3", "F4", "F4_alt", "F5", "F6"),
  manuscript_component = c(
    "Table 2 year-specific daily age-adjusted models",
    "Table 3/S1 age-year correlations",
    "Table 4 pre-pandemic weekly broad strata (all seasons)",
    "Table 5 pre-pandemic daily severity MA7 (all seasons)",
    "Table 5 pre-pandemic daily severity MA3 + MA7 (all seasons)",
    "Table 6a/6b pandemic phase x age x outcome",
    "Table 7 case-crossover"
  ),
  tests = c(
    nrow(t2_rr) + nrow(t2_or),
    nrow(t3) * 2,
    count_t4,
    count_t5_ma7,
    count_t5_ma3_ma7,
    nrow(t6_rr) + nrow(t6_or),
    nrow(t7)
  ),
  inferential_status = c(
    "Inferential",
    "Descriptive check (non-inferential)",
    "Inferential",
    "Inferential",
    "Inferential sensitivity framing",
    "Inferential",
    "Inferential"
  )
)

write_csv(s3, file.path(out_dir, "Table_S3_multiplicity_inventory_draft.csv"))

s4 <- t5 |>
  filter(exposure_var_wlag_age == "spores_10k_ma7_wlag_age", season_wlag_age == "All seasons") |>
  mutate(
    outcome = "acute_fraction",
    estimate_type = "OR",
    estimate = spores_or_wlag_age,
    ci_low = spores_or_low_wlag_age,
    ci_high = spores_or_high_wlag_age,
    p_nominal = spores_p_wlag_age,
    p_bh_within_family = p.adjust(spores_p_wlag_age, method = "BH")
  ) |>
  transmute(
    year = year_wlag_age,
    age_group = age_group_wlag_age,
    outcome,
    exposure = exposure_var_wlag_age,
    n = n_days_wlag_age,
    estimate_type,
    estimate,
    ci_low,
    ci_high,
    p_nominal,
    p_bh_within_family,
    significant_nominal = p_nominal < 0.05,
    significant_bh = p_bh_within_family < 0.05
  ) |>
  arrange(year, age_group)

write_csv(s4, file.path(out_dir, "Table_S4_prepandemic_daily_severity_MA7_all_strata_draft.csv"))

s5_base <- t6_rr |>
  mutate(
    outcome = "volume",
    estimate_type = "RR",
    estimate = RR_pphase_age,
    p_nominal = wald_p(log_est_pphase_age, log_se_pphase_age)
  )

s6_base <- t6_or |>
  mutate(
    outcome = "acute_fraction",
    estimate_type = "OR",
    estimate = OR_pphase_age,
    p_nominal = wald_p(log_est_pphase_age, log_se_pphase_age)
  )

joint_p <- c(s5_base$p_nominal, s6_base$p_nominal)
joint_bh <- p.adjust(joint_p, method = "BH")

s5 <- s5_base |>
  mutate(p_bh_48test_family = joint_bh[seq_len(n())]) |>
  transmute(
    pandemic_phase = pandemic_phase_pphase_age,
    age_group = age_group_pphase_age,
    outcome,
    estimate_type,
    estimate,
    ci_low = low_pphase_age,
    ci_high = high_pphase_age,
    p_nominal,
    p_bh_48test_family,
    significant_nominal = p_nominal < 0.05,
    significant_bh = p_bh_48test_family < 0.05
  ) |>
  arrange(pandemic_phase, age_group)

s6 <- s6_base |>
  mutate(p_bh_48test_family = joint_bh[nrow(s5_base) + seq_len(n())]) |>
  transmute(
    pandemic_phase = pandemic_phase_pphase_age,
    age_group = age_group_pphase_age,
    outcome,
    estimate_type,
    estimate,
    ci_low = low_pphase_age,
    ci_high = high_pphase_age,
    p_nominal,
    p_bh_48test_family,
    significant_nominal = p_nominal < 0.05,
    significant_bh = p_bh_48test_family < 0.05
  ) |>
  arrange(pandemic_phase, age_group)

write_csv(s5, file.path(out_dir, "Table_S5_pandemic_phase_age_volume_all_strata_draft.csv"))
write_csv(s6, file.path(out_dir, "Table_S6_pandemic_phase_age_acute_fraction_all_strata_draft.csv"))

message("Draft tables written to: ", out_dir)
