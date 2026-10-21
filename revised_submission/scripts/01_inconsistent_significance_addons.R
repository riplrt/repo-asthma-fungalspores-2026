#!/usr/bin/env Rscript
############################################################
## Inconsistent-Significance Add-On Tables
## Focus:  Year-specific all-age estimates and pre-pandemic pediatric reversal
## Assembles reviewer-response tables from saved model output (no refitting):
##   - Table 2-style yearly RR/OR with FDR-adjusted p-values
##   - Pediatric (0-17) 2018 vs 2019 MA7 reversal
##   - 2018/2019 spore distributions (Table 1 add-on row)
##   - Key clustered interaction tests used in the revised framing
##
## Outputs (CSV) written to: R2_inconsistent_significance/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_inconsistent_significance/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Table 2 style yearly all-age estimates (age-adjusted series)
rr <- read_csv("statistical_results/SJ_region_spores_asthma_regression_by_year_age_v2_a_adj.csv", show_col_types = FALSE) |>
  transmute(
    year = year_v2_a_adj,
    outcome = "Encounter volume",
    effect = spores_rr_v2_a_adj,
    ci_low = spores_rr_low_v2_a_adj,
    ci_high = spores_rr_high_v2_a_adj,
    p_value = spores_p_v2_a_adj
  )

or <- read_csv("statistical_results/SJ_region_spores_asthma_logit_by_year_age_v2_a_adj.csv", show_col_types = FALSE) |>
  transmute(
    year = year_v2_a_adj,
    outcome = "ED/inpatient encounter fraction",
    effect = spores_or_v2_a_adj,
    ci_low = spores_or_low_v2_a_adj,
    ci_high = spores_or_high_v2_a_adj,
    p_value = spores_p_v2_a_adj
  )

tab2_like <- bind_rows(rr, or) |>
  arrange(outcome, year) |>
  mutate(
    q_bh_8 = p.adjust(p_value, method = "BH"),
    sig_nominal = p_value < 0.05,
    sig_bh = q_bh_8 < 0.05
  )

write_csv(tab2_like, file.path(out_dir, "Table_R2_table2_like_yearly_all_age_with_fdr_draft.csv"))

# Pediatric (0-17) reversal in pre-pandemic yearly MA7 model
pediatric <- read_csv("statistical_results/SJ_region_spores_asthma_logit_DAILY_lags_by_year_season_age_wlag_age.csv", show_col_types = FALSE) |>
  filter(
    year_wlag_age %in% c(2018, 2019),
    season_wlag_age == "All seasons",
    age_group_wlag_age == "0-17",
    exposure_var_wlag_age == "spores_10k_ma7_wlag_age"
  ) |>
  transmute(
    year = year_wlag_age,
    age_group = age_group_wlag_age,
    exposure = exposure_var_wlag_age,
    n_days = n_days_wlag_age,
    or = spores_or_wlag_age,
    ci_low = spores_or_low_wlag_age,
    ci_high = spores_or_high_wlag_age,
    p_value = spores_p_wlag_age
  )

write_csv(pediatric, file.path(out_dir, "Table_R2_pediatric_2018_2019_reversal_MA7_draft.csv"))

# 2018/2019 spore distributions for Table 1 add-on row
sp <- read_excel("data/Datos de esporas San Juan Jan 2018-2022.xlsx") |>
  transmute(
    date = as.Date(date),
    year = as.integer(format(date, "%Y")),
    spores_raw = as.numeric(spores_total),
    spores_10k = spores_raw / 10000
  ) |>
  filter(year %in% c(2018, 2019), is.finite(spores_10k))

sp_summary <- sp |>
  group_by(year) |>
  summarise(
    n_days = n(),
    min_spores_10k = min(spores_10k),
    p25_spores_10k = quantile(spores_10k, 0.25),
    median_spores_10k = median(spores_10k),
    mean_spores_10k = mean(spores_10k),
    p75_spores_10k = quantile(spores_10k, 0.75),
    p90_spores_10k = quantile(spores_10k, 0.90),
    max_spores_10k = max(spores_10k),
    .groups = "drop"
  )

write_csv(sp_summary, file.path(out_dir, "Table_R2_spore_distribution_2018_2019_draft.csv"))

# Pull interaction-test rows used in revised framing
int_tbl <- read_csv("Journal_Aerobiology_Submission_response/R2_interaction_tests/tables_drafts/interaction_tests_spore_age_period_draft.csv", show_col_types = FALSE) |>
  filter(
    spec == "age-specific seasonality/trend",
    outcome %in% c("Encounter volume (NB)", "ED/inpatient encounter fraction (logit)"),
    test %in% c("spore x age (all phases)", "spore x period (all ages)", "spore x age x period")
  ) |>
  dplyr::select(spec, outcome, test, df, chisq_clustered, p_clustered)

write_csv(int_tbl, file.path(out_dir, "Table_R2_key_interaction_tests_clustered_draft.csv"))

message("Inconsistent-significance add-on drafts written to: ", out_dir)
