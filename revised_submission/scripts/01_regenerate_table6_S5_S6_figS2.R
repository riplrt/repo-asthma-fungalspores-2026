#!/usr/bin/env Rscript
############################################################
## Regenerate Tables 6a/6b, S5/S6 and Figure S2
## Focus:  Pandemic phase x age MA7 spore slopes, San Juan spore series
## Age-specific seasonality/trend; date-clustered SEs; FDR column
##
## Input: slopes from 01_spore_age_period_interaction_tests.R (run that first).
##   The original Tables 6a/6b used the Caguas spore column of the corrected
##   spore file; see that script's header.
## Also: comparison with the original (Caguas-series) Table 6 and
##   significance counts for the Table 6 family
##
## Outputs (CSV, PNG) written to: R2_table6_regenerated/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_file <- "Journal_Aerobiology_Submission_response/R2_interaction_tests/tables_drafts/phase_age_slopes_clustered_by_date_draft.csv"
out_dir <- "Journal_Aerobiology_Submission_response/R2_table6_regenerated/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

age_labels <- c("0-4", "5-11", "12-17", "18-39", "40-64", "65-74", "75-84", "85+")
phase_labels <- c(pre2020 = "Pre-2020 (2018-2019)", y2020 = "2020", y2021plus = "2021")
outcome_labels <- c(
  "Encounter volume (NB)" = "Encounter volume (RR)",
  "ED/inpatient encounter fraction (logit)" = "ED/inpatient encounter fraction (OR)"
)

slopes <- read_csv(in_file, show_col_types = FALSE) |>
  mutate(
    age_group = factor(age_group, levels = age_labels),
    phase_lab = factor(phase_labels[phase], levels = phase_labels),
    outcome_lab = factor(outcome_labels[outcome], levels = outcome_labels)
  )

primary <- slopes |> filter(spec == "age-specific seasonality/trend")
sensitivity <- slopes |> filter(spec == "common seasonality/trend")

fmt_ci <- function(est, lo, hi) sprintf("%.3f (%.3f-%.3f)", est, lo, hi)
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

# ---- Main-text Tables 6a/6b: all strata, clustered CIs, FDR column --------
main_table <- function(outcome_name) {
  primary |>
    filter(outcome == outcome_name) |>
    arrange(phase_lab, age_group) |>
    transmute(
      `Pandemic phase` = phase_lab,
      `Age group (years)` = age_group,
      `Estimate (95% CI)` = fmt_ci(estimate, ci_low_clustered, ci_high_clustered),
      `p` = fmt_p(p_clustered),
      `FDR q` = fmt_p(p_bh_48_clustered)
    )
}

write_csv(main_table("Encounter volume (NB)"), file.path(out_dir, "Table_6a_encounter_volume_RR.csv"))
write_csv(main_table("ED/inpatient encounter fraction (logit)"), file.path(out_dir, "Table_6b_ed_inpatient_fraction_OR.csv"))

# ---- Supplementary Tables S5/S6: full numeric detail + sensitivity -------
supp_table <- function(outcome_name) {
  sens <- sensitivity |>
    filter(outcome == outcome_name) |>
    dplyr::select(phase, age_group,
                  sens_estimate = estimate, sens_ci_low = ci_low_clustered,
                  sens_ci_high = ci_high_clustered, sens_q = p_bh_48_clustered)
  primary |>
    filter(outcome == outcome_name) |>
    left_join(sens, by = c("phase", "age_group")) |>
    arrange(phase_lab, age_group) |>
    transmute(
      phase = phase_lab, age_group, estimate,
      ci_low_clustered, ci_high_clustered, p_clustered, q_bh48_clustered = p_bh_48_clustered,
      ci_low_model, ci_high_model, p_model, q_bh48_model = p_bh_48_model,
      sens_common_season_estimate = sens_estimate,
      sens_common_season_ci_low = sens_ci_low,
      sens_common_season_ci_high = sens_ci_high,
      sens_common_season_q = sens_q
    )
}

write_csv(supp_table("Encounter volume (NB)"), file.path(out_dir, "Table_S5_encounter_volume_all_strata.csv"))
write_csv(supp_table("ED/inpatient encounter fraction (logit)"), file.path(out_dir, "Table_S6_ed_inpatient_fraction_all_strata.csv"))

# ---- Figure S2 ------------------------------------------------------------
fig_dat <- primary |>
  mutate(
    log_effect = log(estimate),
    marker = case_when(p_bh_48_clustered < 0.05 ~ "*", p_clustered < 0.05 ~ "\u2020", TRUE ~ ""),
    label = paste0(sprintf("%.2f", estimate), marker)
  )

lim <- max(abs(fig_dat$log_effect))

fig <- ggplot(fig_dat, aes(x = phase_lab, y = age_group, fill = log_effect)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3.2) +
  facet_wrap(~outcome_lab) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim), name = "log(effect)") +
  labs(
    x = "Pandemic phase", y = "Age group (years)",
    title = "Phase- and age-specific MA7 spore associations",
    caption = paste0("Effect per 10,000 spores (7-day moving average). Date-clustered SEs; ",
                     "age-specific seasonality and trend.\n* FDR q < 0.05 across 48 estimates; ",
                     "\u2020 nominal p < 0.05 only.")
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.caption = element_text(hjust = 0))

ggsave(file.path(out_dir, "Figure_S2_phase_age_MA7.png"), fig, width = 10, height = 5.5, dpi = 300, bg = "white")

# ---- Comparison with the original (Caguas-series) Table 6 -----------------
wald_p <- function(b, se) 2 * pnorm(-abs(b / se))
orig <- bind_rows(
  read_csv("statistical_results/pooled_daily_NB_spores_ma7_by_pandemic_phase_age_pphase_age.csv", show_col_types = FALSE) |>
    transmute(outcome = "Encounter volume (NB)", phase = pandemic_phase_pphase_age, age_group = age_group_pphase_age,
              orig_estimate = RR_pphase_age, orig_p = wald_p(log_est_pphase_age, log_se_pphase_age)),
  read_csv("statistical_results/pooled_daily_LOGIT_spores_ma7_by_pandemic_phase_age_pphase_age.csv", show_col_types = FALSE) |>
    transmute(outcome = "ED/inpatient encounter fraction (logit)", phase = pandemic_phase_pphase_age, age_group = age_group_pphase_age,
              orig_estimate = OR_pphase_age, orig_p = wald_p(log_est_pphase_age, log_se_pphase_age))
) |>
  mutate(orig_q = p.adjust(orig_p, method = "BH"))

comparison <- primary |>
  mutate(age_group = as.character(age_group)) |>
  dplyr::select(outcome, phase, age_group, new_estimate = estimate, new_p = p_clustered, new_q = p_bh_48_clustered) |>
  left_join(orig, by = c("outcome", "phase", "age_group")) |>
  mutate(
    orig_status = case_when(orig_q < 0.05 ~ "FDR", orig_p < 0.05 ~ "nominal", TRUE ~ "ns"),
    new_status = case_when(new_q < 0.05 ~ "FDR", new_p < 0.05 ~ "nominal", TRUE ~ "ns")
  )

write_csv(comparison, file.path(out_dir, "comparison_original_caguas_vs_regenerated.csv"))

counts <- slopes |>
  group_by(spec) |>
  summarise(
    n_estimates = n(),
    nominal_clustered = sum(p_clustered < 0.05),
    fdr_clustered = sum(p_bh_48_clustered < 0.05),
    nominal_model_based = sum(p_model < 0.05),
    fdr_model_based = sum(p_bh_48_model < 0.05),
    .groups = "drop"
  ) |>
  bind_rows(tibble(
    spec = "original Table 6 (Caguas series, model-based SEs)",
    n_estimates = nrow(orig),
    nominal_model_based = sum(orig$orig_p < 0.05),
    fdr_model_based = sum(orig$orig_q < 0.05)
  ))

write_csv(counts, file.path(out_dir, "table6_family_significance_counts.csv"))

message("Regenerated Table 6 family written to: ", out_dir)
