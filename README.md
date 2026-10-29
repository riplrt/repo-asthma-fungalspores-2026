# Fungal Spores and Asthma Healthcare Utilization in San Juan, Puerto Rico (2018–2021)

Analysis code for the manuscript:

> **Age- and Period-Stratified Associations Between Ambient Fungal Spore Exposure and Asthma Healthcare Utilization in the San Juan Metropolitan Area, Puerto Rico, 2018–2021**
> Benjamín Bolaños-Rosero, Claudia I. Pereira Moreira, Catalina S. Robles-Rivera, Félix E. Rivera-Mariani


The R scripts in this repository produced the tables, figures and supplementary material reported in that manuscript.

## Study overview

This ecological time-series study linked daily Burkard spore counts from the University of Puerto Rico Medical Sciences Campus rooftop sampler to 168,578 asthma encounters in the San Juan health region (Loíza, San Juan, Trujillo Alto, Guaynabo, Canóvanas, Carolina), with Caguas analyzed as a sensitivity case. It asked whether spore–asthma associations held steady across two contextual eras: post-Hurricane María recovery (2018–2019) and the COVID-19 pandemic (2020–2021). It also asked whether they held steady across age groups.

- **Outcomes:** encounter volume (negative binomial rate ratios) and the ED/inpatient encounter fraction (logistic odds ratios)
- **Exposure:** spores per 10,000 spores·m⁻³ as weekly means and 3- and 7-day moving averages (MA3, MA7), after a weekly missing-data QC rule
- **Methods:** age-adjusted and age-stratified models, period × age interaction tests with date-clustered standard errors, Benjamini–Hochberg FDR control, distributed lag non-linear models (DLNMs), time-stratified case-crossover analysis of ED/inpatient events, and geographic sensitivity analyses
- **Covariates:** daily NOAA weather from Luis Muñoz Marín International Airport (station RQW00011641)

The associations describe population-level co-occurrence, not individual-level risk.

## Structure

- `prior_submission/scripts/`: the original analysis pipeline (yearly, weekly/lag, pandemic-phase and case-crossover models).
- `revised_submission/scripts/`: reanalyses written in response to reviewer comments (multiplicity, interaction tests, corrected Table 6 family, exposure window and case-crossover, ED/inpatient fraction, geographic sensitivity, Figure 1).

Each script opens with a header block that describes its purpose, models and output location.

## Mapping scripts to the manuscript

| Manuscript component | Script(s) |
|---|---|
| Figure 1 (study-area map) | `revised_submission/scripts/Figure_1_study_area_map_UPR_sampler.R` |
| Table 1, age and period descriptives | `01_utilization_interpretation_tables.R`, `01_ed_inpatient_fraction_reanalysis.R`, `02_age_stratum_descriptives.R` |
| Table 2 / Table 3 / Table S1 (yearly age-adjusted models, partial Spearman) | `prior_submission/scripts/asthma_spores_yearly_v2_age_adjusted.R` |
| Tables 4–5 (pre-pandemic weekly and lag models) | `prior_submission/scripts/asthma_spores_yearly_weekly_lags*.R` |
| Tables 6a/6b, S5/S6, Figure S2 (period × age, corrected San Juan series) | `01_spore_age_period_interaction_tests.R` → `01_regenerate_table6_S5_S6_figS2.R` |
| Table 7, Table S2, Figure 2 (case-crossover, DLNM lag-response) | `01_window_casecrossover_reanalysis.R` |
| Table 8, Table S7 (geographic sensitivity) | `01_geographic_sensitivity_reanalysis.R`, `01_representativeness_map_tables.R` |
| Table S3 (multiplicity inventory), S4 | `01_build_draft_tables_S3_S6.R` |

Unless a full path is given, scripts are in `revised_submission/scripts/`.

## Notes

- The data used in these analyses were obtained with permission and are not included in this repository:
  - **Claims data:** de-identified administrative records from the Puerto Rico Health Insurance Administration.
  - **Fungal spore data:** the National Allergy Bureau (NAB)-certified San Juan and Caguas, PR aeroallergen stations.
- Access to either data source requires permission from the data provider.
