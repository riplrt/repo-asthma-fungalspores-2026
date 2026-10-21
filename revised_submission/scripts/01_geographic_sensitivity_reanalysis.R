#!/usr/bin/env Rscript
############################################################
## Geographic Sensitivity of the Spore Exposure Assignment
## Focus:  Pooled municipality panels (4 proximal, 6 study, 6 + Caguas)
## Tests whether spore-asthma DLNM estimates depend on the municipalities included
##   - Daily NB DLNM: spores lag 0-7 + temperature lag 0-3, municipality fixed effects
##   - Weekly NB DLNM: spores lag 0-2, season + year + municipality fixed effects
##   - Cumulative RR at the 75th/90th spore percentiles vs the median
## Also: municipality-resolved DLNMs (Table S7), station-to-centroid distances,
##   Caguas exposure-assignment summary
##
## Outputs (CSV) written to: R2_geographic_sensitivity/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(MASS)
  library(splines)
  library(dlnm)
  library(maps)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_geographic_sensitivity/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_daily <- read_csv("asthma_aeroallergen_municipality/data/dlnm_outputs/analytic_daily_panel.csv", show_col_types = FALSE) |>
  mutate(date = as.Date(date), municipality_desc = as.character(municipality_desc))

panel_weekly <- read_csv("asthma_aeroallergen_municipality/data/dlnm_outputs/analytic_weekly_panel.csv", show_col_types = FALSE) |>
  mutate(week = as.Date(week), municipality_desc = as.character(municipality_desc))

muni_rr_full <- read_csv("asthma_aeroallergen_municipality/data/dlnm_outputs_by_municipality_no_weather/overall_rr_summaries_by_municipality.csv", show_col_types = FALSE)

muni_full6 <- c("San Juan", "Guaynabo", "Trujillo Alto", "Carolina", "Canovanas", "Loiza")
muni_prox4 <- c("San Juan", "Guaynabo", "Trujillo Alto", "Carolina")
muni_full7 <- c(muni_full6, "Caguas")

fit_daily_pool <- function(df) {
  d <- df |>
    arrange(date, municipality_desc) |>
    mutate(
      dow = wday(date, label = TRUE, abbr = TRUE, week_start = 1),
      time_num = as.numeric(date - min(date)) + 1,
      municipality_desc = factor(municipality_desc)
    )

  cb_spores <- crossbasis(d$spores_10k, lag = 7,
                          argvar = list(fun = "ns", df = 3),
                          arglag = list(fun = "ns", df = 3))

  cb_temp <- crossbasis(d$temp_mean_c, lag = 3,
                        argvar = list(fun = "ns", df = 3),
                        arglag = list(fun = "ns", df = 2))

  n_years <- max(1, length(unique(year(d$date))))
  time_df <- 8 * n_years

  fit <- glm.nb(asthma_anydx_n ~ cb_spores + cb_temp + ns(rh_mean, 3) + ns(prcp_mm, 3) +
                  dow + ns(time_num, time_df) + municipality_desc,
                data = d, na.action = na.exclude)

  cen_val <- median(d$spores_10k, na.rm = TRUE)
  pred_at <- quantile(d$spores_10k, probs = c(0.75, 0.90), na.rm = TRUE)
  cp <- crosspred(cb_spores, fit, cen = cen_val, at = unname(pred_at), bylag = 1)

  tibble(
    exposure_level = c("75%", "90%"),
    exposure_value = as.numeric(pred_at),
    rr = unname(cp$allRRfit),
    low = unname(cp$allRRlow),
    high = unname(cp$allRRhigh)
  )
}

fit_weekly_pool <- function(df) {
  d <- df |>
    arrange(week, municipality_desc) |>
    mutate(
      municipality_desc = factor(municipality_desc),
      year_f = factor(year(week)),
      season = factor(case_when(
        month(week) %in% c(12, 1, 2) ~ "winter",
        month(week) %in% c(3, 4, 5) ~ "spring",
        month(week) %in% c(6, 7, 8) ~ "summer",
        TRUE ~ "fall"
      ), levels = c("winter", "spring", "summer", "fall"))
    )

  cb_spores <- crossbasis(d$spores_10k, lag = 2,
                          argvar = list(fun = "ns", df = 3),
                          arglag = list(fun = "integer"))

  fit <- glm.nb(asthma_anydx_n ~ cb_spores + ns(temp_mean_c, 3) + ns(rh_mean, 3) + ns(prcp_mm, 3) +
                  season + year_f + municipality_desc,
                data = d, na.action = na.exclude)

  cen_val <- median(d$spores_10k, na.rm = TRUE)
  pred_at <- quantile(d$spores_10k, probs = c(0.75, 0.90), na.rm = TRUE)
  cp <- crosspred(cb_spores, fit, cen = cen_val, at = unname(pred_at), bylag = 1)

  tibble(
    exposure_level = c("75%", "90%"),
    exposure_value = as.numeric(pred_at),
    rr = unname(cp$allRRfit),
    low = unname(cp$allRRlow),
    high = unname(cp$allRRhigh)
  )
}

tab8 <- bind_rows(
  fit_daily_pool(panel_daily |> filter(municipality_desc %in% muni_full6)) |>
    mutate(comparison = "Full six municipalities (excluding Caguas)", model = "Daily pooled lag 0-7"),
  fit_daily_pool(panel_daily |> filter(municipality_desc %in% muni_prox4)) |>
    mutate(comparison = "Proximity-restricted four (~10 km)", model = "Daily pooled lag 0-7"),
  fit_daily_pool(panel_daily |> filter(municipality_desc %in% muni_full7)) |>
    mutate(comparison = "Full six + Caguas", model = "Daily pooled lag 0-7"),
  fit_weekly_pool(panel_weekly |> filter(municipality_desc %in% muni_full6)) |>
    mutate(comparison = "Full six municipalities (excluding Caguas)", model = "Weekly pooled lag 0-2"),
  fit_weekly_pool(panel_weekly |> filter(municipality_desc %in% muni_prox4)) |>
    mutate(comparison = "Proximity-restricted four (~10 km)", model = "Weekly pooled lag 0-2"),
  fit_weekly_pool(panel_weekly |> filter(municipality_desc %in% muni_full7)) |>
    mutate(comparison = "Full six + Caguas", model = "Weekly pooled lag 0-2")
) |>
  arrange(model, exposure_level, comparison)

write_csv(tab8, file.path(out_dir, "Table_8_geographic_sensitivity_draft.csv"))

s7 <- muni_rr_full |>
  filter(model %in% c("daily_any_lag0_7", "weekly_any_lag0_2"),
         exposure_level %in% c("75%", "90%")) |>
  transmute(
    municipality = municipality_desc,
    model,
    exposure_level,
    exposure_value,
    rr = rr_overall,
    low = rr_low,
    high = rr_high,
    includes_caguas_run = TRUE
  ) |>
  arrange(model, municipality, exposure_level)

write_csv(s7, file.path(out_dir, "Table_S7_municipality_resolved_DLNMs_draft.csv"))

stn <- tibble(
  station_id = "RQW00011641",
  station_name = "SAN JUAN L M MARIN INTL AP",
  station_lat = 18.4325,
  station_lon = -66.0106,
  station_elev_m = 3.0
)

cities <- world.cities |>
  filter(country.etc == "Puerto Rico",
         name %in% c("San Juan", "Guaynabo", "Trujillo Alto", "Carolina", "Canovanas", "Loiza", "Caguas")) |>
  transmute(municipality = name, centroid_lat = lat, centroid_lon = long)

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- 6371
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * r * asin(sqrt(a))
}

dist_tbl <- cities |>
  rowwise() |>
  mutate(
    station_id = stn$station_id,
    station_name = stn$station_name,
    station_lat = stn$station_lat,
    station_lon = stn$station_lon,
    station_elev_m = stn$station_elev_m,
    distance_km = haversine_km(station_lat, station_lon, centroid_lat, centroid_lon),
    in_prox4 = municipality %in% muni_prox4
  ) |>
  ungroup() |>
  arrange(distance_km)

write_csv(dist_tbl, file.path(out_dir, "station_to_municipality_centroid_distances_draft.csv"))

caguas_sum <- panel_daily |>
  filter(municipality_desc == "Caguas") |>
  summarise(
    n_days = n(),
    total_asthma_anydx = sum(asthma_anydx_n, na.rm = TRUE),
    mean_daily_asthma_anydx = mean(asthma_anydx_n, na.rm = TRUE),
    min_spores_10k = min(spores_10k, na.rm = TRUE),
    median_spores_10k = median(spores_10k, na.rm = TRUE),
    max_spores_10k = max(spores_10k, na.rm = TRUE)
  ) |>
  mutate(exposure_assignment = "San Juan station spores/pollen/weather series")

write_csv(caguas_sum, file.path(out_dir, "caguas_exposure_assignment_summary_draft.csv"))

message("Geographic sensitivity draft tables written to: ", out_dir)
