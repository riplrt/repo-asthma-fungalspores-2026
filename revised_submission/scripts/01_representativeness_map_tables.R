#!/usr/bin/env Rscript
############################################################
## Representativeness of Spore Monitoring Across Municipalities
## Focus:  Six study municipalities + Caguas (2020 Census / TIGER)
## Descriptive only (no models):
##   - Encounter share by municipality
##   - Station-to-municipality centroid distances
##   - 2020 Census population by tract and by municipality
##   - Draft Figure 1 map (station point = NOAA RQW00011641; superseded by
##     Figure_1_study_area_map_UPR_sampler.R)
##   - Re-exports Table 8 from the geographic sensitivity script if present
##
## Outputs (CSV, PNG) written to: R2_representativeness_monitoring/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(tigris)
  library(tidycensus)
  library(ggplot2)
  library(scales)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_representativeness_monitoring/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

norm_key <- function(x) {
  x |>
    iconv(to = "ASCII//TRANSLIT", sub = "") |>
    toupper() |>
    gsub("[^A-Z ]", "", x = _) |>
    trimws()
}

muni_study6 <- c("San Juan", "Guaynabo", "Carolina", "Trujillo Alto", "Canovanas", "Loiza")
muni_with_caguas <- c(muni_study6, "Caguas")

station_tbl <- tibble::tibble(
  station_id = "RQW00011641",
  station_name = "SAN JUAN L M MARIN INTL AP",
  station_lat = 18.4325,
  station_lon = -66.0106,
  station_elev_m = 3
)

panel_daily <- read_csv("asthma_aeroallergen_municipality/data/dlnm_outputs/analytic_daily_panel.csv", show_col_types = FALSE)

encounter_share <- panel_daily |>
  mutate(municipality_raw = as.character(municipality_desc), municipality_key = norm_key(municipality_raw)) |>
  filter(municipality_key %in% norm_key(muni_with_caguas)) |>
  group_by(municipality_key) |>
  summarise(encounters = sum(asthma_anydx_n, na.rm = TRUE), .groups = "drop") |>
  mutate(municipality = muni_with_caguas[match(municipality_key, norm_key(muni_with_caguas))]) |>
  mutate(
    in_primary_study6 = municipality %in% muni_study6,
    share_of_total_6plus_caguas = encounters / sum(encounters),
    share_of_primary_6 = if_else(in_primary_study6, encounters / sum(encounters[in_primary_study6]), NA_real_)
  ) |>
  dplyr::select(municipality, encounters, in_primary_study6, share_of_total_6plus_caguas, share_of_primary_6) |>
  arrange(desc(encounters))

write_csv(encounter_share, file.path(out_dir, "Table_R2_municipality_encounter_share_draft.csv"))

options(tigris_use_cache = TRUE)

pr_counties <- counties(state = "72", year = 2020, cb = TRUE, class = "sf")

pr_counties <- pr_counties |>
  mutate(name_key = norm_key(NAME))

target_key <- norm_key(muni_with_caguas)

muni_shapes <- pr_counties |>
  filter(name_key %in% target_key)

if (nrow(muni_shapes) != length(muni_with_caguas)) {
  stop("Could not match all target municipalities in TIGER county geometries.")
}

station_sf <- st_as_sf(
  station_tbl,
  coords = c("station_lon", "station_lat"),
  crs = 4326,
  remove = FALSE
)

centroid_tbl <- muni_shapes |>
  st_transform(32161) |>
  st_centroid() |>
  st_transform(4326) |>
  mutate(
    municipality = NAME,
    centroid_lon = st_coordinates(geometry)[, 1],
    centroid_lat = st_coordinates(geometry)[, 2],
    distance_km = as.numeric(st_distance(geometry, st_geometry(station_sf))[, 1]) / 1000,
    station_id = station_tbl$station_id,
    station_name = station_tbl$station_name,
    station_lat = station_tbl$station_lat,
    station_lon = station_tbl$station_lon,
    station_elev_m = station_tbl$station_elev_m,
    centroid_method = "Centroid of 2020 TIGER/Line county polygon (Puerto Rico municipio equivalent)",
    distance_method = "Great-circle distance from municipio polygon centroid to station point"
  ) |>
  st_drop_geometry() |>
  dplyr::select(
    municipality, centroid_lat, centroid_lon,
    station_id, station_name, station_lat, station_lon, station_elev_m,
    distance_km, centroid_method, distance_method
  ) |>
  arrange(distance_km)

write_csv(centroid_tbl, file.path(out_dir, "Table_R2_station_to_municipality_centroid_distances_2020_tiger_draft.csv"))

api_key <- Sys.getenv("CENSUS_API_KEY")

tract_pop <- NULL
if (api_key != "") {
  tract_pop <- tryCatch(
    get_decennial(
      geography = "tract",
      state = "72",
      county = sort(unique(muni_shapes$COUNTYFP)),
      variables = "P1_001N",
      year = 2020,
      sumfile = "pl",
      geometry = TRUE,
      key = api_key
    ) |>
      mutate(COUNTYFP = substr(GEOID, 3, 5)) |>
      left_join(st_drop_geometry(muni_shapes) |> dplyr::select(COUNTYFP, municipio = NAME), by = "COUNTYFP"),
    error = function(e) e
  )
}

if (inherits(tract_pop, "error") || is.null(tract_pop)) {
  issue_msg <- if (inherits(tract_pop, "error")) conditionMessage(tract_pop) else "CENSUS_API_KEY not set"
  readr::write_lines(
    c(
      "Population-layer retrieval did not complete.",
      paste0("Reason: ", issue_msg),
      "Map was generated without tract population shading."
    ),
    file.path(out_dir, "population_data_status.txt")
  )
} else {
  status_file <- file.path(out_dir, "population_data_status.txt")
  if (file.exists(status_file)) file.remove(status_file)

  tract_pop_tbl <- tract_pop |>
    st_drop_geometry() |>
    transmute(
      geoid = GEOID,
      municipio,
      tract_name = NAME,
      variable,
      population_2020 = value
    )

  write_csv(tract_pop_tbl, file.path(out_dir, "Table_R2_population_by_tract_2020_census_draft.csv"))

  muni_pop <- tract_pop_tbl |>
    group_by(municipio) |>
    summarise(population_2020 = sum(population_2020, na.rm = TRUE), .groups = "drop")

  write_csv(muni_pop, file.path(out_dir, "Table_R2_population_by_municipality_2020_census_draft.csv"))
}

bbox <- st_bbox(st_transform(muni_shapes, 32161))
xpad <- (bbox$xmax - bbox$xmin) * 0.12
ypad <- (bbox$ymax - bbox$ymin) * 0.12

map_plot <- ggplot() +
  geom_sf(data = st_transform(pr_counties, 32161), fill = "grey97", color = "grey85", linewidth = 0.2)

if (!inherits(tract_pop, "error") && !is.null(tract_pop)) {
  map_plot <- map_plot +
    geom_sf(data = st_transform(tract_pop, 32161), aes(fill = value), color = NA, alpha = 0.95) +
    scale_fill_viridis_c(option = "C", trans = "sqrt", labels = label_comma(), name = "2020 Census\ntract population")
}

map_plot <- map_plot +
  geom_sf(data = st_transform(muni_shapes |> filter(NAME != "Caguas"), 32161), fill = NA, color = "grey10", linewidth = 0.5) +
  geom_sf(data = st_transform(muni_shapes |> filter(NAME == "Caguas"), 32161), fill = NA, color = "#D95F02", linewidth = 1) +
  geom_sf(data = st_transform(station_sf, 32161), shape = 21, fill = "black", color = "white", size = 2.8, stroke = 0.4) +
  coord_sf(
    xlim = c(bbox$xmin - xpad, bbox$xmax + xpad),
    ylim = c(bbox$ymin - ypad, bbox$ymax + ypad),
    expand = FALSE
  ) +
  labs(
    title = "Figure 1 draft: Monitoring station and study municipalities",
    subtitle = "Six-municipality primary study area with Caguas highlighted for sensitivity analysis",
    caption = paste(
      "Boundaries: 2020 U.S. Census TIGER/Line county polygons (Puerto Rico municipios).",
      "Point: NOAA station RQW00011641 (18.4325, -66.0106; elevation 3 m)."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "right", plot.caption = element_text(hjust = 0))

ggsave(
  filename = file.path(out_dir, "Figure_1_monitoring_station_municipalities_population_2020_draft.png"),
  plot = map_plot,
  width = 9.5,
  height = 6.4,
  dpi = 320,
  bg = "white"
)

tab8_in <- "Journal_Aerobiology_Submission_response/R2_geographic_sensitivity/tables_drafts/Table_8_geographic_sensitivity_draft.csv"
if (file.exists(tab8_in)) {
  tab8 <- read_csv(tab8_in, show_col_types = FALSE)
  write_csv(tab8, file.path(out_dir, "Table_8_geographic_sensitivity_reused_draft.csv"))
}

message("Representativeness map and tables written to: ", out_dir)
