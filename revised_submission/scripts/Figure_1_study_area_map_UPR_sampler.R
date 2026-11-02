#!/usr/bin/env Rscript
############################################################
## Figure 1: Study-Area Map with UPR Burkard Sampler
## Focus:  Six study municipalities + Caguas, 2020 Census tract population
## Station point: UPR Medical Sciences Campus rooftop Burkard sampler
##   (18.3986 N, -66.0747 W; coordinates confirmed 2026-10-19), replacing the NOAA
##   airport point (RQW00011641) used in 01_representativeness_map_tables.R
##
## Data-only re-export (boundaries + census population + station point);
##   no statistical model is run.
## Requirements: sf, tigris, tidycensus, ggplot2, scales, dplyr, readr; a Census
##   API key available to tidycensus. Run from any working directory.
##
## Outputs (PNG, PDF, CSV) written to: OUT_DIR below
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(sf); library(tigris)
  library(tidycensus); library(ggplot2); library(scales)
})

OUT_DIR <- "/Users/felixrivera-mariani/Downloads/Asthma-and-Aeroallergens/Journal_Aerobiology_Submission_response/Resubmission_Package_2026-10-20_run01/02_Figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

norm_key <- function(x) x |> iconv(to = "ASCII//TRANSLIT", sub = "") |> toupper() |> gsub("[^A-Z ]", "", x = _) |> trimws()

muni_study6 <- c("San Juan", "Guaynabo", "Carolina", "Trujillo Alto", "Canovanas", "Loiza")
muni_with_caguas <- c(muni_study6, "Caguas")

# Burkard sampler, rooftop, UPR Medical Sciences Campus (coordinates confirmed 2026-10-19)
station_tbl <- tibble::tibble(station_name = "UPR MSC rooftop Burkard sampler (NAB San Juan)",
                              station_lat = 18.3986, station_lon = -66.0747)

options(tigris_use_cache = TRUE)
pr_counties <- counties(state = "72", year = 2020, cb = TRUE, class = "sf") |> mutate(name_key = norm_key(NAME))
muni_shapes <- pr_counties |> filter(name_key %in% norm_key(muni_with_caguas))
stopifnot(nrow(muni_shapes) == length(muni_with_caguas))

station_sf <- st_as_sf(station_tbl, coords = c("station_lon", "station_lat"), crs = 4326, remove = FALSE)

# Centroid distances (cross-check against 03_Internal_Audit/station_to_municipality_centroid_distances_UPR_sampler.csv)
centroid_tbl <- muni_shapes |> st_transform(32161) |> st_centroid() |> st_transform(4326) |>
  mutate(municipality = NAME,
         centroid_lon = st_coordinates(geometry)[, 1], centroid_lat = st_coordinates(geometry)[, 2],
         distance_km = as.numeric(st_distance(geometry, st_geometry(station_sf))[, 1]) / 1000) |>
  st_drop_geometry() |> dplyr::select(municipality, centroid_lat, centroid_lon, distance_km) |> arrange(distance_km)
write_csv(centroid_tbl, file.path(OUT_DIR, "Figure_1_centroid_distances_UPR_sampler_check.csv"))

tract_pop <- get_decennial(geography = "tract", state = "72", county = sort(unique(muni_shapes$COUNTYFP)),
                           variables = "P1_001N", year = 2020, sumfile = "pl", geometry = TRUE)

bbox <- st_bbox(st_transform(muni_shapes, 32161))
xpad <- (bbox$xmax - bbox$xmin) * 0.12; ypad <- (bbox$ymax - bbox$ymin) * 0.12

map_plot <- ggplot() +
  geom_sf(data = st_transform(pr_counties, 32161), fill = "grey97", color = "grey85", linewidth = 0.2) +
  geom_sf(data = st_transform(tract_pop, 32161), aes(fill = value), color = NA, alpha = 0.95) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", labels = label_comma(), name = "2020 Census\ntract population") +
  geom_sf(data = st_transform(muni_shapes |> filter(NAME != "Caguas"), 32161), fill = NA, color = "grey10", linewidth = 0.5) +
  geom_sf(data = st_transform(muni_shapes |> filter(NAME == "Caguas"), 32161), fill = NA, color = "#D95F02", linewidth = 1) +
  geom_sf(data = st_transform(station_sf, 32161), shape = 21, fill = "black", color = "white", size = 3.2, stroke = 0.5) +
  geom_sf_text(data = st_transform(muni_shapes, 32161), aes(label = NAME), size = 3, color = "grey15") +
  coord_sf(xlim = c(bbox$xmin - xpad, bbox$xmax + xpad), ylim = c(bbox$ymin - ypad, bbox$ymax + ypad), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "right")

# 600 dpi PNG for MDPI; sources are stated in the manuscript caption (not inside the image)
ggsave(file.path(OUT_DIR, "Figure_1.png"), map_plot, width = 7.5, height = 5.2, dpi = 600, bg = "white")
ggsave(file.path(OUT_DIR, "Figure_1.pdf"), map_plot, width = 7.5, height = 5.2, bg = "white")
message("Figure 1 written to ", OUT_DIR)

