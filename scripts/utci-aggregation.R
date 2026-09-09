# 100-m aggregation of UTCI

# mean, median, majority category
# full, non-building area, pedestrian area

library(here)
library(tidyverse)
library(terra)
library(sf)
library(glue)
library(patchwork)
library(tidyterra)

source("~/Documents/github/cities-heat-resilient-infrastructure/tiling-scripts/utils.R")

# Login once (browser auth)
system("aws sso login --profile cities-data-dev")

# City
city <- "ZAF-Cape_Town"
aoi_name <- "business_district"
city_folder <- file.path("city_projects", city, aoi_name)
baseline_folder <- file.path(city_folder, "scenarios", "baseline", "baseline")

# AWS
bucket   <- "wri-cities-tcm"
aws_http <- "https://wri-cities-tcm.s3.us-east-1.amazonaws.com"

# AOI
aoi <- st_read(glue("{aws_http}/{baseline_folder}/aoi__baseline__baseline.geojson"))

# UTCI
tiles <- list_tiles(paste0("s3://", bucket, "/", baseline_folder))

utci_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/utci-1500__baseline__baseline.tif")
utci <- load_and_merge(utci_paths) |> 
  crop(aoi) |> 
  mask(aoi)

utci_cat_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/utci-cat-1500__baseline__baseline.tif")
utci_cat <- load_and_merge(utci_cat_paths) |> 
  crop(aoi) |> 
  mask(aoi)

# Shade
shade_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/shade-1500__baseline__baseline.tif")
shade <- (load_and_merge(shade_paths) |> 
  crop(aoi) |> 
  mask(aoi)) > 0

# WorldPop
world_pop <- rast(glue("https://wri-cities-indicators.s3.us-east-1.amazonaws.com/data/published/layers/WorldPop/tif/{city}__urban_extent__urban_extent__WorldPop__StartYear_2020_EndYear_2020.tif")) |> 
  crop(aoi) |> 
  mask(aoi)

# Water
open_urban_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/raster_files/cif_open_urban.tif")
open_urban <- load_and_merge(open_urban_paths) |> 
  crop(aoi) |> 
  mask(aoi)
water <- open_urban == 300

# Non-building area
non_build_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/non-building-areas__baseline__baseline.tif")
non_build <- load_and_merge(non_build_paths) |> 
  crop(aoi) |> 
  mask(aoi)

# Pedestrian area
ped_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/pedestrian-areas__baseline__baseline.tif")
ped <- load_and_merge(ped_paths) |> 
  crop(aoi) |> 
  mask(aoi)

# Functions

align_to_worldpop <- function(r, wp, mask_raster = NULL, fun = mean, categorical = FALSE) {
  
  # water mask: aggregate open urban water to raster grid to get the percent
  # water per cell, then choose cells with >= 50% water
  # fraction of water in each coarse cell
  # water_frac <- aggregate(
  #   water,
  #   fact = round(res(r) / res(water)),
  #   fun  = mean,
  #   na.rm = TRUE
  # )
  # 
  # # align to LST grid
  # water_frac_r <- resample(water_frac, r, method = "bilinear")
  # 
  # # majority-water mask (TRUE where >50% water)
  # water_mask <- water_frac_r >= 0.5
  # r <- mask(r, water_mask, maskvalues = 1)
  
  # Mask raster
  if(!is.null(mask_raster)) {
    r <- mask(r, mask_raster, maskvalues = 0)
  }
  
  # aggregate to roughly WorldPop resolution
  fact <- round(res(wp) / res(r))
  r_agg <- aggregate(r, fact = fact, fun = fun, na.rm = TRUE)
  
  # resample to match WorldPop grid
  resample_method <- if (categorical) "near" else "bilinear"
  r_wp <- resample(r_agg, wp, method = resample_method)
  
  return(r_wp)
}

normalize_percentile <- function(r, probs = seq(0, 1, by = 0.01)) {
  
  qs <- quantile(values(r), probs = probs, na.rm = TRUE)
  
  classify(
    r,
    rcl = cbind(qs[-length(qs)], qs[-1], probs[-1]),
    include.lowest = TRUE
  )
}

plot_raster <- function(r, title, fill_label = "100-m UTCI") {
  ggplot() +
    geom_spatraster(data = r, na.rm = TRUE) +
    scale_fill_viridis_c(option = "magma", na.value = NA, name = fill_label) +
    ggtitle(title) +
    theme_void() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
}

# UTCI --------------------------------------------------------------------


## Mean --------------------------------------------------------------------

# Full
full_mean <- align_to_worldpop(utci, world_pop) |> 
  normalize_percentile()

# Non-building area
non_build_mean <- align_to_worldpop(utci, world_pop, non_build) |> 
  normalize_percentile()

# Pedestrian area
ped_mean <- align_to_worldpop(utci, world_pop, ped) |> 
  normalize_percentile()

plot_raster(utci, "1-m UTCI") + 
  plot_raster(full_mean, "Full AOI") + 
  plot_raster(non_build_mean, "Non-building area") + 
  plot_raster(ped_mean, "Pedestrian area") +
  plot_annotation(title = "Mean")

## Median --------------------------------------------------------------------

# Full
full_median <- align_to_worldpop(utci, world_pop, fun = median) |> 
  normalize_percentile()

# Non-building area
non_build_median <- align_to_worldpop(utci, world_pop, non_build, fun = median) |> 
  normalize_percentile()

# Pedestrian area
ped_median <- align_to_worldpop(utci, world_pop, ped, fun = median) |> 
  normalize_percentile()

plot_raster(utci, "1-m UTCI") + 
  plot_raster(full_median, "Full AOI") + 
  plot_raster(non_build_median, "Non-building area") + 
  plot_raster(ped_median, "Pedestrian area") +
  plot_annotation(title = "Median")


# Mean vs. median ---------------------------------------------------------


diff_r <- non_build_mean - non_build_median

change_r <- as.factor(classify(diff_r, rbind(c(-Inf, 0, -1), c(0, 0, 0), c(0, Inf, 1))))
levels(change_r) <- data.frame(value = c(-1, 0, 1), label = c("Median higher", "No change", "Mean higher"))


diff_plot <- ggplot() +
  geom_spatraster(data = change_r, na.rm = TRUE) +
  scale_fill_manual(values = c("Median higher" = "#4575B4", "No change" = "#D9D9D9", "Mean higher" = "#D73027"), na.value = NA) +
  theme_void()

plot_raster(utci, "1-m UTCI") + 
  plot_raster(non_build_mean, "Mean") + 
  plot_raster(non_build_median, "Median") +
  diff_plot +
  plot_annotation(title = "Non-building area")

# UTCI density
median_values <- values(align_to_worldpop(utci, world_pop, non_build, fun = median))
mean_values <- values(align_to_worldpop(utci, world_pop, non_build, fun = mean))
utci_values <- values(utci)

df <- data.frame(
  value = c(utci_values, mean_values, median_values),
  source = rep(c("UTCI (original)", "Mean-aggregated", "Median-aggregated"),
               times = c(length(utci_values), length(mean_values), length(median_values)))
)

ggplot(df, aes(x = value, fill = source, color = source)) +
  geom_density(alpha = 0.3, na.rm = TRUE) +
  labs(x = "UTCI", y = "Density", fill = NULL, color = NULL) +
  theme_minimal()


## Sun only ------------------------------

utci_sun_only <- mask(utci, shade, maskvalues = 1)
utci_sun_agg <- align_to_worldpop(utci_sun_only, world_pop, non_build, fun = mean)
utci_sun_agg_norm <- utci_sun_agg |> 
  normalize_percentile()

plot_raster(utci, "1-m UTCI") + 
  plot_raster(non_build_mean, "Mean") + 
  plot_raster(utci_sun_agg_norm, "Unshaded mean") +
  plot_annotation(title = "Non-building area")

# UTCI density
sun_values <- values(utci_sun_agg)
mean_values <- values(align_to_worldpop(utci, world_pop, non_build, fun = mean))
utci_values <- values(utci)

df <- data.frame(
  value = c(utci_values, mean_values, sun_values),
  source = rep(c("UTCI (original)", "Non-building mean", "Unshaded non-building mean"),
               times = c(length(utci_values), length(mean_values), length(sun_values)))
)

ggplot(df, aes(x = value, fill = source, color = source)) +
  geom_density(alpha = 0.3, na.rm = TRUE) +
  labs(x = "UTCI", y = "Density", fill = NULL, color = NULL) +
  theme_minimal()


# Shade- percent shade in non-building areas ------------------------------

shade_nonbuild <- align_to_worldpop(shade, world_pop, non_build, fun = mean)
shade_nonbuild_norm <- shade_nonbuild |> 
  normalize_percentile()

plot_raster(shade, title = "1-m Shade", fill_label = "Shade %") + 
  plot_raster(align_to_worldpop(shade, world_pop, fun = mean), "% Shade", fill_label = "Shade %") +
  plot_raster(shade_nonbuild_norm, "% Shade (Non-building)", fill_label = "Shade %") 


# Testing -----------------------------------------------------------------

cor_check <- tibble(
  mean_utci   = values(non_build_mean),
  median_utci = values(non_build_median),
  pct_shade   = values(shade_nonbuild)
) |> na.omit()

cor(cor_check$mean_utci, cor_check$pct_shade)
cor(cor_check$median_utci, cor_check$pct_shade)

ggplot(cor_check, aes(pct_shade, mean_utci)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  theme_minimal()

# Unshaded non-building mean and non-building mean
cor(values(utci_sun_agg_norm), values(non_build_mean), method = "spearman", use = "complete.obs")

# Unshaded non-building mean and non-building mean
cor(values(shade_nonbuild_norm), values(non_build_mean), method = "spearman", use = "complete.obs")

## Majority category --------------------------------------------------------------------

# Full
full_cat <- align_to_worldpop(utci_cat, world_pop, fun = "modal", categorical = TRUE) |> 
  normalize_percentile()

# Doesn't make sense to do this!!

