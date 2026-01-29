library(here)
library(tidyverse)
library(terra)
library(sf)
library(glue)
library(patchwork)


# City
city_name <- "ZAF-Cape_Town"
city_folder <- here("data", city_name)

# World Pop from CIF
# use as target raster for aggregation and resampling
world_pop <- rast(glue("https://wri-cities-indicators.s3.us-east-1.amazonaws.com/data/published/layers/WorldPop/tif/{city_name}__urban_extent__urban_extent__WorldPop__StartYear_2020_EndYear_2020.tif"))
utm <- st_crs(world_pop)

# Should be the urban extent. Was using an aoi for testing. We should buffer the aoi by
# 100 m to make sure we have some wiggle room for shifting the grids
aoi <- st_bbox(
  c(xmin = 260748, ymin = 6242306,
    xmax = 262897, ymax = 6244720),
  crs = 32734
) %>%
  st_as_sfc() %>%
  st_as_sf()

world_pop <- world_pop %>% 
  crop(aoi)

# Surface characteristics ---------------------------------------------------------------------
# From CIF. Currently LST is not being calculated according to the tech note, 
# and I'm not sure about albedo and fractional vegetation. These datasets are just
# used for proof of concept. They should not need reprojecting because they are
# all coming from CIF.

# Example for LST calculation: https://code.earthengine.google.com/?scriptPath=users%2Felizabethjanewesley%2Fheat-resilient-infrastructure%3Asurface-characteristics%2Fhot-season-LST-function
# Hot season defined as the three month window centered on the hottest day from
# ERA5. Use Landsat surface temperature (B10) for the most recent three hot seasons
# and calculate the pixel-wise 95th percentiles.

# Example for albedo calculation: https://code.earthengine.google.com/?scriptPath=users%2Felizabethjanewesley%2Fheat-resilient-infrastructure%3Asurface-characteristics%2Fsummer-albedo-function
# Summer of the most recent year. Sentinel-2 with the Cloud Score + cloudmask.

# Example for Fr calculation: https://code.earthengine.google.com/?scriptPath=users%2Felizabethjanewesley%2Fheat-resilient-infrastructure%3Asurface-characteristics%2Fsummer-Fr-function
# Summer of the most recent year. Sentinel-2 with the Cloud Score + cloudmask.


lst <- rast("https://wri-cities-heat.s3.us-east-1.amazonaws.com/index/lst_p95_2022-01-22_1.tif") %>% 
  project(utm$wkt) %>% 
  crop(aoi) 
alb <- rast("https://wri-cities-heat.s3.us-east-1.amazonaws.com/index/albedo_2024_1.tif") %>% 
  project(utm$wkt) %>% 
  crop(aoi)
fr <- rast("https://wri-cities-heat.s3.us-east-1.amazonaws.com/index/Fr-ZAF-Cape_Town-aoi.tif") %>% 
  project(utm$wkt) %>% 
  crop(aoi)
tree <- rast("https://wri-cities-tcm.s3.us-east-1.amazonaws.com/city_projects/ZAF-Cape_Town/OLD-business_district/scenarios/baseline/baseline/tile_00001/ccl_layers/tree-cover__baseline__baseline.tif") %>% 
  project(utm$wkt) %>% 
  crop(aoi)

open_urban <- rast("https://wri-cities-heat.s3.us-east-1.amazonaws.com/ZAF-Cape_Town/scenarios/street-trees/rasters/lulc.tif") %>% 
  project(utm$wkt) %>% 
  crop(aoi) 
water <- open_urban == 300


# Processing --------------------------------------------------------------

# Aggregate & resample
# LST, albedo, Fr, and tree cover
# aggregate to the mean per world pop cell then use bilinear resampling
align_to_worldpop <- function(r, wp, fun = mean) {
  
  # water mask: aggregate open urban water to raster grid to get the percent
  # water per cell, then choose cells with >= 50% water
  # fraction of water in each coarse cell
  water_frac <- aggregate(
    water,
    fact = round(res(r) / res(water)),
    fun  = mean,
    na.rm = TRUE
  )
  
  # align to LST grid
  water_frac_r <- resample(water_frac, r, method = "bilinear")
  
  # majority-water mask (TRUE where >50% water)
  water_mask <- water_frac_r >= 0.5
  r <- mask(r, water_mask, maskvalues = 1)
  
  # aggregate to roughly WorldPop resolution
  fact <- round(res(wp) / res(r))
  r_agg <- aggregate(r, fact = fact, fun = fun, na.rm = TRUE)
  
  # resample to match WorldPop grid
  r_wp <- resample(r_agg, wp, method = "bilinear")
  
  return(r_wp)
}

lst <- align_to_worldpop(lst, world_pop)
alb <- align_to_worldpop(alb, world_pop)
fr <- align_to_worldpop(fr, world_pop)
tree <- align_to_worldpop(tree, world_pop)

# Most input layers are non-normally distributed. We want to be able to combine
# highly skewed inputs in a meaningful way without one layer dominating. We are 
# also interested in relative prioritization rather than magnitude.

# plot(density(values(lst$`lst_p95_2022-01-22_1`), na.rm = T))
# plot(density(values(alb$albedo_2024_1), na.rm = T))
# plot(density(values(fr$Fr), na.rm = T))
# plot(density(values(tree$height), na.rm = T))
# plot(density(values(world_pop$population), na.rm = T))

# Normalize by percentile rank
normalize_percentile <- function(r, probs = seq(0, 1, by = 0.01)) {
  
  qs <- quantile(values(r), probs = probs, na.rm = TRUE)
  
  classify(
    r,
    rcl = cbind(qs[-length(qs)], qs[-1], probs[-1]),
    include.lowest = TRUE
  )
}

# Categorize into 5 levels
cat5_from_01 <- function(x) {
  
  # convert 0–1 to 1–5
  y <- ceiling(x * 5)
  
  # handle edge case where x == 0
  y[y == 0] <- 1
  
  y
}

# Normalize to percentiles (0–1)
lst_p  <- normalize_percentile(lst)
pop_p  <- normalize_percentile(world_pop)

alb_p  <- normalize_percentile(alb)
tree_p <- normalize_percentile(tree)
fr_p   <- normalize_percentile(fr)

# Invert protective infrastructure layers so higher = more vulnerable
alb_p_v  <- 1 - alb_p
tree_p_v <- 1 - tree_p
fr_p_v   <- 1 - fr_p

# Combine infrastructure into a single, evenly weighted component
infra_p <- (alb_p_v + tree_p_v + fr_p_v) / 3

# Sum index components
index_sum <- lst_p + pop_p + infra_p

# Normalize to percentiles
index_p <- normalize_percentile(index_sum)

# Final index
HVI_p_cat <- cat5_from_01(index_p)



# Plots -------------------------------------------------------------------


p2_lst <- ggplot() +
  tidyterra::geom_spatraster(data = lst_p_cat) +
  scale_fill_viridis_c(
    option = "magma",
    name = "Percentile class"
  ) +
  ggtitle("LST")
p2_infra <- ggplot() +
  tidyterra::geom_spatraster(data = infra_p_cat) +
  scale_fill_viridis_c(
    option = "magma",
    name = "Percentile class"
  ) +
  ggtitle("Infrastructure")
p2_pop <- ggplot() +
  tidyterra::geom_spatraster(data = pop_p_cat) +
  scale_fill_viridis_c(
    option = "magma",
    name = "Percentile class"
  ) +
  ggtitle("Population")
p2_hvi <- ggplot() +
  tidyterra::geom_spatraster(data = HVI_p) +
  scale_fill_viridis_c(
    option = "magma",
    name = "Percentile class"
  ) +
  ggtitle("HVI")

p2_lst + p2_infra + p2_pop + p2_hvi
