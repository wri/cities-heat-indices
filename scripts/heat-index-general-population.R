library(here)
library(tidyverse)
library(terra)
library(sf)
library(glue)

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

# data s3 path
s3_path = "https://wri-cities-data-api.s3.us-east-1.amazonaws.com/data/dev"

# city id
city_name <- "CHN-Chengdu"

############ GEt data 
# world_pop
world_pop <- rast(glue("{s3_path}/WorldPop/cog/{city_name}__urban_extent__WorldPop__Version_2__StartYear_2020_EndYear_2020.tif"))
# lst
lst <- rast(glue("https://wri-cities-data-api.s3.us-east-1.amazonaws.com/data/dev/HighLandSurfaceTemperatureIndex/cog/{city_name}__urban_extent__HighLandSurfaceTemperatureIndex__StartYear_2023_EndYear_2025.tif"))
# alb
alb <- rast(glue("{s3_path}/AlbedoCloudMaskedIndex/cog/{city_name}__urban_extent__AlbedoCloudMaskedIndex__2026.tif"))
# fr
fr <- rast(glue("{s3_path}/FractionalVegetationPercentIndex/cog/{city_name}__urban_extent__FractionalVegetationPercentIndex__StartYear_2025_EndYear_2025.tif"))
# tree
tree <- rast(glue("{s3_path}/TreeCanopyCoverMaskIndex/cog/{city_name}__urban_extent__TreeCanopyCoverMaskIndex__2020.tif"))

############# process data
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
index_p = index_p *100

############## Store data

# define path 
filename = glue("./{city_name}__urban_extent__HeatRiskIndexGeneral__2026.tif")

# write
terra::writeRaster(index_p,
                   filename,
                   filetype = "COG",
                   overwrite = TRUE)
