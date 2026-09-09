cities_to_run <- tribble(
  ~c, ~aoi_name, ~shade,
  "BRA-Rio_de_Janeiro", "low_emission_zone", FALSE,
  "MEX-Monterrey", "mitras_centro", FALSE,
  "BRA-Teresina", "accelerator_area_big", FALSE,
  "ARG-Buenos_Aires", "barrio_20", TRUE,
  "ZAF-Johannesburg", "jukskei-river", FALSE,
  "ZAF-Cape_Town", "business_district", TRUE,
  "IND-Bhopal", "tt_nagar", FALSE,
  "BRA-Campinas", "accelerator_area", FALSE,
  "BRA-Florianopolis", "accelerator_area", FALSE,
  "BRA-Fortaleza", "accelerator_area", FALSE,
  "BRA-Recife", "accelerator_area", FALSE
)

align_to_worldpop <- function(r, wp, mask_raster = NULL, fun = mean, categorical = FALSE) {
  
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

calc_cors <- function(city, aoi_name) {
  
  city_folder <- file.path("city_projects", city, aoi_name)
  baseline_folder <- file.path(city_folder, "scenarios", "baseline", "baseline")
  
  result <- tibble(
    city = city,
    aoi_name = aoi_name,
    cor_sun_vs_pctshade = NA_real_,
    cor_mean_vs_pctshade = NA_real_,
    error = NA_character_
  )
  
  # AOI
  aoi <- st_read(glue("{aws_http}/{baseline_folder}/aoi__baseline__baseline.geojson"))
  
  # UTCI
  tiles <- list_tiles(paste0("s3://", bucket, "/", baseline_folder))
  
  utci_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/utci-1500__baseline__baseline.tif")
  utci <- load_and_merge(utci_paths) |> 
    crop(aoi) |> 
    mask(aoi)
  
  # Shade
  shade_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/shade-1500__baseline__baseline.tif")
  shade <- (load_and_merge(shade_paths) |> 
              crop(aoi) |> 
              mask(aoi)) > 0
  
  # WorldPop
  world_pop <- rast(glue("https://wri-cities-indicators.s3.us-east-1.amazonaws.com/data/published/layers/WorldPop/tif/{city}__urban_extent__WorldPop__StartYear_2020_EndYear_2020.tif")) |> 
    crop(aoi) |> 
    mask(aoi)
  
  # Non-building area
  non_build_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/non-building-areas__baseline__baseline.tif")
  non_build <- load_and_merge(non_build_paths) |> 
    crop(aoi) |> 
    mask(aoi)
  
  # Non-building area mean UTCI
  non_build_mean <- align_to_worldpop(utci, world_pop, non_build) 
  non_build_mean_norm <- non_build_mean |> 
    normalize_percentile()
  
  # Unshaded non-building area mean UTCI
  utci_sun_only <- mask(utci, shade, maskvalues = 1)
  utci_sun_agg <- align_to_worldpop(utci_sun_only, world_pop, non_build, fun = mean)
  utci_sun_agg_norm <- utci_sun_agg |> 
    normalize_percentile()
  
  # Non-building shade
  shade_nonbuild <- align_to_worldpop(shade, world_pop, non_build, fun = mean)
  shade_nonbuild_norm <- shade_nonbuild |> 
    normalize_percentile()
  
  # Unshaded non-building mean and non-building mean
  cor1 <- cor(values(utci_sun_agg_norm), values(shade_nonbuild_norm), method = "spearman", use = "complete.obs")
  
  # Unshaded non-building mean and non-building mean
  cor2 <- cor(values(non_build_mean_norm), values(shade_nonbuild_norm), method = "spearman", use = "complete.obs")
  
  result$cor_sun_vs_pctshade  <- cor1
  result$cor_mean_vs_pctshade <- cor2

  result
}




# AWS
bucket   <- "wri-cities-tcm"
aws_http <- "https://wri-cities-tcm.s3.us-east-1.amazonaws.com"

cor_results <- pmap_dfr(
  list(cities_to_run$c, cities_to_run$aoi_name),
  calc_cors
)


cor_check <- tibble(
  mean_utci   = values(non_build_mean),
  pct_shade   = values(shade_nonbuild)
) |> na.omit()

cor(cor_check$mean_utci, cor_check$pct_shade)

ggplot(cor_check, aes(pct_shade, mean_utci)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  theme_minimal()

ggplot(cor_results) +
  geom_col(aes(y = fct_reorder(city, cor_sun_vs_pctshade, .desc = TRUE), x = cor_sun_vs_pctshade))

extract_utci_values <- function(city, aoi_name) {
  
  message(glue("Extracting {city} / {aoi_name}..."))
  
  city_folder <- file.path("city_projects", city, aoi_name)
  baseline_folder <- file.path(city_folder, "scenarios", "baseline", "baseline")
  
  out <- tibble(city = character(), 
                aoi_name = character(), 
                layer = character(), 
                value = double())
  
  tryCatch({
    
    # AOI
    aoi <- st_read(glue("{aws_http}/{baseline_folder}/aoi__baseline__baseline.geojson"),
                   quiet = TRUE)
    
    # UTCI
    tiles <- list_tiles(paste0("s3://", bucket, "/", baseline_folder))
    
    utci_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/utci-1500__baseline__baseline.tif")
    utci <- load_and_merge(utci_paths) |>
      crop(aoi) |>
      mask(aoi)
    
    # Shade
    shade_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/shade-1500__baseline__baseline.tif")
    shade <- (load_and_merge(shade_paths) |>
                crop(aoi) |>
                mask(aoi)) > 0
    
    # WorldPop
    world_pop <- rast(glue("https://wri-cities-indicators.s3.us-east-1.amazonaws.com/data/published/layers/WorldPop/tif/{city}__urban_extent__WorldPop__StartYear_2020_EndYear_2020.tif")) |>
      crop(aoi) |>
      mask(aoi)
    
    # Non-building area
    non_build_paths <- glue("{aws_http}/{baseline_folder}/{tiles}/ccl_layers/non-building-areas__baseline__baseline.tif")
    non_build <- load_and_merge(non_build_paths) |>
      crop(aoi) |>
      mask(aoi)
    
    # 1-m non-building UTCI
    utci_non_build <- mask(utci, non_build, maskvalues = 0)
    nb_values <- as.vector(values(utci_non_build))
    nb_values <- nb_values[!is.na(nb_values)]
    
    out <- bind_rows(out, tibble(
      city = city, aoi_name = aoi_name, layer = "Non-building UTCI (1-m)", value = nb_values
    ))
    
    # 1-m non-building unshaded UTCI
    utci_non_build_unshaded <- mask(utci_non_build, shade, maskvalues = 1)
    nb_unshaded_values <- as.vector(values(utci_non_build_unshaded))
    nb_unshaded_values <- nb_unshaded_values[!is.na(nb_unshaded_values)]
    
    out <- bind_rows(out, tibble(
      city = city, aoi_name = aoi_name, layer = "Non-building unshaded UTCI (1-m)", value = nb_unshaded_values
    ))
    
    # 100-m mean-aggregated non-building UTCI
    non_build_mean <- align_to_worldpop(utci_non_build, world_pop, fun = mean)
    nb_mean_values <- as.vector(values(non_build_mean))
    nb_mean_values <- nb_mean_values[!is.na(nb_mean_values)]
    
    out <- bind_rows(out, tibble(
      city = city, aoi_name = aoi_name, layer = "Non-building UTCI (100-m mean)", value = nb_mean_values
    ))
    
    # 100-m mean-aggregated non-building unshaded UTCI
    non_build_unshaded_mean <- align_to_worldpop(utci_non_build_unshaded, world_pop, fun = mean)
    nb_unshaded_mean_values <- as.vector(values(non_build_unshaded_mean))
    nb_unshaded_mean_values <- nb_unshaded_mean_values[!is.na(nb_unshaded_mean_values)]
    
    out <- bind_rows(out, tibble(
      city = city, 
      aoi_name = aoi_name, 
      layer = "Non-building unshaded UTCI (100-m mean)", 
      value = nb_unshaded_mean_values
    ))
    
  }, error = function(e) {
    message(glue("  Error in {city}: {conditionMessage(e)}"))
  })
  
  out
}

all_utci_values <- pmap_dfr(
  list(cities_to_run$c, cities_to_run$aoi_name),
  extract_utci_values
)

ggplot(all_utci_values, aes(x = value, fill = city, color = city)) +
  geom_density(alpha = 0.3, na.rm = TRUE) +
  facet_wrap(~ layer, scales = "free") +
  labs(x = "UTCI", y = "Density", fill = NULL, color = NULL) +
  theme_minimal()

ggplot(all_utci_values, aes(x = value, color = layer)) +
  geom_density(alpha = 0.3, na.rm = TRUE) +
  facet_wrap(~ city, scales = "free") +
  labs(x = "UTCI", y = "Density", fill = NULL, color = NULL) +
  theme_minimal()
