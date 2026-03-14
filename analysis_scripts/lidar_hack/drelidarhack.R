
# -----------------------------
# CHESS Workshop 2026 LiDAR Analysis
# AOP CHM vs Drone CHM vs Field Height
# -----------------------------

library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(googledrive)
library(lidR)
library(tidyr)
library(progress)

terraOptions(tempdir = "E:/terra_tmp")

# -----------------------------
# Workspace setup
# -----------------------------
setwd("C:/Users/andre/OneDrive/Documents/Repositories/chess_workshop/analysis_scripts/lidar_hack")

config <- config::get(file=file.path('config.yml'))

lapply(list.files(file.path('.', 'R'), full.names=TRUE), source)
load.pkgs(config$pkgs)

drive_auth(email='andres03678@gmail.com', cache='.')

# -----------------------------
# Field metadata
# -----------------------------
trees_meta <- read.csv(ingest.drive(config$extdata$treedemoid)$local_path)
trees_meta <- trees_meta |> mutate(Site_Number = as.character(Site_Number))

# -----------------------------
# Drone crown dataset (LIVE TREES ONLY)
# -----------------------------
drone_crowns <- st_read("E:/CHESS Campaign Data/CrownShp/TreeCrowns.shp") |>
  filter(Alive == "Y") |>  # keep live trees only
  select(Site_ID, Site_Number, Sampling_Area)

# -----------------------------
# Drone CHM file list
# -----------------------------
drone_files <- list.files(
  "E:/CHESS Campaign Data/CrownShp/Drone_CHMs",
  pattern = "\\.tif$",
  full.names = TRUE
)

# =====================================================
# DRONE EXTRACTION FUNCTION (USING Sampling_Area)
# =====================================================
extract_drone_heights <- function(site_code){
  
  # Filter crowns for domain using Sampling_Area prefix
  site_crowns <- drone_crowns |> 
    filter(grepl(paste0("^", site_code, ":"), Sampling_Area))
  
  if(nrow(site_crowns) == 0){
    message("No crowns found for domain ", site_code)
    return(NULL)
  }
  
  site_results <- list()
  pb <- progress_bar$new(
    format = "Extracting Drone CHM [:bar] :current/:total (:percent) in :elapsed",
    total = length(drone_files),
    clear = FALSE, width = 60
  )
  
  for(i in seq_along(drone_files)){
    pb$tick()
    
    chm <- rast(drone_files[i])
    
    crowns_proj <- st_transform(site_crowns, crs(chm))
    
    # Extract max CHM value within each polygon
    vals <- terra::extract(
      chm,
      vect(crowns_proj),
      fun = max,
      na.rm = TRUE
    )
    
    vals <- as.data.frame(vals)
    
    if(ncol(vals) > 1){
      colnames(vals)[2] <- "drone_ch"
      
      tmp <- crowns_proj |>
        st_drop_geometry() |>
        bind_cols(vals |> select(drone_ch))
      
      site_results[[i]] <- tmp
    }
  }
  
  if(length(site_results) == 0){
    message("No overlapping CHM tiles found for ", site_code)
    return(NULL)
  }
  
  # Combine all results and aggregate in case of overlapping tiles
  drone_df <- bind_rows(site_results) |>
    group_by(Site_ID, Site_Number, Sampling_Area) |>
    summarise(drone_ch = max(drone_ch, na.rm = TRUE), .groups = "drop")
  
  message("Drone extraction complete for ", nrow(drone_df), " crowns in ", site_code)
  
  return(drone_df)
}

# =====================================================
# DOMAIN PROCESSING FUNCTION
# =====================================================
process_domain <- function(domain, chm_path, crowns_path){
  
  cat("Processing domain:", domain, "\n")
  
  # Read crown polygons
  crowns <- st_read(crowns_path) |>
    rename(Site_Number = site_number)
  
  # Read AOP CHM raster
  chm <- rast(chm_path)
  
  # Join with field metadata
  trees <- left_join(crowns, trees_meta, by="Site_Number")
  trees_proj <- st_transform(trees, crs=crs(chm))
  
  # Crop CHM to convex hull
  hull <- trees_proj |> st_union() |> st_convex_hull()
  chm_crop <- crop(chm, vect(hull))
  
  # Rasterize polygons for zonal extraction
  trees_proj$cid <- seq_len(nrow(trees_proj))
  zones <- rasterize(vect(trees_proj), chm_crop, field="cid", background=NA, touches=TRUE)
  
  # Extract max AOP CHM
  aop_vals <- zonal(chm_crop, zones, fun="max", na.rm=TRUE) |> as.data.frame()
  colnames(aop_vals) <- c("cid","aop_ch")
  
  # Extract Drone CHM for the domain
  drone_vals <- extract_drone_heights(domain)
  
  # Combine into final dataframe
  df <- trees_proj |>
    st_drop_geometry() |>
    select(cid, Site_Number, Stem_Height) |>
    left_join(aop_vals, by="cid") |>
    left_join(drone_vals |> select(Site_Number, drone_ch), by="Site_Number")
  
  # Pivot for plotting
  plot_df <- df |>
    pivot_longer(
      cols = c(aop_ch, drone_ch),
      names_to = "CHM_Type",
      values_to = "CHM_Height"
    )
  
  # Plot
  p <- ggplot(plot_df, aes(x=Stem_Height, y=CHM_Height, color=CHM_Type)) +
    geom_point(alpha=0.7) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    geom_smooth(method="lm", se=FALSE) +
    labs(
      title = paste0(domain, ": Field vs AOP + Drone CHM"),
      x = "Field Stem Height (m)",
      y = "CHM Height (m)"
    ) +
    coord_equal() +
    theme_minimal()
  
  print(p)  # ensures plot shows
  
  return(df)
}

# =====================================================
# RUN FOR ALL DOMAINS
# =====================================================
almo_results <- process_domain(
  "ALMO",
  chm_path = "G:/.shortcut-targets-by-id/1CfLlW5BOhsk8a855sbGHpivEvMhBGaaY/CHESS_workshop_data/lidar_workshop_data/CHESS25_ALMO_CHM_1m_v1.tif",
  crowns_path = "G:/.shortcut-targets-by-id/1CfLlW5BOhsk8a855sbGHpivEvMhBGaaY/CHESS_workshop_data/lidar_workshop_data/CHESS_ALMO_Crowns_CHM-Aligned.geojson"
)

crbu_results <- process_domain(
  "CRBU",
  chm_path = "G:/.shortcut-targets-by-id/1CfLlW5BOhsk8a855sbGHpivEvMhBGaaY/CHESS_workshop_data/lidar_workshop_data/CHESS25_CRBU_CHM_1m_v2.tif",
  crowns_path = "G:/.shortcut-targets-by-id/1CfLlW5BOhsk8a855sbGHpivEvMhBGaaY/CHESS_workshop_data/lidar_workshop_data/CHESS_CRBU_Crowns_CHM-Aligned.geojson"
)

upta_results <- process_domain(
  "UPTA",
  chm_path = "G:/.shortcut-targets-by-id/1CfLlW5BOhsk8a855sbGHpivEvMhBGaaY/CHESS_workshop_data/lidar_workshop_data/CHESS25_UPTA_CHM_1m_v1.tif",
  crowns_path = "G:/.shortcut-targets-by-id/1CfLlW5BOhsk8a855sbGHpivEvMhBGaaY/CHESS_workshop_data/lidar_workshop_data/CHESS_UPTA_Crowns_CHM-Aligned.geojson"
)