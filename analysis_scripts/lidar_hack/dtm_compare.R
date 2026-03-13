# CHESS Workshop 2026
# Startup script for CHESS LiDAR data analysis
# Sets up workspace and loads some critical data

### --- Workspace setup --- ###

# Source config file
setwd(file.path('~', 'Repos', 'rmbl-chess', 'chess_workshop', 'analysis_scripts', 'lidar_hack'))
config <- config::get(file=file.path('config.yml'))

# Source helpers and other package functions
lapply(list.files(file.path('.', 'R'), full.names=T), source)

# Load packages
load.pkgs(config$pkgs)

# Config drive auth
drive_auth(email='hmworsham@lbl.gov')

# NEON grid
neon.grid <- st_read(file.path(config$data, 'intermediate', 'neon_2025_lidar_utm_grid.geojson'))

# Ingest tree stuff
trees.demog <- read.csv(ingest.drive(config$extdata$treedemoid)$local_path) # yields a dataframe
trees.poly1 <- st_read(ingest.drive(config$extdata$treegeoid1)$local_path) # yields an sf object with geometry and attributes
trees.poly2 <- st_read(ingest.drive(config$extdata$treegeoid2)$local_path) # yields an sf object with geometry and attributes

shrub.poly <-

# Join tree geospatial and demographic data
trees <- left_join(trees.poly1 %>% dplyr::select(-Sampling_Area), trees.demog, by='Site_Number')

# LAS catalog
lascat <- readLAScatalog(file.path(config$extdata$pclocalpath))
plot(st_geometry(neon.grid))
plot(lascat, col='dodgerblue', add=T)
plot(chess.areas, col='red', add=T)

# CHESS sampling area polygons
chess.areas <- drive_download(as_id(config$extdata$chessareasid), path=file.path(tempdir(), 'areas.geojson'), overwrite=T)
chess.areas <- st_read(chess.areas$local_path)
chess.areas <- st_transform(chess.areas, 'EPSG:32613')

trees.poly1 %>%
  group_by(Sampling_Area) %>%
  summarise(n=n())

trees.poly2 %>%
  group_by(Sampling_A) %>%
  summarise(n=n())

# NEON-delivered DTM and CHM
neon_crbu_chm <- rast('/Users/hmworsham/Library/CloudStorage/GoogleDrive-hmworsham@lbl.gov/.shortcut-targets-by-id/1I9nIoDCDNpJ3KQLKR9lCGFKdg533HaDQ/CHESS/Data/CHESS_workshop_data/lidar_workshop_data/CHESS25_CRBU_CHM_1m_v2.tif')
neon_crbu_dtm <- rast('/Users/hmworsham/Library/CloudStorage/GoogleDrive-hmworsham@lbl.gov/.shortcut-targets-by-id/1I9nIoDCDNpJ3KQLKR9lCGFKdg533HaDQ/CHESS/Data/CHESS_workshop_data/lidar_workshop_data/CHESS25_CRBU_DTM_1m_v2.tif')

neon_almo_chm <- rast('/Users/hmworsham/Library/CloudStorage/GoogleDrive-hmworsham@lbl.gov/.shortcut-targets-by-id/1I9nIoDCDNpJ3KQLKR9lCGFKdg533HaDQ/CHESS/Data/CHESS_workshop_data/lidar_workshop_data/CHESS25_ALMO_CHM_1m_v1.tif')
neon_almo_dtm <- rast('/Users/hmworsham/Library/CloudStorage/GoogleDrive-hmworsham@lbl.gov/.shortcut-targets-by-id/1I9nIoDCDNpJ3KQLKR9lCGFKdg533HaDQ/CHESS/Data/CHESS_workshop_data/lidar_workshop_data/CHESS25_ALMO_DTM_1m_v1.tif')

neon_upta_chm <- rast('/Users/hmworsham/Library/CloudStorage/GoogleDrive-hmworsham@lbl.gov/.shortcut-targets-by-id/1I9nIoDCDNpJ3KQLKR9lCGFKdg533HaDQ/CHESS/Data/CHESS_workshop_data/lidar_workshop_data/CHESS25_UPTA_CHM_1m_v1.tif')
neon_upta_dtm <- rast('/Users/hmworsham/Library/CloudStorage/GoogleDrive-hmworsham@lbl.gov/.shortcut-targets-by-id/1I9nIoDCDNpJ3KQLKR9lCGFKdg533HaDQ/CHESS/Data/CHESS_workshop_data/lidar_workshop_data/CHESS25_UPTA_DTM_1m_V1.tif')

# Clip LAS and write temp
las_clipped <- clip_roi(lascat, st_buffer(chess.areas[14,],50))
writeLAS(las_clipped, '~/Downloads/classified_point_cloud_NB.laz')

# CHM clipped
chm_clipped <- crop(neon_crbu_chm, st_buffer(chess.areas[14,], 50))
names(chm_clipped) <- 'Z'

# Check clipped las
plot(las_clipped)
hist(las_clipped$Z)
length(which(las_clipped$Z > 3485))

# Remove noise
las_clipped_denoised <- las_clipped[las_clipped$Z < 3485]
las_clipped_denoised <- las_clipped_denoised[!las_clipped_denoised$Classification==7]

# Check
class.table1 <- table(las_clipped$Classification)

# Add zI
las_clipped_denoised <- add_attribute(las_clipped_denoised, scale(las_clipped_denoised@data$Intensity), name='zI')

# Normalize
las_norm <- normalize_height(las_clipped, algorithm=tin())

# Check
hist(las_norm$Z)
summary(las_norm$Z)
length(which(las_norm$Z < 0))/length(las_norm@data$Z)
las_norm$Classification <- ifelse(las_norm$Z < 0, 7L, las_norm$Classification)

# Locate trees
# f <- function(x, y, z) {x * 0.07 + y * 0.01 + 2}
f <- function(x) {x * 0.05 + 3}
# ws_args <- list(x = "Z", y = "zI")
ws_args   <- list(x = "Z")
las_trees <- locate_trees(las_norm, algorithm=lmf(f, ws_args = ws_args, hmin=2))
las_trees <- locate_trees(las_norm, algorithm=lmf(ws=5, hmin=2))

# Check: how many trees?
print(las_trees)
print(length(unique(las_trees$treeID)))

# Plot with CHM
ggplot() +
  geom_spatraster(data=chm_clipped, aes(fill=Z)) +
  scale_fill_distiller(palette='Spectral', na.value=NA) +
  geom_sf(data=las_trees, shape=3, size=0.8) +
  theme_classic()

# Plot height histogram
ggplot(las_trees) +
  geom_histogram(aes(x=Z), fill='dodgerblue', color='blue3') +
  theme_classic()

# Segment trees
las_segment <- lidR::segment_trees(las_norm[!las_norm$Classification == 7], algorithm=dalponte2016(chm=neon_crbu_chm, treetops=las_trees), attribute='treeID')
#plot(las_segment, bg = "white", size = 4, color = "treeID") # visualize trees

writeLAS(las_segment, '~/Downloads/NB_pointcloud_segmented.laz')

# Convex hull
las_crowns <- delineate_crowns(las_segment, 'convex')

plot(las_crowns)
plot(trees.poly1, add=T, col='red')

st_write(st_as_sf(las_crowns), '~/Downloads/NB_delineated_crowns.geojson')

#################

make.new.chm <- function(sampling.area, lascat) {

  cat('clipping')
  l <- clip_roi(lascat, st_buffer(sampling.area, 200))
  cat('clipped')
  l@data$Classification <- 0

  ws <- seq(3, 15, 3)
  th <- seq(0.1, 3, length.out = length(ws))
  cat('classifying')
  l <- classify_ground(l, algorithm = pmf(ws = ws, th = th))
  cat('classified')
  cat('normalizing')
  l <- normalize_height(l, algorithm=knnidw(k=21L, p=2))
  cat('normalized')
  l <- l[l$Z >= 0]
  cat('rasterizing canopy')
  chm <- rasterize_canopy(las_norm_dtm, res=0.25, algorithm=p2r(0.2), na.fill=tin())
  cat('complete')
  chm
  }

jj <- make.new.chm(chess.areas[1,], x.lascat)

# Remove classification
las_clipped_noclass <- las_clipped_denoised
las_clipped_noclass@data$Classification <- 0


class.table2 <- table(las_clipped_noclass$Classification)
class.table2







#################

library(lidR)
library(terra)
library(sf)

#' Process a single CHESS sampling area: clip, denoise, normalize,
#' segment trees, and delineate crowns.
#'
#' @param area_index    Integer. Row index of the target polygon in chess.areas.
#' @param chess.areas   sf object with 65 sampling area polygons.
#' @param lascat        LAS catalog covering all domains.
#' @param chm_list      Named list of SpatRaster CHMs, one per domain.
#'                      Names must match values in domain_col (e.g.,
#'                      list(CRBU = neon_crbu_chm, ALMO = neon_almo_chm,
#'                           UPTA = neon_upta_chm)).
#' @param area_id_col   Character or NULL. Column used as the area label in
#'                      output filenames. If NULL, the row index is used.
#' @param buffer_dist   Numeric. Buffer radius (m) applied before clipping.
#' @param z_noise_max   Numeric. Absolute Z ceiling for noise removal.
#' @param output_dir    Character. Directory for .laz and .geojson outputs.
#' @param ws            Numeric. Fixed window size for lmf() treetop detection.
#' @param hmin          Numeric. Minimum tree height for lmf().
#' @param write_outputs Logical. Write LAZ and GeoJSON files to disk.
#'
#' @return Invisibly, a list: las_segment, crowns (sf), trees (sf).
process_chess_area <- function(
    area_index,
    chess.areas,
    lascat,
    chm_list,
    area_id_col   = "SamplingAreaCode",
    buffer_dist   = 50,
    z_noise_max   = 3485,
    output_dir    = ".",
    ws            = 5,
    hmin          = 2,
    write_outputs = TRUE
) {

  # ── Identifiers & paths ─────────────────────────────────────────────────────
  area    <- chess.areas[area_index, ]
  domain  <- toupper(trimws(area[[domain_col]]))
  area_id <- if (is.null(area_id_col)) area_index else area[[area_id_col]]

  # base_name    <- paste0(domain, "_", area_id)
  base_name    <- paste0(domain, '_', area_id)
  laz_path     <- file.path(output_dir, paste0(base_name, "_segmented.laz"))
  geojson_path <- file.path(output_dir, paste0(base_name, "_crowns.geojson"))

  message(sprintf("[%d / %d] %s", area_index, nrow(chess.areas), base_name))

  # ── Validate domain ──────────────────────────────────────────────────────────
  if (!domain %in% names(chm_list)) {
    stop(sprintf(
      "Domain '%s' (area %d) not in chm_list. Available: %s",
      domain, area_index, paste(names(chm_list), collapse = ", ")
    ))
  }
  chm <- chm_list[[domain]]

  # ── Clip LAS and CHM to buffered area ────────────────────────────────────────
  area_buf    <- st_buffer(area, buffer_dist)
  las_clipped <- clip_roi(lascat, area_buf)
  chm_clipped <- crop(chm, area_buf)   # retained if needed for visualisation
  names(chm_clipped) <- "Z"

  # ── Denoise ──────────────────────────────────────────────────────────────────
  las_denoised <- las_clipped[las_clipped$Z < z_noise_max]
  las_denoised <- las_denoised[las_denoised$Classification != 7L]
  las_denoised <- add_attribute(
    las_denoised,
    scale(las_denoised@data$Intensity),
    name = "zI"
  )

  # ── Normalize height ──────────────────────────────────────────────────────────
  # NOTE: normalization is applied to the denoised cloud (original code
  # applied it to las_clipped; adjust here if that behaviour is needed).
  las_norm <- normalize_height(las_denoised, algorithm = tin())
  las_norm$Classification <- ifelse(las_norm$Z < 0, 7L, las_norm$Classification)

  # ── Detect treetops ───────────────────────────────────────────────────────────
  las_trees <- locate_trees(las_norm, algorithm = lmf(ws = ws, hmin = hmin))
  message(sprintf("  treetops detected: %d", length(unique(las_trees$treeID))))

  # ── Segment trees ─────────────────────────────────────────────────────────────
  las_segment <- lidR::segment_trees(
    las_norm[las_norm$Classification != 7L],
    algorithm = dalponte2016(chm = chm, treetops = las_trees),
    attribute = "treeID"
  )

  # ── Crown delineation ─────────────────────────────────────────────────────────
  las_crowns <- delineate_crowns(las_segment, type = "convex")

  # ── Write outputs ─────────────────────────────────────────────────────────────
  if (write_outputs) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    writeLAS(las_segment, laz_path)
    sf::st_write(las_crowns, geojson_path, driver = "GeoJSON", delete_dsn = TRUE)
    message(sprintf("  -> %s", laz_path))
    message(sprintf("  -> %s", geojson_path))
  }

  invisible(list(
    las_segment = las_segment,
    crowns      = las_crowns,
    trees       = las_trees
  ))
}


#' Batch-process CHESS sampling areas
#'
#' Calls process_chess_area() for every row in chess.areas (or a chosen
#' subset).  Errors on individual areas are caught and reported without
#' stopping the loop.
#'
#' @param chess.areas sf object with 65 sampling area polygons.
#' @param lascat      LAS catalog covering all domains.
#' @param chm_list    Named list of SpatRaster CHMs (see process_chess_area).
#' @param indices     Integer vector of row indices to process.
#'                    Defaults to all rows.
#' @param ...         Additional arguments forwarded to process_chess_area().
#'
#' @return Invisibly, a named list of per-area results (NULL on error).
process_all_chess_areas <- function(
    chess.areas,
    lascat,
    chm_list,
    indices = NULL,
    ...
) {
  if (is.null(indices)) indices <- seq_len(nrow(chess.areas))

  results <- lapply(indices, function(i) {
    tryCatch(
      process_chess_area(
        area_index  = i,
        chess.areas = chess.areas,
        lascat      = lascat,
        chm_list    = chm_list,
        ...
      ),
      error = function(e) {
        message(sprintf("  ERROR area %d: %s", i, conditionMessage(e)))
        NULL
      }
    )
  })

  names(results) <- as.character(indices)
  invisible(results)
}




# Down the line...


