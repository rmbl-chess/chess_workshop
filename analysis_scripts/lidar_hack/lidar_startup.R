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

# Join tree geospatial and demographic data
trees <- left_join(trees.poly1 %>% dplyr::select(-Sampling_Area), trees.demog, by='Site_Number')

# LAS catalog
lascat <- readLAScatalog(file.path(config$extdata$pclocalpath))
plot(st_geometry(neon.grid))
plot(lascat, col='dodgerblue', add=T)

# CHESS sampling area polygons
chess.areas <- drive_download(as_id(config$extdata$chessareasid), path=file.path(tempdir(), 'areas.geojson'), overwrite=T)
chess.areas <- st_read(chess.areas$local_path)
chess.areas <- st_transform(chess.areas, 'EPSG:32613')

# NEON-delivered DTM and CHM
neon_chm <- ingest.drive.tif(config$extdata$chmid)
neon_dtm <- ingest.drive.tif(config$extdata$dtmid)

