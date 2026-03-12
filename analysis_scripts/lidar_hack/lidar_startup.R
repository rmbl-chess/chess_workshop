# Source config file
config <- config::get(file=file.path('config', 'config.yml'))

# Source helpers and other package functions
lapply(list.files(file.path('.', 'R'), full.names=T), source)

# Load packages
load.pkgs(config$pkgs)

# Config drive auth
drive_auth(email='hmworsham@lbl.gov')



neon.grid <- st_read(file.path(config$data, 'intermediate', 'neon_2025_lidar_utm_grid.geojson'))


# Ingest tree stuff
trees.demog <- read.csv(ingest.drive(config$extdata$treedemoid)$local_path) # yields a dataframe
trees.geo <- st_read(ingest.drive(config$extdata$treegeoid)$local_path) # yields an sf object with geometry and attributes

# Join geospatial and demographic data
trees <- left_join(trees.geo %>% dplyr::select(-Sampling_Area), trees.demog, by='Site_Number')


x.las <- readLAS(ingest.drive(config$extdata$xlasid))