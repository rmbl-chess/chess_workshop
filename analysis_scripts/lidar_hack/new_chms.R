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

### --- Ingest --- ###

# LAS catalog
lascat <- readLAScatalog(file.path(config$extdata$pclocalpath))

# CHESS sampling area polygons
chess.areas <- drive_download(as_id(config$extdata$chessareasid), path=file.path(tempdir(), 'areas.geojson'), overwrite=T)
chess.areas <- st_read(chess.areas$local_path)
chess.areas <- st_transform(chess.areas, 'EPSG:32613')


make.new.chm <- function(chunk) {

  l = readLAS(chunk)
  if (is.empty(l)) return(NULL)        # check if it actually contain points
  cat('clipping')

  #l <- clip_roi(lascat, st_buffer(sampling.area, 100))
  cat('clipped')

  # l@data$Classification <- 0
  # ws <- seq(3, 15, 3)
  # th <- seq(0.1, 3, length.out = length(ws))
  # cat('classifying')
  # l <- classify_ground(l, algorithm = pmf(ws = ws, th = th))
  # cat('classified')

  cat('normalizing')
  l <- normalize_height(l, algorithm=knnidw(k=21L, p=2))
  cat('normalized')
  l <- l[l$Z >= 0]
  cat('rasterizing canopy')
  chm <- rasterize_canopy(l, res=0.25, algorithm=p2r(0.2), na.fill=tin())
  rm(l)
  gc()
  cat('complete')
  #terra::writeRaster(chm, file.path('~', 'Desktop', 'CHESS_CHM', paste0(sampling.area$SamplingAreaCode, '_CHM.tif')))
  return(chm)
  }

future::plan(future::multisession, workers=6L)
opt_chunk_buffer(lascat) <- 30
opt_stop_early(lascat) <- T
chms <- catalog_apply(lascat, make.new.chm)

chms <- mclapply(1:nrow(chess.areas), \(i) make.new.chm(chess.areas[i,], lascat), mc.cores=getOption('mc.cores', 6L))
