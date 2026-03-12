# Helper functions

#' Load packages
#' @description Loads new packages, installing if they're not already installed
#' @param pkg Character string. Package name
#' @return NULL. Loads packages in global environment
#' @export load.pkgs
#'
load.pkgs <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#' Ingest generic from Drive
#' @param id Character string. Unique identifier of file
#' @export ingest.drive.tif
#' 
ingest.drive <- function(id) {
  
  f.tmp <- drive_download(as_id(id),
                          path=file.path(tempdir(), id),
                          overwrite=T)
  
  return(f.tmp)
  
}

#' Ingest tif from Drive
#' @param id Character string. Unique identifier of file
#' @export ingest.drive.tif
#' 
ingest.drive.tif <- function(id) {
  
  f.tmp <- drive_download(as_id(id),
                          path=file.path(tempdir(), id),
                          overwrite=T)
  
  f <- rast(f.tmp$local_path)
  names(f) <- f.tmp$name
  
  return(f)
  
}

#' Ingest shapefile from Drive
#' @export ingest.drive.shp
#'
ingest.drive.shp <- function(id, pattern) {
  
  shpobj <- drive_ls(
    as_id(id)
  )
  
  shpobj <- shpobj[grepl(pattern, shpobj$name),]
  
  tmpfiles <- apply(shpobj, 1, function(x) {
    if(!file_ext(x[['name']]) %in% c('xml', 'kmz')) {
      tmpfile <- drive_download(
        as_id(x[['id']]),
        path=file.path(
          tempdir(),
          x[['name']]),
        overwrite=T)$local_path}
    else tmpfile <- NULL
    return(tmpfile)
  })
  
  shpfile <- tmpfiles[file_ext(tmpfiles)=='shp']
  shp <- st_read(shpfile)
  shp <- st_transform(shp, 'EPSG:32613')
  
  return(shp)
}

# Safe read csv with error handling

#'@export saferead.csv
saferead.csv <- function(infile){
  tryCatch(fread(infile), 
           error = function(cond) {
             message(paste('Reading csv failed'))
           })
}

# Custom ggplot stuff

## Geom defaults
ggplot2::update_geom_defaults("line", list(linewidth = 1))
ggplot2::update_geom_defaults("point", list(shape=16, size=4, stroke=2))

## Matlab similar theme
theme_matlab <- function() {
  theme_base(base_size=18) + 
    theme(legend.position='inside', 
          legend.position.inside=c(.1,.9), 
          axis.line=element_line(linewidth=1), 
          axis.ticks=element_line(linewidth=1) #, 
          # axis.title = element_text(size = 12, face = "bold"),
          # axis.text = element_text(size = 10),
          )
}
