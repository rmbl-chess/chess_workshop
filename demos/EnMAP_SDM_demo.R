## Script using EnMAP imaging spectroscopy data for species distribution modeling of RMBL wildflowers.

# install.packages("remotes","terra")
# remotes::install_github("rmbl-sdp/rSDP")
# library(rSDP)

library(terra) # Package for geospatial data manipulation.
library(sdm) # Package for species distribution modeling.
library(stringr) # Package for string manipulation.
library(sf) # Package for handling spatial data frames.

## Set working directory.
setwd("/Users/ian/Library/CloudStorage/GoogleDrive-ibreckhe@gmail.com/My\ Drive/BreckheimerLab2025/Projects/CHESS/")

#### ------ Prepare data, this was done ahead of time ------ ####

## Load catalog of Spatial Data Platform raster datasets.
#cat <- sdp_get_catalog()

## Grabs a landcover map for the study area.
#dem <- sdp_get_raster("R3D009", download_files=TRUE, download_path="~/Downloads")
#lc <- sdp_get_raster("R3D018", download_files=TRUE, download_path="~/Downloads")

## Computes proportion of meadow vegetation for every 30m pixel.
#lc_meadow <- lc == 3
#lc_paved <- lc == 7 | lc == 8
#lc_meadow_prop <- project(lc_meadow,en_mask_crop,method="average")
#lc_paved_prop <- project(lc_paved,en_mask_crop,method="average")
#lc_meadow_focal <- focal(lc_meadow_prop, w=matrix(1,5,5), fun="mean", na.rm=TRUE)
#lc_focal_thresh <- (lc_meadow_focal > 0.3 & lc_meadow_prop > 0.2 & lc_paved_prop < 0.05) * en_mask_crop

## Computes elevation map at 30m.
#dem_30m <- project(dem, en_mask_crop, method="average") * en_mask_crop

## Thresholds dem.
#dem_thresh <- dem_30m < 3700
#lc_focal_thresh2 <- lc_focal_thresh * dem_thresh
#lc_focal_thresh3 <- lc_focal_thresh2 / lc_focal_thresh2
#writeRaster(lc_focal_thresh3, 
#            filename="Data/CHESS_workshop_data/EnMAP_tutorial/Meadow_mask_30m.tif", 
#            overwrite=TRUE)

#### ------ Prepare EnMAP data for SDM ------ ####

## Loads the EnMAP mosaic for the study area.
enmap <- rast("Data/CHESS_workshop_data/EnMAP_tutorial/EnMAP_reflectance_clipped_2025_06_23_24_v3.tif")

## Loads the DEM for the study area.
dem_30m <- rast("Data/CHESS_workshop_data/EnMAP_tutorial/UG_DEM_30m.tif") * en_mask_crop

## Loads the meadow mask for the study area.
meadow_mask <- rast("Data/CHESS_workshop_data/EnMAP_tutorial/Meadow_mask_30m.tif")

## Creates a mask based on EnMAP grid.
en_mask <- enmap$EnMAP_reflectance_clipped_2025_06_23_24_v3_1 > -9999
mask_poly <- as.polygons(en_mask,round=FALSE)
mask_buff <- buffer(mask_poly,width=60)
mask_ext <- ext(mask_buff)
en_mask_crop <- crop(en_mask,mask_ext)

## Drops water bands with NA values.
enmap_drop <- enmap[[-c(131:135)]]

## Separates VNIR and SWIR bands.
enmap_vnir <- enmap_drop[[1:89]]
enmap_swir <- enmap_drop[[90:219]]

## Normalizes VNIR and SWIR bands by average VNIR and SWIR brightness.
enmap_brightness_vnir <- terra::mean(enmap_vnir, na.rm=TRUE)
enmap_norm_vnir <- enmap_vnir - enmap_brightness_vnir
names(enmap_norm_vnir) <- paste0("band",1:nlyr(enmap_vnir))

enmap_brightness_swir <- terra::mean(enmap_swir, na.rm=TRUE)
enmap_norm_swir <- enmap_swir - enmap_brightness_swir
names(enmap_norm_swir) <- paste0("band",1:nlyr(enmap_swir))

## Masks EnMAP brightness to meadow areas.
enmap_vnir_meadow <- mask(crop(enmap_norm_vnir, lc_focal_thresh3), lc_focal_thresh3)
names(enmap_vnir_meadow) <- paste0("band",1:nlyr(enmap_vnir_meadow))

enmap_swir_meadow <- mask(crop(enmap_norm_swir, lc_focal_thresh3), lc_focal_thresh3)
names(enmap_swir_meadow) <- paste0("band",1:nlyr(enmap_swir_meadow))

## Runs a principal components analysis on all meadow pixels (takes a few minutes).
meadow_pca_vnir <- prcomp(enmap_vnir_meadow, center=TRUE, scale.=TRUE)
meadow_pca_swir <- prcomp(enmap_swir_meadow, center=TRUE, scale.=TRUE)

enmap_rotated_vnir <- predict(enmap_norm_vnir, meadow_pca_vnir)
enmap_rotated_swir <- predict(enmap_norm_swir, meadow_pca_swir)

writeRaster(enmap_rotated_vnir[[1:20]], 
            filename="Data/CHESS_workshop_data/EnMAP_tutorial/EnMAP_VNIR_PCA_meadow_20bands.tif", 
            overwrite=TRUE)
writeRaster(enmap_rotated_swir[[1:20]], 
            filename="Data/CHESS_workshop_data/EnMAP_tutorial/EnMAP_SWIR_PCA_meadow_20bands.tif", 
            overwrite=TRUE)

#### ------ Fit SDMs and make predictions ------ ####

## Re-loads the rasters
preds_enmap_vnir <- rast("Data/CHESS_workshop_data/EnMAP_tutorial/EnMAP_VNIR_PCA_meadow_20bands.tif")
preds_enmap_swir <- rast("Data/CHESS_workshop_data/EnMAP_tutorial/EnMAP_SWIR_PCA_meadow_20bands.tif")
elev <- rast("Data/CHESS_workshop_data/EnMAP_tutorial/UG_DEM_30m.tif")

## Loads the observation data for the target species.
obs <- read.csv("Data/CHESS_workshop_data/EnMAP_tutorial/iNat_RMBL_wildflower_observations_2014_2025.csv")

## Creates a predictor stack for SDM fitting and prediction.
preds <- c(crop(preds_enmap_vnir[[1:3]], elev), crop(preds_enmap_swir[[1:3]], elev), elev)
names(preds) <- c("VNIR_PC1", "VNIR_PC2", "VNIR_PC3",
                   "SWIR_PC1", "SWIR_PC2", "SWIR_PC3",
                   "elev")

## Filters iNat observations by geographic accuracy.
obs <- obs[obs$positional_accuracy <= 30, ]
obs$genus <- str_split_fixed(obs$scientific_name, " ", 3)[, 1]
obs$species <- str_split_fixed(obs$scientific_name, " ", 3)[, 2]
obs$binomial <- paste(obs$genus, obs$species, sep = " ")
obs <- obs[!is.na(obs$genus), ]

## Create spatial points and reproject to predictor CRS.
obs_sf <- st_as_sf(obs, coords = c("longitude", "latitude"), crs = 4326)
obs_sf <- st_transform(obs_sf, crs(preds))

## Filter to observations within the predictor extent.
obs_sf$elev <- terra::extract(elev, obs_sf)[, 2]
obs_sf <- obs_sf[!is.na(obs_sf$elev), ]

## Target species for modeling.
target_species <- c("Delphinium barbeyi",
                    "Erythronium grandiflorum",
                    "Eriogonum umbellatum")

## Split into target presences and background absences.
obs_target <- obs_sf[obs_sf$binomial %in% target_species, ]
obs_target$presence <- 1
obs_sf$presence <- 0

## SDM methods to fit.
sdm_methods <- c("gbm", "glm")

## Create output directory if needed.
dir.create("output", showWarnings = FALSE)

## Store fitted model objects for evaluation/plotting.
all_models <- list()

## Loop through each species and fit models.
for (sp in target_species) {

  sp_label <- gsub(" ", "_", sp)
  message("Fitting models for: ", sp)

  ## Subset presence and absence data for this species.
  pres_sp <- obs_target[obs_target$binomial == sp, ]
  abs_sp  <- obs_sf[obs_sf$binomial != sp, ]

  ## Combine presence and absence data.
  data_sp <- vect(rbind(pres_sp, abs_sp))
  extract_sp <- terra::extract(preds, data_sp)
  data_sp <- cbind(data_sp$presence, extract_sp[, -1])
  names(data_sp)[1] <- "presence"

  ## Create sdmData object and fit models.
  model_data_sp <- sdmData(presence ~ ., train = data_sp)

  med_list <- list()
  for (method in sdm_methods) {
    message("  Fitting ", method, "...")
    model_sp <- sdm(data = model_data_sp,
                    methods = method, test.percent = 20,
                    replication = "cv", cv.folds = 5)

    pred_sp <- predict(model_sp, newdata = preds)
    med_sp  <- app(pred_sp, median)
    med_list[[method]] <- med_sp

    writeRaster(med_sp,
                file.path("output", paste0(sp_label, "_", method, "_median.tif")),
                overwrite = TRUE)
  }

  ## Compute ensemble mean across methods.
  med_stack <- rast(med_list)
  med_ensemble <- app(med_stack, mean)
  writeRaster(med_ensemble,
              file.path("output", paste0(sp_label, "_ensemble_mean.tif")),
              overwrite = TRUE)
  message("  Wrote ensemble prediction for ", sp)

  ## Store the last fitted model object for evaluation/plotting below.
  all_models[[sp]] <- model_sp
}

#### ------ Visualize SDM outputs and performance ------ ####

## Plot ensemble prediction maps side by side.
par(mfrow = c(1, length(target_species)), mar = c(2, 2, 3, 4))
for (sp in target_species) {
  sp_label <- gsub(" ", "_", sp)
  r <- rast(file.path("output", paste0(sp_label, "_ensemble_mean.tif")))
  plot(r,
       main = sp,
       col = hcl.colors(50, "Lajolla"),
       range = c(0, 1))
}

## Collect evaluation statistics across all models and species.
eval_list <- list()

for (sp in target_species) {
  if (is.null(all_models[[sp]])) next
  ev <- getEvaluation(all_models[[sp]],
                      stat = c("AUC", "TSS", "Kappa", "COR"),
                      opt = 2)
  ev$species <- sp
  eval_list[[sp]] <- ev
}

eval_df <- do.call(rbind, eval_list)
rownames(eval_df) <- NULL

eval_df

## Plot ROC curves for each species.
par(mfrow = c(1, length(target_species)), mar = c(4, 4, 3, 1))
for (sp in target_species) {
  if (is.null(all_models[[sp]])) next
  roc(all_models[[sp]], main = sp)
}

