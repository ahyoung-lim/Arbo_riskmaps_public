# ==========================================================
# 
# Masking output maps and converting them into binary maps
# 
# ==========================================================

source("script/00_setup.R")

# Masking ==================================================
# Load in unmaksed raster maps
Ras_DEN <- raster("outputs/Rasters/DEN_riskmap_wmean_unmasked.tif")
Ras_ZIK <- raster("outputs/Rasters/ZIK_riskmap_wmean_unmasked.tif")
Ras_CHIK <- raster("outputs/Rasters/CHIK_riskmap_wmean_unmasked.tif")
Ras_YF <- raster("outputs/Rasters/YF_riskmap_wmean_unmasked.tif")

# load in masking layers (temperature suitability and yf risk classification map)
tsuit <- raster("data/covariate_rasters/Dengue_temperature_suitaiblity_masked_.tif")
yf_eye <- raster("data/covariate_rasters/YF_EYE_mask.tif")

# convert temperature suitability to 0-1 index
# following: https://www.nature.com/articles/s41564-019-0476-8 
tsuit = tsuit / max(as.vector(tsuit), na.rm = T)

# temperature suitaiblity masking - 0.04 = 2 weeks of suitability = 1 serial interval
mask_ts <- (tsuit > 0.04)
mask_yf = yf_eye > 0 


# masking
DEN_range_mask <- Ras_DEN * mask_ts
CHIK_range_mask <- Ras_CHIK * mask_ts
ZIK_range_mask <- Ras_ZIK * mask_ts
YF_range_mask <- Ras_YF * mask_yf


# Collapsing dengue, chikungunya, and zika maps
# DCZ maps
dcz <- stack(DEN_range_mask, CHIK_range_mask, ZIK_range_mask)

dcz_range <- terra::mean(dcz)

# saving outputs
writeRaster(DEN_range_mask, filename=paste0("outputs/Rasters/DEN_riskmap_wmean_masked.tif"), overwrite=T)
writeRaster(CHIK_range_mask, filename=paste0("outputs/Rasters/CHIK_riskmap_wmean_masked.tif"), overwrite=T)
writeRaster(ZIK_range_mask, filename=paste0("outputs/Rasters/ZIK_riskmap_wmean_masked.tif"), overwrite=T)
writeRaster(YF_range_mask, filename=paste0("outputs/Rasters/YF_riskmap_wmean_masked.tif"), overwrite=T)
writeRaster(dcz_range,  filename="outputs/Rasters/DCZ_riskmap_wmean_masked.tif", overwrite=T)




# Calculating 95% confidence interval ======================

path = "C:/Users/user/Dropbox/WORK/WHO_GAI/ENMs/Multi_arbo_mapping/"

source("script/01_get_covariates.R")
# load in covariate rasters 
arbo_cov_list <- loadRasters("arbo")[[1]]
# make prediction data frame
# a standardised set of covariates (5 km x 5 km global scale)
pred.data <- as.data.frame(arbo_cov_list, xy=T) 
rm("arbo_cov_list"); gc()


ncpus=ncpus-2
disease_names = c("dengue", "chikungunya", "zika", "yf")

lapply(1:4, function(i){
  disease = disease_names[i]
  dname = ifelse(disease =="dengue", "DEN", ifelse(disease == "zika", "ZIK", ifelse(disease == "chikungunya", "CHIK", "YF")))
  
  # load 100 predictions and weights (AUC values)
  # saved in separate private repository 
  preds <- readRDS(paste0(path, "cross_validation/", dname, "_100fold_pred.rds"))
  aucs <- readRDS(paste0(path, "cross_validation/", dname, "_100fold_AUC.rds"))
  
  # preds <- readRDS(paste0("outputs/cross_validation/", dname, "_100fold_pred.rds"))
  # aucs <- readRDS(paste0("outputs/cross_validation/", dname, "_100fold_AUC.rds"))
    
  auc = unlist(aucs)
  weights = auc / sum(auc)
  
  # reduced size preds
  NAindex = which(rowSums(is.na(preds[, !colnames(preds) %in% c("wmean", "median", "IQR")])) == 100)
  pred_small <- preds[-NAindex,] %>% select(X1:X100)
  
  
  # register parallel backend
  cl <- makeCluster(ncpus) 
  registerDoParallel(cl) 
  set.seed(1259)
  
  # calculate row-wise confidence intervals
  tic(); wpreds <- rowCIs(pred_small, weights, 1000); toc()
  stopCluster(cl); gc()
  
  preds$lwr <- NA
  preds$upr <- NA
  
  preds[-NAindex, "lwr"] <- wpreds[,1]
  preds[-NAindex, "upr"] <- wpreds[,2]
  
  # save outputs
  writeRaster(bootsRas(preds, "lwr"), filename=paste0("outputs/Rasters/", dname, "_riskmap_lwr_masked.tif"), overwrite=T)
  writeRaster(bootsRas(preds, "upr"), filename=paste0("outputs/Rasters/", dname, "_riskmap_upr_masked.tif"), overwrite=T)
  rm(preds, aucs, NAindex, pred_small, wpreds); gc()
})


# Binary maps ==============================================
# load in rasters
DEN_ex <- raster(paste0(path, "published_datasets/DENV_Messina/dengue_2015_mean_mask.gri"))
CHIK_ex <- raster(paste0(path, "published_datasets/CHIKV_Nsoesie/CHIKV_prediction_map_Kraemer.tif"))
ZIK_ex <- raster(paste0(path, "published_datasets/ZIKV_Messina/ZIKV_prediction_map_Messina.tif"))
YF_ex <- raster(paste0(path, "published_datasets/YFV_Shearer/prediction_who_extent_.tif"))

DEN_range_mask <- raster("outputs/Rasters/DEN_riskmap_wmean_masked.tif")
CHIK_range_mask <- raster("outputs/Rasters/CHIK_riskmap_wmean_masked.tif")
ZIK_range_mask <- raster("outputs/Rasters/ZIK_riskmap_wmean_masked.tif")
YF_range_mask <- raster("outputs/Rasters/YF_riskmap_wmean_masked.tif")

crs(ZIK_ex) <- "+proj=longlat +datum=WGS84 +no_defs"

DEN_ex <- crop(DEN_ex, DEN_range_mask)
CHIK_ex <- crop(CHIK_ex, DEN_range_mask)
ZIK_ex <- crop(ZIK_ex, DEN_range_mask)

# adding back in GAUL code to the original arbovirus occurrence points
# thinning separately on each arboviral disease
source("script/01_get_data_arbo_model.R")

cnames <- c("Longitude", "Latitude", "Admin", "disease")
arbo_occ <- rbind(den[, cnames],
                  chik[, cnames],
                  zik[, cnames],
                  yf[, cnames])


# random background point generation
set.seed(1221)
arbo_bg <- thin_bg_arbo(arbo_occ)
nrow(arbo_bg)#[1] 13127

# combining presence and background
arbo_pres <- rbind(arbo_occ,
                   data.frame(arbo_bg[, c("Longitude", "Latitude", "disease")],
                              Admin = -999))
arbo_pres$PA <- as.factor(c(rep(1, nrow(arbo_occ)), rep(0, nrow(arbo_bg))))

points_locs <- data.frame(a0 = extract(admin$ad0, arbo_pres[, c("Longitude", "Latitude")]),
                          a1 = extract(admin$ad1, arbo_pres[, c("Longitude", "Latitude")]),
                          a2 = extract(admin$ad2, arbo_pres[, c("Longitude", "Latitude")]))
points_locs = points_locs[(arbo_pres$Admin != -999), ]
arbo_pres$GAUL = NA
arbo_pres$Admin <- as.numeric(arbo_pres$Admin)
v = cbind(1:nrow(points_locs),
          (arbo_pres$Admin[(arbo_pres$Admin != -999)] + 1))
arbo_pres$GAUL[(arbo_pres$Admin != -999)] = points_locs[v]

# Extract map values based occ+bg points
den_pres <- mapExtract(arbo_pres, "dengue", "DEN" )
chik_pres <- mapExtract(arbo_pres, "chik", "CHIK")
zik_pres <- mapExtract(arbo_pres, "zik", "ZIK")
yf_pres <- mapExtract(arbo_pres, "yf", "YF")

saveRDS(den_pres, paste0("data/intermediate_datasets/DEN_occ_map_extract.rds"))
saveRDS(chik_pres, paste0("data/intermediate_datasets/CHIK_occ_map_extract.rds"))
saveRDS(zik_pres, paste0("data/intermediate_datasets/ZIK_occ_map_extract.rds"))
saveRDS(yf_pres, paste0("data/intermediate_datasets/YF_occ_map_extract.rds"))

# load in datasets
# previously published map and new map values were extracted based occ+bg points 
den_pres <- readRDS("data/intermediate_datasets/DEN_occ_map_extract.rds")
chik_pres <- readRDS("data/intermediate_datasets/CHIK_occ_map_extract.rds")
zik_pres <- readRDS("data/intermediate_datasets/ZIK_occ_map_extract.rds")
yf_pres <- readRDS("data/intermediate_datasets/YF_occ_map_extract.rds")

# Calculate threshold and return binary maps
den_bin <- mapROC(den_pres, "DEN", "new_only")$bin_map 
chik_bin <- mapROC(chik_pres, "CHIK", "new_only")$bin_map 
zik_bin <- mapROC(zik_pres, "ZIK",  "new_only")$bin_map 

# collapsing dengue, chikungunya, and zika maps
dcz_bin <- calc(stack(den_bin, chik_bin, zik_bin), sum)
dcz_bin[dcz_bin > 0] <- 1

yf_bin <- mapROC(yf_pres, "YF",  "new_only")$bin_map 

writeRaster(dcz_bin, "outputs/Rasters/DCZ_binmap_mean.tif", overwrite=T)
writeRaster(yf_bin, "outputs/Rasters/YF_binmap_mean.tif", overwrite=T)


summary(values(den_bin))
summary(values(chik_bin))
summary(values(zik_bin))
summary(values(dcz_bin))

# save binary rasters for previously published maps - for public uses
den_past <- mapROC(den_pres, "DEN", "past_only")
zik_past <- mapROC(zik_pres, "ZIK", "past_only")
chik_past <- mapROC(chik_pres, "CHIK", "past_only")
yf_past <- mapROC(yf_pres, "YF", "past_only")

writeRaster(den_past, "data/intermediate_datasets/DEN_previous_binrast.tif", overwrite=T)
writeRaster(zik_past, "data/intermediate_datasets/ZIK_previous_binrast.tif", overwrite=T)
writeRaster(chik_past, "data/intermediate_datasets/CHIK_previous_binrast.tif", overwrite=T)
writeRaster(yf_past, "data/intermediate_datasets/YF_previous_binrast.tif", overwrite=T)


# convert lower and upper bound rasters into binary maps
disease_names <- c("DEN", "CHIK", "ZIK", "YF")
masked_lwr <- stack(); masked_upr <- stack(); dcz_lwr_bin <- stack(); dcz_upr_bin <- stack()

# Loop through each disease
for (disease in disease_names) {
  # Read in the raster files
  lwr <- raster(paste0("outputs/Rasters/", disease, "_riskmap_lwr_masked.tif"))
  upr <- raster(paste0("outputs/Rasters/", disease, "_riskmap_upr_masked.tif"))
  
  if(disease == "YF") { mask <-  mask_yf } else { mask <- mask_ts }
  pres <- get(paste0(tolower(disease), "_pres"))
  cut_new <- mapROC(pres, disease, "new_only")$cut_new
  
  # Mask rasters
  masked_lwr <- stack(masked_lwr, lwr * mask)
  masked_upr <- stack(masked_upr, upr * mask)
  
  # Rename layers
  names(masked_lwr) <- c(names(masked_lwr)[-length(names(masked_lwr))], tolower(disease))
  names(masked_upr) <- c(names(masked_upr)[-length(names(masked_upr))], tolower(disease))
  
  if(disease == "YF") { 
   
    yf_lwr_bin <- masked_lwr[[tolower(disease)]] > cut_new ; 
    yf_upr_bin <- masked_upr[[tolower(disease)]] > cut_new 
    
    } else { 
    
  dcz_lwr_bin <- stack(dcz_lwr_bin, masked_lwr[[tolower(disease)]] > cut_new)
  dcz_upr_bin <- stack(dcz_upr_bin, masked_upr[[tolower(disease)]] > cut_new)
  
  }
}

# collapsing dengue, chikungunya, and zika maps
dcz_lwr_bin <- calc(dcz_lwr_bin, sum)
dcz_upr_bin <- calc(dcz_upr_bin, sum)
dcz_lwr_bin[dcz_lwr_bin > 0] <- 1
dcz_upr_bin[dcz_upr_bin > 0] <- 1


# save rasters
writeRaster(dcz_lwr_bin, "outputs/Rasters/DCZ_binmap_lwr.tif", overwrite=T)
writeRaster(dcz_upr_bin, "outputs/Rasters/DCZ_binmap_upr.tif", overwrite=T)

writeRaster(yf_lwr_bin, "outputs/Rasters/YF_binmap_lwr.tif", overwrite=T)
writeRaster(yf_upr_bin, "outputs/Rasters/YF_binmap_upr.tif", overwrite=T)