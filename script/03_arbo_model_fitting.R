# ==========================================================
# 
# Arbovirus model fitting and cross-validation
# 
# ==========================================================

source("script/00_setup.R")
source("script/01_get_covariates.R")
# load in covariate rasters 
arbo_cov_list <- loadRasters("arbo")[[1]]

# load in intermediate dataset
# this includes predictions from the surveillance capability model
arbo_dat <- readRDS("data/intermediate_datasets/Arbo_model_fit_data.rds")



# Arbovirus model building =================================
# building random forest model 
arbo_term_names <- c("Tsuit", "Tcold", "PRCP", "NDVI", "DHI", "Albo", "Aegypti",
                    "GDP_national", "Pop", "Urban",
                    "disease", "offset(Surv)"
                    )

yf_term_names <- c("Tsuit", "Tcold", "PRCP", "NDVI", "DHI", "Aegypti", 
                  "GDP_national", "Pop", "Urban",
                  "NHP", "Hg", 
                  "disease", "offset(Surv*Vaccine_offset)"
                  )

# plain vanilla RF models
# mod_arbo <- RF_mod(data = arbo_dat, 
#                    cov_names = arbo_term_names, 
#                    importance = TRUE)
# 
# mod_yf <- RF_mod(data = yf_dat,
#                  cov_names = yf_term_names, 
#                  importance = TRUE)


# make prediction data frame
# a standardised set of covariates (5 km x 5 km global scale)
pred.data <- as.data.frame(arbo_cov_list, xy=T) 
rm("arbo_cov_list"); gc()

# create an NA index
all_cols = names(pred.data)
col_to_exclude = c("x", "y", "GDP", "disease", "Surv")
col_to_include <- all_cols[!(all_cols %in% col_to_exclude)]

NAindex = which(!complete.cases(pred.data[, col_to_include]))

# reduced size pred.data
pred.data.small <- pred.data[-NAindex, ]
pred.data.small_yf <- subset(pred.data.small, select = -c(Albo, GDP)) # for yellow fever model
pred.data.small_arbo <- subset(pred.data.small, select = -c(NHP, Vaccine, Vaccine_offset, Hg, GDP)) # for dcz


# Spatial block CV and prediction (arbo) ===================
nfolds = 100

# create 500 km x 500 km blocks
# create blocks and assign data to different blocks/folds (returns dataframe, cv_spatial object, and block geometry)
set.seed(13456)
arbo_bList <- create_blocks(arbo_dat, nfolds = nfolds) 

# prepare datasets for training/test
set.seed(1259)
arboCV_splits <- CV_data(arbo_bList, nfolds = nfolds)

# CV and prediction for dengue, zika, and chikungunya
pred_rast <- list() # to store prediction rasters
disease <- "dengue"

tic(); Sys.time() # get start time
pred.data.small_arbo$disease <- factor("dengue", levels = levels(arbo_dat$disease))

cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

pred = TRUE
oob = TRUE

CV <- foreach(CV_train = arboCV_splits$trainSet, CV_test = arboCV_splits$testSet,
              .packages = c("randomForest", "Metrics", "data.table", "dplyr")) %dopar% {
                fit_mod(CV_train, CV_test,
                        cov_names = arbo_term_names,
                        pred_data_parallel = pred.data.small_arbo,
                        disease_name = disease,
                        pred= pred, 
                        oob = oob)
              }; toc()
stopCluster(cl); gc()

if (oob){ 
  OOB_out <- OOBout(CV)
  
  saveRDS(OOB_out, file = "outputs/cross_validation/Arbo_OOB_data.rds")
  
 }


if (pred) { 
  
  cl <- makeCluster(ncpus-2) 
  registerDoParallel(cl) 
  # to aggregate results across different folds
  Pred_out <- Predout(CV, pred=pred) 
  
  # save plots and rasters in a local directory and return weighted mean raster only
  dname = ifelse(disease =="dengue", "DEN", 
                 ifelse(disease == "zika", "ZIK", 
                        ifelse(disease == "chikungunya", "CHIK", "YF")))
  
  # save 100 predictions and weights
  # saveRDS(Pred_out$pred, file = paste0("outputs/cross_validation/", dname, "_", nfolds, "fold_pred.rds"))
  # saveRDS(Pred_out$AUC, file = paste0("outputs/cross_validation/", dname, "_", nfolds, "fold_AUC.rds"))
  
  pred_rast[[disease]] <- saveRasPlot(Pred_out$pred, dname, color_opt = "rocket", color_direction = -1) 
}


stopCluster(cl)
gc()

# Spatial block CV and prediction (yellow fever) ===========
# create 500 km x 500 km blocks
# create blocks and assign data to different blocks/folds (returns dataframe, cv_spatial object, and block geometry)
set.seed(13456)
yf_bList <- create_blocks(yf_dat, nfolds = nfolds)

# prepare datasets for training/test
set.seed(1259)
yfCV_splits <- CV_data(yf_bList, nfolds = nfolds)

# CV and prediction for yellow fever
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 
pred.data.small_yf$disease = factor("yf", levels = levels(arbo_dat$disease))

yf_CV <-  foreach(CV_train = yfCV_splits$trainSet, CV_test = yfCV_splits$testSet,
                    .packages=c("randomForest", "Metrics")) %dopar%{ 
                    fit_mod(CV_train, CV_test, 
                            cov_names = yf_term_names, 
                            pred_data_parallel = pred.data.small_yf, 
                            disease_name = "yf",
                            pred=TRUE,
                            oob=TRUE
                    )}
stopCluster(cl)

cl <- makeCluster(ncpus-2) 
registerDoParallel(cl) 

# to aggregate results across different folds
yf_pred_out <-  Predout(yf_CV, pred=TRUE) 

# save 100 predictions and weights
# saveRDS(yf_pred_out$pred, file = paste0("outputs/cross_validation/YF_, nfolds, "fold_pred.rds"))
# saveRDS(yf_pred_out$AUC, file = paste0("outputs/cross_validation/YF_, nfolds, "fold_AUC.rds"))

# save plots and rasters in a local directory and return weighted mean raster only
saveRasPlot(yf_pred_out$pred, "YF", color_opt = "rocket", color_direction = -1)

# spatially stratified AUCs
yf_OOB_out <- OOBout(yf_CV)
saveRDS(yf_OOB_out, file = "outputs/cross_validation/YF_OOB_data.rds")



# calculate the variable importance and partial dependence values ======================
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 
tic()
vi_arbo <- foreach(df = arboCV_splits$trainSet,  
                   .packages=c('randomForest', 'dplyr', 'tidyr', 'data.table', 'pdp')) %dopar% {
                     vImp(df, 
                          term_names = arbo_term_names, 
                          feature_names =  arbo_term_names[!grepl("offset", arbo_term_names)], 
                          pdp = TRUE)
                   }; toc()
vi_yf <- foreach(df = yfCV_splits$trainSet,  
                   .packages=c('randomForest', 'dplyr', 'tidyr', 'data.table', 'pdp')) %dopar% {
                     vImp(df, 
                          term_names = yf_term_names, 
                          feature_names = yf_term_names[!grepl("offset", yf_term_names)], 
                                pdp = TRUE)
                   }

viPlot_arbo <- viPlot(vi_arbo)
viPlot_yf <- viPlot(vi_yf)



