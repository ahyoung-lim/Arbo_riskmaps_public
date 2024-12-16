# ==========================================================
# 
# Surveillance model fitting and cross-validation
# 
# ==========================================================

path = "C:/Users/AhyoungLim/Dropbox/WORK/WHO_GAI/ENMs/Multi_arbo_mapping/"

source("script/00_setup.R")
source("script/01_get_covariates.R")

# load in covariate rasters 
surv_cov_list <- loadRasters("surv")[[1]]

# load in intermediate dataset
dat <- readRDS("data/intermediate_datasets/Surv_model_fit_data.rds")

# names of covariates
surv_cov_names <- c("GDP", "GDP_national", "Urban", "Acc_walk", "Acc_city", "Trmt", "U5M", "GovEff", "Physician")


# plain vanilla RF model
# mod_surv <- RF_mod(data = dat, 
#                    cov_names = surv_cov_names) # returns an RF model object


# make prediction data frame
# a standardised set of covariates (5 km x 5 km global scale)
pred.data <- as.data.frame(surv_cov_list, xy=T) 
rm("surv_cov_list"); gc()

# create an NA index
NAindex = which(!complete.cases(pred.data[, surv_cov_names]))

# reduced size pred.data
pred.data.small <- pred.data[-NAindex, ]



# Spatial block N-fold CV and prediction ==================
set.seed(123) 
nfolds = 100 

# create 500 km x 500 km blocks
# create blocks and assign data to different blocks/folds (returns dataframe, cv_spatial object, and block geometry)
vList <- create_blocks(dat, nfolds = nfolds) 
vdat <- vList$dat 

# prepare datasets for training/test
set.seed(123) 
vCV_splits <- CV_data(vList, nfolds = nfolds) # returns a trianSet and testSet

# register the parallel backend
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 

tic() # get start time
v_CV <-  foreach(CV_train = vCV_splits$trainSet, 
                 CV_test = vCV_splits$testSet,
                 .packages=c("randomForest", "Metrics", "data.table","dplyr")
                 ) %dopar%{ 
                   fit_mod_surv(CV_train, CV_test, 
                           cov_names = surv_cov_names, 
                           pred_data_parallel = pred.data.small,
                           pred = TRUE, # for predictions
                           oob = TRUE # for calculating spatially stratified AUC
                           )
                                }; toc() 
stopCluster(cl)

# Processing the CV and prediction results
# Spatially stratified AUC
v_OOB_out <- OOBout(v_CV)  # get the out-of-bag samples
saveRDS(v_OOB_out, file = "outputs/cross_validation/Surv_OOB_data.rds")

# AUC_Surv <- AUC_strata(v_OOB_out, "Surv") # returns a table and map of spatially stratified AUCs


# Weighted average of predictions
cl <- makeCluster(ncpus-2) 
registerDoParallel(cl)

v_pred_out <- Predout(v_CV, pred=TRUE) # returns a dataframe with N predictions and their weighted average and IQR
stopCluster(cl)

# saveRDS(v_pred_out$pred, file = paste0("outputs/cross_validation/Surv_", nfolds, "fold_pred.rds"))
# saveRDS(v_pred_out$AUC, file = paste0("outputs/cross_validation/Surv_", nfolds, "fold_AUC.rds"))

# save figures and rasters in a local directory and return weighted mean raster only
surv_rast <- saveRasPlot(v_pred_out$pred, "Surv", color_opt = "viridis", color_direction = 1)


# calculate the variable importance and partial dependence values
cl <- makeCluster(ncpus) 
registerDoParallel(cl) 
tic()
vi_Surv <- foreach(df = vCV_splits$trainSet,  
                    .packages=c('randomForest', 'dplyr', 'tidyr', 'data.table', 'pdp')) %dopar% {
                      vImp(df, 
                           term_names = surv_cov_names, 
                           feature_names = surv_cov_names,
                           pdp = TRUE)
                    }; toc()
stopCluster(cl)
viPlot_surv <- viPlot(vi_Surv)

saveRDS(vi_Surv, "outputs/cross_validation/Surv_VI.rds")

