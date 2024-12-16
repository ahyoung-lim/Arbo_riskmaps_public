# ==========================================================
# 
# Arbovirus model fitting and cross-validation
# 
# ==========================================================
path = "C:/Users/user/Dropbox/WORK/WHO_GAI/ENMs/Multi_arbo_mapping/"

source("script/00_setup.R")
source("script/01_get_covariates.R")
# load in covariate rasters 
arbo_cov_list <- loadRasters("arbo")[[1]]

# load in intermediate dataset
# this includes predictions from the surveillance capability model
arbo_dat <- readRDS("data/intermediate_datasets/Arbo_model_fit_data.rds")
yf_dat <- readRDS("data/intermediate_datasets/YF_model_fit_data.rds")



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
pred.data.small$DHI <- as.factor(round(pred.data.small$DHI, 0))
pred.data.small_yf <- subset(pred.data.small, select = -c(Albo, GDP)) # for yellow fever model
pred.data.small_arbo <- subset(pred.data.small, select = -c(NHP, Vaccine, Vaccine_offset, Hg, GDP)) # for dcz

rm("pred.data.small"); gc()

# Spatial block CV and prediction (arbo) ===================

spatial_CV <- function(disease_name, oob, pred) { 
  rm(CV)
  nfolds = 100
  
  if (disease_name == "yf") {
    data = get(paste0(disease_name, "_dat"))
    pred.dat = get(paste0("pred.data.small_yf"))
    term_names = get("yf_term_names")
    file_head = "YF"
  } else { 
    data = get("arbo_dat")
    pred.dat = get(paste0("pred.data.small_arbo"))
    term_names = get("arbo_term_names")
    file_head = "Arbo"
  }
  
  
  if (oob & disease_name %in% c("dengue", "yf")) { 
    oob_d = TRUE
    } else { oob_d = FALSE}

  # create 500 km x 500 km blocks
  # create blocks and assign data to different blocks/folds (returns dataframe, cv_spatial object, and block geometry)
  set.seed(13456)
  bList <- create_blocks(data, nfolds = nfolds) 

  # prepare datasets for training/test
  set.seed(1259)
  CV_splits <- CV_data(bList, nfolds = nfolds)

  # CV and prediction for dengue, zika, and chikungunya

  tic(); Sys.time() # get start time
  pred.dat$disease <- factor(disease_name, levels = levels(data$disease))

  cl <- makeCluster(ncpus) 
  registerDoParallel(cl) 
  
  CV <- foreach(CV_train = CV_splits$trainSet, CV_test = CV_splits$testSet,
                .packages = c("randomForest", "Metrics", "data.table", "dplyr"), 
                .export = c("fit_mod", "RF_mod")) %dopar% {
                  fit_mod(CV_train, CV_test,
                          cov_names = term_names,
                          pred_data_parallel = pred.dat,
                          pred= pred, 
                          oob = oob_d)
                }; toc()
  stopCluster(cl); gc()

  if (oob_d){ 
    OOB_out <- OOBout(CV)
    
    saveRDS(OOB_out, file = paste0("outputs/cross_validation/", file_head, "_OOB_data.rds"))
    
   }


  if (pred) {

    cl <- makeCluster(ncpus-2)
    registerDoParallel(cl)
    # to aggregate results across different folds
    Pred_out <- Predout(CV, pred=pred)

    # save plots and rasters in a local directory and return weighted mean raster only
    dname = ifelse(disease_name =="dengue", "DEN",
                   ifelse(disease_name == "zika", "ZIK",
                          ifelse(disease_name == "chikungunya", "CHIK", "YF")))

    # save 100 predictions and weights
    saveRDS(Pred_out$pred, file = paste0("outputs/cross_validation/", dname, "_", nfolds, "fold_pred.rds"))
    saveRDS(Pred_out$AUC, file = paste0("outputs/cross_validation/", dname, "_", nfolds, "fold_AUC.rds"))

    saveRasPlot(Pred_out$pred, dname, color_opt = "rocket", color_direction = -1)

  }

  stopCluster(cl)
  gc()
  return(CV)
}

disease_names <- c("dengue", "chikungunya", "zika", "yf")
list <- lapply(disease_names, function(name) {spatial_CV(name, oob = F, pred = T)})
       





# calculate the variable importance and partial dependence values ======================
get_vi <- function(disease_name) { 
  data = get(paste0(disease_name, "_dat"))
  nfolds = 100
  
  if (disease_name == "yf") {
    data = get(paste0(disease_name, "_dat"))
    term_names = get("yf_term_names")
    file_head = "YF"
  } else { 
    data = get("arbo_dat")
    term_names = get("arbo_term_names")
    file_head = "Arbo"
  }
  
  set.seed(13456)
  bList <- create_blocks(data, nfolds = nfolds) 

  # prepare datasets for training/test
  set.seed(1259)
  CV_splits <- CV_data(bList, nfolds = nfolds)

    
  cl <- makeCluster(ncpus) 
  registerDoParallel(cl) 
  tic()
  vi <- foreach(df = CV_splits$trainSet,  
                     .packages=c('randomForest', 'dplyr', 'tidyr', 'data.table', 'pdp'), 
                     .export=c('vImp')) %dopar% {
                       vImp(df, 
                            term_names = term_names, 
                            feature_names =  term_names[!grepl("offset", term_names)], 
                            pdp = TRUE)
                     }; toc()
  
  
  # viPlot <- viPlot(vi)

  saveRDS(vi, paste0("outputs/cross_validation/", file_head, "_VI.rds"))
}
