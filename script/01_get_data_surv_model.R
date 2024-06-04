# Load occurrence data for surveillance model ==============

# Load in data in separate private repository: 
# https://github.com/ahyoung-lim/Multi_arbo_mapping

path = "C:/Users/AhyoungLim/Dropbox/WORK/WHO_GAI/ENMs/Multi_arbo_mapping/"

# published datasets
viral <- read.csv(paste0(path, "published_datasets/Shearer_background/pathogen_class_5_viral.csv"))

# assume all background data are points (not polygons)
# reformatting
viral <- data.frame(Longitude = viral$longitude,
                    Latitude = viral$latitude,
                    Admin = -999,
                    disease = viral$disease_name)

# Updated HealthMap data 2006-2019
viral2 <- read.csv(paste0(path, "New_occurrence_datasets/HealthMap_2006_2019_cleaned.csv"))%>%
  # remove some errors
  filter(!(disease == "YELLOW FEVER" & 
           country %in% c("Haiti", "Jamaica", "Mexico", "Belize", "Armenia",
                          "Republic of Georgia", "Saudi Arabia", "Jordan", "United Arab Emirates",
                          "Sri Lanka", "Maldives", "Australia", "New Zealand",
                          "Singapore", "Brunei", "Bhutan", "Vietnam")))
  
viral2 <- data.frame(Longitude = viral2$Longitude,
                    Latitude = viral2$Latitude,
                    Admin = viral2$Admin,
                    disease = viral2$disease)

# recombine an all viral disease dataset
viral_all <- rbind(viral[, c("Longitude",  "Latitude", "Admin")], 
                   viral2[, c("Longitude",  "Latitude", "Admin")]
)


# Data processing ==========================================

# thinning of all viral diseases occurrence points
nrow(viral_all) #[1] 455370

viral_occ <- thin_occ_viral(viral_all)
nrow(viral_occ) #[1] 23632


# random background point generation
set.seed(1030)

viral_bg <- thin_bg_viral(viral_occ)
nrow(viral_bg) #[1] 23632
viral_bg$Admin = -999


# dat <- rbind(data.frame(viral_occ[, c("Longitude", "Latitude", "Admin")],
#                         disease = "occurrence"),
#              viral_bg[, c("Longitude", "Latitude", "Admin", "disease")])


# # add admin relevant admin unit GAUL codes
# points_locs <- data.frame(a0 = raster::extract(admin$ad0, dat[, c("Longitude", "Latitude")]),
#                           a1 = raster::extract(admin$ad1, dat[, c("Longitude", "Latitude")]),
#                           a2 = raster::extract(admin$ad2, dat[, c("Longitude", "Latitude")]))
# points_locs = points_locs[(dat$Admin != -999), ]
# dat$GAUL = NA
# dat$Admin <- as.numeric(dat$Admin)
# 
# v = cbind(1:nrow(points_locs),
#                  (dat$Admin[(dat$Admin != -999)] + 1))
# dat$GAUL[(dat$Admin != -999)] = points_locs[v]
# 
# # extract covariate values (points and polygons)
# c_vals <- adminExtract(dat, surv_cov_list, admin, fun = "mean")
# 
# dat = data.frame(dat, c_vals)
# 
# # create a numeric response variable
# dat$PA <- as.factor(as.numeric(!(dat$disease == "background")))


# source("script/01_get_covariates.R")
# surv_cov_list <- loadRasters("surv")[[1]]


# fix for as many NAs as possible, then drop records that are still NA 
# dat <- fixNAs(dat, surv_cov_list)
# summary(dat)
# dat = dat[!is.na(rowSums(dat[, names(surv_cov_list)])), ]

# saveRDS(dat, "data/indermediate_datasets/Surv_model_fit_data.rds")

