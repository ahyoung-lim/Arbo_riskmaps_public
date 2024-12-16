# Load occurrence data for surveillance model ==============
source("script/00_setup.R")

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

# standardise disease names
viral$disease[viral$disease == "Avian influenza virus subtype H5N1"] = "Avian Influenza H5N1"
viral$disease[grepl("Herpes simplex", viral$disease)] = "Herpes"
viral$disease[grepl("Herpes zoster", viral$disease)] = "Shingles"

viral$disease[grepl("Marburg", viral$disease)] = "Marburg fever"
viral$disease[grepl("mononucleosis", viral$disease)] = "Mononucleosis"
viral$disease[grepl("Monkey", viral$disease)] = "Mpox"
viral$disease[grepl("Parvovirus", viral$disease)] = "Parvovirus"
viral$disease[grepl("Polio", viral$disease)] = "Polio"
viral$disease[grepl("Respiratory syncytial virus", viral$disease)] = "Respiratory syncytial virus"
viral$disease[grepl("Roseola", viral$disease)] = "Roseola"
viral$disease[grepl("Rota", viral$disease)] = "Rotavirus"
viral$disease[grepl("Meningitis", viral$disease)] = "Viral meningitis"
viral$disease[grepl("Machupo", viral$disease)] = "Machupo virus"
viral$disease[grepl("Varicella", viral$disease)] = "Chicken pox"
viral$disease[grepl("California", viral$disease)] = "La crosse encephalitis"


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
viral_all <- rbind(viral[, c("Longitude",  "Latitude", "Admin", "disease")], 
                   viral2[, c("Longitude",  "Latitude", "Admin", "disease")]
)

viral_all$disease <- toupper(viral_all$disease)

rm(viral, viral2)


# Data processing ==========================================

viral_all <- viral_all %>%
  filter(!disease %in% c(
                        "EBOLA", "HIV/AIDS", "POLIO" , "MEASLES"
  ))
                         # "MUMPS", "HEPATITIS A", "ROTAVIRUS"
                        # "VIRAL MENINGITIS", no points available in Africa
  #                        ))%>%
  # select(-disease)
nrow(viral_all) # 338005
viral_occ <- thin_occ_viral(viral_all)
nrow(viral_occ) #[1] 21700

# Supplementary table
viral_occ %>%
  group_by(disease)%>%
  tally()

# random background point generation
# set.seed(1030)
# 
# viral_bg <- thin_bg_viral(viral_occ)
# nrow(viral_bg) #[1] 21700
# viral_bg$Admin = -999
# 
# 
# dat <- rbind(data.frame(viral_occ[, c("Longitude", "Latitude", "Admin")],
#                         disease = "occurrence"),
#              viral_bg[, c("Longitude", "Latitude", "Admin", "disease")])
# 
# 
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
# source("script/01_get_covariates.R")
# surv_cov_list <- loadRasters("surv")[[1]]
# 
# # extract covariate values (points and polygons)
# c_vals <- adminExtract(dat, surv_cov_list, admin, fun = "mean")
# 
# dat = data.frame(dat, c_vals)
# 
# # create a numeric response variable
# dat$PA <- as.factor(as.numeric(!(dat$disease == "background")))
# 
# source("functions/fixNAs.R")
# # fix for as many NAs as possible, then drop records that are still NA
# dat <- fixNAs(dat, surv_cov_list)
# summary(dat)
# dat = dat[!is.na(rowSums(dat[, names(surv_cov_list)])), ]

# saveRDS(dat, "data/intermediate_datasets/Surv_model_fit_data.rds")

