# Load occurrence data for arbovirus model ================

# Load in data in separate private repository: 
# https://github.com/ahyoung-lim/Multi_arbo_mapping
path = "C:/Users/AhyoungLim/Dropbox/WORK/Dengue Project 2022/ENMs/Multi_arbo_mapping/"

# Occurrence (for different arboviral diseases) ===========
# Dengue --------------------------------------------------

# published datasets
# NB dengue occurrence dataset is make up of both points and polygons in separate files
den_pt <- read.csv(paste0(path, "published_datasets/DENV_Messina/Occurrence_points_standard_checkedEC.csv"))
den_py <- read.csv(paste0(path, "published_datasets/DENV_Messina/Occurrence_poly_standard_checkedEC.csv"))
den1 <- rbind(den_pt, den_py)

# Updated dengue occ data 
den_who <- read.csv(paste0(path, "New_occurrence_datasets/DENV_oc_who_update.csv"))
den_hm <- read.csv(paste0(path, "New_occurrence_datasets/DENV_oc_HM.csv")) %>% rename(Year = year)
den_supp <- read.csv(paste0(path, "New_occurrence_datasets/DENV_oc_supp.csv")) # from ECDC, WHO DONS

den_new <- rbind(
             den_who[, c("Longitude",  "Latitude", "Admin", "Year")], 
             den_hm[, c("Longitude",  "Latitude", "Admin", "Year")], 
             den_supp[, c("Longitude",  "Latitude", "Admin", "Year")]
             )
den <- rbind(den1[, c("Longitude",  "Latitude", "Admin", "Year")], 
             den_new[, c("Longitude",  "Latitude", "Admin", "Year")]
             )


# Chikungunya-----------------------------------------------
# published datasets
chik1 <- read.csv(paste0(path,"published_datasets/CHIKV_Nsoesie/CHIKV_oc_Kraemer_V2.csv"))%>%
  filter(!occurrence.id %in% c('PM_391', 'PM_392', 'PM_458'))%>% #wrong points removed (USA, China)
  rename(Admin = precision)%>%
  mutate(Admin = gsub("ADMIN", "", Admin))%>%
  mutate(Admin = ifelse(Admin %in% c("COUNTRY"), "0", 
                        ifelse(Admin %in% c("PRECISE"), "-999", Admin)))

# new occurrence datasets
chik_dmmg <- read.csv(paste0(path, "New_occurrence_datasets/CHIKV_oc_dmmg.csv"))
chik_who <- read.csv(paste0(path,"New_occurrence_datasets/CHIKV_oc_who_update.csv"))
chik_bra <- read.csv(paste0(path,"New_occurrence_datasets/CHIKV_oc_BRA_V2.csv"))
chik_hm <- read.csv(paste0(path,"New_occurrence_datasets/CHIKV_oc_HM.csv")) %>% rename(Year = year)
chik_supp <- read.csv(paste0(path, "New_occurrence_datasets/CHIKV_oc_supp.csv"))
chik_new <- rbind(
              chik_dmmg[, c("Longitude",  "Latitude", "Admin", "Year")],
              chik_who[, c("Longitude",  "Latitude", "Admin", "Year")], 
              chik_bra[, c("Longitude",  "Latitude", "Admin", "Year")], 
              chik_hm[, c("Longitude",  "Latitude", "Admin", "Year")], 
              chik_supp[, c("Longitude",  "Latitude", "Admin", "Year")])

chik <- rbind(chik1[, c("Longitude",  "Latitude", "Admin", "Year")], 
              chik_new[, c("Longitude",  "Latitude", "Admin", "Year")])

# Zika ----------------------------------------------------
# published datasets
zik1 <- read.csv(paste0(path, "published_datasets/ZIKV_Messina/ZIKV_oc_Messina.csv"))

# new occurrence datasets
zik_dmmg <- read.csv(paste0(path,"New_occurrence_datasets/ZIKV_oc_dmmg.csv"))%>% 
  filter(!Country %in% c("ISRAEL", "CHINA"))
zik_hm <- read.csv(paste0(path, "New_occurrence_datasets/ZIKV_oc_HM.csv")) %>% rename(Year = year)

zik_new <- rbind(
             zik_dmmg[, c("Longitude",  "Latitude", "Admin", "Year")], 
             zik_hm[, c("Longitude",  "Latitude", "Admin", "Year")] )

zik <- rbind(zik1[, c("Longitude",  "Latitude", "Admin", "Year")], 
             zik_new[, c("Longitude",  "Latitude", "Admin", "Year")])

# Yellow fever----------------------------------------------
# published datasets
yf1 <- read.csv(paste0(path, "published_datasets/YFV_Shearer/YFV_occurrence_clean_updated.csv"))
yf1 <- yf1 %>%
  rename(Year = year, 
         Longitude = longitude, 
         Latitude = latitude,
         Admin = admin_level)

# new occurrence datasets
yf_who <- read.csv(paste0(path, "New_occurrence_datasets/YFV_oc_who_update.csv"))
yf_hm <- read.csv(paste0(path, "New_occurrence_datasets/YFV_oc_HM.csv")) %>% rename(Year = year) %>%
  filter(!country %in% c("Haiti", "Jamaica", "Mexico", "Belize"))
yf_supp <- read.csv(paste0(path, "New_occurrence_datasets/YFV_oc_supp.csv"))

yf_new <- rbind(
            yf_who[, c("Longitude",  "Latitude", "Admin", "Year")],
            yf_hm[, c("Longitude",  "Latitude", "Admin", "Year")],

            yf_supp[, c("Longitude",  "Latitude", "Admin", "Year")]
            )

yf <- rbind(yf1[, c("Longitude",  "Latitude", "Admin", "Year")], 
            yf_new[, c("Longitude",  "Latitude", "Admin", "Year")]
            )


# admin 3 and 4 are considered points
den$Admin[den$Admin > 2] = -999
zik$Admin[zik$Admin > 2] = -999
chik$Admin[chik$Admin > 2] = -999

yf$Admin[is.na(yf$Admin)] = -999 # NA or higher than admin2 = point
yf$Admin[yf$Admin > 2] = -999


den$disease <- "dengue"
yf$disease <- "yf"
chik$disease <- "chik"
zik$disease <- "zik"

rm(chik_bra, chik_dmmg, chik_who, chik_hm, chik_supp, den_supp, den_hm, den_pt, den_py, den_who, yf_hm,yf_who, yf_supp, zik_dmmg, zik_hm)


# Arbovirus data preparation ======================================
# thinning separately on each arboviral disease
den <- thin_occ_arbo(den); nrow(den) # 5867
chik <- thin_occ_arbo(chik); nrow(chik) # 4727
zik <- thin_occ_arbo(zik); nrow(zik) # 1138
yf <- thin_occ_arbo(yf); nrow(yf) # 1395


# combine into one arbovirus dataset
arbo_occ<- rbind(den[, c("Longitude", "Latitude", "Admin")],
              chik[, c("Longitude", "Latitude", "Admin")],
              zik[, c("Longitude", "Latitude", "Admin")],
              yf[, c("Longitude", "Latitude", "Admin")])
arbo_occ$disease = as.factor(c(rep("dengue", nrow(den)),
                            rep("chikungunya", nrow(chik)),
                            rep("zika", nrow(zik)),
                            rep("yf", nrow(yf))))
nrow(arbo_occ) 
# write.csv(arbo_occ, "data/intermediate_datasets/arbo_occ_thinned.csv", row.names=F)

# random background point generation
# set.seed(1221)
# arbo_bg <- thin_bg_arbo(arbo_occ) 
# nrow(arbo_bg) 
# 
# # combining presence and background
# arbo_dat <- rbind(arbo_occ,
#                   data.frame(arbo_bg[, c("Longitude", "Latitude", "disease")],
#                              Admin = -999))
# 
# 
# # add admin relevant admin unit GAUL codes
# points_locs <- data.frame(a0 = raster::extract(admin$ad0, arbo_dat[, c("Longitude", "Latitude")]),
#                           a1 = raster::extract(admin$ad1, arbo_dat[, c("Longitude", "Latitude")]),
#                           a2 = raster::extract(admin$ad2, arbo_dat[, c("Longitude", "Latitude")]))
# points_locs = points_locs[(arbo_dat$Admin != -999), ]
# arbo_dat$GAUL = NA
# arbo_dat$Admin <- as.numeric(arbo_dat$Admin)
# v = cbind(1:nrow(points_locs),
#           (arbo_dat$Admin[(arbo_dat$Admin != -999)] + 1))
# arbo_dat$GAUL[(arbo_dat$Admin != -999)] = points_locs[v]
# 
# # extract covariate values (points and polygons)
# c_vals <- adminExtract(arbo_dat, arbo_cov_list, admin, fun = "mean")
# arbo_dat = data.frame(arbo_dat, c_vals)
# 
# # DHI should be treated as a factor
# arbo_dat$DHI = as.factor(round(arbo_dat$DHI, 0))
# 
# # add presence absence column
# arbo_dat$PA <- as.factor(c(rep(1, nrow(arbo_occ)), rep(0, nrow(arbo_bg))))
# 
# yf_dat = subset(arbo_dat, select = -c(GDP))
# arbo_dat <- subset(yf_dat, select = -c(NHP, Vaccine, Vaccine_offset, Hg))
# 
# arbo_cov <- names(arbo_dat %>% select(Tsuit:Surv))
# yf_cov <- names(yf_dat %>% select(Tsuit:Surv))
# 
# source("script/01_get_covariates.R")
# arbo_cov_list <- loadRasters("arbo")[[1]]
#
# # fix NAs
# arbo_dat <- fixNAs(arbo_dat, subset(arbo_cov_list, arbo_cov))
# yf_dat <- fixNAs(yf_dat, subset(arbo_cov_list, yf_cov))
# 
# summary(arbo_dat)
# summary(yf_dat)
# 
# # trim to datapoints will no NAs
# arbo_dat = arbo_dat[!apply(arbo_dat[, c(arbo_cov)], 1, function(x) any(is.na(x))), ]
# yf_dat = yf_dat[!apply(yf_dat[, c(yf_cov)], 1, function(x) any(is.na(x))), ]
# 
# gc()
# 
# saveRDS(arbo_dat, file = paste0("data/intermediate_datasets/Arbo_model_fit_dat.rds"))
# saveRDS(yf_dat, file = paste0("data/intermediate_datasets/YF_model_fit_dat.rds"))


