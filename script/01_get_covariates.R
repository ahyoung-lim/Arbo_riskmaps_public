# Load in covariate rasters ================================
# in separate private repository: 
# https://github.com/ahyoung-lim/Multi_arbo_mapping
path = "C:/Users/AhyoungLim/Dropbox/WORK/Dengue Project 2022/ENMs/Multi_arbo_mapping/"

# function for scaling and centring covariate raster data
# includes (default) option for first zero inflating and log transforming
raster.scale <- function(ras, logT = TRUE){
  vec = as.vector(ras)
  if(logT){
    zinf <- min(vec[vec > 0], na.rm = T)
    vec = log(vec + 0.5 * zinf)
  }
  vec = scale(vec)
  values(ras) = vec
  return(ras)
}

loadRasters <- function(type) { 
  
  # Covariates for surveillance capability model
  surv_cov_list <- stack(
    raster.scale(raster(paste0(path, "Covariates/final/GDP_2009_2019_1km_masked_.tif"))), 
    raster.scale(raster(paste0(path, "Covariates/final/GDP_2009_2019_National_masked_.tif"))),
    raster(paste0(path, "Covariates/final/urban_2010_2020_1km_masked_fillNAs.tif")),
    raster.scale(raster(paste0(path, "Covariates/final/2020_walking_only_travel_time_to_healthcare_masked_.tif"))),
    raster.scale(raster(paste0(path, "Covariates/final/2015_accessibility_to_cities_v1.0_masked_.tif")), logT=FALSE),
    raster.scale(raster(paste0(path, "Covariates/final/treatment_seeking_2023_masked_.tif")), logT = FALSE),
    raster.scale(raster(paste0(path, "Covariates/final/Child_Mortality_masked_.tif"))),
    raster(paste0(path, "Covariates/final/Government_Effectiveness_masked_.tif")),
    raster.scale(raster(paste0(path, "Covariates/final/Physicians_Density_masked_.tif")))
  )
  
  names(surv_cov_list) = c("GDP", "GDP_national", "Urban", "Acc_walk", "Acc_city", "Trmt", "U5M", "GovEff", "Physician")
  
  # Covariates for arbovirus risk model
  
  arbo_cov_list <- stack(

    raster.scale(raster(paste0(path, "Covariates/final/Dengue_temperature_suitaiblity_masked_.tif"))),
    raster.scale(raster(paste0(path, "Covariates/final/Tmean_coldest_TerraClim_2010_2020_005dg_masked_.tif")), logT=FALSE),
    raster.scale(raster(paste0(path, "Covariates/final/PRCP_TerraClim_2010_2020_005dg_masked_.tif"))),
    raster.scale(raster(paste0(path, "Covariates/final/NDVI_2010_2020_005dg_masked_.tif")), logT = FALSE),
    raster(paste0(path, "Covariates/final/DHI_global_clusters_14c_sv20_masked_.tif")),
    raster(paste0(path, "Covariates/final/Albopictus_mean_2020_rcp60_spreadXsuit_masked_.tif")),
    raster(paste0(path, "Covariates/final/Aegypti_mean_2020_rcp60_spreadXsuit_masked_.tif")),
    raster.scale(raster(paste0(path, "Covariates/final/GDP_2009_2019_1km_masked_.tif"))),
    raster.scale(raster(paste0(path, "Covariates/final/GDP_2009_2019_National_masked_.tif"))),
    raster(paste0(path, "Covariates/final/urban_2010_2020_1km_masked_fillNAs.tif")),
    raster.scale(raster(paste0(path, "Covariates/final/landscan_global_2022_masked_.tif"))),
    raster.scale(raster(paste0(path, "Covariates/final/YF_combined_monkeys_updated_.tif")), logT=FALSE),
    raster.scale(raster(paste0(path, "Covariates/final/YF_vaccine_adm1_2020_imperial_masked_.tif")), logT = FALSE) ,
    (100-raster(paste0(path, "Covariates/final/YF_vaccine_adm1_2020_imperial_masked_.tif")))*0.01, # (1-vaccine)
    raster.scale(raster(paste0(path, "Covariates/final/YF_Haemagogus_masked_.tif")), logT=FALSE)
    # raster.scale(raster("Covariates/final/Tcur_DEN_Ae_masked_.tif")),
    # raster.scale(raster("Covariates/final/Tcur_ZIK_Ae_masked_.tif"))
    
  )
  
  names(arbo_cov_list) = c("Tsuit", "Tcold", "PRCP", "NDVI", "DHI",
                            "Albo", "Aegypti", 
                           "GDP",  "GDP_national", "Urban", "Pop",
                           "NHP", "Vaccine","Vaccine_offset", "Hg"
                           #"Tcur_DEN", "Tcur_ZIK"
                           )
  
  # always load the latest surveillance capability raster together
  surv_rast_name <- list.files(paste0(path, "Maps_and_plots/Rasters"), pattern = "Surveillance_map_100pred_wmean.*.tif", full.names = FALSE)[which.max(file.info(list.files(paste0(path, "Maps_and_plots/Rasters"), pattern = "Surveillance_map_100pred_wmean.*.tif", full.names = TRUE))$mtime)]
  
  bootsRas_Surv <- raster(paste0(path, "Maps_and_plots/Rasters/", surv_rast_name))
  
  arbo_cov_list <- setNames(addLayer(arbo_cov_list, bootsRas_Surv), c(names(arbo_cov_list), "Surv"))
  
  # print messages
  if (type == "all") { 
    print(paste0("loading the latest surveillance raster: ", surv_rast_name))
    List <- list(surv_cov_list, arbo_cov_list)
  } else if (type == "surv") { 
      List <- list(surv_cov_list)
  } else { 
      print(paste0("loading the latest surveillance raster: ", surv_rast_name))
      List <- list(arbo_cov_list)    
    } 
  
  
  return(List)
  
  }
