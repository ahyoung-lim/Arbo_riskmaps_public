
# Function to calculate quick AUC
auc.f <- function(pred, obs){

  # Calculate the AUC
  
  ROC_perf <- performance(prediction(pred,obs),"tpr","fpr")
  ROC_sens <- performance(prediction(pred,obs),"sens","spec")
  ROC_auc <- performance(prediction(pred,obs),"auc")
  AUC <- ROC_auc@y.values[[1]] # AUC
  
  # Mean sensitivity across all thresholds
  x.Sens <- mean(as.data.frame(ROC_sens@y.values)[,1], na.rm=T)
  # Mean specificity across all thresholds
  x.Spec <- mean(as.data.frame(ROC_sens@x.values)[,1], na.rm=T)
  
  # Create output table
  cbind(AUC, x.Sens, x.Spec)
}


AUC_map_buffer <- function(P_results, buffer, stratified=TRUE) { 
   
   points_sf <- P_results %>%
    st_as_sf(coords= c("Longitude", "Latitude"), crs=4326) 
    
    points_sf <- points_sf %>%
    mutate(within_buffer = st_is_within_distance(x = .,
                                              y = st_as_sf(points_sf),
                                              dist = buffer*1000), # m to km
           within_buffer_n = lengths(within_buffer))
  
  AUC_all <- list()
  for (i in 1:nrow(P_results)) { 
    ind <- points_sf$within_buffer[i][[1]] # get the index for all points that fall within 250 km buffer
    x <- P_results[ind, ] # subset the dataframe
    x$AUC <- Metrics::auc(x$PA, x$pred) # calculate the buffer-wise AUC
    
    if(!is.na(mean(x$AUC))){
    x$Sens <- auc.f(x$pred, x$PA)[1,2]
    x$Spec <- auc.f(x$pred, x$PA)[1,3]
    } else { 
      x$Sens <- NA
      x$Spec <- NA
    }
    
    AUC_all[[i]] <- x # save AUC values into a list
    cat(paste("\rcompleted", i))
  }
  
  AUC_dat <- rbindlist(AUC_all) %>%
    
    group_by(pred, PA, Longitude, Latitude)%>%
    summarise(AUC_mean = mean(AUC, na.rm=T), 
              Sens_mean = mean(Sens, na.rm=T), 
              Spec_mean = mean(Spec, na.rm=T))%>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs=4326)


  world <- read_sf("data/admin_rasters/world_no_lakes.shp")
  
  AUC_map <- 
    
    ggplot()+
      geom_sf(data=AUC_dat[!is.na(AUC_dat$AUC_mean),], aes(color=AUC_mean))+
      geom_sf(data=world, fill="transparent", color="grey40")+
      scale_color_distiller(
        name = "mean AUC",
        palette = "Spectral",
        values = c(0, 0.5,  0.7, 0.8, 1),
        breaks = c(0, 0.5,  0.7, 0.8, 0.9, 1),
        limits = c(0, 1),
        na.value = "grey80",
        direction = 1
  
      )+
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(),
            legend.key.height = unit(10, "mm"),
            plot.background = element_rect(fill="white"), 
            panel.background = element_rect(fill="white"))+
      coord_sf(xlim = c(-180, 180), ylim = c(-59, 75), expand = FALSE)
  
  if (stratified){ 
    cat("calculating spatially stratified AUCs")
    tab <- AUC_region(AUC_dat)
    List <- list("dat" = AUC_dat, 
                 "map" = AUC_map, 
                 "tab" = tab)

  } else { 
    List <- list("dat" = AUC_dat, 
                 "map" = AUC_map)
    }

  return(List)
}


AUC_region <- function(AUC_dat) {
  # load standardised UN states namaes
  source("functions/RNE_country_names.R")

  # occurrence/bg point
  x <- AUC_dat %>%  st_as_sf(coords = c("Longitude", "Latitude"), crs=4326)
  x$continent <- NA
  x$endemic <- NA
  
  # get the index for the nearest feature
  sf_use_s2(FALSE)
  nearest_continent <- st_nearest_feature(x, world)
  sf_use_s2(TRUE)
  
  # Extract continent/region/country names
  for (i in 1:nrow(x)){ 
    index <- nearest_continent[i]
    x$continent[i] <- world$continent[as.numeric(index)]
    x$endemic[i] <- world$endemic[as.numeric(index)]
  }
  

  tab <- x %>%
    st_drop_geometry()%>%
    group_by(continent)%>%
    summarise(mAUC = mean(AUC_mean, na.rm=T), 
              mSens = mean(Sens_mean, na.rm=T), 
              mSpec = mean(Spec_mean, na.rm=T))
  
  names(tab) <- c("Region", "meanAUC", "meanSens", "meanSpec")
  tab <- tab[!tab$Region %in% c("Seven seas (open ocean)"), ]
  
  tab$Region <- factor(tab$Region, levels = c("North America", "South America", 
                                              "Europe", "Africa", 
                                              "Asia", "Oceania"))
  tab2 <- x %>%
    st_drop_geometry()%>%
    group_by(endemic)%>%
    summarise(mAUC = mean(AUC_mean, na.rm=T), 
              mSens = mean(Sens_mean, na.rm=T), 
              mSpec = mean(Spec_mean, na.rm=T))
  
  names(tab2) <- c("Region", "meanAUC", "meanSens", "meanSpec")
  tab_all <- rbind(tab, tab2)
  
return(tab_all)
}

AUC_strata <- function(OOB, disease_abb) { 
  
  disease_abb2 <- ifelse(disease_abb == "YF", tolower("YF"), disease_abb)
  # Overall AUC (across region)
  overall_auc <- data.frame(AUC = Metrics::auc(OOB$PA, OOB$pred),
                            Sens = auc.f(OOB$pred, OOB$PA)[1,2],
                            Spec = auc.f(OOB$pred, OOB$PA)[1,3])
  
  # regionally stratified AUCs
  AUC <- AUC_map_buffer(OOB, buffer = 250, stratified = TRUE)
  
  AUC$tab <- rbind(AUC$tab, 
                   data.frame(Region = "Overall",
                                        meanAUC = overall_auc$AUC,
                                         meanSens = overall_auc$Sens,
                                         meanSpec = overall_auc$Spec)) %>%
                    mutate(disease = disease_abb2)
  
  write.csv(AUC$tab, paste0("outputs/cross_validation/", disease_abb, "_AUC_table.csv"), row.names=F)

  # ggsave(filename = paste0("outputs/cross_validation/", disease_abb, "_AUC", gsub("-", "_", Sys.Date()), ".png"),
  # AUC$map$data, height=6, width=12, dpi=300)
  
  
  List <- list("tab" = AUC$tab, 
               "map" = AUC$map)
  return(List)  
} 

AUC_strata_arbo <- function(OOB) { 
  # Overall AUC (across region)
  overall_auc <- OOB %>% group_by(disease) %>% summarise(AUC = Metrics::auc(PA, pred),
                                                         Sens = auc.f(pred, PA)[1,2],
                                                         Spec = auc.f(pred, PA)[1,3])

  # regionally stratified AUCs
  den_AUC <- AUC_map_buffer(subset(OOB, OOB$disease == "dengue"),  buffer = 250, stratified = TRUE)
  chik_AUC <- AUC_map_buffer(subset(OOB, OOB$disease == "chikungunya"),  buffer = 250, stratified = TRUE)
  zik_AUC <- AUC_map_buffer(subset(OOB, OOB$disease == "zika"),  buffer = 250, stratified = TRUE)
  
  den_AUC$tab <- rbind(den_AUC$tab, data.frame(Region = "Overall",
                                               meanAUC = overall_auc$AUC[overall_auc$disease == "dengue"],
                                               meanSens = overall_auc$Sens[overall_auc$disease == "dengue"],
                                               meanSpec = overall_auc$Spec[overall_auc$disease == "dengue"])) %>%
    mutate(disease = "dengue")
  chik_AUC$tab <- rbind(chik_AUC$tab, data.frame(Region = "Overall",
                                                 meanAUC = overall_auc$AUC[overall_auc$disease == "chikungunya"],
                                                 meanSens = overall_auc$Sens[overall_auc$disease == "chikungunya"],
                                                 meanSpec = overall_auc$Spec[overall_auc$disease == "chikungunya"]))%>%
   mutate(disease = "chikungunya")
  zik_AUC$tab <- rbind(zik_AUC$tab, data.frame(Region = "Overall",
                                               meanAUC = overall_auc$AUC[overall_auc$disease == "zika"],
                                               meanSens = overall_auc$Sens[overall_auc$disease == "zika"],
                                               meanSpec = overall_auc$Spec[overall_auc$disease == "zika"]))%>%
    mutate(disease = "zika")
  all_tab <- rbindlist(list(den_AUC$tab, chik_AUC$tab, zik_AUC$tab))
  all_map <- list(den_AUC$map, chik_AUC$map, zik_AUC$map)
  
  write.csv(all_tab, paste0("outputs/cross_validation/Arbo_AUC_table.csv"), row.names=F)
  

  
  List <- list("tab" = all_tab, 
               "map" = all_map)
  return(List)  
} 


