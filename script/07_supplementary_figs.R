# ==========================================================
# 
# Producing figures in the Supplementary Information
# 
# ==========================================================

source("script/00_setup.R")

# SFig 1 ===================================================

compare_datasets <- function(dat1, dat2) { 
  dat1$Admin[dat1$Admin > 2] = -999
  dat2$Admin[dat2$Admin > 2] = -999
  
  dat1 <- dat1[, c("Longitude", "Latitude", "Admin", "Year")]
  
  v1 <- nrow(dat1) 
  dat_occ <- thin_occ_arbo(dat1) 
  v1_thinned <- nrow(dat_occ) 
  
  v2 <- nrow(dat2)
  dat_occ2 <- thin_occ_arbo(dat2)
  v2_thinned <- nrow(dat_occ2) 
  
  all_dat <- rbind(dat_occ %>%  mutate(id = "v1") , 
                   dat_occ2 %>%  mutate(id = "v2"))
  
  all_dat <- all_dat %>% arrange(id)
  all_dat <- all_dat[!duplicated(all_dat$uID), ]
  
  new_occ <- nrow(all_dat)-nrow(dat_occ) # no. of unique locations added
  
  print(paste0("v1: ", v1_thinned, ", v2: ", v2_thinned, ", no. of unique locations added: ", new_occ ))
  
  all_dat$Type = c("Polygon", "Point")[((all_dat$Admin == -999) + 1)]
  
  occ_pt <- geom_point(data=all_dat, 
                        aes(x=Longitude, y=Latitude, color = id, shape = Type),
                        size=1.5)
  
  p <- base_map(world) + occ_pt + 
    scale_color_manual(values = c("#6388b4", "#Ef6f6a"), name = "Datasets") + 
    scale_shape_manual(values=c(1, 2))
  
  return(p)
}

# load in occurrence datasets
source("script/01_get_data_arbo_model.R")

p1 <- compare_datasets(den1, den_new) + ggtitle("A  Dengue (+495)")
p2 <- compare_datasets(chik1, chik_new) + ggtitle("B  Chikungunya (+3994)")
p3 <- compare_datasets(zik1, zik_new) + ggtitle("C  Zika (+948)")
p4 <- compare_datasets(yf1, yf_new) + ggtitle("D  Yellow fever (+372)")

SI_fig1 <- ggarrange(p1, p2, p3, p4, 
          nrow=2, ncol = 2, 
          common.legend = T)

ggsave(filename = paste0("outputs/Figures/SI_Fig1.png"), SI_fig1, bg = "white", height=9, width=18, dpi=300)


# SFig 2 ===================================================
# load in occurrence datasets
source("script/01_get_data_surv_model.R")

# visualize and save
viral_for_map <- rbind(viral_bg[, c("Longitude", "Latitude", "Admin")],
                       viral_occ[, c("Longitude", "Latitude", "Admin")])
viral_for_map$Type = c("Polygon", "Point")[((viral_for_map$Admin == -999) + 1)]
viral_for_map$Type[1:nrow(viral_bg)] = "Background point"

viral_for_map$Type <- factor(viral_for_map$Type, levels = c("Point", "Background point"))
levels(viral_for_map$Type)[levels(viral_for_map$Type) == "Point"] <- "Occurrence"

occ_pt <- geom_point(data=viral_for_map[viral_for_map$Type == "Occurrence",], 
                     aes(x=Longitude, y=Latitude, color = Type), alpha=0.7,
                     shape = 1, size=1)

# global map
(SI_fig2 <- base_map(world) + occ_pt +
    scale_color_manual(values = c("#Ef6f6a","#6388b4"))+
    guides(color = guide_legend(override.aes = list(size = 4)))+
    theme(legend.position = "none")
)

ggsave(filename = paste0("outputs/Figures/SI_Fig2.png"), SI_fig2, height = 6, width = 12, dpi=300)


# SFig 3 ===================================================
# variable importance and partial dependence plots for surv model
# load in variable importance and partial dependence values
vi_Surv <- readRDS("outputs/cross_validation/Surv_VI.rds")

viPlot_surv <- viPlot(vi_Surv)

ggsave(filename = "outputs/Figures/SI_Fig3a.png", viPlot_surv[[1]], height=6, width=9, dpi=300)
ggsave(filename = "outputs/Figures/SI_Fig3b.png", viPlot_surv[[2]], height=6, width=10, dpi=300)


# SFig 4 ===================================================
# IQR maps for surv model
# load in IQR raster

Ras <- raster("outputs/Rasters/Surveillance_map_IQR.tif")

SI_fig4 <- plotRaster(Ras) + 
    theme(legend.title = element_text(family = "Arial", vjust = 1))+
    scale_fill_viridis(option = "viridis", direction = 1, na.value="white", 
                       name= "Uncertainty")

ggsave(filename = "outputs/Figures/SI_Fig4.png", SI_fig4, height=6, width=12, dpi=300)

# SFig5 ====================================================
# Overall and regionally-stratified model performance (surv model)
# load in source data
v_OOB_out <- readRDS("outputs/cross_validation/Surv_OOB_data.rds")
AUC_Surv <- AUC_strata(v_OOB_out, "Surv") # calculate spatially stratified AUCs
surv_tab <- AUC_Surv$tab

surv_tab$Region <- factor(surv_tab$Region, levels = c("Overall", "Endemic", "Non-endemic", 
                                                      "South America", "North America", "Africa", 
                                                      "Europe", "Asia", "Oceania"))

# visualise the AUC tables as heatmap
SI_fig5 <- surv_tab %>%
  gather(., metric, value, meanAUC:meanSpec, factor_key=TRUE)%>%
  filter(!Region %in% c("Endemic", "Non-endemic"))%>%
  ggplot()+
  geom_tile(aes(x = metric, y=Region, fill=value), color="white")+
  geom_text(aes(x = metric, y=Region,
                label = paste(sprintf("%.2f", round(value, 2)))), color = "black", size= 4)+
  scale_y_discrete(limits=rev)+
  scale_x_discrete(position = "top", labels = c("AUC", "Sensitivity", "Specificity"))+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1.5),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_blank(), 
        axis.line = element_blank(),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = 'white'), 
        panel.grid.major = element_blank(),
        text= element_text(family = "Arial")
  )+  
  ylab(NULL)+xlab(NULL)+
  scale_fill_distiller(
    name = "Metric",
    palette = "Spectral",
    values = c(0, 0.5,  0.7, 0.8, 1),
    breaks = c(0, 0.5,  0.7, 0.8, 0.9, 1),
    na.value = "grey80",
    direction = 1,
    limits = c(0, 1)
  ); surv_p


ggsave(filename = "outputs/Figures/SI_Fig5.png", SI_fig5, height=6, width=6, dpi=300, background = "white")

# SFig6 ====================================================
# Spatial map of AUC (surv model)
# load in source data
v_OOB_out <- readRDS("outputs/cross_validation/Surv_OOB_data.rds")
AUC_Surv <- AUC_strata(v_OOB_out, "Surv") # calculate spatially stratified AUCs

ggsave(filename = "outputs/Figures/SI_Fig6.png", AUC_Surv$map, height=6, width=12, dpi=300, bg = "white")


# SFig7 ====================================================
# variable importance and partial dependence plots for arbo model
# load in variable importance and partial dependence values
vi_arbo <- readRDS("outputs/cross_validation/Arbo_VI.rds")

viPlot_arbo <- viPlot(vi_arbo)

ggsave(filename = "outputs/Figures/SI_Fig7a.png", viPlot_arbo[[1]], height=6, width=9, dpi=300)
ggsave(filename = "outputs/Figures/SI_Fig7b.png", viPlot_arbo[[2]], height=6, width=10, dpi=300)


# SFig8 ====================================================
# Model predicted environmental suitability for dengue, chikungunya, and zika seperately
# load in rasters

DEN_range_mask <- raster("outputs/Rasters/DEN_riskmap_wmean_masked.tif")
CHIK_range_mask <- raster("outputs/Rasters/CHIK_riskmap_wmean_masked.tif")
ZIK_range_mask <- raster("outputs/Rasters/ZIK_riskmap_wmean_masked.tif")


# plotting
r_DEN <- plotRaster(DEN_range_mask)+
  scale_fill_viridis(option = "rocket", direction = -1, na.value = "transparent",
                     name= "Environmental\nsuitability",
                     limits = c(0,1))

r_CHIK <- plotRaster(CHIK_range_mask)+
  scale_fill_viridis(option = "rocket", direction = -1, na.value = "transparent",
                     name= "Environmental\nsuitability",
                     limits = c(0,1))


r_ZIK <- plotRaster(ZIK_range_mask)+
  scale_fill_viridis(option = "rocket", direction = -1, na.value = "transparent",
                     name= "Environmental\nsuitability",
                     limits = c(0,1))

SI_fig8 <- r_DEN + ggtitle("a  Dengue") +
  r_CHIK + ggtitle("b  Chikungunya") +
  r_ZIK + ggtitle("c  Zika") +
  plot_layout(ncol=1, guides = "collect")&
  theme(text = element_text(family = "Arial"))&
  scale_fill_viridis(option = "rocket", direction = -1, na.value = "transparent",
                     name= "Environmental\nsuitability",
                     limits = c(0,1)); SI_fig8

ggsave(filename = "outputs/Figures/SI_Fig8.png", SI_fig8, bg="white", height=10, width=10, dpi=300)


# SFig9 ====================================================
# Sensitive analysis
# comparison of model performance with and without disease-specific thermal suitability included as a covariate

# SI_fig9a -------------------------------------------------
# visualise the AUC tables as heatmap

base <- read.csv("outputs/cross_validation/SA_TCur_base_AUC_table.csv") # base model
mod2 <- read.csv("outputs/cross_validation/SA_TCur_DenAe_AUC_table.csv") # base model + TCur_den_ae
mod3 <- read.csv("outputs/cross_validation/SA_TCur_DenZikAe_AUC_table.csv") # base model + TCur_den_ae + TCur_Zik_ae

all_tab <- rbind(base %>% mutate(id = "base"), 
                 mod2 %>% mutate(id = "DEN_ae") 
                 )
all_tab <- rbind(all_tab, 
                 mod3 %>% mutate(id = "DEN_ZIK_ae"))
all_tab$Region <- factor(all_tab$Region, levels = c("Overall", "Endemic", "Non-endemic", 
                                                    "South America", "North America", "Africa", 
                                                    "Europe", "Asia", "Oceania"))
all_tab$disease <- factor(all_tab$disease, levels = c("dengue", "chikungunya", "zika"))

metrics <- names(all_tab)[2:4]

AUC_hm <- list()
for (metric in metrics) {
  m_title <- ifelse(metric == "meanAUC", "AUC", 
                    ifelse(metric == "meanSens", "Sensitivity", "Specificity"))
  
   p <- ggplot(data=all_tab[all_tab$Region == "Overall",])+
      geom_tile(aes(x=id, y=disease, fill=!!sym(metric)), color="white")+
      
      geom_text(aes(x = id, y=disease,
                    label = paste(sprintf("%.3f", round(!!sym(metric), 3)))), color = "black", size= 5)+
      scale_y_discrete(limits=rev)+
      scale_x_discrete(position = "top")+
      theme(plot.title = element_text(hjust = 0.5, vjust = 1.5),
            axis.text.x.top = element_text(size=11, margin = margin(t = 0.1, unit = "in")),
            axis.text.y = element_text(size=12),
            axis.ticks.x = element_blank(),        
            axis.ticks.y = element_blank(), 
            axis.line = element_blank(),
            panel.background = element_rect(fill = 'white'), 
            panel.grid.major = element_blank(),
            text= element_text(family = "Arial")
      )+  
      ylab(NULL)+xlab(NULL)+
      ggtitle(m_title)+
      scale_fill_distiller(
        name = "Metric",
        palette = "Spectral",
        values = c(0, 0.5,  0.7, 0.8, 1),
        breaks = c(0, 0.5,  0.7, 0.8, 0.9, 1),
        na.value = "grey80",
        direction = 1,
        limits = c(0, 1)
      )
      AUC_hm[[metric]] <- p

} 

SI_fig9a <- ggarrange(AUC_hm[["meanAUC"]], 
                   AUC_hm[["meanSens"]], 
                   AUC_hm[["meanSpec"]], 
                   ncol=3, nrow=1, 
                   common.legend = TRUE, 
                   legend = "right"); SI_fig9a

# SI_fig9b -------------------------------------------------


AUC_hm_d <- list()
diseases = c("dengue", "chikungunya", "zika")
for (dis in diseases) { 
  d_title <- ifelse(dis == "dengue", "Dengue", 
                    ifelse(dis == "zika", "Zika", "Chikungunya"))
  p <- ggplot(data=all_tab[all_tab$disease == dis & !all_tab$Region %in% c("Endemic", "Non-endemic"),])+
    geom_tile(aes(x=id, y=Region, fill=meanAUC), color="white")+
    
    geom_text(aes(x = id, y=Region,
                  label = paste(sprintf("%.3f", round(meanAUC, 3)))), color = "black", size= 5)+
    scale_y_discrete(limits=rev)+
    scale_x_discrete(position = "top")+
    theme(plot.title = element_text(hjust = 0.5, vjust = 2),
          axis.text.x.top = element_text(size=11, margin = margin(t = 0.1, unit = "in")),
          axis.text.y = element_text(size=12, margin = margin(r = 0, unit = "in")),
          axis.ticks.x = element_blank(),        
          axis.ticks.y = element_blank(), 
          axis.line = element_blank(),
          panel.background = element_rect(fill = 'white'), 
          panel.grid.major = element_blank(),
          text= element_text(family = "Arial")
    )+  
    ylab(NULL)+xlab(NULL)+
    ggtitle(d_title)+
    scale_fill_distiller(
      name = "Metric",
      palette = "Spectral",
      values = c(0, 0.5,  0.7, 0.8, 1),
      breaks = c(0, 0.5,  0.7, 0.8, 0.9, 1),
      na.value = "grey80",
      direction = 1,
      limits = c(0, 1)
    )
  
  AUC_hm_d[[dis]] <- p
}

SI_fig9b <- ggarrange(AUC_hm_d[["dengue"]], 
                   AUC_hm_d[["chikungunya"]], 
                   AUC_hm_d[["zika"]], 
                   ncol=3, nrow=1, 
                   common.legend = TRUE, 
                   legend = "right"); SI_fig9b

# Merge everything
SI_fig9  <- ggarrange(SI_fig9a, SI_fig9b,
                  ncol=1, nrow=2,
                  heights = c(0.4, 0.6),
                  common.legend = TRUE, 
                  legend = "right",
                  align = "hv",
                  labels = c("a", "b")); all

ggsave(filename = "outputs/Figures/SI_Fig9.png", SI_fig9, height=9, width=12, dpi=300, bg = "white")


# SFig10 ====================================================
# Uncertainty maps for dengue, zika, chikungunya, yellow fever

IQR_list <- list.files("outputs/Rasters", pattern = "*IQR_unmasked.tif", full.names = FALSE); IQR_list
IQR_plot <- list()
for (i in 1:4) { 
  Ras <- raster(paste0("outputs/Rasters/", IQR_list[i]))
  
  disease_abb <- unlist(strsplit(names(Ras), "_"))[[1]] 
  plot_dtitle <- ifelse(disease_abb=="DEN", "a  Dengue", 
                        ifelse(disease_abb=="CHIK", 'b  Chikungunya', 
                               ifelse(disease_abb=="ZIK", "c  Zika", "d  Yellow fever")))
  
  IQR_plot[[disease_abb]] <- plotRaster(Ras) + 
    theme(legend.title = element_text(family = "Arial", vjust = 1))+

    scale_fill_viridis(option = "rocket", direction = -1, na.value="white", 
                       name= "Uncertainty")+
    labs(title = plot_dtitle)
  
}

SI_fig10 <- ggarrange(IQR_plot$DEN + theme(legend.key.width=unit(1,"cm")), 
                     IQR_plot$CHIK, 
                     IQR_plot$ZIK, 
                     IQR_plot$YF, 
                     ncol = 2, nrow = 2, 
                     common.legend = TRUE, 
                     legend = "bottom")

ggsave(filename = "outputs/Figures/SI_Fig10.png", SI_fig10, height=9, width=18, dpi=300, bg = "white")


# SFig11 ====================================================
# Overall and regionally-stratified model performance for arbo and yf model

# load in tables
OOB_out <- readRDS("outputs/cross_validation/Arbo_OOB_data.rds")
AUC_arbo <- AUC_strata_arbo(OOB_out) 

yf_OOB_out <- readRDS("outputs/cross_validation/YF_OOB_data.rds")
AUC_yf <- AUC_strata(subset(yf_OOB_out, yf_OOB_out$disease == "yf"), "YF")

arbo_tab <- AUC_arbo$tab
yf_tab <- AUC_yf$tab

all_tab <- rbind(arbo_tab, yf_tab) 
all_tab$Region <- factor(all_tab$Region, levels = c("Overall", "Endemic", "Non-endemic", 
                                                    "South America", "North America", "Africa", 
                                                    "Europe", "Asia", "Oceania"))
all_tab$disease <- factor(all_tab$disease, levels = c("dengue", "chikungunya", "zika", "yf"))

metrics <- names(all_tab)[2:4]

AUC_hm <- list()
for (metric in metrics) { 
  m_title <- ifelse(metric == "meanAUC", "AUC", 
                    ifelse(metric == "meanSens", "Sensitivity", "Specificity"))
  p <- ggplot(data=all_tab[!all_tab$Region %in% c("Endemic", "Non-endemic"),])+
    geom_tile(aes(x=disease, y=Region, fill=!!sym(metric)), color="white")+
    
    geom_text(aes(x = disease, y=Region,
                  label = paste(sprintf("%.2f", round(!!sym(metric), 2)))), color = "black", size= 3)+
    scale_y_discrete(limits=rev)+
    scale_x_discrete(position = "top", labels = c("Dengue", "Chikungunya", "Zika", "Yellow fever"))+
    theme(plot.title = element_text(hjust = 0.5, vjust = 1.5),
          axis.ticks.x = element_blank(),        
          axis.ticks.y = element_blank(), 
          axis.line = element_blank(),
          panel.background = element_rect(fill = 'white'), 
          panel.grid.major = element_blank(),
          text= element_text(family = "Arial")
    )+  
    ylab(NULL)+xlab(NULL)+
    ggtitle(m_title)+
    scale_fill_distiller(
      name = "Metric",
      palette = "Spectral",
      values = c(0, 0.5,  0.7, 0.8, 1),
      breaks = c(0, 0.5,  0.7, 0.8, 0.9, 1),
      na.value = "grey80",
      direction = 1,
      limits = c(0, 1)
    )
  
  AUC_hm[[metric]] <- p
}

SI_fig11 <- ggarrange(AUC_hm[["meanAUC"]], 
                   AUC_hm[["meanSens"]], 
                   AUC_hm[["meanSpec"]], 
                   ncol=3, nrow=1, 
                   common.legend = TRUE, 
                   legend = "right")

ggsave(filename = "outputs/Figures/SI_Fig11.png", SI_fig11, height=6, width=12, dpi=300, bg = "white")


# SFig12 ====================================================
# Overall and regionally-stratified model performance for arbo and yf model

OOB_out <- readRDS("outputs/cross_validation/Arbo_OOB_data.rds")
AUC_arbo <- AUC_strata_arbo(OOB_out) 

yf_OOB_out <- readRDS("outputs/cross_validation/YF_OOB_data.rds")
AUC_yf <- AUC_strata(subset(yf_OOB_out, yf_OOB_out$disease == "yf"), "YF")

# save AUC maps 
(SI_fig12 <- ggarrange(AUC_arbo$map[[1]] + ggtitle("A  Dengue") +
                        theme(legend.key.width=unit(1,"cm")), 
                      AUC_arbo$map[[2]] + ggtitle("B  Chikungunya") ,
                      AUC_arbo$map[[3]] + ggtitle("C  Zika"), 
                      AUC_yf$map +  ggtitle("D  Yellow fever"),
                      ncol=2, nrow=2,
                      common.legend = TRUE, 
                      legend = "bottom"))

ggsave(filename = "outputs/Figures/SI_Fig12.png", SI_fig12, height=9, width=18, dpi=300, bg = "white")

# SFig13 ====================================================
# variable importance and partial dependence plots for yf model
# load in variable importance and partial dependence values
vi_yf <- readRDS("outputs/cross_validation/YF_VI.rds")

viPlot_yf <- viPlot(vi_yf)

ggsave(filename = "outputs/Figures/SI_Fig13a.png", viPlot_yf[[1]], height=6, width=9, dpi=300)
ggsave(filename = "outputs/Figures/SI_Fig13b.png", viPlot_yf[[2]], height=6, width=10, dpi=300)


# SFig14 ====================================================
# Sensitive analysis
# comparison of model performance with and without disease-specific thermal suitability included as a covariate

# SI_fig14a -------------------------------------------------  
ras <- raster("outputs/cross_validation/SA_YF_noAe_riskmap_wmean_unmasked.tif")
yf_eye <- raster("data/covariate_rasters/YF_EYE_mask.tif")
mask_yf = yf_eye > 0 

YF_range_mask <- ras * mask_yf
SI_fig14a <- plotRaster(YF_range_mask)+
  scale_fill_viridis(option = "rocket", direction = -1, na.value = "transparent", 
                     name= "Environmental\nsuitability", 
                     limits = c(0,1))+
  coord_sf(xlim = c(-115, 65), ylim = c(-59, 38), expand = FALSE) 


# SI_fig14b -------------------------------------------------  
base <- read.csv("outputs/cross_validation/SA_YF_base_AUC_table.csv") # base model
mod2 <- read.csv("outputs/cross_validation/SA_YF_noAe_AUC_table.csv") # without Aegypti included as a covariate


# visualise the AUC tables as heatmap
all_tab <- rbind(base %>% mutate(id = "base"), 
                 mod2 %>% mutate(id = "wo_aegypti") 
                 )

all_tab$Region <- factor(all_tab$Region, levels = c("Overall", "Endemic", "Non-endemic", 
                                                    "South America", "North America", "Africa", 
                                                    "Europe", "Asia", "Oceania"))

all_tab_long <- gather(all_tab, metrics, value, meanAUC:meanSpec, factor_key=TRUE)
all_tab_long$metrics <- ifelse(all_tab_long$metrics == "meanAUC", "AUC", 
                    ifelse(all_tab_long$metrics == "meanSens", "Sensitivity", "Specificity"))
  
SI_fig14b <- ggplot(data=all_tab_long[all_tab_long$Region %in% c("Overall", "South America", "Africa"),])+
    geom_tile(aes(x=id, y=Region, fill=value), color="white")+
    
    geom_text(aes(x = id, y=Region,
                  label = paste(sprintf("%.3f", round(value, 3)))), color = "black", size= 4)+
    scale_y_discrete(limits=rev)+
    scale_x_discrete(position = "top")+
    theme(plot.title = element_text(hjust = 0.5, vjust = 1.5),
          axis.text.x.top = element_text(size=11),
          axis.text.y = element_text(size=11),
          axis.ticks.x = element_blank(),        
          axis.ticks.y = element_blank(), 
          axis.line = element_blank(),
          panel.background = element_rect(fill = 'white'), 
          panel.grid.major = element_blank(),
          text= element_text(family = "Arial")
    )+  
    ylab(NULL)+xlab(NULL)+
    
    scale_fill_distiller(
      name = "Metric",
      palette = "Spectral",
      values = c(0, 0.5,  0.7, 0.8, 1),
      breaks = c(0, 0.5,  0.7, 0.8, 0.9, 1),
      na.value = "grey80",
      direction = 1,
      limits = c(0, 1)
    )+
   theme(strip.placement = 'outside', 
         strip.background = element_blank(), 
         strip.text = element_text(size = 12, face = "bold"))+
   facet_wrap(metrics~., nrow = 1, ncol = 3); SI_fig14b
 

SI_fig14  <- ggarrange(SI_fig14a, SI_fig14b, 
                  ncol=1, nrow=2,
                  # heights = c(0.4, 0.6),
                  # common.legend = TRUE, 
                  legend = "right",
                  align = "hv",
                  labels = c("a", "b")); SI_fig14


ggsave(filename = "outputs/Figures/SI_Fig14.png", SI_fig14, height=8, width=8, dpi=300, bg = "white")


# SFig15 ====================================================

# covariates for surveillance model
source("script/01_get_covariates.R")
surv_cov_list <- loadRasters("surv")[[1]]
arbo_cov_list <- loadRasters("arbo")[[1]]

brews <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Oranges", 
           "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "YlGn", "YlGnBu", 
           "YlOrBr")

title_lookup <- c(
    "GDP" = "GDP",
    "GDP_national" = "GDP (national)",
    "Urban" = "Urbanisation", 
    "Acc_walk" = "Travel time to health facilities (by walk)", 
    "Acc_city" = "Travel time to cities", 
    "Trmt" = "Fever treatment seeking rate", 
    "U5M" = "Child mortality", 
    "GovEff" = "Government effectiveness",
    "Physician" = "Physicians density", 
    "Tcold" = "Temperature of the coldest month", 
    "Tsuit" = "Temperature suitability", 
    "PRCP" = "Precipitation", 
    "NDVI" = "Normalized difference vegetation index", 
    "DHI" = "Dynamic habitat indices", 
    "Aegypti" = "Aedes aegypti", 
    "Albo" = "Aedes albopictus", 
    "Pop" = "Population", 
    "Hg" = "Haemagogus janthinomys", 
    "NHP" = "Non-human primates" 
    
)


all_stack <- stack(surv_cov_list, 
                   arbo_cov_list[[-which(names(arbo_cov_list) %in% c("GDP", "GDP_national", "Urban", "Vaccine", "Vaccine_offset", "Surv"))]])
all_names <- names(all_stack)

plist <- list()
for (name in all_names) { 
  full_title = title_lookup[[name]]
  
  # Choose a palette from brews
  pal.index <- sample(length(brews), 1) # Choose a random index
  pal.name <- brews[pal.index]
  
  p <-  plotRaster(all_stack[[name]]) +
   scale_fill_distiller(type= "seq", 
                        palette=pal.name, 
                        direction =1, 
                        na.value="transparent"
                        )+
   theme(legend.position = "bottom")+
   labs(fill = paste0(full_title))
  
  plist[[name]] <- p

}

names(plist)
SI_fig15 <- wrap_plots(plist[1:9], ncol=2)
SI_fig16 <- wrap_plots(plist[10:19], ncol=2)

ggsave(filename = "outputs/Figures/SI_Fig15.png", SI_fig15, height = 18, width = 15, dpi=300)
ggsave(filename = "outputs/Figures/SI_Fig16.png", SI_fig16, height = 18, width = 15, dpi=300)

