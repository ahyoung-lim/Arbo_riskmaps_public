# ==========================================================
# 
# Producing figures in the manuscript
# 
# ==========================================================

# Figure 1 =================================================
# The temporal and spatial distribution of Aedes-borne arbovirus occurrence points
# Fig 1a stacked bar ---------------------------------------

arbo_occ_thinned <- read.csv("data/intermediate_datasets/arbo_occ_thinned.csv")

arbo_occ_thinned$disease <- factor(arbo_occ_thinned$disease, levels = c("dengue", "zika", "chikungunya", "yf"))
arbo_occ_thinned$Type = c("Polygon", "Point")[((arbo_occ_thinned$Admin == -999) + 1)]


n_by_year <- arbo_occ_thinned %>%
  dplyr::filter(!Year == "0")%>%
  group_by(disease, Year)%>%
  tally()%>% 
  group_by(Year)%>%
  summarise(total = sum(n))

(arbo_bar <- arbo_occ_thinned %>%
    filter(!Year == "0")%>%
    filter(as.numeric(Year)>1984)%>%
    group_by(disease, Year)%>%
    tally()%>% 
    ggplot()+
    geom_bar(aes(x=as.numeric(Year), y=n, fill=disease), color = "grey20", 
             stat = "identity", 
             position = position_stack(reverse = TRUE))+
    scale_fill_manual(values = c("#EE6677", "#228833","#4477AA", "#AA3377"),
                      name = "Disease", 
                      labels = c("Dengue",  "Chikungunya", "Zika", "Yellow fever"))+
    scale_x_continuous(position = "bottom", 
                       breaks=seq(1985, 2020, by=5), 
                       expand=c(0,0))+
    scale_y_continuous(breaks = seq(0, 2000, by=500),
                       limits=c(0, 2000),
                        expand=c(0,0))+
    xlab(NULL)+ylab("Number of new unique occurrence points")+
    labs(tag = "a")+
    theme(text= element_text(family = "Arial"), 
          axis.title.y = element_text(size=15, margin = margin(r = 15)),
          axis.text.x = element_text(size=14, margin = margin(t = 8)), 
          axis.text.y = element_text(size=12, margin = margin(r = 8)), 
          axis.ticks.x = element_blank(),        
          axis.ticks.y = element_blank(), 
          plot.tag = element_text(size = 20),
          panel.background = element_rect(fill = 'white'), 
          plot.margin = unit(c(1,1,1,1), "cm"),
          axis.line.x = element_line(color="black", linewidth=1), 
          axis.line.y = element_blank(),
          panel.grid.major.y = element_line(color="#D9D9D9", size=0.8, linetype = "solid"), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.title = element_text(size=16, family = "Arial"),
          legend.text = element_text(size=16, family = "Arial", vjust=0.5),
          legend.text.align = 0, 
          legend.margin = ggplot2::margin(11,11,11,11), 
          legend.key = element_rect(fill="white")
    )+
    theme(legend.position="top")   
)

ggsave(filename = paste0("outputs/Figures/Fig1a.png"), arbo_bar, height = 6, width = 19.2, dpi=300)


# Fig 1b-e disease separately --------------------------------

(dengue_data_plot <- base_map(world) + 
    geom_point(data=subset(arbo_occ_thinned, disease == "dengue"),
               aes(x=Longitude, y=Latitude,  shape = Type), color = "#EE6677",
               alpha = 0.5, size = 2)+
    scale_shape_manual(values=c(1, 2))+
    guides(color="none")+
    ggtitle("b  Dengue")+
    theme(text = element_text(family = "Arial"), 
          plot.title = element_text(size=16, family = "Arial"),
          legend.title = element_text(size=10, family = "Arial"),
          legend.text = element_text(size=10, family = "Arial", vjust=0.5),
          legend.text.align = 0, 
          legend.position="top",
          legend.margin = margin(11,11,11,11), 
          legend.key = element_rect(fill="white")
    ))

(zika_data_plot <- base_map(world) + 
    geom_point(data=subset(arbo_for_map, disease == "zika"), color = "#4477AA",
               aes(x=Longitude, y=Latitude,  shape = Type), 
               alpha = 0.5, size = 2
    )+
    scale_shape_manual(values=c(1, 2))+
    guides(color="none")+
    ggtitle("d  Zika")+
    theme(text = element_text(family = "Arial"), 
          plot.title = element_text(size=16, family = "Arial"),
          legend.title = element_text(size=10, family = "Arial"),
          legend.text = element_text(size=10, family = "Arial", vjust=0.5),
          legend.text.align = 0, 
          legend.position="top",
          legend.margin = margin(11,11,11,11), 
          legend.key = element_rect(fill="white")
    ))

(chik_data_plot <- base_map(world) + 
    geom_point(data=subset(arbo_for_map, disease == "chikungunya"), color = "#228833",
               aes(x=Longitude, y=Latitude,  shape = Type), 
               alpha = 0.5, size = 2
    )+
    scale_shape_manual(values=c(1, 2))+
    guides(color="none")+
    ggtitle("c  Chikungunya")+
    theme(text = element_text(family = "Arial"), 
          plot.title = element_text(size=16, family = "Arial"),
          legend.title = element_text(size=10, family = "Arial"),
          legend.text = element_text(size=10, family = "Arial", vjust=0.5),
          legend.text.align = 0, 
          legend.position="top",
          legend.margin = margin(11,11,11,11), 
          legend.key = element_rect(fill="white")
    ))

(yf_data_plot <- base_map(world) + 
    geom_point(data=subset(arbo_for_map, disease == "yf"), color = "#AA3377",
               aes(x=Longitude, y=Latitude,  shape = Type), 
               alpha = 0.5, size = 2
    )+
    scale_shape_manual(values=c(1, 2))+
    guides(color="none")+
    ggtitle("e  Yellow fever")+
    theme(text = element_text(family = "Arial"), 
          plot.title = element_text(size=16, family = "Arial"),
          legend.title = element_text(size=10, family = "Arial"),
          legend.text = element_text(size=10, family = "Arial", vjust=0.5),
          legend.text.align = 0, 
          legend.position="top",
          legend.margin = margin(11,11,11,11), 
          legend.key = element_rect(fill="white")
    ))

fig1b_e <- dengue_data_plot + chik_data_plot + 
  zika_data_plot + yf_data_plot + plot_layout(ncol=2)

ggsave(filename = paste0("outputs/Figures/Fig1b_e.png"), fig1b_e, height = 9.6, width = 19.2, dpi=300)


# Figure 2 =================================================

Ras_Surv <- raster("outputs/Rasters/Surveillance_map_wmean.tif")

fig2 <- plotRaster(Ras_Surv) + 
  scale_fill_viridis(option = "viridis", direction = 1, na.value="white", 
                         name= "Surveillance\ncapability")+
  theme(legend.title = element_text(family = "Arial", vjust = 2.5)); fig2

ggsave("outputs/Figures/Fig2.png", fig2, height=6, width=12, dpi=300)


# Figure 3 =================================================

dcz_range <- raster("outputs/Rasters/DCZ_riskmap_wmean_masked.tif")
YF_range_mask <- raster("outputs/Rasters/YF_riskmap_wmean_masked.tif")

(r_DCZ <- plotRaster(dcz_range)+
    scale_fill_viridis(option = "rocket", direction = -1, na.value = "transparent",
                     name= "Environmental\nsuitability",
                     limits = c(0,1)))

(r_YF <- plotRaster(YF_range_mask)+
  scale_fill_viridis(option = "rocket", direction = -1, na.value = "transparent",
                     name= "Environmental\nsuitability",
                     limits = c(0,1))+
  coord_sf(xlim = c(-115, 65), ylim = c(-59, 38), expand = FALSE))

fig3 <- r_DCZ + ggtitle("a  Dengue, chikungunya, and Zika") +
  r_YF + ggtitle("b  Yellow fever") +
  plot_layout(ncol=1, widths=c(1,2), guides = "collect") &
  theme(text = element_text(family = "Arial"))&
  scale_fill_viridis(option = "rocket", direction = -1, na.value = "transparent",
                     name= "Environmental\nsuitability",
                     limits = c(0,1));fig3


ggsave(filename = "outputs/Figures/Fig3.png", fig3, bg="white", height=10, width=12, dpi=300)


# Figure 4 =================================================

# load in rasters
DEN_ex <- raster("data/intermediate_datasets/DEN_previous_binrast.tif")
CHIK_ex <- raster("data/intermediate_datasets/DEN_previous_binrast.tif")
ZIK_ex <- raster("data/intermediate_datasets/DEN_previous_binrast.tif")
YF_ex <- raster("data/intermediate_datasets/DEN_previous_binrast.tif")

DEN_range_mask <- raster("outputs/Rasters/DEN_riskmap_wmean_masked.tif")
CHIK_range_mask <- raster("outputs/Rasters/CHIK_riskmap_wmean_masked.tif")
ZIK_range_mask <- raster("outputs/Rasters/ZIK_riskmap_wmean_masked.tif")
YF_range_mask <- raster("outputs/Rasters/YF_riskmap_wmean_masked.tif")


# comparison maps
den_com <- mapROC(den_pres, "DEN", "compare")
zik_com <- mapROC(zik_pres, "ZIK", "compare")
chik_com <- mapROC(chik_pres, "CHIK", "compare")
yf_com <- mapROC(yf_pres, "YF", "compare")

fig4 <- ggarrange(den_com,
          chik_com,
          zik_com,
          yf_com,
          ncol=2, nrow=2,
          common.legend = TRUE,
          legend = "bottom")+
    theme(text = element_text(family = "Arial"));fig4


ggsave(filename = "outputs/Figures/Fig4.png",fig4, bg="white", height=9, width=18, dpi=300)

