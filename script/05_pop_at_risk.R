# ==========================================================
# 
# Calculating population at risk across regions and countries
# 
# ==========================================================

source("script/00_setup.R")

# load in binary maps 
dcz_bin <- raster("outputs/Rasters/DCZ_binmap_mean.tif")
dcz_lwr_bin <- raster("outputs/Rasters/DCZ_binmap_lwr.tif")
dcz_upr_bin <- raster("outputs/Rasters/DCZ_binmap_upr.tif")

yf_bin <- raster("outputs/Rasters/YF_binmap_mean.tif")
yf_lwr_bin <- raster("outputs/Rasters/YF_binmap_lwr.tif")
yf_upr_bin <- raster("outputs/Rasters/YF_binmap_upr.tif")

# load in population raster (landscan 2022)
pop <- raster("data/covariate_rasters/landscan_global_2022_masked_.tif")

# load in surveillance capability map
Surv <- raster("outputs/Rasters/Surveillance_map_wmean.tif")
Surv_w <- Surv*pop  # for population-weighted mean

# calculating population at risk per pixel
dcz_risk_pop <- dcz_bin*pop
dcz_lwr_pop <- dcz_lwr_bin*pop
dcz_upr_pop <- dcz_upr_bin*pop

yf_risk_pop <- yf_bin*pop
yf_lwr_pop <- yf_lwr_bin*pop
yf_upr_pop <- yf_upr_bin*pop


# load standardised UN states namaes
source("functions/RNE_country_names.R")
world$continent <- ifelse(world$continent %in% c("North America", "South America"), "Americas", world$continent)

world$pop <- exact_extract(pop, world, "sum") #  the sum of fractions of raster cells with non-NA values covered by the polygon

# extract pixel values based on country boundaries
world$Surv_w <- exact_extract(Surv_w, world, "sum")

world$dcz_risk_pop <- exact_extract(dcz_risk_pop, world, "sum")
world$dcz_lwr_pop <- exact_extract(dcz_lwr_pop, world, "sum") 
world$dcz_upr_pop <- exact_extract(dcz_upr_pop, world, "sum") 

world$yf_risk_pop <- exact_extract(yf_risk_pop, world, "sum") 
world$yf_lwr_pop <- exact_extract(yf_lwr_pop, world, "sum") 
world$yf_upr_pop <- exact_extract(yf_upr_pop, world, "sum") 

# percentage of population at risk of dengue, chikungunya, and zika 
sum(values(dcz_risk_pop), na.rm = TRUE)/ sum(values(pop), na.rm=T) # 0.734476

# Pop at risk by continent
global <- 
  data.frame(
  continent = "Global",
  dcz = round(sum(values(dcz_risk_pop), na.rm = TRUE) * 0.000000001, 2),
  dcz_lwr = round(sum(values(dcz_lwr_pop), na.rm = TRUE) * 0.000000001, 2),
  dcz_upr = round(sum(values(dcz_upr_pop), na.rm = TRUE) * 0.000000001, 2),
  yf = round(sum(values(yf_risk_pop), na.rm = TRUE) * 0.000000001, 2),
  yf_lwr = round(sum(values(yf_lwr_pop), na.rm = TRUE) * 0.000000001, 2),
  yf_upr = round(sum(values(yf_upr_pop), na.rm = TRUE) * 0.000000001, 2)
)

continent <- world %>%
  st_drop_geometry()%>%
  group_by(continent)%>%
  summarise(dcz = round(sum(dcz_risk_pop, na.rm=T)*0.000000001,2), 
            dcz_lwr = round(sum(dcz_lwr_pop, na.rm=T)*0.000000001,2),
            dcz_upr = round(sum(dcz_upr_pop, na.rm=T)*0.000000001,2),
            yf = round(sum(yf_risk_pop, na.rm=T)*0.000000001,2),
            yf_lwr = round(sum(yf_lwr_pop, na.rm=T)*0.000000001,2),
            yf_upr = round(sum(yf_upr_pop, na.rm=T)*0.000000001,2)
            
            )

# for table formatting
format_with_range <- function(value, lower, upper) {
  formatted_value <- format(value, nsmall = 2)
  formatted_lower <- format(lower, nsmall = 2)
  formatted_upper <- format(upper, nsmall = 2)
  paste0(formatted_value, " (", formatted_lower, "-", formatted_upper, ")")
}

table <- bind_rows(global, continent)%>%
  mutate(
    DCZ = format_with_range(dcz, dcz_lwr, dcz_upr),
    YF = format_with_range(yf, yf_lwr, yf_upr)
  )%>%
  select(-(dcz:yf_upr))%>%
  filter(!continent %in% c("Antarctica", "Seven seas (open ocean)")); table

write.csv(table, paste0("outputs/Tables/Table1_Pop_at_risk_continent.csv"), row.names=F)

# Pop at risk by country
unique(world$standard_sovereignt) # 193 UN states

country_at_risk <- world %>%
  st_drop_geometry()%>%
  group_by(standard_sovereignt)%>%

  summarise(tot_pop = sum(pop, na.rm=T),
            tot_Surv = sum(Surv_w, na.rm=T),
            sum_dcz_risk_pop = sum(dcz_risk_pop, na.rm=T),
            sum_dcz_lwr_pop = sum(dcz_lwr_pop, na.rm=T),
            sum_dcz_upr_pop = sum(dcz_upr_pop, na.rm=T),
            sum_yf_risk_pop = sum(yf_risk_pop, na.rm=T), 
            sum_yf_lwr_pop = sum(yf_lwr_pop, na.rm=T), 
            sum_yf_upr_pop = sum(yf_upr_pop, na.rm=T), 
            
            prop_dcz = sum_dcz_risk_pop/tot_pop, 
            prop_yf = sum_yf_risk_pop/tot_pop, 
            Surv_wm = tot_Surv/tot_pop)


world %>%
  st_drop_geometry() %>%
  group_by(continent, standard_sovereignt)%>% 
  tally() %>% arrange(desc(n))

# Number of countries at risk
# def: >10% of total population at-risk
country_at_risk %>% filter(prop_dcz >0.1) #169
country_at_risk %>% filter(prop_yf > 0.1) # 54

# formatting table
table2 <- country_at_risk %>%
  mutate(
    surv_wm = as.numeric(format(round(Surv_wm, 3), nsmall=3)),
    tot_pop = as.numeric(format(round(tot_pop*0.000001, 3), nsmall=3)),
    pop_dcz = as.numeric(format(round(sum_dcz_risk_pop*0.0000001, 3), nsmall=3)),
    pop_dcz_lwr = as.numeric(format(round(sum_dcz_lwr_pop*0.0000001, 3), nsmall=3)),
    pop_dcz_upr = as.numeric(format(round(sum_dcz_upr_pop*0.0000001, 3), nsmall=3)),
    
    pop_yf = as.numeric(format(round(sum_yf_risk_pop*0.0000001, 3), nsmall=3)), 
    pop_yf_lwr = as.numeric(format(round(sum_yf_lwr_pop*0.0000001, 3), nsmall=3)), 
    pop_yf_upr = as.numeric(format(round(sum_yf_upr_pop*0.0000001, 3), nsmall=3)), 
    
    prop_dcz = as.numeric(format(round(prop_dcz,3), nsmall=3)),
    prop_yf = as.numeric(format(round(prop_yf, 3), nsmall=3))
  )%>%
  select( standard_sovereignt, surv_wm, tot_pop, pop_dcz:pop_yf_upr)

names(table2) <- c("Sovereignt", 
                   "Surveillance_score", 
                   "Total_population", 
                   "DCZ_pop_risk_mean", 
                   "DCZ_pop_risk_lower",
                   "DCZ_pop_risk_upper", 
                   "YF_pop_risk_mean", 
                   "YF_pop_risk_lower", 
                   "YF_pop_risk_upper")

# remove latin accent
table2$Sovereignt <- stri_trans_general(str = table2$Sovereignt, id = "Latin-ASCII")

write.csv(table2, paste0("outputs/Tables/STable_Countries_at_risk.csv"), row.names=F, fileEncoding = "UTF-8")
