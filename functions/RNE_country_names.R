library(countrycode)
dt <- codelist %>% select(region, region23, country.name.en, iso3c, un.name.en)

dt_un <- dt %>% 
  filter(!is.na(un.name.en))

rne <- rnaturalearth::ne_countries(scale =10 , type= "countries", returnclass = "sf") %>%
  select(continent, sovereignt, geounit, adm0_a3, iso_a3)

rne_countries <- rne %>%
  st_drop_geometry()%>%
  group_by(sovereignt)%>%
  tally()

# un states that are not in rne list
dt_un %>% 
  select(country.name.en, un.name.en) %>% 
  filter(!country.name.en %in% rne_countries$sovereignt)%>%
  filter(!un.name.en %in% rne_countries$sovereignt)

# countries that are not standardised
not_standard <- rne_countries %>%
  filter(!sovereignt %in% dt_un$un.name.en)

standard <- dt_un %>%
  select(country.name.en, un.name.en) %>% 
  filter(country.name.en %in% not_standard$sovereignt)

rne_countries <- rne_countries %>%
  merge(., standard, by.x="sovereignt", by.y="country.name.en", all=T)%>%
  mutate(standard_sovereignt = ifelse(!is.na(un.name.en), paste0(un.name.en), paste0(sovereignt)))%>%
  select(-n, -un.name.en)

# manual changes
rne_countries$standard_sovereignt[rne_countries$sovereignt == "The Bahamas"] <- "Bahamas"
rne_countries$standard_sovereignt[rne_countries$sovereignt == "Republic of the Congo"] <- "Congo"
rne_countries$standard_sovereignt[rne_countries$sovereignt == "Ivory Coast"] <- "Côte d’Ivoire"
rne_countries$standard_sovereignt[rne_countries$sovereignt == "eSwatini"] <- "Eswatini"
rne_countries$standard_sovereignt[rne_countries$sovereignt == "Federated States of Micronesia"] <- "Micronesia (Federated States of)"
rne_countries$standard_sovereignt[rne_countries$sovereignt == "Republic of Serbia"] <- "Serbia"
rne_countries$standard_sovereignt[rne_countries$sovereignt == "São Tomé and Principe"] <- "Sao Tome and Principe"
rne_countries$standard_sovereignt[rne_countries$sovereignt == "East Timor"] <- "Timor-Leste"
rne_countries$standard_sovereignt[rne_countries$sovereignt %in% c("Cyprus No Mans Area", "Northern Cyprus")] <- "Cyprus"

# attach region 
rne_countries <- rne_countries %>%
  merge(., dt_un[, c("region", "un.name.en")], by.x="standard_sovereignt", by.y="un.name.en", all.x=T)

rm(dt, rne, not_standard, standard)

world <- rnaturalearth::ne_countries(scale =10 , type= "countries", returnclass = "sf") %>%
  select(continent, sovereignt, geounit, adm0_a3, iso_a3)%>%
  merge(., rne_countries, by=c("sovereignt"), all.x=T, all.y=T)%>%
  filter(!sovereignt == "Antartica")

world$standard_sovereignt[world$standard_sovereignt == "Taiwan"] <- "China"
world$standard_sovereignt[world$standard_sovereignt == "Kosovo"] <- "Serbia"
world$standard_sovereignt[world$standard_sovereignt == "Somaliland"] <- "Somalia"
world$standard_sovereignt[world$standard_sovereignt == "Kashmir"] <- "India"
world$standard_sovereignt[world$standard_sovereignt == "Vatican"] <- "Italy"

world <- world %>%
  filter(standard_sovereignt %in% dt_un$un.name.en)

world$continent[world$geounit == "Reunion"] <- "Africa"
world$continent[world$geounit == "Seychelles"] <- "Africa"
world$continent[world$geounit == "Heard Island and McDonald Islands"] <- "Oceania"
world$continent[world$geounit == "Saint Helena"] <- "Africa"
world$continent[world$geounit == "Mauritius"] <- "Africa"
world$continent[world$geounit == "Azores"] <- "Europe"
world$continent[world$geounit == "Maldives"] <- "Asia"
world$continent[world$geounit == "Clipperton Island"] <- "Oceania"


# endemic countries (based on cdc yellow book 2024: https://wwwnc.cdc.gov/travel/yellowbook/2024/infections-diseases/dengue)
world <- world %>%
  mutate(endemic = ifelse(continent %in% c("Europe", "Antarctica") | 
                          standard_sovereignt %in% c("Greenland" , "United States", "Canada", "Morocco", "Algeria", "Tunisia", 
                                                       "Libya", "Botswana", "South Africa", "Lesotho", "Eswatini", 
                                                       "Saudi Arabia", "Iran (Islamic Republic of)", "Iraq", "Jordan","Israel",
                                                       "Lebanon", "Cyprus", "Syrian Arab Republic", "Türkiye", "Georgia", 
                                                       "Armenia", "Azerbaijan", "Kuwait", "United Arab Emirates",
                                                       "Kazakhstan", "Turkmenistan", "Uzbekistan", "Kyrgyzstan",
                                                       "Tajikistan","China", "Mongolia", "Democratic People’s Republic of Korea", 
                                                       "Republic of Korea", "Japan", "Australia"), "Non-endemic",  "Endemic"))



# c <- ggplot(world %>%
#               filter(continent != "Seven seas (open ocean)"))+
#   geom_sf(aes(fill = continent))+
#   theme_bw()
# 
# ggsave(filename = paste0("Maps_and_plots/Figures/world_map", gsub("-", "_", Sys.Date()), ".png"), 
#        c, height=6, width=12, dpi=300, bg = "white")
