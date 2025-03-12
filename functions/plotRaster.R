base_map <- function(world, lake_mask = T) {
  # world map
  if (lake_mask) {
    world <- st_read("data/admin_rasters/world_no_lakes.shp", quiet = TRUE)
  } else {
    world <- ne_countries(scale = 50, type = "countries", returnclass = "sf")
  }

  ggplot() +
    geom_sf(data = world, fill = "white", color = "grey70", linewidth = 0.1) +
    theme_bw() +
    coord_sf(xlim = c(-180, 180), ylim = c(-59, 75), expand = FALSE) +
    theme(
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank()
    )
}

# to convert dataframe into a raster
mapRas <- function(pred_df) {
  output <- template
  values(output) = as.numeric(pred_df[, 2])
  return(output)
}

# to convert summary data (weighted average/IQR) into a raster
bootsRas <- function(all.p, measure) {
  df <- pred.data[, c("x", "y")]
  df$z <- all.p[[measure]]
  bootsRas <- df %>% rasterFromXYZ()
  return(bootsRas)
}

# to plot a raster
plotRaster <- function(raster) {
  world <- st_read("data/admin_rasters/world_no_lakes.shp", quiet = TRUE)

  ggplot() +
    geom_spatraster(data = rast(raster)) +
    geom_sf(data = world, fill = "transparent", color = "grey70", lwd = 0.1) +
    theme_bw() +
    coord_sf(xlim = c(-180, 180), ylim = c(-59, 75), expand = FALSE) +
    theme(
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(family = "sans", vjust = 2.5),
      legend.text = element_text(family = "sans"),
      panel.grid.major = element_blank()
    ) +
    scale_fill_viridis(na.value = "white", limits = c(0, 1))
}

# to compare binary maps and convert the results into a plot (for Fig4)
plotFRaster <- function(raster) {
  # brewer.pal(4, "YlGnBu")
  cols <- c(
    "4" = "#FC8D59",
    "3" = "#74c476", #"#FFDF7D",
    "2" = "#2B83BA",
    "1" = "#F0F0F0" #white"
  )
  world <- st_read("data/admin_rasters/world_no_lakes.shp", quiet = TRUE)

  ggplot() +
    geom_spatraster(data = rast(raster)) +
    geom_sf(data = world, fill = "transparent", color = "grey70", lwd = 0.1) +
    theme_bw() +
    coord_sf(xlim = c(-180, 180), ylim = c(-59, 75), expand = FALSE) +
    theme(
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank()
    ) +

    scale_fill_manual(
      name = "",
      values = cols,
      breaks = c("4", "3", "2", "1"),
      labels = c("Both high", "New high", "Previous high", "Both low"),
      na.value = "white"
    )
}

today = gsub("-", "_", Sys.Date())
# to save plots and rasters
saveRasPlot <- function(df, disease, color_opt, color_direction) {
  bootsRas2 <- bootsRas(df, "IQR")
  # p1 <- plotRaster(bootsRas2) + scale_fill_viridis(option = color_opt, direction = color_direction, na.value="white")

  bootsRas3 <- bootsRas(df, "wmean")
  # p2 <- plotRaster(bootsRas3) + scale_fill_viridis(option = color_opt, direction = color_direction, na.value="white")

  if (!disease == "Surv") {
    name = "_riskmap_"
    writeRaster(
      bootsRas2,
      filename = paste0(
        "outputs/Rasters/",
        disease,
        name,
        "IQR_unmasked",
        today,
        ".tif"
      ),
      overwrite = T
    )
    writeRaster(
      bootsRas3,
      filename = paste0(
        "outputs/Rasters/",
        disease,
        name,
        "wmean_unmasked",
        today,
        ".tif"
      ),
      overwrite = T
    )
  } else {
    name = "eillance_map_"
    writeRaster(
      bootsRas2,
      filename = paste0(
        "outputs/Rasters/",
        disease,
        name,
        "IQR",
        today,
        ".tif"
      ),
      overwrite = T
    )
    writeRaster(
      bootsRas3,
      filename = paste0(
        "outputs/Rasters/",
        disease,
        name,
        "wmean",
        today,
        ".tif"
      ),
      overwrite = T
    )
  }

  return(bootsRas3)
}
