# function that extracts both point and polygon data from rasters

adminExtract <- function (occurrence, covariates, admin, fun = "mean") 
{
  # point extraction
  point_idx <- which(occurrence$Admin == -999)
  point_coords <- st_as_sf(occurrence[point_idx, ], coords = c("Longitude", "Latitude"))
  st_crs(point_coords) = 4326
  vals  = raster::extract(covariates, point_coords, method = "simple")
  
  # polygon extraction
  poly_idx <- which(occurrence$Admin != -999)
  ex <- matrix(NA, nrow = length(poly_idx), ncol = nlayers(covariates))
  colnames(ex) <- names(covariates)
  all_admins <- occurrence$Admin[poly_idx]
  all_GAULs <- occurrence$GAUL[poly_idx]
  for (level in 0:2) {
    level_idx <- all_admins == level
    if (any(level_idx)) {
      level_GAULs <- all_GAULs[level_idx]
      ad <- admin[[level + 1]]
      keep <- function(cells, ...) {
        ifelse(cells %in% level_GAULs, cells, NA)
      }
      ad <- calc(ad, keep)
      zones <- zonal(covariates, ad, fun = fun)
      which_zones <- match(level_GAULs, zones[, 1])
      ex[level_idx, ] <- zones[which_zones, -1]
    }
  }
  
  # recombining poitn and polygon results
  rex <- matrix(NA, nrow = nrow(occurrence), ncol = nlayers(covariates))
  colnames(rex) <- names(covariates)
  rex[point_idx, ] = vals
  rex[poly_idx, ] = ex
  
  return(rex)
}