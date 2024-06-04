
fixNAs <- function(dataInput, covariates){
  #fix NAs in points 
  point_idx <- which(dataInput$Admin == -999)
  xy <- c("Longitude", "Latitude")
  cov_names <- names(covariates)

  # Loop through additional columns
  for (col in cov_names) { 
    selected_columns <- c(xy, col)
  
    dt <- dataInput[point_idx, selected_columns][is.na(dataInput[point_idx, col]), ]
  
    dat_pt <- terra::vect(dt, geom=c("Longitude", "Latitude"))
    crs(dat_pt) <- "epsg:4326"
  
    terra_buffer <- dat_pt |>
      buffer(width = 500000) %>% # Unit is meter if x has a longitude/latitude CRS
      st_as_sf()
  
    dt[, col] <- exactextractr::exact_extract(covariates[[col]], terra_buffer, 'mean')
    dataInput[point_idx, selected_columns][is.na(dataInput[point_idx, col]), ] <- dt
    cat('Processing NAs in ', col, '\n')
  }
 
  return(
    dataInput
  )
}


