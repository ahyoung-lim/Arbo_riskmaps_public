thin_occ_viral <- function(data){
  trv <- 1:length(as.vector(template)) # unique pixel values
  trv = trv * !is.na(as.vector(template))
  pixelID = template
  values(pixelID) = trv

  # extract unique ID for each occurrence point
  data$uID <- raster::extract(pixelID, data[, c("Longitude",  "Latitude")])
 
  # remove duplicates
  data <- data %>%
    distinct(uID, .keep_all = TRUE)%>%
    group_by(uID)%>%
    slice_head(n=1)%>%
    ungroup()
  
  return(data)
}

thin_occ_arbo <- function(data){
  trv <- 1:length(as.vector(template)) # unique pixel values
  trv = trv * !is.na(as.vector(template))
  pixelID = template
  values(pixelID) = trv

  # extract unique ID for each occurrence point
  data$uID <- raster::extract(pixelID, data[, c("Longitude",  "Latitude")])
 
  # remove duplicates
  data <- data %>%
    distinct(uID, Year, .keep_all = TRUE)%>%
    group_by(uID)%>%
    arrange(Year)%>%
    slice_head(n=1)%>%
    ungroup()
  
  return(data)
}

thin_bg_viral <- function(occ){
  trv <- 1:length(as.vector(template)) # unique pixel values
  trv = trv * !is.na(as.vector(template))
  pixelID = template
  values(pixelID) = trv

  bg <- data.frame(Longitude = runif(nrow(occ) * 100, min = -180, max = 180),
                 Latitude = runif(nrow(occ) * 100, min = -59, max = 75),
                 disease = "background")

  # trim to those on land and trim to size of occurrence data
  onland <- !is.na(raster::extract(template, bg[, 1:2]))
  bg = bg[onland, ]

  # thinning
  bg$uID <- raster::extract(pixelID, bg[, 1:2])
  # remove duplicates
  bg = bg[!duplicated(bg$uID), ]
  
  # trim to size of occurrence dataset
  bg = bg[1:nrow(occ), ]
  return(bg)
} 

thin_bg_arbo <- function(occ){
  trv <- 1:length(as.vector(template)) # unique pixel values
  trv = trv * !is.na(as.vector(template))
  pixelID = template
  values(pixelID) = trv
  
  bg <- data.frame(Longitude = runif(nrow(occ) * 100, min = -180, max = 180),
                   Latitude = runif(nrow(occ) * 100, min = -59, max = 75),
                   disease = sample(unique(occ$disease), 
                                    nrow(occ) * 100, 
                                    prob = as.numeric(table(occ$disease))[order(unique(occ$disease))] / nrow(occ),
                                    replace = T))
  
  # trim to those on land and trim to size of occurrence data
  onland <- !is.na(raster::extract(template, bg[, 1:2]))
  bg = bg[onland, ]
  
  # thinning
  bg$uID <- raster::extract(pixelID, bg[, 1:2])
  # remove duplicates
  bg = bg[!duplicated(bg$uID), ]
  
  # trim to size of occurrence dataset
  bg = bg[1:nrow(occ), ]
  return(bg)
} 
