# Setup: packages formatting and directories ===============

# Packages
if (!require("pacman")) install.packages("pacman")
if (!require("conflicted")) install.packages("conflicted")
# fmt: skip
pkgs =
  c(
    # data manipulation
    "tidyr","data.table","dplyr","countrycode","conflicted",
    "stringi","lubridate",
    # visualizations
    "ggplot2","patchwork","ggpubr","viridis","ggtext",
    "glue","extrafont","tidyterra",
    # spatial data
    "sf","mapview","exactextractr","raster","rnaturalearth",
    "terra",
    # modelling
    "randomForest","rsample","blockCV","boot",
    "cutpointr","Metrics","pROC","ROCR","matrixStats",
    "pdp",
    # programming
    "here","tictoc","doParallel"
  )

conflicts_prefer(
  dplyr::filter,
  dplyr::mutate,
  dplyr::select,
  raster::extract,
  ggplot2::margin,
  pROC::roc,
  data.table::shift,
  base::intersect,
  base::union,
  .quiet = TRUE
)

pacman::p_load(pkgs, character.only = T)

# Setting working directory
setwd(here())

# load functions
source("functions/fixNAs.R")
source("functions/plotRaster.R")
source("functions/thinning.R")
source("functions/adminExtract.R")
source("functions/buildRF.R")
source("functions/spatialAUC.R")
source("functions/spatialCV.R")
source("functions/calc_95_CIs.R")
source("functions/compare_maps.R")

# for parallel
ncpus = 18

# template raster
template = raster("data/template_raster.tif")

# load in Administrative unit rasters
admin <- stack(
  ad0 <- raster("data/admin_rasters/Admin0_5k_raster.tif"),
  ad1 <- raster("data/admin_rasters/Admin1_5k_raster.tif"),
  ad2 <- raster("data/admin_rasters/Admin2_5k_raster.tif")
)
names(admin) <- c("ad0", "ad1", "ad2")
rm("ad0", "ad1", "ad2")
