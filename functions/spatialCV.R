
create_blocks <- function(dat, nfolds) { 
  sf <- sf::st_as_sf(dat, coords = c("Longitude", "Latitude"), crs = 4326)
  sf$PA <- as.integer(as.character(sf$PA)) 

  sb <- cv_spatial(x = sf,
                  column = "PA", # the response column (binary or multi-class)
                  k =  nfolds,  # number of folds
                  size = 500000, # size of the blocks in metres
                  # rows_cols = c(5, 100),
                  # extend = 1,
                  seed = 13456,
                  hexagon = FALSE, # use square blocks
                  selection = "random", # random blocks-to-fold
                  iteration = 50, # find evenly dispersed folds
                  biomod2 = FALSE) # also create folds for biomod2
  
  block <- sb$blocks[, c("folds", "block_id", "geometry")]
  
  sf$folds_ids <- sb$folds_ids
  # intersect and extract block id (only one overlapping block)
  sf$block_id <- apply(st_intersects(block, sf, sparse = FALSE), 2, 
   function(col) { 
      intersecting_blocks <-
        block[which(col), ]$block_id
      if (length(intersecting_blocks) == 1) {
         return(intersecting_blocks)  
      } else {
         return(NA)  
      }
   }
  )
   
   # intersect and extract block id (more than two overlapping blocks)
   choose_block <- sf[is.na(sf$block_id), ] %>% st_drop_geometry()
   
   if (nrow(choose_block)>0) { 
    print("choosing block ids")
     intersect <- apply(st_intersects(block, sf[is.na(sf$block_id),], sparse = FALSE), 2, 
   function(col) { 
       intersecting_blocks <-  block[which(col), ]$block_id

   return(intersecting_blocks)
   }
  ) 
   intersect_df <- as.data.frame(matrix(unlist(intersect), nrow = nrow(choose_block), byrow = TRUE))

  # Combine the original data frame with the intersections data frame
  choose_block <- cbind(choose_block, intersect_df)
  choose_block <- choose_block %>%
    left_join(st_drop_geometry(block), by = c("V1" = "block_id")) %>%
    left_join(st_drop_geometry(block), by = c("V2" = "block_id")) %>%
    mutate(same_fold = folds.x == folds.y)%>%
    mutate(block_id = ifelse(!is.na(same_fold) & folds.x == folds_ids, paste0(V1), 
                      ifelse(!is.na(same_fold) & folds.y == folds_ids, paste0(V2), paste0(V1))))

  sf$block_id[is.na(sf$block_id)] <- choose_block$block_id
  sf$block_id <- as.factor(as.character(sf$block_id))
  dat$folds_id <- sf$folds_id
  dat$block_id <- sf$block_id
  bList <- list("dat" = dat, # dataframe with fold id and block id
                  "sb" = sb, # cv_spatial object 
                  "block" = block # block geometry
  )
   } else { 
     sf$block_id <- as.factor(as.character(sf$block_id))
     dat$folds_id <- sf$folds_id
     dat$block_id <- sf$block_id
     bList <- list("dat" = dat, # dataframe
                     "sb" = sb, # cv_spatial object 
                     "block" = block # block geometry
     )
     
     }

  return(bList)
}

block_CV <- function(bList, nfolds, cov_names){ 
    
    dat <- bList[[1]]
    sb <- bList[[2]]
  
    OOB_prediction <- data.frame()
    testSetCopy <- data.frame()
    AUC=c()

  for (k in 1:nfolds) { 

    train_ind <- sb$folds_list[[k]][[1]] # [[k]] fold id [[1]] 80% [[2]] 20%
    trainSet <- dat[train_ind, ]
    testSet <- dat[-train_ind,]

    # train the model (80%)
    mod <- RF_mod(trainSet, cov_names, importance=FALSE)

    # test the model on held-out dataset (20%)
    OOB_pred <- as.data.frame(predict(mod, testSet, type="prob")[,2])
    names(OOB_pred)[1] <- "OOB_prediction"
    OOB_prediction <- rbind(OOB_prediction, OOB_pred) #save OOB prediction outputs
    
    testSetCopy <- rbind(testSetCopy, as.data.frame(testSet[,c("Longitude", "Latitude", "PA", "block_id", "folds_id")])) 
    AUC[k]=Metrics::auc(testSetCopy$PA, OOB_prediction[,1])
  
    cat(k, "\rcompleted")
    }
  
  OOB <- cbind(OOB_prediction, testSetCopy[, c(1:5)])
  
  CV_list <- list("AUC" = AUC, 
                  "OOB" = OOB)
  return(CV_list)
}  

CV_data <- function(bList, nfolds){

  trainSet_ls <- list()
  testSet_ls <- list()

  dat <- bList[[1]]
  sb <- bList[[2]]

  for (k in 1:nfolds) { 
  
  train_ind <- sb$folds_list[[k]][[1]] # [[k]] fold id [[1]] 80% [[2]] 20%
  trainSet <- dat[train_ind, ]
  testSet <- dat[-train_ind,]
  
  trainSet_ls[[k]] <- trainSet
  testSet_ls[[k]] <- testSet
  }

  return(List = list("trainSet" = trainSet_ls, 
                     "testSet" = testSet_ls))
}

# for CV and weighted predictions
fit_mod <- function(CV_train, CV_test, cov_names, pred_data_parallel, disease_name, pred=TRUE, oob=TRUE) {
  
  if (pred) { 
  set.seed(1239023)  
  # train the model on nfolds-1/nfolds
  mod <- RF_mod(CV_train, cov_names, importance = FALSE)
  
  testSetCopy <- as.data.frame(CV_test[,c("Longitude", "Latitude", "PA", "disease", "block_id", "folds_id")]) 
  
  # test the model on held-out dataset 
  OOB_pred <- as.data.frame(predict(mod, CV_test, type="prob")[,2])
  names(OOB_pred)[1] <- "pred"
  
  OOB <- cbind(OOB_pred, testSetCopy)
  
  # save the model performance 
  AUC = Metrics::auc(OOB$PA, OOB$pred)
 
  # make predictions to new data
  mod_prediction <- predict(mod, pred_data_parallel, type = "prob")
  
  if (oob){
  List <- list("OOB" = OOB,
               "pred" = mod_prediction[,2],
               "AUC" = AUC)
  } else { 
  List <- list("pred" = mod_prediction[,2],
                "AUC" = AUC)
           }
           
  return(List)
  
  } else { 
    # train the model on nfolds-1/nfolds
    mod <- RF_mod(CV_train, cov_names)
    
    testSetCopy <- as.data.frame(CV_test[,c("Longitude", "Latitude", "PA", "disease", "block_id", "folds_id")]) 
    # test the model on held-out dataset 
    OOB_pred <- as.data.frame(predict(mod, CV_test, type="prob")[,2])
    names(OOB_pred)[1] <- "pred"
    
    OOB <- cbind(OOB_pred, testSetCopy)
    
    # save the model performance 
    AUC = Metrics::auc(OOB$PA, OOB$pred)
     
    List <- list("OOB" = OOB,
                 "AUC" = AUC)
   return(List)
  
  } 
 
}

fit_mod_surv <- function(CV_train, CV_test, cov_names, pred_data_parallel, pred=TRUE, oob=TRUE) {
  if (pred) { 
    # train the model on nfolds-1/nfolds
    mod <- RF_mod(CV_train, cov_names)
    
    testSetCopy <- as.data.frame(CV_test[,c("Longitude", "Latitude", "PA", "block_id", "folds_id")]) 
    
    # test the model on held-out dataset 
    OOB_pred <- as.data.frame(predict(mod, CV_test, type="prob")[,2])
    names(OOB_pred)[1] <- "pred"
    
    OOB <- cbind(OOB_pred, testSetCopy)
    
    # save the model performance 
    AUC = Metrics::auc(OOB$PA, OOB$pred)
    
    mod_prediction <- predict(mod, pred_data_parallel, type = "prob")
    
    if (oob){
      List <- list("OOB" = OOB,
                   "pred" = mod_prediction[,2],
                   "AUC" = AUC)
    } else { 
      
      List <- list("pred" = mod_prediction[,2],
                   "AUC" = AUC)
    }
    
    return(List)
    
  } else { 
    # train the model on nfolds-1/nfolds
    mod <- RF_mod(CV_train, cov_names)
    
    testSetCopy <- as.data.frame(CV_test[,c("Longitude", "Latitude", "PA", "block_id", "folds_id")]) 
    # test the model on held-out dataset 
    OOB_pred <- as.data.frame(predict(mod, CV_test, type="prob")[,2])
    names(OOB_pred)[1] <- "pred"
    
    OOB <- cbind(OOB_pred, testSetCopy)
    
    # save the model performance 
    AUC = Metrics::auc(OOB$PA, OOB$pred)
    
    List <- list("OOB" = OOB,
                 "AUC" = AUC)
    return(List)
    
  } 
  
}


OOBout <- function(CV_result) { 
  # Access the results
  OOB_list <- do.call(rbind, lapply(CV_result, function(result) result$OOB))

  }

Predout <- function(CV_result, pred=TRUE) {
  
  # to calculate the weighted mean of predictions
  AUC_list <- sapply(CV_result, function(result) result$AUC)
  
  if(pred) {
  pred_list <- do.call(cbind, lapply(CV_result, function(result) result$pred))
  # when returns predictions, bring back NAs and estimate the weighted average, median and IQR
  pred.data.large = pred_data_large(pred_list, AUC_list, sum = TRUE)
    # Merge the lists into a single list
    List <- list(
      pred = pred.data.large,
      AUC = AUC_list
    )
  } else { 
    
    List <- list(
      AUC = AUC_list
    )
    }
  
  return(List)
}


pred_data_large <- function(boots_df, AUC_list, sum = FALSE){ 
  template_df <- data.frame(matrix(NA, nrow = nrow(pred.data), 
                                       ncol = ncol(boots_df)))  
  template_df[-NAindex, ] <- boots_df
  
  
  if (sum) { 
    template_df$median <- NA
    template_df$IQR <- NA
    template_df$wmean <- NA
  
    # weighted average of predictions according to AUCs
    auc = unlist(AUC_list)
    weights <- auc / sum(auc)
  
    template_df[-NAindex, "median"] <- rowMedians(as.matrix(boots_df), na.rm = TRUE)
    template_df[-NAindex, "IQR"] <- rowIQR(boots_df)$IQR
    template_df[-NAindex, "wmean"] <- rowSums(as.matrix(boots_df) %*% diag(weights))
    
    
  }
  
  return(template_df)
  
}


rowIQR <- function(matrix) { 
  
  num_subsets <- 100
  # Calculate the number of rows in each subset
  subset_size <- ceiling(nrow(matrix) / num_subsets)
  
  # Create a list to store the results
  subset_iqr <- vector("list", length = num_subsets)
  
  # Function to calculate rowIQRs for a subset
  calculate_iqr_subset <- function(start_row, end_row, df) {
    subset_df <- df[start_row:end_row, ]
    IQR <- matrixStats::rowIQRs(as.matrix(subset_df), na.rm=TRUE)
    return(IQR)
  }
  
  subset_iqr <- foreach(i = 1:num_subsets, .combine = c) %dopar% { 
    
    start_row <- (i - 1) * subset_size + 1
    end_row <- min(i * subset_size, nrow(matrix))
    
    calculate_iqr_subset(start_row, end_row, matrix)
  }
  
  
  # Merge the results
  allIQRs <- data.frame(IQR = unlist(subset_iqr))
  
  return(allIQRs)
} 


