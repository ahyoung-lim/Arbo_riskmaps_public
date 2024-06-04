rowCIs <- function(matrix, weights, R) { 
  
  num_subsets <- 500
  # Calculate the number of rows in each subset
  subset_size <- ceiling(nrow(matrix) / num_subsets)
  
  weighted_bootstrap <- function(data, weights, R) { 
    # returns an output matrix of weighted average where each row represents each pixel
    estimates <- matrix(NA, ncol = R, nrow = nrow(data))  
    set.seed(1239873)
    for (i in seq_len(R)) {
      # columns are different models X1 to X100
      chosen_columns <- sample(names(data), size = ncol(data), replace = TRUE, prob = weights)
      estimates[,i] <- rowWeightedMeans(as.matrix(data[, chosen_columns]))  # weighted average of each pixel

          }
    
    # find confidence interval for each pixel
    CI <- rowQuantiles(estimates, probs=c(.025,  .975) )
    
    return(list( CI = CI))
  }
  
  # Function for applying weighted_bootstrap for a subset
  weighted_bootstrap_subset <- function(start_row, end_row, df) {
    subset_df <- df[start_row:end_row, ]
    wboots <- weighted_bootstrap(subset_df, weights, R)
    return(wboots)
  }
  
  subset_bootstrap <- foreach(i = 1:num_subsets,  
                              .packages = c("matrixStats")) %dopar% { 
                                
                                start_row <- (i - 1) * subset_size + 1
                                end_row <- min(i * subset_size, nrow(matrix))
                                
                                out <-  weighted_bootstrap_subset(start_row, end_row, matrix)$CI
                                
                                return(out)
                              }
  
  # Merge the results
  allCI <- do.call(rbind, subset_bootstrap)

  return(allCI)
} 