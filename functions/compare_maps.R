mapExtract <- function(disease_data, disease_name, disease_abb) {
  new_map <- get(paste0(disease_abb, "_range_mask"))
  past_map <- get(paste0(disease_abb, "_ex"))

  # subset occ points by disease
  pres_subset <- disease_data[disease_data$disease == disease_name, ]

  names(new_map) <- paste0(disease_abb, "_oc_preds")
  names(past_map) <- paste0(disease_abb, "_oc_exs")

  # extract the map values
  new_vals <- adminExtract(pres_subset, new_map, admin)
  past_vals <- adminExtract(pres_subset, past_map, admin)

  # append extracted values to occ dataset
  pres_subset <- cbind(pres_subset, new_vals)
  pres_subset <- cbind(pres_subset, past_vals)

  # Fix NA values as much as possible (using buffers)
  pres_subset <- fixNAs(pres_subset, stack(new_map, past_map))

  return(pres_subset)
}


mapROC <- function(disease_data, disease_abb, output = "compare") {
  data <- disease_data
  cols <- c(paste0(disease_abb, "_oc_preds"), paste0(disease_abb, "_oc_exs"))

  dat_new <- data %>%
    dplyr::select(PA, cols[1]) %>%
    na.omit()
  dat_old <- data %>%
    dplyr::select(PA, cols[2]) %>%
    na.omit()

  cut_new <- cutpointr(
    dat_new,
    !!cols[1],
    PA,
    pos_class = "1",
    neg_class = "0",
    method = maximize_metric,
    metric = sum_sens_spec
  )$optimal_cutpoint

  cut_old <- cutpointr(
    dat_old,
    !!cols[2],
    PA,
    pos_class = "1",
    neg_class = "0",
    method = maximize_metric,
    metric = sum_sens_spec
  )$optimal_cutpoint
  cat("Optimal Threshold:", cut_new, " and ", cut_old, "\n")

  auc_new <- sprintf(
    "%.3f",
    round(Metrics::auc(dat_new$PA, dat_new[, cols[1]]), 3)
  )
  auc_old <- sprintf(
    "%.3f",
    round(Metrics::auc(dat_old$PA, dat_old[, cols[2]]), 3)
  )

  # auc_new <- round(Metrics::auc(ifelse(dat_new[,cols[1]] >= cut_new, 1, 0), dat_new$PA),3)
  # auc_old <- round(Metrics::auc(ifelse(dat_old[,cols[2]] >= cut_old, 1, 0), dat_old$PA),3)

  cat("AUC:", auc_new, " and ", auc_old, "\n")

  new_risk <- get(paste0(disease_abb, "_range_mask")) > cut_new
  past_risk <- get(paste0(disease_abb, "_ex")) > cut_old

  if (output == "new_only") {
    return(list(bin_map = new_risk, cut_new = cut_new))
  } else if (output == "compare") {
    plot_dtitle <- ifelse(
      disease_abb == "DEN",
      "a  Dengue",
      ifelse(
        disease_abb == "CHIK",
        'b  Chikungunya',
        ifelse(disease_abb == "ZIK", "c  Zika", "d  Yellow fever")
      )
    )
    plot_stitle <- paste0("New (", auc_new, "), Previous (", auc_old, ")")

    comparison_map <- plotFRaster(compare_maps(new_risk, past_risk)) +
      labs(title = plot_dtitle, subtitle = plot_stitle) +
      theme(text = element_text(family = "sans"))

    return(comparison_map)
  } else {
    return(past_risk)
  }
}

compare_maps <- function(new_map, ex_map) {
  compare_map <- new_map
  compare_map[new_map == 1 & ex_map == 0] <- "3" # New high
  compare_map[new_map == 0 & ex_map == 1] <- "2" # Ex high
  compare_map[new_map == 1 & ex_map == 1] <- "4" # Both high
  compare_map[new_map == 0 & ex_map == 0] <- "1" # Both low
  compare_map[is.na(new_map) & is.na(ex_map)] <- NA

  return(compare_map)
}
