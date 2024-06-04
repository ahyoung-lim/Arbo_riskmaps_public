RF_mod <- function(data, cov_names, importance = TRUE) {
  formula <- as.formula(paste("PA ~", paste(cov_names, collapse = " + ")))
  model <- randomForest(formula, data = data, importance = importance)
  return(model)
}

# extract and save variable importance and partial dependence values for each RF model fit
vImp <- function(bootstraps, term_names, feature_names, pdp = TRUE) {
  mod <- RF_mod(bootstraps, term_names, importance = TRUE)
  vi <- importance(mod)
  
  if (pdp) { 
    
    pd <- list()
    for (i in seq_along(feature_names)) {
      tmp <- partial(mod, 
                     pred.var = feature_names[i], 
                     train = bootstraps,
                     prob = TRUE, 
                     which.class = 2,
                     plot = FALSE,
                     rug = TRUE,
                     parallel = TRUE, 
                     paropts = list(.packages = c("randomForest")))
      
      tmp <- tmp %>% 
        gather(., xvar, xvalue, feature_names[i], factor_key=TRUE)
      
      pd[[i]] <- tmp 
    }
  # stopCluster(cl) 
  
  pd_all <- rbindlist(pd)
  
  List <- list("vi" = vi,
               "pdp" = pd_all)  
  } else { 
    List <- vi
    
    }

  return(List)
}


# Variable importance plot
viPlot <- function(vi)  { 
  
  vi_list <-  lapply(vi, function(list) {
    data.frame(list[1])
  })
  

  vi_list <- lapply(vi_list, function(df) { 
    total_a <- sum(df$vi.MeanDecreaseAccuracy)
    total_g <- sum(df$vi.MeanDecreaseGini)
    df %>% 
      # normalization
      mutate(n_Acc = vi.MeanDecreaseAccuracy / total_a) %>%
      mutate(n_Gini = vi.MeanDecreaseGini / total_g)
    
  })
  
  vi_avg <- data.frame(
    Var = rownames(vi_list[[1]]),
    MeanDecreaseAccuracy = rowMeans(sapply(vi_list, function(df) df$n_Acc)), 
    MeanDecreaseGini = rowMeans(sapply(vi_list, function(df) df$n_Gini))
  )
  
  cov_mapping <- c(
    "GDP" = "GDP",
    "GDP_national" = "GDP (national)",
    "Urban" = "Urbanisation", 
    "Acc_walk" = "Travel time (health)", 
    "Acc_city" = "Travel time (cities)", 
    "Trmt" = "Treatment seeking", 
    "U5M" = "Child mortality", 
    "GovEff" = "Government effectiveness",
    "Physician" = "Physicians density", 
    "Tcold" = "Temperature of the coldest month", 
    "Tsuit" = "Temperature suitability", 
    "PRCP" = "Precipitation", 
    "NDVI" = "NDVI", 
    "DHI" = "DHI", 
    "Aegypti" = "Aedes aegypti", 
    "Albo" = "Aedes albopictus", 
    "Pop" = "Population", 
    "Hg" = "Haemagogus janthinomys", 
    "NHP" = "NHP", 
    "disease" = "Arbovirus disease"
)
  
  vip <- function(dat, measure) { 
    # 1_variable importance plot
    
    measure_title <- 
      ifelse(measure == "MeanDecreaseAccuracy", "Mean decrease in accuracy", "Mean decrease in node impurity"
      )
    
    for (i in 1:nrow(dat)) { 
      new_name <- cov_mapping[dat$Var[i]]
      dat$labels[i] <- ifelse(dat$Var[i] %in% names(cov_mapping), paste0(new_name), NA)
    }
    
    dat <- dat %>%
      mutate(labels2 = ifelse(grepl("Aedes|Haemagogus", labels), glue("*{labels}*"), paste0(labels)))
    
    ggplot(data=dat, aes(x=!!sym(measure), y = reorder(Var, !!sym(measure))))+
      geom_col(fill = "#5fa2ce", width = 0.6)+
      geom_richtext(aes(0, y = Var, label = labels2), 
                hjust = 0, nudge_x = 0.001, 
                fill= NA, label.colour = NA, colour = "grey15", size=5) + 
      ylab(NULL)+
      xlab(NULL)+
      labs(title = paste0(measure_title)) + 
      scale_x_continuous( expand = c(0, 0), position = "top" )+
      scale_y_discrete(expand = expansion(c(0, 0))) +
      theme(
            plot.title = element_text( face = "bold", size = 18), 
            plot.margin = margin(0.05, 0.1, 0.1, 0.05, "npc"),
            axis.title.x = element_text( size=15, face="bold", vjust=-1),
            axis.title.y = element_text( size=15, face="bold", vjust=2), 
            axis.text.x = element_text( size = 12),
            axis.text.y = element_blank(),
            axis.line.y.left = element_line(color = "black"),
            axis.ticks.length = unit(0, "mm"))+
      
      theme(panel.border = element_blank(),
            panel.spacing = unit(0.5, "lines"), 
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3)
      )
    
  }
  
  p1 <- vip(vi_avg, "MeanDecreaseAccuracy") 
    
  p2 <- vip(vi_avg, "MeanDecreaseGini") 
  
  # 2_Partial dependence plots =======================
  
  pdp_list <-  lapply(vi, function(list) {
    data.frame(list[2])%>%
      mutate(pdp.xvalue = ifelse(pdp.xvar == "DHI", paste0("DHI", pdp.xvalue), pdp.xvalue))
  })
  
  xvar_n <- pdp_list[[1]]%>% group_by(pdp.xvar)%>% tally()
  # Function to calculate row means with factor handling
  calculate_row_means <- function(df) {
    values <- df %>% 
      as.data.frame()%>%
      mutate_all(., ~coalesce(as.numeric(as.character(.)), NA))
    
    if (all(!is.na(values))) {
      return(rowMeans(df, na.rm = TRUE))
    } else {
      ind <- which(!complete.cases(values))
      result <- rep(NA, times = nrow(df))
      result[ind] <- df[ind]
      result[-ind] <- rowMeans(values[-ind, ], na.rm=TRUE)
      return(result)
      
    }
  }
  
  pdp_avg <- data.frame(
    Var = array(unlist(mapply(rep, xvar_n$pdp.xvar, xvar_n$n))),
    pdp.yhat = rowMeans(sapply(pdp_list, function(df) df$pdp.yhat)), 
    pdp.xvalue = calculate_row_means(as.matrix(sapply(pdp_list, function(df) df$pdp.xvalue))),
    pdp.ysd = rowSds(sapply(pdp_list, function(df) df$pdp.yhat))
  )
   pdp_avg <- pdp_avg %>%
     mutate(pdp.upper = pdp_avg$pdp.yhat+1.96*pdp_avg$pdp.ysd,
            pdp.lower = pdp_avg$pdp.yhat-1.96*pdp_avg$pdp.ysd)
  
   for (i in 1:nrow(pdp_avg)) { 
      new_name <- cov_mapping[pdp_avg$Var[i]]
      pdp_avg$labels[i] <- ifelse(pdp_avg$Var[i] %in% names(cov_mapping), paste0(new_name), NA)
    }
 
  if (all(is.numeric(pdp_avg$pdp.xvalue))) { 
    
    p3 <- ggplot(pdp_avg, aes(x=pdp.xvalue, y=pdp.yhat, group=Var))+
      geom_line(color = "red", aes(group=Var))+
    
      geom_ribbon(
                  aes(ymin = pdp.lower, ymax = pdp.upper, group = Var),
                  fill="grey70", alpha=0.5)+
      
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
      xlab(NULL)+
      ylab("Partial dependence values")+
      theme(axis.title.x = element_text( size=15, face="bold", vjust=-1))+ 
      theme(axis.title.y = element_text( size=15, face="bold", vjust=2))+  
      theme(axis.text = element_text( size=10), 
            axis.ticks = element_line(color = "grey80"))+
      
      theme(strip.background = element_blank(),
            strip.text = element_text( size=11, face="bold")
            # axis.line = element_line(colour = "grey70"),
      )+
      theme(panel.border = element_rect(colour = "grey70", fill = NA),
            panel.spacing = unit(0.5, "lines"), 
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "grey80", linetype = 3),
            panel.grid.minor = element_blank())+
      facet_wrap(labels~., scales="free", 
                 labeller = labeller(facet = function(labels) {
    ifelse(grepl("Aedes|Haemagogus", labels), parse(text = italic(labels)), labels)
  })
                 )
    
    
  } else { 
  dat_numeric <- pdp_avg %>%
    filter(!Var %in% c("DHI", "disease"))%>%
    mutate(pdp.xvalue = as.numeric(pdp.xvalue))%>%
    group_split(Var)
  
  dat_cat <- pdp_avg %>%
    filter(Var %in% c("DHI", "disease"))%>%
    mutate(pdp.xvalue = ifelse(Var == "DHI", as.factor(grep("DHI", pdp.xvalue)), pdp.xvalue))%>%
    group_split(Var)
  
  dat_cat[[1]]$pdp.xvalue <- factor(dat_cat[[1]]$pdp.xvalue)
  dat_cat[[2]]$pdp.xvalue <- factor(dat_cat[[2]]$pdp.xvalue, 
                                    levels = c("dengue", "zika", "chikungunya", "yf"))
  
  p_num <- function(dat_numeric){ 
    
    ggplot()+
    
    # Line plot for numeric values
    geom_line(data = dat_numeric,
              aes(x=as.numeric(pdp.xvalue), y=pdp.yhat, group =1), color = "red") +
    geom_ribbon(data = dat_numeric,
                aes(x=as.numeric(pdp.xvalue), y=pdp.yhat, ymin = pdp.lower, ymax = pdp.upper, group = Var),
                fill="grey70", alpha=0.5)+
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    xlab(NULL)+
    ylab(NULL)+
    ggtitle(ifelse(grepl("Aedes|Haemagogus", dat_numeric$Var),
                  expression(italic(dat_numeric$labels)),
                  dat_numeric$labels))+

    theme(plot.title = element_text( size=12, face = "bold"),
          axis.text = element_text( size=10), 
          axis.ticks = element_line(color = "grey80"))+
    theme(panel.border = element_rect(colour = "grey70", fill = NA),
          panel.spacing = unit(0.5, "lines"), 
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey80", linetype = 3),
          panel.grid.minor = element_blank())
  } 
  
  p_cat <- function(dat_cat) { 
    if(nrow(dat_cat) == 15) { 
      level_order <- as.character(1:15)
    } else { 
      level_order <- c("dengue", "chikungunya", "zika", "yf")  
      }
    
    ggplot()+  
     # Bar chart for categorical variables
      geom_bar(data = dat_cat,
               aes(x= factor(pdp.xvalue, level = level_order), y=pdp.yhat, group = 1),
               stat = "identity", position = "dodge") +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
      xlab(NULL)+
      ylab(NULL)+
      ggtitle(dat_cat$labels)+
      theme(plot.title = element_text(size=12, face="bold"),
            axis.text = element_text(size=8), 
            axis.ticks = element_line(color = "grey80"))+
      theme(panel.border = element_rect(colour = "grey70", fill = NA),
            panel.spacing = unit(0.5, "lines"), 
            panel.background = element_blank(),
            panel.grid.major = element_line(colour = "grey80", linetype = 3),
            panel.grid.minor = element_blank())
    
  
  } 
  
  all_plots <- c(lapply(dat_numeric, p_num), lapply(dat_cat, p_cat))
  
  
  p3 <- wrap_plots(all_plots, nrow = 4, ncol = 3) 
  }
   vip <- wrap_plots(list(p1, p2), ncol = 2, nrow = 1) +
     plot_annotation(title = "a", 
                     theme = theme(plot.title = element_text(size = 16, face= "bold", vjust=-1.5, hjust = -0.06)))
   pdp <- p3+
     plot_annotation(title = "b", 
                     theme = theme(plot.title = element_text(size = 16, face= "bold")))
   List <- list(vip, pdp)

    
  return(List)
}
  
