# Cross-sectional work on the classification of responses to adversity.

## Packages ####
library(dplyr) # dataframe managment
library(ggplot2) # ploting
library(ggrepel)
library(car)  # For VIF
library(caret)  # for createFolds
library(randomForest)
library(tidyr) # dataframe management
library(gridExtra) # plotinh
library(haven) # read sav data
library(missForest) # to impute data
library(olsrr) # influence statistic
library(rpart) # classification trees
library(rstanarm) # bayesian lm

## Set-up ####
setwd("~/Ecole/M1/Stage/Internship_repo/Classification/Cross-sectional/RYSE data")
load("~/Ecole/M1/Stage/Internship_repo/Classification/Cross-sectional/clean_rdata.RData")

RYSE_master_dataset <- read_sav("RYSE_master_dataset_08082022.sav")
df_CA <- RYSE_master_dataset[RYSE_master_dataset$Country==1,]
df_SA <- RYSE_master_dataset[RYSE_master_dataset$Country==2 
                             & RYSE_master_dataset$Site != 4,] # we exclude Zamdela


## Adjusted fit function ####
# Function for adjusted linear regression using influencial statistics
# Return : list containing :
# plot : plot with data points, standard regression, adjusted regression and bayesian regression
# influencers_indices : index of points from the original df that are considered influencial and removed to get the adjusted linear regression
# lm_adjusted : adjusted linear regression
# lm_ajusted_cred : adjusted linear regression using bayesian regression
# residuals_adjusted : residuals of the adjusted linear regression
adjusted_fit <- function(df,adversity,outcome,main="Adjusted and unadjusted linear regression",xlab="Adversity",ylab="Outcome"){
  rownames(df) <- NULL
  
  # Unadjusted linear model
  lm_unadjusted <- lm(as.formula(paste(outcome, "~", adversity)), data = df)
  
  # Identification of influencial points using Cook's D.
  used_data <- model.frame(lm_unadjusted)
  influencial_points <- which(cooks.distance(lm_unadjusted) > 4 / nrow(used_data))
  used_rows <- as.numeric(rownames(used_data))
  original_indices <- used_rows[influencial_points]
  
  # Cleaned data for the LM
  df_clean <- df[-original_indices, ]
  
  # Adjusted linear model
  lm_adjusted <- lm(as.formula(paste(outcome, "~", adversity)), data = df_clean)
  
  # Residuals of the adjusted linear model
  predicted_all <- predict(lm_adjusted, newdata = df)
  residuals_all <- df[[outcome]] - predicted_all
  
  # Bayesian glm for credibility intervals
  lm_adjusted_cred <- stan_glm(as.formula(paste(outcome, "~", adversity)), data = df,refresh=0)
  
  # Plot
  # Influencer binary value to plot
  influencer_flags <- rep(FALSE, nrow(df))
  influencer_flags[original_indices] <- TRUE
  
  df_plot <- df %>%
    mutate(influencer = influencer_flags)
  
  
  # Lines
  lines_df <- data.frame(
    intercept = c(coef(lm_unadjusted)[1],
                  coef(lm_adjusted)[1]),
    slope = c(coef(lm_unadjusted)[2],
              coef(lm_adjusted)[2]),
    model = c("Unadjusted", 
              "Adjusted")
  )
  
  
  # ggplot
  plot <- ggplot(df_plot, aes_string(x = adversity, y = outcome)) +
    
    # Influencers
    geom_point(aes(shape = influencer), size = 3) +
    
    # Regression lines
    geom_abline(data = lines_df,
                aes(intercept = intercept, slope = slope, color = model),
                size = 1)+
    
    # Legend
    scale_color_manual(values = c("Unadjusted" = "deepskyblue4",
                                  "Adjusted" = "aquamarine2"),
                       name="Model") +
    
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                       labels = c("Normal", "Influencer"),
                       name = "Point type") +
    
    # Title, xlab, ylab, theme
    labs(
      x = xlab,
      y = ylab,
      title = main
    ) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(size = 18)
    )
  
  # Resilient/Vulnerable anotations
  if(lines_df$slope[[1]] < 0){
    plot <-plot + geom_text(x=max(na.omit(df[[adversity]]))-5,y=max(na.omit(df[[outcome]]))-5,label="Resilient",alpha=0.2,color="grey",size=7) + geom_text(x=min(na.omit(df[[adversity]]))+5,y=min(na.omit(df[[outcome]]))+5,label="Non-resilient",alpha=0.2,color="grey",size=7)
  }
  else{
    plot <- plot + geom_text(x=min(na.omit(df[[adversity]]))+5,y=max(na.omit(df[[outcome]]))-5,label="Non-resilient",alpha=0.2,color="grey",size=7) + geom_text(x=max(na.omit(df[[adversity]]))-5,y=min(na.omit(df[[outcome]]))+5,label="Resilient",alpha=0.2,color="grey",size=7)
  }
  
  return(list(plot=plot,
              influencers_indices=original_indices,
              lm_adjusted=lm_adjusted,
              lm_adjusted_cred=lm_adjusted_cred,
              residuals_adjusted=residuals_all))
}


## Naive approach functions ####
# Functions to get the groups from the raw residuals
get_groups_raw_residuals <- function(residuals,is_resilience_positive=FALSE){
  res_list <- list()
  for(i in 1:length(residuals)){
    
    # Result depends of the slope sign
    if(is_resilience_positive){
      res_list[i] <- if(is.na(residuals[i])){NA}
      else if(residuals[i]>0){"resilient"}
      else if(residuals[i]==0){"average"}
      else{"non_resilient"}
    }
    else{
      res_list[i] <- if(is.na(residuals[i])){NA}
      else if(residuals[i]<0){"resilient"}
      else if(residuals[i]==0){"average"}
      else{"non_resilient"}
    }
  }
  return(res_list)
}

# Visualization function for raw residuals
visualization_raw_residuals <- function(df, adversity, outcome, adjusted_lm, groups, main = "Groups using naive method",xlab="",ylab="") {
  
  # Adjusted linear regression coefficient for the plot
  intercept <- coef(adjusted_lm)[1]
  slope     <- coef(adjusted_lm)[2]
  
  # Add groups to the temporary df to color the points
  df$group <- factor(groups, levels = c("resilient", "average", "non_resilient"))
  
  # Viz
  plot <- ggplot(df, aes(x = .data[[adversity]], y = .data[[outcome]], color = group)) +
    geom_point(shape=16,size=1.5) +
    geom_abline(intercept = intercept, slope = slope, color = "grey", linetype = "solid") +
    labs(
      x = if(xlab==""){adversity}else{xlab},
      y = if(xlab==""){outcome}else{ylab},
      title = main,
      color = "Group",
    ) +
    theme_minimal(base_size = 18) +
    theme(plot.title = element_text(size = 18))+
    scale_color_manual(values=c("resilient"="cadetblue4","non_resilient"="coral","average"="grey60"),labels = c("resilient" = "Resilient", "average" = "Average", "non_resilient" = "Non-resilient"))
  
  # Add resilient and non_resilient text
  if(slope < 0){
    plot <-plot + geom_text(x=max(na.omit(df[[adversity]]))-5,y=max(na.omit(df[[outcome]]))-5,label="Resilient",alpha=0.2,color="grey",size=7) + geom_text(x=min(na.omit(df[[adversity]]))+5,y=min(na.omit(df[[outcome]]))+5,label="Non-resilient",alpha=0.2,color="grey",size=7)
  }
  else{
    plot <- plot + geom_text(x=min(na.omit(df[[adversity]]))+5,y=max(na.omit(df[[outcome]]))-5,label="Non-resilient",alpha=0.2,color="grey",size=7) + geom_text(x=max(na.omit(df[[adversity]]))-5,y=min(na.omit(df[[outcome]]))+5,label="Resilient",alpha=0.2,color="grey",size=7)
  }
  
  return(plot)
}


## Confidence / Prediction / Credibility intervals approach  ####
get_credibility_intervals <- function(lm_adjusted_cred,newdata,adversity_string,lwr=0.025,upr=0.975){
  
  # Posterior linear prediction
  preds <- posterior_linpred(lm_adjusted_cred,
                             newdata = newdata[!is.na(newdata[[adversity_string]]), , drop = FALSE],
                             draws = 1000,
                             transform = TRUE)
  
  # Ensure preds is a matrix with rows = draws, cols = observations
  if (is.null(dim(preds))) {
    preds <- matrix(preds, nrow = 1000)
  }
  
  # Compute credible intervals per observation (apply over columns)
  intervals <- t(apply(preds, 2, quantile, probs = c(lwr, upr)))
  
  # Formating the result
  res <- data.frame(lwr=c(), upr=c())
  counter <- 1
  
  # Adding NA or computed value for every point of the dataset
  for(i in 1:nrow(newdata)){
    if(!is.na(newdata[i,1])){
      res[i,"lwr"] <- intervals[counter,1]
      res[i,"upr"] <- intervals[counter,2]
      counter <- counter + 1
    }
    else{
      res[i,"lwr"] <- NA
      res[i,"upr"] <- NA
    }
  }
  return(res)
}


# Function to return groups based on confidence/prediction/credibility intervals
get_groups_intervals <- function(actual, pred, is_resilience_positive = FALSE) {
  # Initialization of the result vector
  groups <- rep(NA_character_, length(actual))
  
  # Valid lines
  valid <- !is.na(actual) & !is.na(pred$lwr) & !is.na(pred$upr)
  
  # Classification with respect to the sign
  if (is_resilience_positive) {
    groups[valid & actual < pred$lwr] <- "non_resilient"
    groups[valid & actual > pred$upr] <- "resilient"
    groups[valid & actual >= pred$lwr & actual <= pred$upr] <- "average"
  } else {
    groups[valid & actual < pred$lwr] <- "resilient"
    groups[valid & actual > pred$upr] <- "non_resilient"
    groups[valid & actual >= pred$lwr & actual <= pred$upr] <- "average"
  }
  
  return(groups)
}

visualization_intervals <- function(df, adversity, outcome, adjusted_lm, preds, labels, 
                                    main = "Intervals", 
                                    colors = c("#deebf7", "#9ecae1","skyblue", "#6baed6", "#3182bd", "#08519c"),xlab="",ylab="") {
  # Coefficients of the linear regression
  intercept <- coef(adjusted_lm)[1]
  slope     <- coef(adjusted_lm)[2]
  
  # Base graph with points and linear regression 
  plot <- ggplot(df, aes(x = .data[[adversity]], y = .data[[outcome]])) +
    geom_point(shape = 16, size = 1.5) +
    geom_abline(intercept = intercept, slope = slope, color = "grey", linetype = "solid") +
    labs(
      x = if(xlab==""){adversity}else{xlab},
      y = if(xlab==""){outcome}else{ylab},
      title = main,
      fill = "Interval"
    ) +
    theme_minimal(base_size = 18) +
    theme(plot.title = element_text(size = 18))
  
  # Combine all intervals parsed in the function
  all_preds <- do.call(rbind, lapply(seq_along(preds), function(i) {
    pred <- preds[[i]]
    pred$label <- labels[i]
    pred
  }))
  all_preds$label <- factor(all_preds$label, levels = labels)
  
  # Add ribbons for every intervals
  plot <- plot +
    geom_ribbon(
      data = all_preds,
      aes(
        x = .data[[adversity]],
        ymin = lwr,
        ymax = upr,
        fill = label,
        group = label
      ),
      alpha = 0.4,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = colors)
  
  # Resilient / Vulnerable annotations
  if (slope < 0) {
    plot <- plot +
      geom_text(x = max(na.omit(df[[adversity]])) - 5,
                y = max(na.omit(df[[outcome]])) - 5,
                label = "Resilient", alpha = 0.2, color = "grey",size=7) +
      geom_text(x = min(na.omit(df[[adversity]])) + 5,
                y = min(na.omit(df[[outcome]])) + 5,
                label = "Non-resilient", alpha = 0.2, color = "grey",size=7)
  } else {
    plot <- plot +
      geom_text(x = min(na.omit(df[[adversity]])) + 5,
                y = max(na.omit(df[[outcome]])) - 5,
                label = "Non-resilient", alpha = 0.2, color = "grey",size=7) +
      geom_text(x = max(na.omit(df[[adversity]])) - 5,
                y = min(na.omit(df[[outcome]])) + 5,
                label = "Resilient", alpha = 0.2, color = "grey",size=7)
  }
  
  return(plot)
} 

## Quantile approach ####
classification_quantiles <- function(data_training, new_data, quantile_sub, quantile_sup, model, outcome_string, adversity_string) {
  residuals <- data_training[[outcome_string]] - (model$coefficients[[1]] + data_training[[adversity_string]] * model$coefficients[[2]])
  
  res_sub <- quantile(residuals, quantile_sub, na.rm = TRUE)
  res_sup <- quantile(residuals, quantile_sup, na.rm = TRUE)
  residuals_new <- new_data[[outcome_string]] - (model$coefficients[[1]] + new_data[[adversity_string]] * model$coefficients[[2]])
  is_resilience_positive <- model$coefficients[[2]] < 0
  
  res <- rep(NA, length(residuals_new))
  valid <- !is.na(residuals_new)
  res[valid] <- "average"
  
  if (is_resilience_positive) {
    res[valid & residuals_new <= res_sub] <- "non_resilient"
    res[valid & residuals_new >= res_sup] <- "resilient"
  } else {
    res[valid & residuals_new <= res_sub] <- "resilient"
    res[valid & residuals_new >= res_sup] <- "non_resilient"
  }
  
  return(res)
}


## SD-based approach ####
get_sd_in_bins <- function(data_training, lm, bins,multiplicator,outcome_string, adversity_string){
  df <- data_training
  residuals <- data_training[[outcome_string]] - (lm$coefficients[[1]] + data_training[[adversity_string]] * lm$coefficients[[2]])
  res_SD <- list()
  
  # Calculate the SD for each bin
  for (i in 1:(length(bins) - 1)) {
    in_bin <- df[[adversity_string]] >= bins[i] & df[[adversity_string]] < bins[i + 1]
    
    residuals_bin <- residuals[in_bin]
    if(length(residuals_bin)>0){
      res_SD[i] <- sd(residuals_bin, na.rm = TRUE)*multiplicator
    }
    else{
      res_SD[i] <- 0
    }
    
  }
  return(res_SD)
}

classification_sd <- function(new_data,bins,res_SD,model,outcome_string,adversity_string,sd_multiplicator){
  bin_labels <- rep(NA, nrow(new_data))  # To stock the bin of each line.
  
  # Find the bin corresponding to each individual
  for (i in 1:(length(bins) - 1)) {
    in_bin <- df[[adversity_string]] >= bins[i] & df[[adversity_string]] < bins[i + 1]
    
    bin_labels[in_bin] <- i #We save the i indice of the bin for each line that's in the bin
  }
  
  # Get the residual of each individual
  residuals <- new_data[[outcome_string]] - (model$coefficients[[1]] + new_data[[adversity_string]] * model$coefficients[[2]])
  
  # Flag each residual as resilient, average or non_resilient
  groups_sd <- rep(NA, nrow(new_data))
  
  #Fin the sign of the relation
  is_resilience_positive <- model$coefficients[[2]]<0
  
  for (i in seq_along(residuals)) {
    bin_i <- bin_labels[i]
    
    # Skip if bin or residual is NA
    if (is.na(bin_i) || is.na(residuals[i])) next
    
    # Recuperate the SD for the bin and the residual of the current point
    sd_i <- res_SD[[bin_i]]*sd_multiplicator
    res <- residuals[i]
    
    # Look at the value of the residuals with respect to the SD
    if (abs(res) <= sd_i) {
      groups_sd[i] <- "average"
    }
    # Bigger residual -> resilient
    else if(is_resilience_positive){
      if (res > sd_i) {
        groups_sd[i] <- "resilient"
      } else if (res < -sd_i) {
        groups_sd[i] <- "non_resilient"
      }
    }
    # Bigger residual -> non_resilient
    else{
      if (res > sd_i) {
        groups_sd[i] <- "non_resilient"
      } 
      else if (res < -sd_i) {
        groups_sd[i] <- "resilient"
      }
    }
  }
  return(groups_sd)
}

# Visualization function for the SD intervals
visualization_sd_intervals <- function(df,adversity,outcome,adjusted_lm,bins,res,names_sd,main="SD Intervals",xlab="",ylab=""){
  # Adjusted linear model coefficients
  intercept <- coef(adjusted_lm)[1]
  slope     <- coef(adjusted_lm)[2]
  
  # Base graph with points and regression line
  plot <- ggplot(df, aes(x = .data[[adversity]], y = .data[[outcome]])) +
    geom_point(shape=16,size=1.5) +
    geom_abline(intercept = intercept, slope = slope, color = "grey", linetype = "solid") +
    labs(
      x = if(xlab==""){adversity}else{xlab},
      y = if(xlab==""){outcome}else{ylab},
      title = main,
      fill = "Interval"
    ) +
    theme_minimal(base_size = 18) +
    theme(plot.title = element_text(size = 18))
  
  # Color the area for each result
  all_polygons <- data.frame()
  
  for (i in seq_along(res)){
    res_SD <- res[[i]]
    for(j in 1:(length(bins)-1)){
      point_j <- bins[[j]]
      point_j1 <- bins[[j+1]]
      fit_j <- intercept + slope * point_j
      fit_j1 <- intercept + slope * point_j1
      sd_j <- res_SD[[j]]
      
      
      df_area <- data.frame(
        x = c(point_j, point_j, point_j1, point_j1),
        y = c(fit_j - sd_j, fit_j + sd_j, fit_j1 + sd_j, fit_j1 - sd_j),
        group = paste0("poly_", i, "_", j),
        fill_factor = names_sd[[i]]
      )
      
      all_polygons <- rbind(all_polygons, df_area)
    }
  }
  
  # Add polygons to the main plot
  plot <- plot +
    geom_polygon(data = all_polygons, aes(x = x, y = y, group = group, fill = fill_factor), alpha = 0.3, color = NA) +
    scale_fill_manual(values = c("2SD" = "lightblue1", "1SD" = "skyblue2", "0.5SD" = "deepskyblue3"))
  
  
  # Add resilient and non_resilient text
  if(slope < 0){
    plot <-plot + geom_text(x=max(na.omit(df[[adversity]]))-5,y=max(na.omit(df[[outcome]]))-5,label="Resilient",alpha=0.2,color="grey",size=7) + geom_text(x=min(na.omit(df[[adversity]]))+5,y=min(na.omit(df[[outcome]]))+5,label="Non-resilient",alpha=0.2,color="grey",size=7)
  }
  else{
    plot <- plot + geom_text(x=min(na.omit(df[[adversity]]))+5,y=max(na.omit(df[[outcome]]))-5,label="Non-resilient",alpha=0.2,color="grey",size=7) + geom_text(x=max(na.omit(df[[adversity]]))-5,y=min(na.omit(df[[outcome]]))+5,label="Resilient",alpha=0.2,color="grey",size=7)
  }
  
  return(plot)
}

## K-means approach ####
res_kmeans <- function(data_training,lm){
  #Get the residuals of all individuals
  residuals <- data_training[[outcome_string]] - (lm$coefficients[[1]] + data_training[[adversity_string]] * lm$coefficients[[2]])
  
  # Create dataframe for clustering
  non_na_indices <- which(!is.na(residuals))
  clean_residuals <- residuals[non_na_indices]
  df <- data.frame(residual = clean_residuals)
  
  # Run K-means clustering with 3 centers
  clust <- kmeans(df, centers = 3)
  df$cluster <- as.factor(clust$cluster)
  
  # Calculate mean residual for each cluster
  cluster_means <- tapply(df$residual, df$cluster, mean)
  
  # Identify cluster corresponding to min, average, and max
  ordered_clusters <- names(sort(cluster_means))
  average_group <- ordered_clusters[2]
  if(lm$coefficients[[2]]<0){
    resilient_group <- ordered_clusters[3]
    non_resilient_group <- ordered_clusters[1]
  }
  else{
    resilient_group <- ordered_clusters[1]
    non_resilient_group <- ordered_clusters[3]
  }
  
  return(list(non_resilient_center=clust$centers[as.numeric(non_resilient_group)],average_center=clust$centers[as.numeric(average_group)],resilient_center=clust$centers[as.numeric(resilient_group)]))
  
}

classification_kmeans <- function(new_data, res_kmeans, model, outcome_string, adversity_string) {
  # Residuals of new data
  residuals <- new_data[[outcome_string]] - (model$coefficients[[1]] + new_data[[adversity_string]] * model$coefficients[[2]])
  
  # centers
  centers <- c(
    non_resilient = res_kmeans$non_resilient_center,
    average = res_kmeans$average_center,
    resilient = res_kmeans$resilient_center
  )
  
  # Initialization of res vector
  classification <- rep(NA_character_, length(residuals))
  
  # Valid data
  valid <- !is.na(residuals)
  
  # Pour chaque individu, on cherche le centre le plus proche
  dists <- sapply(centers, function(center) abs(residuals[valid] - center))
  closest <- apply(dists, 1, function(row) names(centers)[which.min(row)])
  
  classification[valid] <- closest
  
  return(classification)
}

visualization_groups <- function(df,adversity,outcome,adjusted_lm,groups,main="Clusturing results",xlab="",ylab=""){
  # Adjusted linear regression coefficient for the plot
  intercept <- coef(adjusted_lm)[1]
  slope     <- coef(adjusted_lm)[2]
  
  # Add groups to the temporary df to color the points
  df$group <- factor(groups, levels = c("resilient", "average", "non_resilient",NA))
  
  
  # Get the limit values of the residuals to visualize the separation of the groups
  df$residuals <- df[[outcome]] - (intercept + slope*df[[adversity]])
  # Resilient = positive residuals, non_resilient = negative residuals
  if(slope<0){
    resilient_limit <- min(df[df$group=="resilient",]$residuals, na.rm = TRUE)
    non_resilient_limit <- max(df[df$group=="non_resilient",]$residuals, na.rm = TRUE)
  }
  # Resilient = negative residuals, non_resilient = positive residuals
  else{
    resilient_limit <- max(df[df$group=="resilient",]$residuals, na.rm = TRUE)
    non_resilient_limit <- min(df[df$group=="non_resilient",]$residuals, na.rm = TRUE)
  }
  
  # Viz
  plot <- ggplot(df, aes(x = .data[[adversity]], y = .data[[outcome]], color = group)) +
    geom_point(shape=16,size=1.5) +
    geom_abline(intercept = intercept, slope = slope, color = "grey", linetype = "solid") +
    labs(
      x = if(xlab==""){adversity}else{xlab},
      y = if(xlab==""){outcome}else{ylab},
      title = main,
      color = "Group"
    ) +
    theme_minimal(base_size = 18) +
    theme(plot.title = element_text(size = 18))+
    scale_color_manual(values = c("resilient" = "cadetblue4", "average" = "grey60", "non_resilient" = "coral"),labels = c("resilient" = "Resilient", "average" = "Average", "non_resilient" = "Non-resilient"))
  
  # Add resilient and non_resilient text + lines for the seperation of the groups
  if(slope < 0){
    plot <-plot + geom_text(x=max(na.omit(df[[adversity]]))-5,y=max(na.omit(df[[outcome]]))-5,label="Resilient",alpha=0.2,color="grey",size=7) + geom_text(x=min(na.omit(df[[adversity]]))+5,y=min(na.omit(df[[outcome]]))+5,label="Non-resilient",alpha=0.2,color="grey",size=7)
    plot <- plot + 
      geom_abline(intercept = intercept+resilient_limit, slope = slope, color = "grey", linetype = "dashed")+
      geom_abline(intercept = intercept+non_resilient_limit, slope = slope, color = "grey", linetype = "dashed")
  }
  else{
    plot <- plot + geom_text(x=min(na.omit(df[[adversity]]))+5,y=max(na.omit(df[[outcome]]))-5,label="Non-resilient",alpha=0.2,color="grey",size=7) + geom_text(x=max(na.omit(df[[adversity]]))-5,y=min(na.omit(df[[outcome]]))+5,label="Resilient",alpha=0.2,color="grey",size=7)
    plot <- plot + 
      geom_abline(intercept = intercept+resilient_limit, slope = slope, color = "grey", linetype = "dashed")+
      geom_abline(intercept = intercept+non_resilient_limit, slope = slope, color = "grey", linetype = "dashed")
  }
  
  return(plot)
}
## Get all groups ####

# Function to get a dataframe with all of the grouping methods result and the dataframe with the sizes of each group for each method
get_all_groups <- function(df,adversity_string,outcome_string,bins,res,visualization=TRUE,xlab="Adversity",ylab="Outcome"){
  
  # Get the info about the lm
  lm_adjusted <- res$lm_adjusted
  lm_adjusted_cred <- res$lm_adjusted_cred
  resilience_sign <- lm_adjusted$coefficients[2]<0
  data_training <- res$lm_adjusted$model
  
  # Initialize the result data_frames : one with the grouping for each person and each method and one with the number of people in each group for each method
  outcome <- df[[outcome_string]]
  adversity <- df[[adversity_string]]
  residuals_new_data <- outcome - (lm_adjusted$coefficients[[1]] + adversity * lm_adjusted$coefficients[[2]])
  
  df_n_groups <- data.frame(resilient=c(),average=c(),non_resilient=c())
  df_result <- data.frame(residuals=residuals_new_data,adversity=adversity)
  
  # Raw residuals
  groups_raw <- get_groups_raw_residuals(residuals_new_data,is_resilience_positive=resilience_sign)
  df_n_groups <- rbind(df_n_groups,data.frame(resilient = sum(groups_raw=="resilient", na.rm=TRUE), average = sum(groups_raw=="average", na.rm=TRUE), non_resilient = sum(groups_raw=="non_resilient", na.rm=TRUE), row.names=c("raw")))
  df_result[["raw"]] <- groups_raw
  
  if(visualization){
    print(visualization_raw_residuals(df,adversity_string,outcome_string,lm_adjusted,groups_raw,xlab=xlab,ylab=ylab))
  }
  
  # Prediction and confidence intervals
  preds_conf <- list(
    as.data.frame(predict(lm_adjusted, newdata = df, interval = "prediction", level = 0.75)),
    as.data.frame(predict(lm_adjusted, newdata = df, interval = "prediction", level = 0.6)),
    as.data.frame(predict(lm_adjusted, newdata = df, interval = "prediction", level = 0.5)),
    as.data.frame(predict(lm_adjusted, newdata = df, interval = "confidence", level = 0.99)),
    as.data.frame(predict(lm_adjusted, newdata = df, interval = "confidence", level = 0.95))
  )
  names_conf <- list(
    "pred. residuals (75%)",
    "pred. residuals (60%)",
    "pred. residuals (50%)",
    "conf. residuals (99%)",
    "conf. residuals (95%)"
  )
  for(i in 1:length(preds_conf)){
    groups <- get_groups_intervals(outcome, preds_conf[[i]],is_resilience_positive=resilience_sign)
    df_n_groups <- rbind(df_n_groups,
                         data.frame(resilient = sum(groups=="resilient", na.rm=TRUE), average = sum(groups=="average", na.rm=TRUE), non_resilient = sum(groups=="non_resilient", na.rm=TRUE), row.names=c(names_conf[[i]])))
    df_result[[names_conf[[i]]]] <- groups
  }
  
  if(visualization){
    # Create grid of adversity scores for smooth prediction
    x_grid <- data.frame(adversity_var = seq(min(df[[adversity_string]], na.rm = TRUE),
                                             max(df[[adversity_string]], na.rm = TRUE),
                                             length.out = 200))

    # Rename column to match what the model expects
    colnames(x_grid) <- adversity_string

    # Generate predictions over the grid
    preds_conf2 <- list(
      as.data.frame(cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "prediction", level = 0.75)))),
      as.data.frame(cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "prediction", level = 0.6)))),
      as.data.frame(cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "prediction", level = 0.5)))),
      as.data.frame(cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "confidence", level = 0.99)))),
      as.data.frame(cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "confidence", level = 0.95))))
    )

    print(visualization_intervals(df=df,adversity=adversity_string,outcome=outcome_string,adjusted_lm =lm_adjusted,preds_conf2,names_conf,main="Confidence and prediction intervals",xlab=xlab,ylab=ylab))
  }
  
  
  # Credibility intervals
  preds_cred <- list(as.data.frame(cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df,adversity_string,lwr=0.0005,upr=0.9995))),
                     as.data.frame(cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df,adversity_string,lwr=0.005,upr=0.995))),
                     as.data.frame(cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df,adversity_string,lwr=0.025,upr=0.975))),
                     as.data.frame(cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df,adversity_string,lwr=0.05,upr=0.95))),
                     as.data.frame(cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df,adversity_string,lwr=0.125,upr=0.875))),
                     as.data.frame(cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df,adversity_string,lwr=0.25,upr=0.75)))
    
)
  names_cred <- list(
    "cred. 99.9%",
    "cred. 99%",
    "cred. 95%",
    "cred. 90%",
    "cred. 75%",
    "cred. 50%"
  )
  for(i in 1:length(preds_cred)){
    groups <- get_groups_intervals(outcome, preds_cred[[i]],is_resilience_positive=resilience_sign)
    df_n_groups <- rbind(df_n_groups,
                         data.frame(resilient = sum(groups=="resilient", na.rm=TRUE), average = sum(groups=="average", na.rm=TRUE), non_resilient = sum(groups=="non_resilient", na.rm=TRUE), row.names=c(names_cred[[i]])))
    df_result[[names_cred[[i]]]] <- groups
  }
  
  if(visualization){
    # Create grid of adversity scores for smooth prediction
    x_grid <- data.frame(adversity_var = seq(min(df[[adversity_string]], na.rm = TRUE),
                                             max(df[[adversity_string]], na.rm = TRUE),
                                             length.out = 200))
    
    # Rename column to match what the model expects
    colnames(x_grid) <- adversity_string
    
    # Generate predictions over the grid
    preds_cred2 <- list(as.data.frame(cbind(x_grid,get_credibility_intervals(lm_adjusted_cred,newdata=x_grid,adversity_string,lwr=0.0005,upr=0.9995))),
                       as.data.frame(cbind(x_grid,get_credibility_intervals(lm_adjusted_cred,newdata=x_grid,adversity_string,lwr=0.005,upr=0.995))),
                       as.data.frame(cbind(x_grid,get_credibility_intervals(lm_adjusted_cred,newdata=x_grid,adversity_string,lwr=0.025,upr=0.975))),
                       as.data.frame(cbind(x_grid,get_credibility_intervals(lm_adjusted_cred,newdata=x_grid,adversity_string,lwr=0.05,upr=0.95))),
                       as.data.frame(cbind(x_grid,get_credibility_intervals(lm_adjusted_cred,newdata=x_grid,adversity_string,lwr=0.125,upr=0.875))),
                       as.data.frame(cbind(x_grid,get_credibility_intervals(lm_adjusted_cred,newdata=x_grid,adversity_string,lwr=0.25,upr=0.75)))
    )
    print(visualization_intervals(df=df,adversity=adversity_string,outcome=outcome_string,adjusted_lm =lm_adjusted_cred,preds_cred2,names_cred,main="Credibility intervals",xlab=xlab,ylab=ylab))
  }
  
  # Quantiles
  list_quantile_sub <- list(0.05,0.1,0.15,0.2,0.25,0.3,0.35)
  list_quantile_sup <- list(0.95,0.9,0.85,0.8,0.75,0.7,0.65)
  names_quant <- list(
    "quantiles (5%)",
    "quantiles (10%)",
    "quantiles (15%)",
    "quantiles (20%)",
    "quantiles (25%)",
    "quantiles (30%)",
    "quantiles (35%)"
  )
  for(i in 1:length(list_quantile_sub)){
    groups <- classification_quantiles(data_training = data_training, new_data = df,quantile_sub = list_quantile_sub[[i]],quantile_sup = list_quantile_sup[[i]],model = lm_adjusted,outcome_string,adversity_string)
    df_n_groups <- rbind(df_n_groups,
                         data.frame(resilient = sum(groups=="resilient", na.rm=TRUE), average = sum(groups=="average", na.rm=TRUE), non_resilient = sum(groups=="non_resilient", na.rm=TRUE), row.names=c(names_quant[[i]])))
    df_result[[names_quant[[i]]]] <- groups
  }
  
  # Standard deviation
  list_sd_multiplicator <- list(2,1,0.5)
  names_sd <- list("2SD","1SD","0.5SD")
  res_sd <- list()
  
  for(i in 1:length(list_sd_multiplicator)){
    res_sd_i <- get_sd_in_bins(data_training,lm_adjusted,bins,list_sd_multiplicator[[i]],outcome_string,adversity_string)
    res_sd[[i]] <- res_sd_i
    groups <- classification_sd(df,bins,res_sd_i,lm_adjusted,outcome_string,adversity_string,list_sd_multiplicator[[i]])
    df_n_groups <- rbind(df_n_groups,
                         data.frame(resilient = sum(groups=="resilient", na.rm=TRUE), average = sum(groups=="average", na.rm=TRUE), non_resilient = sum(groups=="non_resilient", na.rm=TRUE), row.names=c(names_sd[[i]])))
    df_result[[names_sd[[i]]]] <- groups
  }
  if(visualization){
    print(visualization_sd_intervals(df,adversity=adversity_string,outcome=outcome_string,adjusted_lm=lm_adjusted,bins=bins,res=res_sd,names_sd=names_sd,main="SD Intervals",xlab=xlab,ylab=ylab))
  }
  
  # Kmeans (only residuals)
  res_kmeans <- res_kmeans(data_training,lm_adjusted)
  groups_kmeans <- classification_kmeans(df, res_kmeans, lm_adjusted, outcome_string, adversity_string)
  df_n_groups <- rbind(df_n_groups,
                       data.frame(resilient = sum(groups_kmeans=="resilient", na.rm=TRUE), average = sum(groups_kmeans=="average", na.rm=TRUE), non_resilient = sum(groups_kmeans=="non_resilient", na.rm=TRUE), row.names=c("Kmeans")))
  df_result[["Kmeans"]] <- groups_kmeans
  if(visualization){
    print(visualization_groups(df,adversity_string,outcome_string,lm_adjusted,groups_kmeans,main="Groups using k-means algorithm",xlab=xlab,ylab=ylab))
  }
  
  return(list(df_result=df_result,df_n_groups=df_n_groups))
}



## Commands for the illustrative example ####
# Variables used
residuals_vars <- c("T1_SES_total_SA","T1_WES_total", "T1_BDI_II","T1_edu_1a")
explication_vars <- c("T1_Sex","T1_Age",paste0("T1_CYRM_", 1:28),paste0("T1_PoNS_",1:8),paste0("T1_SF15_",1:15),paste0("T1_CPTS_",1:20),paste0("T1_FAS_",1:9),paste0("T1_BCE_",1:10))
#explication_vars <- c("T1_Sex","T1_Age","T1_CYRM28_total","T1_PoNS","T1_SF_14_PHC","T1_CPTS","T1_FAS","T1_BCE")

# Dataframe with SES or WES + pertinent variables
df_SAr <- df_SA[(!is.na(df_SA$T1_SES_total_SA) | !is.na(df_SA$T1_WES_total)) & !is.na(df_SA$T1_BDI_II),
                c(residuals_vars, explication_vars)]
dim(df_SAr) # 423 individuals

# Suppr individuals with over 30% of missing data
indices <- which(rowMeans(is.na(df_SAr)) > 0.3)
df_SAr <- df_SAr[-indices,]
dim(df_SAr)

# Impute data 
df_SAr_wtWESSES <- df_SAr[explication_vars]
df_SAr_wtWESSES <- as.data.frame(lapply(df_SAr_wtWESSES, function(x) as.numeric(as.character(x))))
df_SAr_wtWESSES.mf <- missForest::missForest(xmis = df_SAr_wtWESSES)
df_SAr[,explication_vars] <- df_SAr_wtWESSES.mf$ximp

# Creation of the engagement variable
df_SAr$T1_Engagement <- NA

# Use SES if we don't have WES
use_SES <- is.na(df_SAr$T1_WES_total)
df_SAr$T1_Engagement[use_SES] <- (df_SAr$T1_SES_total_SA[use_SES] - 33)/(165 - 33) * 100

# Use WES if we don't have SES
use_WES <- is.na(df_SAr$T1_SES_total_SA)
df_SAr$T1_Engagement[use_WES] <- (df_SAr$T1_WES_total[use_WES] - 9)/(63 - 9) * 100

# Then the remaning have both, we look is they are in grade 9->12
use_SES <- is.na(df_SAr$T1_Engagement) & df_SAr$T1_edu_1a %in% 9:12
df_SAr$T1_Engagement[use_SES] <- (df_SAr$T1_SES_total_SA[use_SES] - 33)/(165 - 33) * 100

# The remaning we take WES
use_WES <- is.na(df_SAr$T1_Engagement)
df_SAr$T1_Engagement[use_WES] <- (df_SAr$T1_WES_total[use_WES] - 9)/(63 - 9) * 100


df_SAr <- df_SAr[,c("T1_Engagement","T1_BDI_II",explication_vars)]

df <- df_SAr
adversity_string <- "T1_BDI_II"
outcome_string <- "T1_Engagement"
bins_BDI_II <- c(0,14,20,29,64)
bins <- bins_BDI_II
outcome <- df$T1_Engagement
res <- adjusted_fit(df_SAr,adversity="T1_BDI_II",outcome="T1_Engagement",main="Adjusted and unadjusted linear regression of Engagement (work or school) with BDI-II score",xlab="BDI-II score",ylab="Engagement")

#lm_adjusted <- res$lm_adjusted
#lm_adjusted_cred <- res$lm_adjusted_cred
#residuals <- res$residuals_adjusted

all_groups_BDI_Engagement <- get_all_groups(df,adversity_string,outcome_string,bins,res,visualization = TRUE)
df_result_BDI_Engagement <- all_groups_BDI_Engagement$df_result
df_n_groups_BDI_Engagement <- all_groups_BDI_Engagement$df_n_groups



## Classification functions ####
classification_metrics <- function(true_labels, predicted_labels,levels=c("resilient", "average", "non_resilient")) {
  # Convert to factors with same levels
  true_labels <- factor(true_labels, levels = levels)
  predicted_labels <- factor(predicted_labels, levels = levels)
  
  # Confusion Matrix
  cm <- caret::confusionMatrix(predicted_labels, true_labels)
  
  # Extract raw table
  cm_table <- cm$table
  n_classes <- length(levels)
  
  # Print
  print(cm_table)
  
  # Initialize results list
  results <- list()
  
  # Overall accuracy
  results$accuracy <- cm$overall["Accuracy"]
  
  # Balanced accuracy = mean(recall per class)
  recalls <- numeric(n_classes)
  precisions <- numeric(n_classes)
  f1s <- numeric(n_classes)
  fns <- numeric(n_classes)
  fps <- numeric(n_classes)
  supports <- numeric(n_classes)
  
  for (i in 1:n_classes) {
    class <- levels[i]
    
    TP <- cm_table[i, i]
    FN <- sum(cm_table[, i]) - TP
    FP <- sum(cm_table[i, ]) - TP
    TN <- sum(cm_table) - TP - FP - FN
    
    recall <- if ((TP + FN) == 0) NA else TP / (TP + FN)
    precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
    f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA else 2 * (precision * recall) / (precision + recall)
    
    recalls[i] <- recall
    precisions[i] <- precision
    f1s[i] <- f1
    fns[i] <- FN
    fps[i] <- FP
    supports[i] <- sum(cm_table[, i])
  }
  
  names(recalls) <- levels
  names(precisions) <- levels
  names(f1s) <- levels
  names(fns) <- paste0("FN_", levels)
  names(fps) <- paste0("FP_", levels)
  names(supports) <- levels
  
  results$macro_recall <- mean(recalls, na.rm = TRUE)
  results$macro_precision <- mean(precisions,na.rm=TRUE)
  results$f1_per_class <- f1s
  results$precision_per_class <- precisions
  results$recall_per_class <- recalls
  results$macro_f1 <- mean(f1s, na.rm = TRUE)
  results$false_negatives <- fns
  results$false_positives <- fps
  results$support <- supports
  
  return(results)
}

remove_high_vif <- function(df, predictors, threshold = 5) {
  # Garder seulement les variables non constantes
  non_constant <- predictors[sapply(df[, predictors, drop = FALSE], function(x) length(unique(x)) > 1)]
  
  # Retirer les variables avec trop peu de niveaux uniques (ex : 1 ou 2 niveaux sur une variable censée être continue)
  non_low_variance <- non_constant[sapply(df[, non_constant, drop = FALSE], function(x) length(unique(x)) > 3)]
  
  if (length(non_low_variance) < 2) {
    warning("Pas assez de variables valides pour calculer le VIF.")
    return(non_low_variance)
  }
  
  df$fake_outcome <- 1
  f <- as.formula(paste("fake_outcome ~", paste(non_low_variance, collapse = " + ")))
  lm_fit <- lm(f, data = df)
  vif_values <- car::vif(lm_fit)
  print("aa")
  print(vif_values)
  print("bbb")
  
  kept_predictors <- names(vif_values)[vif_values < threshold]
  return(kept_predictors)
}


estimation_classification <- function(df_train,df_test,adversity_string,outcome_string,bins,list_group_names,predictors){
  res <- data.frame()
  
  # Get the right predictors
  print(length(predictors))
  predictors_filtered <- remove_high_vif(df_train,predictors,threshold = 5)
  print(length(predictors_filtered))
  
  # Get the groupings for test and train
  res_data_train <- adjusted_fit(df_train,adversity_string,outcome_string)
  df_train_result <- get_all_groups(df_train, adversity_string, outcome_string,bins,res_data_train,visualization = FALSE)$df_result
  df_test_result <- get_all_groups(df_test, adversity_string, outcome_string,bins,res_data_train,visualization = FALSE)$df_result

  for(i in 1:length(list_group_names)){
    # Get the grouping method name
    group_name <- list_group_names[[i]]
    print(group_name)
    
    # We get the groups for the chosen grouping for test and train
    # Assign group labels
    df_train[["groups"]] <- as.factor(df_train_result[[group_name]])
    df_test[["groups"]] <- as.factor(df_test_result[[group_name]])
    
    
    # Train random forest
    formula <- as.formula(paste("groups ~", paste(predictors_filtered, collapse = " + ")))
    rf_model <- randomForest(formula, data = df_train, ntree = 500, mtry = floor(sqrt(length(predictors_filtered))), importance = TRUE)
    
    # Predict
    predictions <- predict(rf_model, newdata = df_test)
    
    # Metrics
    metrics <- classification_metrics(df_test[["groups"]],predictions)
    null_model <- max(metrics$support[["average"]],metrics$support[["resilient"]],metrics$support[["non_resilient"]]) / (metrics$support[["resilient"]]+metrics$support[["average"]]+metrics$support[["non_resilient"]])
    
    res <- rbind(res, data.frame(
      group_name = group_name,
      accuracy = metrics$accuracy[[1]],
      null_model = null_model,
      difference = metrics$accuracy[[1]]-null_model,
      macro_precision = metrics$macro_precision[[1]],
      macro_recall = metrics$macro_recall[[1]],
      macro_f1 = metrics$macro_f1[[1]],
      precision_resilient = metrics$precision_per_class[["resilient"]],
      recall_resilient = metrics$recall_per_class[["resilient"]],
      f1score_resilient = metrics$f1_per_class[["resilient"]],
      precision_average = metrics$precision_per_class[["average"]],
      recall_average = metrics$recall_per_class[["average"]],
      f1score_average = metrics$f1_per_class[["average"]],
      precision_non_resilient = metrics$precision_per_class[["non_resilient"]],
      recall_non_resilient = metrics$recall_per_class[["non_resilient"]],
      f1score_non_resilient = metrics$f1_per_class[["non_resilient"]],
      support_resilient = metrics$support[["resilient"]],
      support_average = metrics$support[["average"]],
      support_non_resilient = metrics$support[["non_resilient"]]
    ))
  }
  return(res)
}

estimation_classification_cv <- function(df, adversity_string, outcome_string, bins, list_group_names, predictors, k = 5, seed = 123){
  set.seed(seed)
  res <- data.frame()
  
  # Get the groupings once, from the full dataset
  res_data_train <- adjusted_fit(df, adversity_string, outcome_string)
  df_result <- get_all_groups(df, adversity_string, outcome_string, bins, res_data_train, visualization = FALSE)$df_result
  # Remove high VIF predictors
  predictors <- remove_high_vif(df, predictors, threshold = 5)
  
  for(i in 1:length(list_group_names)){
    group_name <- list_group_names[[i]]
    print(paste("Running CV for group:", group_name))
    
    df[["groups"]] <- as.factor(df_result[[group_name]])
    
    # Create folds (stratified)
    folds <- createFolds(df$groups, k = k, list = TRUE, returnTrain = FALSE)
    
    # Initialize metrics accumulators
    all_metrics <- list()
    
    for (fold_idx in seq_along(folds)) {
      test_indices <- folds[[fold_idx]]
      train_data <- df[-test_indices, ]
      test_data  <- df[test_indices, ]
      
      # Avoid crashing
      train_data$groups <- factor(train_data$groups, levels = levels(df$groups))

      # Fit random forest
      formula <- as.formula(paste("groups ~", paste(predictors, collapse = " + ")))
      rf_model <- randomForest(formula, data = train_data, ntree = 500, mtry = floor(sqrt(length(predictors))), importance = FALSE)
      
      # Predict
      preds <- predict(rf_model, newdata = test_data)
      
      # Evaluate
      metrics <- classification_metrics(test_data[["groups"]], preds)
      all_metrics[[fold_idx]] <- metrics
    }
    
    # Average metrics across folds
    average_metric <- function(metric_name, accessor = NULL) {
      values <- sapply(all_metrics, function(m) if (is.null(accessor)) m[[metric_name]][[1]] else m[[metric_name]][[accessor]])
      return(mean(values, na.rm = TRUE))
    }
    
    support_total <- average_metric("support", "resilient") + average_metric("support", "average") + average_metric("support", "non_resilient")
    null_model <- max(average_metric("support", "resilient"), average_metric("support", "average"), average_metric("support", "non_resilient")) / support_total
    
    res <- rbind(res, data.frame(
      group_name = group_name,
      accuracy = average_metric("accuracy"),
      null_model = null_model,
      difference = average_metric("accuracy") - null_model,
      macro_precision = average_metric("macro_precision"),
      macro_recall = average_metric("macro_recall"),
      macro_f1 = average_metric("macro_f1"),
      precision_resilient = average_metric("precision_per_class", "resilient"),
      recall_resilient = average_metric("recall_per_class", "resilient"),
      f1score_resilient = average_metric("f1_per_class", "resilient"),
      precision_average = average_metric("precision_per_class", "average"),
      recall_average = average_metric("recall_per_class", "average"),
      f1score_average = average_metric("f1_per_class", "average"),
      precision_non_resilient = average_metric("precision_per_class", "non_resilient"),
      recall_non_resilient = average_metric("recall_per_class", "non_resilient"),
      f1score_non_resilient = average_metric("f1_per_class", "non_resilient"),
      support_resilient = average_metric("support", "resilient"),
      support_average = average_metric("support", "average"),
      support_non_resilient = average_metric("support", "non_resilient")
    ))
  }
  
  return(res)
}

estimation_classification_binary_cv <- function(df, adversity_string, outcome_string, bins, list_group_names, predictors, k = 5, seed = 123){
  set.seed(seed)
  res <- data.frame()
  
  # Get groupings once from full data
  res_data_train <- adjusted_fit(df, adversity_string, outcome_string)
  df_result <- get_all_groups(df, adversity_string, outcome_string, bins, res_data_train, visualization = FALSE)$df_result
  
  # Remove high VIF predictors
  predictors <- remove_high_vif(df, predictors, threshold = 5)
  
  for(i in 1:length(list_group_names)){
    group_name <- list_group_names[[i]]
    print(paste("Running binary CV for group:", group_name))
    
    # Recode groups: average vs non-average
    df[["groups"]] <- ifelse(df_result[[group_name]] == "average", "average", "non_average")
    df[["groups"]] <- as.factor(df[["groups"]])
    
    # Stratified k-fold
    folds <- createFolds(df$groups, k = k, list = TRUE, returnTrain = FALSE)
    all_metrics <- list()
    
    for (fold_idx in seq_along(folds)) {
      test_indices <- folds[[fold_idx]]
      train_data <- df[-test_indices, ]
      test_data  <- df[test_indices, ]
      
      # Avoid crash <3
      train_data$groups <- factor(train_data$groups, levels = levels(df$groups))
      
      
      # Train RF
      formula <- as.formula(paste("groups ~", paste(predictors, collapse = " + ")))
      rf_model <- randomForest(formula, data = train_data, ntree = 1000, mtry = floor(sqrt(length(predictors))), importance = FALSE)
      
      # Predict
      preds <- predict(rf_model, newdata = test_data)
      
      # Evaluate
      metrics <- classification_metrics(test_data[["groups"]], preds,levels=c("average","non_average"))
      all_metrics[[fold_idx]] <- metrics
    }
    
    # Average metrics across folds
    average_metric <- function(metric_name, accessor = NULL) {
      values <- sapply(all_metrics, function(m) if (is.null(accessor)) m[[metric_name]][[1]] else m[[metric_name]][[accessor]])
      return(mean(values, na.rm = TRUE))
    }
    
    support_total <- average_metric("support", "average") + average_metric("support", "non_average")
    null_model <- max(average_metric("support", "average"), average_metric("support", "non_average")) / support_total
    
    res <- rbind(res, data.frame(
      group_name = group_name,
      accuracy = average_metric("accuracy"),
      null_model = null_model,
      difference = average_metric("accuracy") - null_model,
      macro_precision = average_metric("macro_precision"),
      macro_recall = average_metric("macro_recall"),
      macro_f1 = average_metric("macro_f1"),
      precision_average = average_metric("precision_per_class", "average"),
      recall_average = average_metric("recall_per_class", "average"),
      f1score_average = average_metric("f1_per_class", "average"),
      precision_non_average = average_metric("precision_per_class", "non_average"),
      recall_non_average = average_metric("recall_per_class", "non_average"),
      f1score_non_average = average_metric("f1_per_class", "non_average"),
      support_average = average_metric("support", "average"),
      support_non_average = average_metric("support", "non_average")
    ))
  }
  
  return(res)
}

## Visualization of classification results functions ####
comparison_accuracy_null_model_classifier <- function(df_perf){
  df_long <- df_perf %>%
    pivot_longer(cols = c(accuracy, null_model),
                 names_to = "metric",
                 values_to = "value") %>%
    mutate(metric = dplyr::recode(metric,
                                  accuracy = "Classifier",
                                  null_model = "Null model"))
  
  plot<- ggplot(df_long, aes(x = support_average, y = value, color = metric)) +
    geom_line(linewidth = 1) +
    geom_point() +
    labs(title = "Accuracy of the classifier vs accuracy of the null model as a function of the average group size",
         x = "Average group size",
         y = "Accuracy",
         color = "Metric",
         size=15) +
    xlim(0, max(df_perf$support_average)) +
    ylim(0, 1) +
    theme_minimal()+
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  return(plot)
}
visualization_recall_precision <- function(df_perf){
  criteria_index <- df_perf$support_average>=df_perf$support_resilient & df_perf$support_average>=df_perf$support_non_resilient & df_perf$recall_resilient>0.0001 & df_perf$precision_resilient>0.0001
  df_perf_class_comparison <- df_perf[criteria_index,]
  
  ggplot(df_perf_class_comparison,aes(x=1-recall_resilient,y=precision_resilient,label=group_name))+
    geom_point(shape=19,size=1.5)+
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2)+
    xlim(0,1.2)+
    ylim(0,1)+
    labs(title="Comparison of the grouping methods",
         x="1- recall of the resilient group",
         y= "precision of the resilient group")+
    theme_minimal()
}

## Commands for the classification on items RYSE ####
groups_to_test <- list("2SD","1SD","0.5SD","quantiles (5%)","quantiles (10%)","quantiles (15%)","quantiles (20%)","quantiles (25%)","quantiles (30%)","quantiles (35%)",
                       "pred. residuals (75%)","pred. residuals (60%)","pred. residuals (50%)","conf. residuals (99%)","conf. residuals (95%)",
                       "cred. 99.9%","cred. 99%","cred. 95%","cred. 90%","cred. 75%","cred. 50%","Kmeans")


# Reproductibility
set.seed(42)

#use 70% of dataset as training set and 30% as test set
sample <- sample(c(TRUE, FALSE), nrow(df_SAr), replace=TRUE, prob=c(0.7,0.3))
df_train <- df_SAr[sample, ]
df_test <- df_SAr[!sample, ]

df_perf_classification_tree42 <- estimation_classification(df_train,df_test,adversity_string,outcome_string,bins,groups_to_test,predictors = explication_vars)

df_perf_cv <- estimation_classification_cv(
  df = df_SAr,
  adversity_string = adversity_string,
  outcome_string = outcome_string,
  bins = bins,
  list_group_names = groups_to_test,
  predictors = explication_vars,
  k = 10
)

df_perf_binary_cv <- estimation_classification_binary_cv(
  df = df_SAr,
  adversity_string = adversity_string,
  outcome_string = outcome_string,
  bins = bins,
  list_group_names = groups_to_test,
  predictors = explication_vars,
  k = 10
)

# Test and train 42
comparison_accuracy_null_model_classifier(df_perf_classification_tree42)
visualization_recall_precision(df_perf_classification_tree42)
# CV
comparison_accuracy_null_model_classifier(df_perf_cv)
visualization_recall_precision(df_perf_cv)
# Binary CV
comparison_accuracy_null_model_classifier(df_perf_binary_cv)


## Commands for the classification on sums RYSE ####
# Variables used
residuals_vars <- c("T1_SES_total_SA","T1_WES_total", "T1_BDI_II","T1_edu_1a")
explication_vars_sum <- c("T1_Sex","T1_Age","T1_CYRM28_total","T1_PoNS","T1_SF_14_PHC","T1_CPTS","T1_FAS","T1_BCE")

# Dataframe with SES or WES + pertinent variables
df_SAr_sum <- df_SA[(!is.na(df_SA$T1_SES_total_SA) | !is.na(df_SA$T1_WES_total)) & !is.na(df_SA$T1_BDI_II),
                    c(residuals_vars, explication_vars_sum)]
dim(df_SAr_sum) # 423 individuals

# Suppr individuals with over 30% of missing data
indices <- which(rowMeans(is.na(df_SAr_sum)) > 0.3)
df_SAr_sum <- df_SAr_sum[-indices,]
dim(df_SAr_sum)

# Impute data 
df_SAr_wtWESSES_sum <- df_SAr_sum[explication_vars_sum]
df_SAr_wtWESSES_sum <- as.data.frame(lapply(df_SAr_wtWESSES_sum, function(x) as.numeric(as.character(x))))
df_SAr_wtWESSES.mf_sum <- missForest::missForest(xmis = df_SAr_wtWESSES_sum)
df_SAr_sum[,explication_vars_sum] <- df_SAr_wtWESSES.mf_sum$ximp

# Creation of the engagement variable
df_SAr_sum$T1_Engagement <- NA

# Use SES if we don't have WES
use_SES <- is.na(df_SAr_sum$T1_WES_total)
df_SAr_sum$T1_Engagement[use_SES] <- (df_SAr_sum$T1_SES_total_SA[use_SES] - 33)/(165 - 33) * 100

# Use WES if we don't have SES
use_WES <- is.na(df_SAr_sum$T1_SES_total_SA)
df_SAr_sum$T1_Engagement[use_WES] <- (df_SAr_sum$T1_WES_total[use_WES] - 9)/(63 - 9) * 100

# Then the remaning have both, we look is they are in grade 9->12
use_SES <- is.na(df_SAr_sum$T1_Engagement) & df_SAr_sum$T1_edu_1a %in% 9:12
df_SAr_sum$T1_Engagement[use_SES] <- (df_SAr_sum$T1_SES_total_SA[use_SES] - 33)/(165 - 33) * 100

# The remaning we take WES
use_WES <- is.na(df_SAr_sum$T1_Engagement)
df_SAr_sum$T1_Engagement[use_WES] <- (df_SAr_sum$T1_WES_total[use_WES] - 9)/(63 - 9) * 100


df_SAr_sum <- df_SAr_sum[,c("T1_Engagement","T1_BDI_II",explication_vars_sum)]

adversity_string <- "T1_BDI_II"
outcome_string <- "T1_Engagement"
bins_BDI_II <- c(0,14,20,29,64)
bins <- bins_BDI_II

groups_to_test <- list("2SD","1SD","0.5SD","quantiles (5%)","quantiles (10%)","quantiles (15%)","quantiles (20%)","quantiles (25%)","quantiles (30%)","quantiles (35%)",
                       "pred. residuals (75%)","pred. residuals (60%)","pred. residuals (50%)","conf. residuals (99%)","conf. residuals (95%)",
                       "cred. 99.9%","cred. 99%","cred. 95%","cred. 90%","cred. 75%","cred. 50%",
                       "Kmeans")


# Reproductibility
set.seed(42)

#use 70% of dataset as training set and 30% as test set
sample <- sample(c(TRUE, FALSE), nrow(df_SAr_sum), replace=TRUE, prob=c(0.7,0.3))
df_train_sum <- df_SAr_sum[sample, ]
df_test_sum <- df_SAr_sum[!sample, ]

df_perf_classification_tree42_sum <- estimation_classification(df_train_sum,df_test_sum,adversity_string,outcome_string,bins,groups_to_test,predictors = explication_vars_sum)

df_perf_cv_sum <- estimation_classification_cv(
  df = df_SAr_sum,
  adversity_string = adversity_string,
  outcome_string = outcome_string,
  bins = bins,
  list_group_names = groups_to_test,
  predictors = explication_vars_sum,
  k = 10
)

df_perf_binary_cv_sum <- estimation_classification_binary_cv(
  df = df_SAr_sum,
  adversity_string = adversity_string,
  outcome_string = outcome_string,
  bins = bins,
  list_group_names = groups_to_test,
  predictors = explication_vars_sum,
  k = 10
)

# Test and train 42
comparison_accuracy_null_model_classifier(df_perf_classification_tree42_sum)
visualization_recall_precision(df_perf_classification_tree42_sum)
# CV
comparison_accuracy_null_model_classifier(df_perf_cv_sum)
visualization_recall_precision(df_perf_cv_sum)
# Binary CV
comparison_accuracy_null_model_classifier(df_perf_binary_cv_sum)





## LORA Data Cleaning ####

# Get the data
df_LORA <- readRDS("C:/Users/garan/Documents/Ecole/M1/Stage/Internship_repo/Longitudinal/LORA/ds_forJan.rds")

# List all the interesting items/variables

# predictive
bfi <- paste0("bfi_",1:10)
cdrisk <- paste0("cdrisk_",1:25)
cerq <- paste0("cerq_",1:25)
cope <- paste0("cope_",1:28)
ctq <- paste0("ctq_",1:25)
fsozu <- paste0("fsozu_",1:14)
gpass <- paste0("gpass_",c(1,2,6,8,10:14,16:18,20,26,27))
gse <- paste0("gse_",1:10)
ielc <- paste0("ielc_",1:28)
le <- paste0("le_",1:27)
lotr <- paste0("lotr_",1:10)

# outcome/adversity
pss <- paste0("pss_",1:10)
dh <- paste0("dh_",c(1:28,44:58))
ghq <- paste0("ghq_", 1:28)

variables <- c("id","age","income","employment_status","gender",bfi,cdrisk,cerq,cope,ctq,fsozu,gpass,gse,ielc,le,lotr,pss,dh,ghq)

# Select variables and t=1
df_LORA$id <- as.character(df_LORA$id)
df_LORA$study <- as.character(df_LORA$study)
df_LORA <- df_LORA[df_LORA$t==1&!is.na(df_LORA$study)&!is.na(df_LORA$id),variables]
df_LORA <- df_LORA[1:1191,]
# Format the dataframe
df_LORA <- as.data.frame(lapply(df_LORA, function(x) as.numeric(as.character(x))))

# Replace negative values by NAs
for(variable in variables){
  df_LORA[[variable]] <- ifelse(df_LORA[[variable]]>=0, df_LORA[[variable]], NA)
}

# Recode the items if needed
recode <- function(df,variables_to_recode,min,max){
  for(variable in variables_to_recode){
    inverse_var <- paste0(variable, "_inverse")
    df[[inverse_var]] <- ifelse(df[[variable]] %in% min:max, max + min - df[[variable]], NA)
  }
  return(df)
}

  # BFI
bfi_to_recode <- paste0("bfi_",c(1,3,4,5,7))
df_LORA <- recode(df_LORA,bfi_to_recode,min=1,max=5)
bfi <- c(paste0("bfi_",c(2,6,8,9,10)),paste0("bfi_",c(1,3,4,5,7),"_inverse"))

  # CTQ
ctq_to_recode <- paste0("ctq_",c(2,5,7,12,17,23,25))
df_LORA <- recode(df_LORA,ctq_to_recode,min=1,max=5)
ctq <- c(paste0("ctq_",c(2,5,7,12,17,23,25),"_inverse"),paste0("ctq_",c(1,3,4,6,8,9,10,11,13,14,15,16,18,19,20,21,22,24)))

  # IELC
ielc_to_recode <- paste0("ielc_",c(3:5,10:13,15,22,26,27))
df_LORA <- recode(df_LORA,ielc_to_recode,min=0,max=1)
ielc <- c(paste0("ielc_",c(3:5,10:13,15,22,26,27),"_inverse"),paste0("ielc_",c(1,2,4:9,14,16:21,23:25,28)))

  # LOTR
lotr_to_recode <- paste0("lotr_",c(3,7,9))
df_LORA <- recode(df_LORA,lotr_to_recode,min=0,max=4)
lotr <- c(paste0("lotr_",c(3,7,9),"_inverse"),paste0("lotr_",c(1,4,10)))

  # PSS
pss_to_recode <- paste0("pss_", c(4, 5, 7, 8))
df_LORA <- recode(df_LORA,pss_to_recode,min=0,max=4)
pss <- c(paste0("pss_",c(1:3,6,9,10)),paste0("pss_",c(4, 5, 7, 8),"_inverse"))


## LORA Sums and subscales ####
# Able coping 1 and 19
cope_able <- paste0("cope_",c(1,19))
# Verl coping 3 and 8
cope_verl <- paste0("cope_",c(3,8))
# emU 5 and 15
cope_emu <- paste0("cope_",c(5,15))
# Ruck 6 and 16
cope_ruck <- paste0("cope_",c(6,16))
# poUm 12 and 17
cope_poum <- paste0("cope_",c(12,17))
# Hum 18 and 28
cope_hum <- paste0("cope_",c(18,28))
# akBe 2 and 7
cope_akbe <- paste0("cope_",c(2,7))
# AlDro 4 and 11
cope_aldro <- paste0("cope_",c(4,11))
# insUn 10 and 23
cope_insun <- paste0("cope_",c(10,23))
# ausE 9 and 21
cope_ause <- paste0("cope_",c(9,21))
# Plan 14 and 25
cope_plan <- paste0("cope_",c(14,25))
# Akze 20 and 24
cope_akze <- paste0("cope_",c(20,24))
# Sebe 13 and 26
cope_sebe <- paste0("cope_",c(13,26))
# Relo 22 and 27
cope_reli <- paste0("cope_",c(22,27))


# Build sums
df_LORA[["bfi_sum"]] <- rowSums(df_LORA[,bfi],na.rm = TRUE)
df_LORA[["cdrisk_sum"]] <- rowSums(df_LORA[,cdrisk],na.rm = TRUE)
df_LORA[["cerq_sum"]] <- rowSums(df_LORA[,cerq],na.rm = TRUE)
df_LORA[["able"]] <- rowSums(df_LORA[,cope_able],na.rm = TRUE)
df_LORA[["verl"]] <- rowSums(df_LORA[,cope_verl],na.rm = TRUE)
df_LORA[["emu"]] <- rowSums(df_LORA[,cope_emu],na.rm = TRUE)
df_LORA[["ruck"]] <- rowSums(df_LORA[,cope_ruck],na.rm = TRUE)
df_LORA[["poum"]] <- rowSums(df_LORA[,cope_poum],na.rm = TRUE)
df_LORA[["hum"]] <- rowSums(df_LORA[,cope_hum],na.rm = TRUE)
df_LORA[["akbe"]] <- rowSums(df_LORA[,cope_akbe],na.rm = TRUE)
df_LORA[["aldro"]] <- rowSums(df_LORA[,cope_aldro],na.rm = TRUE)
df_LORA[["insun"]] <- rowSums(df_LORA[,cope_insun],na.rm = TRUE)
df_LORA[["ause"]] <- rowSums(df_LORA[,cope_ause],na.rm = TRUE)
df_LORA[["plan"]] <- rowSums(df_LORA[,cope_plan],na.rm = TRUE)
df_LORA[["akze"]] <- rowSums(df_LORA[,cope_akze],na.rm = TRUE)
df_LORA[["sebe"]] <- rowSums(df_LORA[,cope_sebe],na.rm = TRUE)
df_LORA[["reli"]] <- rowSums(df_LORA[,cope_reli],na.rm = TRUE)
df_LORA[["ctq_sum"]] <- rowSums(df_LORA[,ctq],na.rm = TRUE)
df_LORA[["dh_sum"]] <- rowSums(df_LORA[,dh],na.rm = TRUE)
df_LORA[["fsozu_sum"]] <- rowSums(df_LORA[,fsozu],na.rm = TRUE)/14
df_LORA[["gpass_sum"]] <- rowSums(df_LORA[,gpass],na.rm = TRUE)
df_LORA[["gse_sum"]] <- rowSums(df_LORA[,gse],na.rm = TRUE)
df_LORA[["ghq_sum"]] <- rowSums(df_LORA[,ghq],,na.rm = TRUE)
df_LORA[["ielc_sum"]] <- rowSums(df_LORA[,ielc],na.rm = TRUE)
df_LORA[["le_sum"]] <- rowSums(df_LORA[,le],na.rm = TRUE)
df_LORA[["lotr_sum"]] <- rowSums(df_LORA[,lotr],na.rm = TRUE)
df_LORA[["pss_sum"]] <- rowSums(df_LORA[,pss],na.rm = TRUE)


# Select final variables
explication_vars_LORA <- c("bfi_sum","cdrisk_sum","cerq_sum",
                           "able","verl","emu","ruck","poum","hum",
                           "akbe","aldro","insun","ause","plan","akze",
                           "sebe","reli","ctq_sum","fsozu_sum",
                           "gpass_sum","gse_sum","ielc_sum","le_sum",
                           "lotr_sum","age","income","employment_status","gender")
adversity_vars_LORA <- c("dh_sum","pss_sum")
outcome_var_LORA <- c("ghq_sum")

df_LORA_final <- df_LORA[,c(explication_vars_LORA,adversity_vars_LORA,outcome_var_LORA)]


# Miss forest
LORAdh <- df_LORA_final[!is.na(df_LORA_final$dh_sum) & !is.na(df_LORA_final$ghq_sum)&!is.na(df_LORA_final$age)&!is.na(df_LORA_final$gender),]
LORApss <- df_LORA_final[!is.na(df_LORA_final$pss_sum) & !is.na(df_LORA_final$ghq_sum)&!is.na(df_LORA_final$age)&!is.na(df_LORA_final$gender),]

LORA_pss.mf <- missForest::missForest(xmis = LORApss)
LORA_pss.r <- LORA_pss.mf$ximp
LORA_dh.mf <- missForest::missForest(xmis = LORAdh)
LORA_dh.r <- LORA_dh.mf$ximp



# LORA - PSS 
# Residualization with PSS
df <- LORA_pss.r
adversity_string <- "pss_sum"
outcome_string <- "ghq_sum"
outcome <- df$ghq_sum
bins_pss <- c(0,14,26,40)
res_pss <- adjusted_fit(df,adversity=adversity_string,outcome=outcome_string,main="Adjusted and unadjusted linear regression for GHQ~PSS",xlab="PSS",ylab="GHQ")


all_groups_pss <- get_all_groups(df,adversity_string,outcome_string,bins,res_pss,visualization = TRUE)
df_result_pss <- all_groups_pss$df_result
df_n_groups_pss <- all_groups_pss$df_n_groups

# Classification test with PSS
groups_to_test <- list("quantiles (5%)","quantiles (10%)","quantiles (15%)","quantiles (20%)","quantiles (25%)","quantiles (30%)","quantiles (35%)",
                       "pred. residuals (75%)","pred. residuals (60%)","pred. residuals (50%)","conf. residuals (99%)","conf. residuals (95%)",
                       "cred. 99.9%","cred. 99%","cred. 95%","cred. 90%","cred. 75%","cred. 50%",
                       "2SD","1SD","0.5SD","Kmeans")


  # Reproductibility
set.seed(42)

#   #use 70% of dataset as training set and 30% as test set
# sample <- sample(c(TRUE, FALSE), nrow(LORA_pss.r), replace=TRUE, prob=c(0.7,0.3))
# df_train_pss <- LORA_pss.r[sample, ]
# df_test_pss <- LORA_pss.r[!sample, ]
# 
# df_perf_classification_tree42_pss <- estimation_classification(df_train_pss,df_test_pss,adversity_string,outcome_string,bins_pss,groups_to_test,predictors = explication_vars_LORA)

df_perf_cv_pss <- estimation_classification_cv(
  df = LORA_pss.r,
  adversity_string = adversity_string,
  outcome_string = outcome_string,
  bins = bins_pss,
  list_group_names = groups_to_test,
  predictors = explication_vars_LORA,
  k = 10
)

# df_perf_binary_cv_pss <- estimation_classification_binary_cv(
#   df = LORA_pss.r,
#   adversity_string = adversity_string,
#   outcome_string = outcome_string,
#   bins = bins_pss,
#   list_group_names = groups_to_test,
#   predictors = explication_vars_LORA,
#   k = 10
# )

#Visualizations
# Test and train 42
# comparison_accuracy_null_model_classifier(df_perf_classification_tree42_pss)
# visualization_recall_precision(df_perf_classification_tree42_pss)
# CV
comparison_accuracy_null_model_classifier(df_perf_cv_pss)
visualization_recall_precision(df_perf_cv_pss)
# Binary CV
# comparison_accuracy_null_model_classifier(df_perf_binary_cv_pss)

# LORA - DH 
# Residualization with DH
df <- LORA_dh.r
adversity_string <- "dh_sum"
outcome_string <- "ghq_sum"
outcome <- df$ghq_sum
bins_dh <- c(0,76,151,226)
res_dh <- adjusted_fit(df,adversity=adversity_string,outcome=outcome_string,main="Adjusted and unadjusted linear regression for GHQ~DH",xlab="DH",ylab="GHQ")


all_groups_dh <- get_all_groups(df,adversity_string,outcome_string,bins_dh,res_dh,visualization = TRUE)
df_result_dh <- all_groups_dh$df_result
df_n_groups_dh <- all_groups_dh$df_n_groups


#Classification test with DH
# Reproductibility
set.seed(42)
# 
# #use 70% of dataset as training set and 30% as test set
# sample <- sample(c(TRUE, FALSE), nrow(LORA_dh.r), replace=TRUE, prob=c(0.7,0.3))
# df_train_dh <- LORA_dh.r[sample, ]
# df_test_dh <- LORA_dh.r[!sample, ]
# 
# df_perf_classification_tree42_dh <- estimation_classification(df_train_dh,df_test_dh,adversity_string,outcome_string,bins_dh,groups_to_test,predictors = explication_vars_LORA)

df_perf_cv_dh <- estimation_classification_cv(
  df = LORA_dh.r,
  adversity_string = adversity_string,
  outcome_string = outcome_string,
  bins = bins_dh,
  list_group_names = groups_to_test,
  predictors = explication_vars_LORA,
  k = 10
)

# df_perf_binary_cv_dh <- estimation_classification_binary_cv(
#   df = LORA_dh.r,
#   adversity_string = adversity_string,
#   outcome_string = outcome_string,
#   bins = bins_dh,
#   list_group_names = groups_to_test,
#   predictors = explication_vars_LORA,
#   k = 10
# )

#Visualization
# Test and train 42
# comparison_accuracy_null_model_classifier(df_perf_classification_tree42_dh)
# visualization_recall_precision(df_perf_classification_tree42_dh)
# CV
comparison_accuracy_null_model_classifier(df_perf_cv_dh)
visualization_recall_precision(df_perf_cv_dh)
# Binary CV
# comparison_accuracy_null_model_classifier(df_perf_binary_cv_dh)
# 


## LORA Items ####
# Sums for PSS, DH and GHQ
df_LORA[["dh_sum"]] <- rowSums(df_LORA[,dh],na.rm = TRUE)
df_LORA[["pss_sum"]] <- rowSums(df_LORA[,pss],na.rm = TRUE)
df_LORA[["ghq_sum"]] <- rowSums(df_LORA[,ghq],na.rm = TRUE)

# Select final variables
explication_vars_items <- c("age","income","employment_status","gender",bfi,cdrisk,cerq,cope,ctq,fsozu,gpass,gse,ielc,le,lotr)
adversity_vars_LORA <- c("dh_sum","pss_sum")
outcome_var_LORA <- c("ghq_sum")

df_LORA_final <- df_LORA[,c(explication_vars_items,adversity_vars_LORA,outcome_var_LORA)]



# Miss forest
LORAdh <- df_LORA_final[!is.na(df_LORA_final$dh_sum) & !is.na(df_LORA_final$ghq_sum)&!is.na(df_LORA_final$age)&!is.na(df_LORA_final$gender),]
LORApss <- df_LORA_final[!is.na(df_LORA_final$pss_sum) & !is.na(df_LORA_final$ghq_sum)&!is.na(df_LORA_final$age)&!is.na(df_LORA_final$gender),]

# Suppr individuals with over 30% of missing data
indices <- which(rowMeans(is.na(LORAdh)) > 0.3)
LORAdh <- LORAdh[-indices,]
dim(LORAdh)

indices <- which(rowMeans(is.na(LORApss)) > 0.3)
LORApss <- LORApss[-indices,]
dim(LORApss)

LORA_pss.mf <- missForest::missForest(xmis = LORApss)
LORA_pss.r <- LORA_pss.mf$ximp
LORA_dh.mf <- missForest::missForest(xmis = LORAdh)
LORA_dh.r <- LORA_dh.mf$ximp


# LORA - PSS 
df <- LORA_pss.r
adversity_string <- "pss_sum"
outcome_string <- "ghq_sum"
outcome <- df$ghq_sum
bins_pss <- c(0,14,26,40)

# Classification test with PSS
groups_to_test <- list("quantiles (5%)","quantiles (10%)","quantiles (15%)","quantiles (20%)","quantiles (25%)","quantiles (30%)","quantiles (35%)",
                       "pred. residuals (75%)","pred. residuals (60%)","pred. residuals (50%)","conf. residuals (99%)","conf. residuals (95%)",
                       "cred. 99.9%","cred. 99%","cred. 95%","cred. 90%","cred. 75%","cred. 50%",
                       "2SD","1SD","0.5SD","Kmeans")


# Reproductibility
set.seed(42)

df_perf_cv_pss_items <- estimation_classification_cv(
  df = LORA_pss.r,
  adversity_string = adversity_string,
  outcome_string = outcome_string,
  bins = bins_pss,
  list_group_names = groups_to_test,
  predictors = explication_vars_items,
  k = 10
)

#Visualizations
comparison_accuracy_null_model_classifier(df_perf_cv_pss_items)
visualization_recall_precision(df_perf_cv_pss_items)

# LORA - DH 
# Residualization with DH
df <- LORA_dh.r
adversity_string <- "dh_sum"
outcome_string <- "ghq_sum"
outcome <- df$ghq_sum
bins_dh <- c(0,76,151,226)

#Classification test with DH
# Reproductibility
set.seed(42)

# CV
df_perf_cv_dh_items <- estimation_classification_cv(
  df = LORA_dh.r,
  adversity_string = adversity_string,
  outcome_string = outcome_string,
  bins = bins_dh,
  list_group_names = groups_to_test,
  predictors = explication_vars_items,
  k = 10
)


#Visualization
comparison_accuracy_null_model_classifier(df_perf_cv_dh_items)
visualization_recall_precision(df_perf_cv_dh_items)


## LORA DH Grouping and visualization ####
df_LORA_grouping_dh <- df_LORA[!is.na(df_LORA$dh_sum)&!is.na(df_LORA$ghq_sum),] 
dim(df_LORA_grouping_dh)
adversity_string <- "dh_sum"
outcome_string <- "ghq_sum"
bins_dh <- c(0,76,151,226)
res_dh <- adjusted_fit(df_LORA_grouping_dh,adversity=adversity_string,outcome=outcome_string,main="Adjusted and unadjusted linear regression for GHQ~DH",xlab="DH",ylab="GHQ")
all_groups_dh <- get_all_groups(df_LORA_grouping_dh,adversity_string,outcome_string,bins_dh,res_dh,visualization = TRUE)
df_result_dh <- all_groups_dh$df_result
df_n_groups_dh <- all_groups_dh$df_n_groups


## LORA PSS Grouping and visualization ####
df_LORA_grouping_pss <- df_LORA[!is.na(df_LORA$pss_sum)&!is.na(df_LORA$ghq_sum),] 
dim(df_LORA_grouping_pss)
adversity_string <- "pss_sum"
outcome_string <- "ghq_sum"
bins_pss <- c(0,14,26,40)
res_pss <- adjusted_fit(df_LORA_grouping_pss,adversity=adversity_string,outcome=outcome_string,main="Adjusted and unadjusted linear regression for GHQ~PSS",xlab="Adversity (PSS)",ylab="Outcome (GHQ)")
res_pss$plot
all_groups_pss <- get_all_groups(df_LORA_grouping_pss,adversity_string,outcome_string,bins,res_pss,visualization = TRUE,xlab="Adversity (PSS)",ylab="Outcome (GHQ)")
df_result_pss <- all_groups_pss$df_result
df_n_groups_pss <- all_groups_pss$df_n_groups

