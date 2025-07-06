# Cross-sectional work on the classification of responses to adversity.

## Packages ####
library(dplyr) # dataframe managment
library(ggplot2) # ploting
library(gridExtra) # plotinh
library(haven) # read sav data
library(missForest) # to impute data
library(olsrr) # influence statistic
library(rpart) # classification trees
library(rstanarm) # bayesian lm

## Set-up ####
setwd("~/Ecole/M1/Stage/Internship_repo/Cross-sectional/RYSE data")
load("~/Ecole/M1/Stage/Internship_repo/Cross-sectional/clean_rdata.RData")

RYSE_master_dataset <- read_sav("RYSE_master_dataset_08082022.sav")
df_CA <- RYSE_master_dataset[RYSE_master_dataset$Country==1,]
df_SA <- RYSE_master_dataset[RYSE_master_dataset$Country==2 
                             & RYSE_master_dataset$Site != 4,] # we exclude Zamdela

bins_BDI_II <- c(0,14,20,29,64)

## Adjusted fit function ####
# Function for adjusted linear regression using influencial statistics
# Return : list containing :
# plot : plot with data points, standard regression, adjusted regression and bayesian regression
# influencers_indices : index of points from the original df that are considered influencial and removed to get the adjusted linear regression
# lm_adjusted : adjusted linear regression
# lm_ajusted_cred : adjusted linear regression using bayesian regression
# residuals_adjusted : residuals of the adjusted linear regression
adjusted_fit <- function(df,adversity,outcome,main="Adjusted and unadjusted linear regression",xlab="Adversity",ylab="Outcome"){
  
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
  df_plot <- df %>%
    mutate(influencer = FALSE)
  df_plot$influencer[original_indices] <- TRUE
  
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
    geom_point(aes(shape = influencer), size = 2) +
    
    # Regression lines
    geom_abline(data = lines_df,
                aes(intercept = intercept, slope = slope, color = model),
                size = 0.8)+
    
    # Legend
    scale_color_manual(values = c("Unadjusted" = "deepskyblue3",
                                  "Adjusted" = "aquamarine3"),
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
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 10)
    )
  
  # Resilient/Vulnerable anotations
  if(lines_df$slope[[1]] < 0){
    plot <-plot + geom_text(x=max(na.omit(df[[adversity]]))-5,y=max(na.omit(df[[outcome]]))-5,label="Resilient",alpha=0.2,color="grey") + geom_text(x=min(na.omit(df[[adversity]]))+5,y=min(na.omit(df[[outcome]]))+5,label="Vulnerable",alpha=0.2,color="grey")
  }
  else{
    plot <- plot + geom_text(x=min(na.omit(df[[adversity]]))+5,y=max(na.omit(df[[outcome]]))-5,label="Vulnerable",alpha=0.2,color="grey") + geom_text(x=max(na.omit(df[[adversity]]))-5,y=min(na.omit(df[[outcome]]))+5,label="Resilient",alpha=0.2,color="grey")
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
      else{"vulnerable"}
    }
    else{
      res_list[i] <- if(is.na(residuals[i])){NA}
      else if(residuals[i]<0){"resilient"}
      else if(residuals[i]==0){"average"}
      else{"vulnerable"}
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
  df$group <- factor(groups, levels = c("resilient", "average", "vulnerable"))
  
  # Viz
  plot <- ggplot(df, aes(x = .data[[adversity]], y = .data[[outcome]], color = group)) +
    geom_point(shape=1,size=0.8) +
    geom_abline(intercept = intercept, slope = slope, color = "grey", linetype = "solid") +
    labs(
      x = if(xlab==""){adversity}else{xlab},
      y = if(xlab==""){outcome}else{ylab},
      title = main,
      color = "Group"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10))+
    scale_color_manual(values=c("resilient"="skyblue","vulnerable"="coral","average"="grey"))
  
  # Add resilient and vulnerable text
  if(slope < 0){
    plot <-plot + geom_text(x=max(na.omit(df[[adversity]]))-5,y=max(na.omit(df[[outcome]]))-5,label="Resilient",alpha=0.2,color="grey") + geom_text(x=min(na.omit(df[[adversity]]))+5,y=min(na.omit(df[[outcome]]))+5,label="Vulnerable",alpha=0.2,color="grey")
  }
  else{
    plot <- plot + geom_text(x=min(na.omit(df[[adversity]]))+5,y=max(na.omit(df[[outcome]]))-5,label="Vulnerable",alpha=0.2,color="grey") + geom_text(x=max(na.omit(df[[adversity]]))-5,y=min(na.omit(df[[outcome]]))+5,label="Resilient",alpha=0.2,color="grey")
  }
  
  return(plot)
}


## Confidence / Prediction / Credibility intervals approach  ####
get_credibility_intervals <- function(lm_adjusted_cred,newdata,lwr=0.025,upr=0.975){
  
  # Posterior linear prediction
  preds <- posterior_linpred(lm_adjusted_cred,
                             newdata = newdata[!is.na(newdata[,c(adversity_string)]),],
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
    groups[valid & actual < pred$lwr] <- "vulnerable"
    groups[valid & actual > pred$upr] <- "resilient"
    groups[valid & actual >= pred$lwr & actual <= pred$upr] <- "average"
  } else {
    groups[valid & actual < pred$lwr] <- "resilient"
    groups[valid & actual > pred$upr] <- "vulnerable"
    groups[valid & actual >= pred$lwr & actual <= pred$upr] <- "average"
  }
  
  return(groups)
}

visualization_intervals <- function(df, adversity, outcome, adjusted_lm, preds, labels, 
                                    main = "Intervals", 
                                    colors = c("#deebf7", "#9ecae1","skyblue2", "#6baed6", "#3182bd", "#08519c"),xlab="",ylab="") {
  # Coefficients of the linear regression
  intercept <- coef(adjusted_lm)[1]
  slope     <- coef(adjusted_lm)[2]
  
  # Base graph with points and linear regression 
  plot <- ggplot(df, aes(x = .data[[adversity]], y = .data[[outcome]])) +
    geom_point(shape = 1, size = 0.8) +
    geom_abline(intercept = intercept, slope = slope, color = "grey", linetype = "solid") +
    labs(
      x = if(xlab==""){adversity}else{xlab},
      y = if(xlab==""){outcome}else{ylab},
      title = main,
      fill = "Interval"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10))
  
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
                label = "Resilient", alpha = 0.2, color = "grey") +
      geom_text(x = min(na.omit(df[[adversity]])) + 5,
                y = min(na.omit(df[[outcome]])) + 5,
                label = "Vulnerable", alpha = 0.2, color = "grey")
  } else {
    plot <- plot +
      geom_text(x = min(na.omit(df[[adversity]])) + 5,
                y = max(na.omit(df[[outcome]])) - 5,
                label = "Vulnerable", alpha = 0.2, color = "grey") +
      geom_text(x = max(na.omit(df[[adversity]])) - 5,
                y = min(na.omit(df[[outcome]])) + 5,
                label = "Resilient", alpha = 0.2, color = "grey")
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
    res[valid & residuals_new <= res_sub] <- "vulnerable"
    res[valid & residuals_new >= res_sup] <- "resilient"
  } else {
    res[valid & residuals_new <= res_sub] <- "resilient"
    res[valid & residuals_new >= res_sup] <- "vulnerable"
  }
  
  return(res)
}


## SD-based approach ####
get_sd_in_bins <- function(data_training, lm, bins,outcome_string, adversity_string){
  df <- data_training
  residuals <- data_training[[outcome_string]] - (lm$coefficients[[1]] + data_training[[adversity_string]] * lm$coefficients[[2]])
  res_SD <- list()
  
  # Calculate the SD for each bin
  for (i in 1:(length(bins) - 1)) {
    in_bin <- df[[adversity_string]] >= bins[i] & df[[adversity_string]] < bins[i + 1]
    
    residuals_bin <- residuals[in_bin]
    res_SD[i] <- sd(residuals_bin, na.rm = TRUE)
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
  
  # Flag each residual as resilient, average or vulnerable
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
        groups_sd[i] <- "vulnerable"
      }
    }
    # Bigger residual -> vulnerable
    else{
      if (res > sd_i) {
        groups_sd[i] <- "vulnerable"
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
    geom_point(shape=1,size=0.8) +
    geom_abline(intercept = intercept, slope = slope, color = "grey", linetype = "solid") +
    labs(
      x = if(xlab==""){adversity}else{xlab},
      y = if(xlab==""){outcome}else{ylab},
      title = main,
      fill = "Interval"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10))
  
  # Color the area for each result
  all_polygons <- data.frame()
  
  for (i in seq_along(res)){
    res_SD <- res[[i]]$res_SD
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
    scale_fill_manual(values = c("2SD" = "lightblue", "1SD" = "skyblue2", "0.5SD" = "deepskyblue3"))
  
  
  # Add resilient and vulnerable text
  if(slope < 0){
    plot <-plot + geom_text(x=max(na.omit(df[[adversity]]))-5,y=max(na.omit(df[[outcome]]))-5,label="Resilient",alpha=0.2,color="grey") + geom_text(x=min(na.omit(df[[adversity]]))+5,y=min(na.omit(df[[outcome]]))+5,label="Vulnerable",alpha=0.2,color="grey")
  }
  else{
    plot <- plot + geom_text(x=min(na.omit(df[[adversity]]))+5,y=max(na.omit(df[[outcome]]))-5,label="Vulnerable",alpha=0.2,color="grey") + geom_text(x=max(na.omit(df[[adversity]]))-5,y=min(na.omit(df[[outcome]]))+5,label="Resilient",alpha=0.2,color="grey")
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
    vulnerable_group <- ordered_clusters[1]
  }
  else{
    resilient_group <- ordered_clusters[1]
    vulnerable_group <- ordered_clusters[3]
  }
  
  return(list(vulnerable_center=clust$centers[as.numeric(vulnerable_group)],average_center=clust$centers[as.numeric(average_group)],resilient_center=clust$centers[as.numeric(resilient_group)]))
  
}

classification_kmeans <- function(new_data, res_kmeans, model, outcome_string, adversity_string) {
  # Residuals of new data
  residuals <- new_data[[outcome_string]] - (model$coefficients[[1]] + new_data[[adversity_string]] * model$coefficients[[2]])
  
  # centers
  centers <- c(
    vulnerable = res_kmeans$vulnerable_center,
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
  df$group <- factor(groups, levels = c("resilient", "average", "vulnerable",NA))
  
  
  # Get the limit values of the residuals to visualize the separation of the groups
  df$residuals <- df[[outcome]] - (intercept + slope*df[[adversity]])
  # Resilient = positive residuals, vulnerable = negative residuals
  if(slope<0){
    resilient_limit <- min(df[df$group=="resilient",]$residuals, na.rm = TRUE)
    vulnerable_limit <- max(df[df$group=="vulnerable",]$residuals, na.rm = TRUE)
  }
  # Resilient = negative residuals, vulnerable = positive residuals
  else{
    resilient_limit <- max(df[df$group=="resilient",]$residuals, na.rm = TRUE)
    vulnerable_limit <- min(df[df$group=="vulnerable",]$residuals, na.rm = TRUE)
  }
  
  # Viz
  plot <- ggplot(df, aes(x = .data[[adversity]], y = .data[[outcome]], color = group)) +
    geom_point(shape=1,size=0.8) +
    geom_abline(intercept = intercept, slope = slope, color = "grey", linetype = "solid") +
    labs(
      x = if(xlab==""){adversity}else{xlab},
      y = if(xlab==""){outcome}else{ylab},
      title = main,
      color = "Group"
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 10))+
    scale_color_manual(values = c("resilient" = "skyblue", "average" = "grey", "vulnerable" = "coral"))
  
  # Add resilient and vulnerable text + lines for the seperation of the groups
  if(slope < 0){
    plot <-plot + geom_text(x=max(na.omit(df[[adversity]]))-5,y=max(na.omit(df[[outcome]]))-5,label="Resilient",alpha=0.2,color="grey") + geom_text(x=min(na.omit(df[[adversity]]))+5,y=min(na.omit(df[[outcome]]))+5,label="Vulnerable",alpha=0.2,color="grey")
    plot <- plot + 
      geom_abline(intercept = intercept+resilient_limit, slope = slope, color = "grey", linetype = "dashed")+
      geom_abline(intercept = intercept+vulnerable_limit, slope = slope, color = "grey", linetype = "dashed")
  }
  else{
    plot <- plot + geom_text(x=min(na.omit(df[[adversity]]))+5,y=max(na.omit(df[[outcome]]))-5,label="Vulnerable",alpha=0.2,color="grey") + geom_text(x=max(na.omit(df[[adversity]]))-5,y=min(na.omit(df[[outcome]]))+5,label="Resilient",alpha=0.2,color="grey")
    plot <- plot + 
      geom_abline(intercept = intercept+resilient_limit, slope = slope, color = "grey", linetype = "dashed")+
      geom_abline(intercept = intercept+vulnerable_limit, slope = slope, color = "grey", linetype = "dashed")
  }
  
  return(plot)
}
## Get all groups ####

# Function to get a dataframe with all of the grouping methods result and the dataframe with the sizes of each group for each method
get_all_groups <- function(df,adversity_string,outcome_string,bins,res,visualization=TRUE){
  
  # Get the info about the lm
  lm_adjusted <- res$lm_adjusted
  lm_adjusted_cred <- res$lm_adjusted_cred
  resilience_sign <- lm_adjusted$coefficients[2]<0
  data_training <- res$lm_adjusted$model
  
  # Initialize the result data_frames : one with the grouping for each person and each method and one with the number of people in each group for each method
  outcome <- df[[outcome_string]]
  adversity <- df[[adversity_string]]
  residuals_new_data <- outcome - (lm_adjusted$coefficients[[1]] + adversity * lm_adjusted$coefficients[[2]])
  
  df_n_groups <- data.frame(resilient=c(),average=c(),vulnerable=c())
  df_result <- data.frame(residuals=residuals_new_data,adversity=adversity)
  
  # Raw residuals
  groups_raw <- get_groups_raw_residuals(residuals_new_data,is_resilience_positive=resilience_sign)
  df_n_groups <- rbind(df_n_groups,data.frame(resilient = sum(groups_raw=="resilient", na.rm=TRUE), average = sum(groups_raw=="average", na.rm=TRUE), vulnerable = sum(groups_raw=="vulnerable", na.rm=TRUE), row.names=c("raw")))
  df_result[["raw"]] <- groups_raw
  
  if(visualization){
    print(visualization_raw_residuals(df,adversity_string,outcome_string,lm_adjusted,groups_raw,xlab="BDI-II score",ylab="Engagement"))
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
                         data.frame(resilient = sum(groups=="resilient", na.rm=TRUE), average = sum(groups=="average", na.rm=TRUE), vulnerable = sum(groups=="vulnerable", na.rm=TRUE), row.names=c(names_conf[[i]])))
    df_result[[names_conf[[i]]]] <- groups
  }
  
  if(visualization){
    # Create grid of adversity scores for smooth prediction
    x_grid <- data.frame(adversity_var = seq(min(adversity, na.rm = TRUE),
                                             max(adversity, na.rm = TRUE),
                                             length.out = 200))
    
    # Rename column to match what the model expects
    colnames(x_grid) <- adversity_string
    
    # Generate predictions over the grid
    preds_conf2 <- list(
      cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "prediction", level = 0.75))),
      cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "prediction", level = 0.6))),
      cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "prediction", level = 0.5))),
      cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "confidence", level = 0.99))),
      cbind(x_grid, as.data.frame(predict(lm_adjusted, newdata = x_grid, interval = "confidence", level = 0.95)))
    )
    
    print(visualization_intervals(df=df,adversity=adversity_string,outcome=outcome_string,adjusted_lm =lm_adjusted,preds_conf,names_conf,main="Confidence and prediction intervals",xlab="BDI-II score",ylab="Engagement"))
  }
  
  
  # Credibility intervals
  preds_cred <- list(
    cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df[c(adversity_string)],lwr=0.0005,upr=0.9995)),
    cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df[c(adversity_string)],lwr=0.005,upr=0.995)),
    cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df[c(adversity_string)],lwr=0.025,upr=0.975)),
    cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df[c(adversity_string)],lwr=0.05,upr=0.95)),
    cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df[c(adversity_string)],lwr=0.125,upr=0.875)),
    cbind(df[c(adversity_string)],get_credibility_intervals(lm_adjusted_cred,newdata=df[c(adversity_string)],lwr=0.25,upr=0.75)),
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
                         data.frame(resilient = sum(groups=="resilient", na.rm=TRUE), average = sum(groups=="average", na.rm=TRUE), vulnerable = sum(groups=="vulnerable", na.rm=TRUE), row.names=c(names_cred[[i]])))
    df_result[[names_cred[[i]]]] <- groups
  }
  
  if(visualization){
    print(visualization_intervals(df=df,adversity=adversity_string,outcome=outcome_string,adjusted_lm =lm_adjusted_cred,preds_cred,names_cred,main="Credibility intervals",xlab="BDI-II score",ylab="Engagement"))
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
                         data.frame(resilient = sum(groups=="resilient", na.rm=TRUE), average = sum(groups=="average", na.rm=TRUE), vulnerable = sum(groups=="vulnerable", na.rm=TRUE), row.names=c(names_quant[[i]])))
    df_result[[names_quant[[i]]]] <- groups
  }
  
  # Standard deviation
  list_sd_multiplicator <- list(2,1,0.5)
  names_sd <- list("2SD","1SD","0.5SD")
  res_sd <- list()
  
  for(i in 1:length(list_sd_multiplicator)){
    res_sd_i <- get_sd_in_bins(data_training,lm_adjusted,bins,outcome_string,adversity_string)
    res_sd[[i]] <- res_sd_i
    groups <- classification_sd(df,bins,res_sd_i,lm_adjusted,outcome_string,adversity_string,list_sd_multiplicator[[i]])
    df_n_groups <- rbind(df_n_groups,
                         data.frame(resilient = sum(groups=="resilient", na.rm=TRUE), average = sum(groups=="average", na.rm=TRUE), vulnerable = sum(groups=="vulnerable", na.rm=TRUE), row.names=c(names_sd[[i]])))
    df_result[[names_sd[[i]]]] <- groups
  }
  #if(visualization){
  #  print(visualization_sd_intervals(df,adversity=adversity_string,outcome=outcome_string,adjusted_lm=lm_adjusted,bins=bins,res=res_sd,names_sd=names_sd,main="SD Intervals",xlab="BDI-II score",ylab="Engagement"))
  #}
  
  # Kmeans (only residuals)
  res_kmeans <- res_kmeans(data_training,lm_adjusted)
  groups_kmeans <- classification_kmeans(df, res_kmeans, lm_adjusted, outcome_string, adversity_string)
  df_n_groups <- rbind(df_n_groups,
                       data.frame(resilient = sum(groups_kmeans=="resilient", na.rm=TRUE), average = sum(groups_kmeans=="average", na.rm=TRUE), vulnerable = sum(groups_kmeans=="vulnerable", na.rm=TRUE), row.names=c("Kmeans")))
  df_result[["Kmeans"]] <- groups_kmeans
  if(visualization){
    print(visualization_groups(df,adversity_string,outcome_string,lm_adjusted,groups_kmeans,main="Groups using k-means algorithm",xlab="BDI-II score",ylab="Engagement"))
  }
  
  return(list(df_result=df_result,df_n_groups=df_n_groups))
}



## Commands for the illustrative example ####
# Variables used
residuals_vars <- c("T1_SES_total_SA","T1_WES_total", "T1_BDI_II","T1_edu_1a")
explication_vars <- c("T1_Sex","T1_Age",paste0("T1_CYRM_", 1:28),paste0("T1_PoNS_",1:8),paste0("T1_SF15_",1:15),paste0("T1_CPTS_",1:20),paste0("T1_FAS_",1:9),paste0("T1_BCE_",1:10))
#explication_vars_sum <- c("T1_Sex","T1_Age","T1_CYRM28_total","T1_PoNS","T1_SF_14_PHC","T1_CPTS","T1_FAS","T1_BCE")

# Dataframe with SES or WES + pertinent variables
df_SAr <- df_SA[(!is.na(df_SA$T1_SES_total_SA)|!is.na(df_SA$T1_WES_total)) & !is.na(df_SA$T1_BDI_II)&!is.na(df_SA$T1_CYRM_10),c(residuals_vars,explication_vars,"Master_ID")]
dim(df_SAr) # 423 individuals

# Impute data 
df_SAr_wtWESSES <- df_SAr[explication_vars]
df_SAr_wtWESSES <- as.data.frame(lapply(df_SAr_wtWESSES, function(x) as.numeric(as.character(x))))
df_SAr_wtWESSES.mf <- missForest::missForest(xmis = df_SAr_wtWESSES)
df_SAr[,explication_vars] <- df_SAr_wtWESSES.mf$ximp

# Creation of the engagement variable
n <- nrow(df_SAr)
for(i in 1:n){
  if(is.na(df_SAr[i,"T1_WES_total"])){
    df_SAr[i,"T1_Engagement"] <- (df_SAr[i,"T1_SES_total_SA"]-33)/(165-33)*100
  }
  else if(is.na(df_SAr[i,"T1_SES_total_SA"])){
    df_SAr[i,"T1_Engagement"] <- (df_SAr[i,"T1_WES_total"]-9)/(63-9)*100
  }
  else{# If both are not NA
    if(df_SAr[i,"T1_edu_1a"] %in% 9:12){
      df_SAr[i,"T1_Engagement"] <- (df_SAr[i,"T1_SES_total_SA"]-33)/(165-33)*100
    }
    else{
      df_SAr[i,"T1_Engagement"] <- (df_SAr[i,"T1_SES_total_SA"]-33)/(165-33)*100
    }
  }
}

df <- df_SAr
adversity_string <- "T1_BDI_II"
outcome_string <- "T1_Engagement"
outcome <- df$T1_Engagement
bins <- bins_BDI_II
res <- adjusted_fit(df_SAr,adversity="T1_BDI_II",outcome="T1_Engagement",main="Adjusted and unadjusted linear regression of Engagement (work or school) with BDI-II score",xlab="BDI-II score",ylab="Engagement")

#lm_adjusted <- res$lm_adjusted
#lm_adjusted_cred <- res$lm_adjusted_cred
#residuals <- res$residuals_adjusted

all_groups_BDI_Engagement <- get_all_groups(df,adversity_string,outcome_string,bins,res,visualization = TRUE)
df_result_BDI_Engagement <- all_groups_BDI_Engagement$df_result
df_n_groups_BDI_Engagement <- all_groups_BDI_Engagement$df_n_groups



## Classification functions ####
classification_metrics <- function(true_labels, predicted_labels) {
  # Convert to factors with same levels
  levels <- c("resilient", "average", "vulnerable")
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

estimation_classification <- function(df_train,df_test,adversity_string,outcome_string,bins,list_group_names,predictors=c("T1_Sex", "T1_Age", paste0("T1_CYRM_", 1:28))){
  set.seed(1) # For reproductibility
  res <- data.frame()
  
  # Get the groupings for test and train
  res_data_train <- adjusted_fit(df_train,adversity_string,outcome_string)
  df_train_result <- get_all_groups(df_train, adversity_string, outcome_string,bins,res_data_train,visualization = FALSE)$df_result
  df_test_result <- get_all_groups(df_test, adversity_string, outcome_string,bins,res_data_train,visualization = FALSE)$df_result

  for(i in 1:length(list_group_names)){
    # Get the grouping method name
    group_name <- list_group_names[[i]]
    print(group_name)
    
    # We get the groups for the chosen grouping for test and train
    df_train[["groups"]] <- df_train_result[[group_name]]
    df_test[["groups"]] <- df_test_result[[group_name]]
    
    # We build the classification tree from the train data
    formula <- as.formula(paste("groups ~", paste(predictors, collapse = " + ")))
    arbre <- rpart(formula, data = df_train, method = "class")
    
    # Get the prediction on the test data
    predictions <- predict(arbre, newdata=df_test, type = "class")
    
    # Metrics
    metrics <- classification_metrics(df_test[["groups"]],predictions)
    null_model <- max(metrics$support[["average"]],metrics$support[["resilient"]],metrics$support[["vulnerable"]]) / (metrics$support[["resilient"]]+metrics$support[["average"]]+metrics$support[["vulnerable"]])
    
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
      precision_vulnerable = metrics$precision_per_class[["vulnerable"]],
      recall_vulnerable = metrics$recall_per_class[["vulnerable"]],
      f1score_vulnerable = metrics$f1_per_class[["vulnerable"]],
      support_resilient = metrics$support[["resilient"]],
      support_average = metrics$support[["average"]],
      support_vulnerable = metrics$support[["vulnerable"]]
    ))
  }
  return(res)
}


## Call of the classification ####
groups_to_test <- list("quantiles (5%)",
                       "quantiles (10%)",
                       "quantiles (15%)",
                       "quantiles (20%)",
                       "quantiles (25%)",
                       "quantiles (30%)",
                       "quantiles (35%)",
                       "pred. residuals (75%)",
                       "pred. residuals (60%)",
                       "pred. residuals (50%)",
                       "conf. residuals (99%)",
                       "conf. residuals (95%)",
                       "cred. 99.9%",
                       "cred. 99%",
                       "cred. 95%",
                       "cred. 90%",
                       "cred. 75%",
                       "cred. 50%",
                       "2SD",
                       "1SD",
                       "0.5SD",
                       "Kmeans")

# Reproductibility
set.seed(1)

#use 70% of dataset as training set and 30% as test set
sample <- sample(c(TRUE, FALSE), nrow(df_SAr), replace=TRUE, prob=c(0.7,0.3))
df_train <- df_SAr[sample, ]
df_test <- df_SAr[!sample, ]

df_perf_classification_tree1 <- estimation_classification(df_train,df_test,adversity_string,outcome_string,bins,groups_to_test)

set.seed(42)

#use 70% of dataset as training set and 30% as test set
sample <- sample(c(TRUE, FALSE), nrow(df_SAr), replace=TRUE, prob=c(0.7,0.3))
df_train <- df_SAr[sample, ]
df_test <- df_SAr[!sample, ]

df_perf_classification_tree42 <- estimation_classification(df_train,df_test,adversity_string,outcome_string,bins,groups_to_test)


# Selection of the best grouping methods 
# Criterias :
# The average group size is at least 1/3 of the whole dataset -> >129/3=
# The model as to predict (rightfully or wrongfully) resilience -> recall >0
# We want good precision and in second a good recall for the resilient group.
criteria_index <- df_perf_classification_tree42$support_average>=43&df_perf_classification_tree42$recall_resilient>0
df_perf_class_comparison <- df_perf_classification_tree[criteria_index,]

ggplot(df_perf_class_comparison,aes(x=1-recall_resilient,y=precision_resilient,label=group_name))+
  geom_point(shape=19,size=1.5)+
  geom_text(hjust=-0.1, vjust=0,size=3)+
  xlim(0,1)+
  ylim(0,1)+
  labs(title="Comparison of the grouping methods",
       x="1- recall of the resilient group",
       y= "precision of the resilient group")+
  theme_minimal()

# Visualization of the evolution of the accuracy and nul-model accuracy
df_long <- df_perf_classification_tree %>%
  pivot_longer(cols = c(accuracy, null_model),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = dplyr::recode(metric,
                         accuracy = "Accuracy du classifieur",
                         null_model = "Accuracy du modèle nul"))

ggplot(df_long, aes(x = support_average, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(title = "Accuracy du modèle nul et du classifieur en fonction de la taille du groupe normal",
       x = "Taille du groupe normal",
       y = "Performance",
       color = "Metrique",
       size=20) +
  xlim(0, 129) +
  ylim(0, 1) +
  theme_minimal()+
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))


