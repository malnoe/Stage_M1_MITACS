classification_naive(df_SAr,adversity_string,outcome_string,lm_adjusted)
classification_prediction_intervals(new_data,outcome_string,confidence,lm)
classification_confidence_intervals(new_data,outcome_string,confidence,lm)
classification_credibility_intervals(new_data,outcome_string,confidence,res$lm_adjusted_cred)
classification_quantiles(df_training, new_df, quantile_sub, quantile_sup, lm, outcome_string, adversity_string)
res_sd <- get_sd_in_bins(data_training, lm,bins,outcome_string,adversity_string)
classification_sd(new_data,bins,res_sd,lm,outcome_string,adversity_string,sd_multiplicator)
res_kmeans <- res_kmeans(data_training,lm)
classification_kmeans(new_data, res_kmeans, lm, outcome_string, adversity_string)

classification_naive <- function(new_data, adversity_string, outcome_string, lm) {
  residuals <- new_data[[outcome_string]] - (lm$coefficients[[1]] + new_data[[adversity_string]] * lm$coefficients[[2]])
  is_resilience_positive <- lm$coefficients[[2]] < 0
  
  res <- rep(NA_character_, length(residuals))
  
  if (is_resilience_positive) {
    res[residuals > 0] <- "resilient"
    res[residuals == 0] <- "average"
    res[residuals < 0] <- "vulnerable"
  } else {
    res[residuals < 0] <- "resilient"
    res[residuals == 0] <- "average"
    res[residuals > 0] <- "vulnerable"
  }
  
  return(res)
}


classification_confidence_intervals <- function(new_data,outcome_string,confidence,lm){
  actual <- new_data[[outcome_string]]
  pred <- as.data.frame(predict(lm_adjusted, newdata = df, interval = "confidence", level = confidence))
  sign_resilience <- lm$coefficients[[2]]<0
  groups <- get_groups_intervals(actual,pred,is_resilience_positive = sign_resilience)
  return(groups)
}

classification_prediction_intervals <- function(new_data,outcome_string,confidence,lm){
  actual <- new_data[[outcome_string]]
  pred <- as.data.frame(predict(lm_adjusted, newdata = df, interval = "prediction", level = confidence))
  sign_resilience <- lm$coefficients[[2]]<0
  groups <- get_groups_intervals(actual,pred,is_resilience_positive = sign_resilience)
  return(groups)
}


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

classification_credibility_intervals <- function(new_data,outcome_string,confidence,lm_cred){
  actual <- new_data[[outcome_string]]
  alpha<- 1-confidence
  pred <- get_credibility_intervals(lm_cred,new_data,lwr=alpha/2,upr=1-alpha/2)
  sign_resilience <- lm$coefficients[[2]]<0
  groups <- get_groups_intervals(actual,pred,is_resilience_positive = sign_resilience)
  return(groups)
}

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
