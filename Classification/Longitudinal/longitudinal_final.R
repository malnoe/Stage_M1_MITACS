# Longitudinal study on LORA

## Packages ####
library(ggplot2) # Ploting
library(dplyr) # Dataframe managment
library(rjags) # Spike and slab
library(SSranef) # Spike and slab
library(tidyr) # Dataframe managment
library(missForest) # MissForest, missing data
library(randomForest) # Random forest
library(car)  # VIF
library(caret)  # CreateFolds, cross-validation


## Set-up and Rdata ####
original_df <- readRDS("C:/Users/garan/Documents/Ecole/M1/Stage/Internship_repo/Classification/Longitudinal/ds_forJan.rds")
load("~/Ecole/M1/Stage/Internship_repo/Classification/Longitudinal/longitudinal_work_data.RData")

## Data preparation ####

# Select main variables
dh_variables <- paste0("dh_",c(1:28,44:58))
ghq_variables <- paste0("ghq_", 1:28)
variables <- c("id","age","gender",ghq_variables,paste0("pss_",1:10),dh_variables)
df <- original_df[1:30966,variables]

# Creat week variable
df[["week"]] <- rep(0:25,1191)

# Format dataframe
df <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))

# Create a variable indicating if the participant was there or not
df$present <- !apply(is.na(df[, ghq_variables]), 1, all)

# Replace negative values by NAs
for(variable in variables){
  df[[variable]] <- ifelse(df[[variable]]>=0, df[[variable]], NA)
}

# Recode PSS
variables_to_recode <- paste0("pss_", c(4, 5, 7, 8))
for (variable in variables_to_recode) {
  inverse_var <- paste0(variable, "_inverse")
  df[[inverse_var]] <- ifelse(df[[variable]] %in% 0:4, 4 - df[[variable]], NA)
}

# Build total scores for main variables
pss_variables <- c(paste0("pss_",c(1:3,6,9,10)),paste0("pss_",c(4, 5, 7, 8),"_inverse"))

df[["pss_sum"]] <- rowSums(df[,pss_variables],na.rm = TRUE)
df[["ghq_sum"]] <- rowSums(df[,ghq_variables],,na.rm = TRUE)
df[["dh_sum"]] <- rowSums(df[,dh_variables],,na.rm = TRUE)

# Select relevant lines and final variables
final_variables <- c("id","week","present","age","gender","pss_sum","ghq_sum","dh_sum")
df <- df[,final_variables]

## Participation diagnostic and selection of individuals ####

# We build a dataframe with the number of participant each weeks and the number of
# participants who partcipated all weeks until the considered week
df_participation <- df %>%
  filter(present) %>%
  count(week, name = "number_participants")

presence_counts <- df %>%
  filter(present) %>%
  group_by(id) %>%
  summarise(weeks_present = list(week), .groups = "drop")

df_full_attendance <- tibble()

for (week_number in 1:25) {
  weeks_required <- 1:week_number
  
  n_fully_present <- sum(sapply(presence_counts$weeks_present, function(weeks) {
    all(weeks_required %in% weeks)
  }))
  
  df_full_attendance <- bind_rows(df_full_attendance, tibble(
    week = week_number,
    fully_present = n_fully_present
  ))
}
df_participation <- df_participation %>%
  left_join(df_full_attendance, by = "week")

# We chose the maximum number week that would allow us to have 500 of more
# participants present every week. According to df_participation it is week 12.
participants_full_12 <- df %>%
  filter(week <= 12) %>%
  group_by(id) %>%
  summarise(n_weeks_present = sum(present), .groups = "drop") %>%
  filter(n_weeks_present == 12)

# d is the data frame with the corresponding individuals.
d <- df %>%
  filter(id %in% participants_full_12$id,week>0,week<=12)

## Residualization / Construction of a measurment of resilience ####

# We want to build two variables : 
# resilience for the GHQ~PSS model -> residuals_ghq_pss
# resilience for the GHQ~DH model -> residuals_ghq_dh

# Initialization with NAs
d[["residuals_ghq_pss"]] <- NA
d[["residuals_ghq_dh"]] <- NA

# Loop to do the residualization every week, applying influence statistics
# to put aside overly-influential individuals when there are enough individuals
for(adversity in c("pss_sum","dh_sum")){
  for(week_number in 1:25){
    # Get the people present that week
    index_present <- d$week==week_number&d$present
    # Check if there are enough people (ie > 2) to do the regression.
    if(sum(index_present)>2){
      data_regression <- d[index_present,c(adversity,"ghq_sum")]
      
      # Unadjusted linear model
      lm_unadjusted <- lm(as.formula(paste("ghq_sum", "~", adversity)), data = data_regression)
      
      # Identification of influencial points using Cook's D.
      used_data <- model.frame(lm_unadjusted)
      influencial_points <- which(cooks.distance(lm_unadjusted) > 4 / nrow(used_data))
      used_rows <- as.numeric(rownames(used_data))
      original_indices <- used_rows[influencial_points]
      
      # Cleaned data for the linear model
      df_clean <- data_regression[-original_indices, ]
      
      # Verify that there are at least two points for the adjusted regression
      if(nrow(df_clean)>2){
        # Adjusted linear model
        lm_adjusted <- lm(as.formula(paste("ghq_sum", "~", adversity)), data = df_clean)
        
        # Residuals of the adjusted linear model
        predicted_all <- predict(lm_adjusted, newdata = data_regression)
        residuals_all <- data_regression[["ghq_sum"]] - predicted_all
        
        # Put the residual in the result df
        if(adversity=="pss_sum"){
          d[index_present,"residuals_ghq_pss"] <- residuals_all
        }
        else{
          d[index_present,"residuals_ghq_dh"] <- residuals_all
        }
      }
    }
  }
}

## Classification : resilient / average / non-resilient -> spike and slab method ####

# Model GHQ~PSS
alpha_pss <- ss_ranef_alpha(y=d$residuals_ghq_pss, unit=d$id)
# Model GHQ~DH
alpha_dh <- ss_ranef_alpha(y=d$residuals_ghq_dh, unit=d$id)

# Function to get the number of people that are categorized as non-average depending on the threshold for the PIP
pct_PIP <- function(alpha_res,percentages){
  res <- data.frame(PIP=c(),pct_non_average=c(),pct_resilient=c(),pct_non_resilient=c())
  for(percentage in percentages){
    ranef_sum <- ranef_summary(alpha_res, ci = 0.95, digits = 2)
    pct_resilient <- sum(ranef_sum$PIP>=percentage&ranef_sum$Post.mean<0,na.rm=TRUE)/515*100
    pct_non_resilient <- sum(ranef_sum$PIP>=percentage&ranef_sum$Post.mean>0,na.rm=TRUE)/515*100
    res <- rbind(res,data.frame(PIP=c(percentage),pct_non_average=c(pct_resilient+pct_non_resilient),pct_resilient=c(pct_resilient),pct_non_resilient=c(pct_non_resilient)))
  }
  return(res)
}


# Results GHQ~PSS
posterior_summary(alpha_pss, ci = 0.90, digits = 2) # Summary of general variable of the model
ranef_summary(alpha_pss, ci = 0.95, digits = 2) # Summary of personal variable of the model
caterpillar_plot(alpha_pss, col_id = FALSE) # Caterpilar plot
pip_plot(alpha_pss, col_id = FALSE) # Values of PIP depending on theta
pct_PIP_alpha_pss<- pct_PIP(alpha_pss,seq(from=0.5,to=1,by=0.05)) # Classification depending on PIP threshold
summary(ranef_summary(alpha_pss, ci = 0.95, digits = 2)$PIP) # Distribution of PIP


# Results GHQ~DH
posterior_summary(alpha_dh, ci = 0.90, digits = 2)# Summary of general variable of the model
ranef_summary(alpha_dh, ci = 0.95, digits = 2)# Summary of personal variable of the model
caterpillar_plot(alpha_dh,col_id = FALSE)# Caterpilar plot
pip_plot(alpha_dh,col_id = FALSE)# Values of PIP depending on theta
pct_PIP_alpha_dh<- pct_PIP (alpha_dh,seq(from=0.5,to=1,by=0.05))# Classification depending on PIP threshold
summary(ranef_summary(alpha_dh, ci = 0.95, digits = 2)$PIP)# Distribution of PIP

## Functions - Visualization of the results ####

# Function to visualize the variation of the size of classes depending on PIP threshold
visualization_variation_size_classes <- function(pct_PIP_alpha,title){
  df_long <- pct_PIP_alpha %>%
    mutate(pct_average = 100 - pct_non_average) %>%
    pivot_longer(cols = c(pct_resilient, pct_non_resilient, pct_average),
                 names_to = "group",
                 values_to = "percentage") %>%
    mutate(group = dplyr::recode(group,
                                 "pct_resilient" = "Resilient",
                                 "pct_non_resilient" = "Non-resilient",
                                 "pct_average" = "Average"))
  
  plot <- ggplot(df_long, aes(x = PIP, y = percentage, color = group)) +
    geom_line(linewidth = 1.8) +
    geom_point(size=2) +
    labs(title = title,
         x = "PIP",
         y = "Percentage",
         color = "Group") +
    ylim(0, 100) +
    scale_color_manual(values = c("Average" = "gray10",
                                  "Resilient" = "gray40",
                                  "Non-resilient" = "gray80")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16)
    )
  return(plot)
}


visualization_longitudinal <- function(df, ranef_summary, posterior_summary, PIP_threshold,type="PSS") {
  
  # Create table for random effects
  id_list <- unique(df$id)
  
  ranef_df <- ranef_summary %>%
    mutate(id = id_list) %>%
    mutate(
      class = case_when(
        PIP < PIP_threshold ~ "Average",
        PIP >= PIP_threshold & Post.mean < 0 ~ "Resilient",
        PIP >= PIP_threshold & Post.mean > 0 ~ "Non-resilient",
        TRUE ~ "Average"
      )
    )
  
  # Join with original data
  d_plot <- df %>%
    left_join(ranef_df, by = "id")
  
  # Get alpha (mean intercept)
  alpha <- posterior_summary$Post.mean[1]
  
  # Get minimal effect sizes selected by spike-and-slab
  min_positive_theta <- min(ranef_df$Post.mean[ranef_df$PIP >= PIP_threshold & ranef_df$Post.mean > 0], na.rm = TRUE)
  max_negative_theta <- max(ranef_df$Post.mean[ranef_df$PIP >= PIP_threshold & ranef_df$Post.mean < 0], na.rm = TRUE)
  
  if(type=="PSS"){
    # Plot
    ggplot(data = d_plot, aes(x = week, y = residuals_ghq_pss, group = id, color = class)) +
      geom_line(alpha = 0.15) +
      geom_point(size = 0.5, alpha = 0.7) +
      
      # Horizontal lines
      geom_hline(yintercept = alpha, linetype = "dashed", color = "black", size = 0.6) +
      geom_hline(yintercept = min_positive_theta, linetype = "dotted", color = "firebrick", size = 0.6) +
      geom_hline(yintercept = max_negative_theta, linetype = "dotted", color = "steelblue", size = 0.6) +
      
      # Annotations for horizontal lines
      annotate("text", x = Inf, y = alpha, label = paste0("Intercept (α = ", round(alpha, 2), ")"),
               hjust = 1.1, vjust = -0.5, color = "black", size = 3) +
      
      annotate("text", x = Inf, y = min_positive_theta,
               label = paste0("Min θ (non_resilient) = ", round(min_positive_theta, 2)),
               hjust = 1.1, vjust = -0.5, color = "firebrick", size = 3) +
      
      annotate("text", x = Inf, y = max_negative_theta,
               label = paste0("Max θ (resilient) = ", round(max_negative_theta, 2)),
               hjust = 1.1, vjust = 1.5, color = "steelblue", size = 3) +
      
      scale_color_manual(values = c(
        "Average" = "gray80",
        "Resilient" = "steelblue",
        "Non-resilient" = "firebrick"
      )) +
      
      labs(
        x = "Week",
        y = "Residuals",
        color = "Class",
        title = "Residuals over time by individual",
        subtitle = paste0("Colored only for PIP ≥ ", PIP_threshold)
      ) +
      theme_minimal() 
  }
  else{
    # Plot
    ggplot(data = d_plot, aes(x = week, y = residuals_ghq_dh, group = id, color = class)) +
      geom_line(alpha = 0.15) +
      geom_point(size = 0.5, alpha = 0.7) +
      
      # Horizontal lines
      geom_hline(yintercept = alpha, linetype = "dashed", color = "black", size = 0.6) +
      geom_hline(yintercept = min_positive_theta, linetype = "dotted", color = "firebrick", size = 0.6) +
      geom_hline(yintercept = max_negative_theta, linetype = "dotted", color = "steelblue", size = 0.6) +
      
      # Annotations for horizontal lines
      annotate("text", x = Inf, y = alpha, label = paste0("Intercept (α = ", round(alpha, 2), ")"),
               hjust = 1.1, vjust = -0.5, color = "black", size = 3) +
      
      annotate("text", x = Inf, y = min_positive_theta,
               label = paste0("Min θ (non_resilient) = ", round(min_positive_theta, 2)),
               hjust = 1.1, vjust = -0.5, color = "firebrick", size = 3) +
      
      annotate("text", x = Inf, y = max_negative_theta,
               label = paste0("Max θ (resilient) = ", round(max_negative_theta, 2)),
               hjust = 1.1, vjust = 1.5, color = "steelblue", size = 3) +
      
      scale_color_manual(values = c(
        "Average" = "gray80",
        "Resilient" = "steelblue",
        "Non-resilient" = "firebrick"
      )) +
      
      labs(
        x = "Week",
        y = "Residuals",
        color = "Class",
        title = "Residuals over time by individual",
        subtitle = paste0("Colored only for PIP ≥ ", PIP_threshold)
      ) +
      theme_minimal()
  }
  # Plot
  ggplot(data = d_plot, aes(x = week, y = residuals_ghq_pss, group = id, color = class)) +
    geom_line(alpha = 0.15) +
    geom_point(size = 0.5, alpha = 0.7) +
    
    # Horizontal lines
    geom_hline(yintercept = alpha, linetype = "dashed", color = "black", size = 0.6) +
    geom_hline(yintercept = min_positive_theta, linetype = "dotted", color = "firebrick", size = 0.6) +
    geom_hline(yintercept = max_negative_theta, linetype = "dotted", color = "steelblue", size = 0.6) +
    
    # Annotations for horizontal lines
    annotate("text", x = Inf, y = alpha, label = paste0("Intercept (α = ", round(alpha, 2), ")"),
             hjust = 1.1, vjust = -0.5, color = "black", size = 3) +
    
    annotate("text", x = Inf, y = min_positive_theta,
             label = paste0("Min θ (non_resilient) = ", round(min_positive_theta, 2)),
             hjust = 1.1, vjust = -0.5, color = "firebrick", size = 3) +
    
    annotate("text", x = Inf, y = max_negative_theta,
             label = paste0("Max θ (resilient) = ", round(max_negative_theta, 2)),
             hjust = 1.1, vjust = 1.5, color = "steelblue", size = 3) +
    
    scale_color_manual(values = c(
      "Average" = "gray80",
      "Resilient" = "steelblue",
      "Non-resilient" = "firebrick"
    )) +
    
    labs(
      x = "Week",
      y = "Residuals",
      color = "Class",
      title = "Residuals over time by individual",
      subtitle = paste0("Colored only for PIP ≥ ", PIP_threshold)
    ) +
    theme_minimal()
}

# Visualization of the main trajectories for one threshold with confidence intervals
visualization_main_trajectories <- function(df, ranef_summary, posterior_summary, PIP_threshold,type="PSS"){
  # Create table for random effects
  id_list <- unique(df$id)
  
  ranef_df <- ranef_summary %>%
    mutate(id = id_list) %>%
    mutate(
      class = case_when(
        PIP < PIP_threshold ~ "Average",
        PIP >= PIP_threshold & Post.mean < 0 ~ "Resilient",
        PIP >= PIP_threshold & Post.mean > 0 ~ "Non-resilient",
        TRUE ~ "Average"
      )
    )
  
  # Join with original data
  d_plot <- df %>%
    left_join(ranef_df, by = "id")
  
  # Get alpha (mean intercept)
  alpha <- posterior_summary$Post.mean[1]
  
  # Mean trajectories with confidence intervals
  if(type=="PSS"){
    df_trajectories <- d_plot %>%
      group_by(week, class) %>%
      summarise(
        value = mean(residuals_ghq_pss, na.rm = TRUE),
        sd = sd(residuals_ghq_pss, na.rm = TRUE),
        n = n(),
        se = sd / sqrt(n),
        lower = value - qt(0.975, df = n-1) * se,
        upper = value + qt(0.975, df = n-1) * se,
        .groups = "drop"
      )
  }
  else{
    df_trajectories <- d_plot %>%
      group_by(week, class) %>%
      summarise(
        value = mean(residuals_ghq_dh, na.rm = TRUE),
        sd = sd(residuals_ghq_pss, na.rm = TRUE),
        n = n(),
        se = sd / sqrt(n),
        lower = value - qt(0.975, df = n-1) * se,
        upper = value + qt(0.975, df = n-1) * se,
        .groups = "drop"
      )
  }
  
  
  # Plot
  ggplot(data = df_trajectories, aes(x = week, y = value, color = class, fill = class)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line(size = 0.8) +
    geom_point(size = 1) +
    
    geom_hline(yintercept = alpha, linetype = "dashed", color = "black", size = 0.6) +
    
    scale_color_manual(values = c(
      "Average" = "gray60",
      "Resilient" = "steelblue",
      "Non-resilient" = "firebrick"
    )) +
    scale_fill_manual(values = c(
      "Average" = "gray80",
      "Resilient" = "steelblue",
      "Non-resilient" = "firebrick"
    )) +
    
    labs(
      x = "Week",
      y = "Mean residuals",
      color = "Class",
      fill = "Class",
      title = "Mean residuals over time by class",
      subtitle = paste0("Class defined with PIP ≥ ", PIP_threshold)
    ) +
    theme_minimal()
}

# Visualization of main trajectories of groups for different thresholds
visualization_trajectories_multiple_thresholds <- function(df, ranef_summary, posterior_summary, thresholds, type ="PSS") {
  
  id_list <- unique(df$id)
  
  all_trajectories <- lapply(thresholds, function(PIP_threshold) {
    # Assign classes based on current threshold
    ranef_df <- ranef_summary %>%
      mutate(id = id_list) %>%
      mutate(
        class = case_when(
          PIP < PIP_threshold ~ "Average",
          PIP >= PIP_threshold & Post.mean < 0 ~ "Resilient",
          PIP >= PIP_threshold & Post.mean > 0 ~ "Non-resilient",
          TRUE ~ "Average"
        )
      )
    
    d_plot <- df %>%
      left_join(ranef_df, by = "id")
    
    # Aggregate trajectories
    if(type=="PSS"){
      d_plot %>%
        group_by(week, class) %>%
        summarise(
          value = mean(residuals_ghq_pss, na.rm = TRUE),
          sd = sd(residuals_ghq_pss, na.rm = TRUE),
          n = n(),
          se = sd / sqrt(n),
          lower = value - qt(0.975, df = n-1) * se,
          upper = value +  qt(0.975, df = n-1) * se,
          .groups = "drop"
        ) %>%
        mutate(PIP_threshold = PIP_threshold)
    }
    else{
      d_plot %>%
        group_by(week, class) %>%
        summarise(
          value = mean(residuals_ghq_dh, na.rm = TRUE),
          sd = sd(residuals_ghq_pss, na.rm = TRUE),
          n = n(),
          se = sd / sqrt(n),
          lower = value - qt(0.975, df = n-1) * se,
          upper = value + qt(0.975, df = n-1) * se,
          .groups = "drop"
        ) %>%
        mutate(PIP_threshold = PIP_threshold)
    }
  }) %>% bind_rows()
  
  alpha <- posterior_summary$Post.mean[1]
  
  # Plot
  # Plot
  ggplot(data = all_trajectories, aes(x = week, y = value, color = class, linetype = as.factor(PIP_threshold), group = interaction(class, PIP_threshold))) +
    geom_line(size = 0.8) +
    scale_x_continuous(breaks = 1:12)+
    scale_color_manual(values = c(
      "Average" = "gray10",
      "Resilient" = "gray40",
      "Non-resilient" = "gray80"
    )) +
    scale_linetype_manual(
      values = c("solid", "dashed", "dotdash", "twodash")[1:length(thresholds)],
      name = "PIP threshold"
    ) +
    
    labs(
      x = "Week",
      y = "Mean residuals",
      color = "Class",
      fill = "Class",
      title = "Mean residual trajectories by class for multiple PIP thresholds",
      linetype = "PIP threshold"
    ) +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16)
    )
}


## Application - Visualization of the results ####

# Variations of class size as a function of the PIP threshold
# DH
visualization_variation_size_classes(pct_PIP_alpha_dh,"Group Proportions by PIP for GHQ~DH")
# PSS
visualization_variation_size_classes(pct_PIP_alpha_pss,"Group Proportions by PIP for GHQ~PSS")

# Trajectories of all individuals colored depending on their group
# DH
ranef_summary_DH <- ranef_summary(alpha_dh, ci = 0.95, digits = 2)
posterior_sumarry_DH <- posterior_summary(alpha_dh, ci = 0.90, digits = 2)
visualization_longitudinal(d,ranef_summary_DH,posterior_sumarry_DH,0.99,type="DH")
# PSS
ranef_summary_PSS <- ranef_summary(alpha_pss, ci = 0.95, digits = 2)
posterior_sumarry_PSS <- posterior_summary(alpha_pss, ci = 0.90, digits = 2)
visualization_longitudinal(d,ranef_summary_PSS,posterior_sumarry_PSS,0.99,type="PSS")

# Visualization of the main trajectories for one threshold with confidence intervals
# DH
visualization_main_trajectories(d,ranef_summary_DH,posterior_sumarry_DH,0.99,type="DH")
# PSS
visualization_main_trajectories(d,ranef_summary_PSS,posterior_sumarry_PSS,0.99,type="PSS")

# Visualization of main trajectories of groups for different thresholds
# DH
visualization_trajectories_multiple_thresholds(d,ranef_summary_DH,posterior_sumarry_DH,list(0.70,0.8,0.9),type="DH")
# PSS
visualization_trajectories_multiple_thresholds(d,ranef_summary_PSS,posterior_sumarry_PSS,list(0.70,0.8,0.9),type="PSS")


## Evaluation of predictive power : Data preparation - Individuals selection, residualization and spike and slab ####

# 1. Get the predictors in the original dataframe and select the individuals who have enough to build a new_d
bfi <- paste0("bfi_",1:10)
cdrisk <- paste0("cdrisk_",1:25)
cerq <- paste0("cerq_",1:25)
cope <- paste0("cope_",1:28)
ctq <- paste0("ctq_",1:25)
fsozu <- paste0("fsozu_",1:14)
gpass <- paste0("gpass_",c(1,2,6,8,10:14,16:18,20,26,27))
pas_content <- paste0("pass_content", 1:14)
gse <- paste0("gse_",1:10)
ielc <- paste0("ielc_",1:28)
le <- paste0("le_",1:27)
lotr <- paste0("lotr_",1:10)
variables <- c("id","age","income","employment_status","gender",bfi,cdrisk,cerq,cope,ctq,fsozu,gpass,pas_content,gse,ielc,le,lotr)
predictive_df <- original_df[original_df$id %in% d$id & original_df$t==1,variables]
# Format the dataframe
predictive_df <- as.data.frame(lapply(predictive_df, function(x) as.numeric(as.character(x))))
# Replace negative values by NAs
for(variable in variables){
  predictive_df[[variable]] <- ifelse(predictive_df[[variable]]>=0, predictive_df[[variable]], NA)
}

# Select people with more than 70% of correct data
dim(predictive_df)
na_ratio <- rowMeans(is.na(predictive_df), na.rm = FALSE)
indices <- which(na_ratio > 0.3 & !is.na(na_ratio))  # on exclut les NA
predictive_df <- predictive_df[-indices,]

# We do missForest
predictive_df.mf <- missForest::missForest(xmis = predictive_df)
predictive_df <- predictive_df.mf$ximp

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
predictive_df <- recode(predictive_df,bfi_to_recode,min=1,max=5)
bfi <- c(paste0("bfi_",c(2,6,8,9,10)),paste0("bfi_",c(1,3,4,5,7),"_inverse"))

# CTQ
ctq_to_recode <- paste0("ctq_",c(2,5,7,12,17,23,25))
predictive_df <- recode(predictive_df,ctq_to_recode,min=1,max=5)
ctq <- c(paste0("ctq_",c(2,5,7,12,17,23,25),"_inverse"),paste0("ctq_",c(1,3,4,6,8,9,10,11,13,14,15,16,18,19,20,21,22,24)))

# IELC
ielc_to_recode <- paste0("ielc_",c(3:5,10:13,15,22,26,27))
predictive_df <- recode(predictive_df,ielc_to_recode,min=0,max=1)
ielc <- c(paste0("ielc_",c(3:5,10:13,15,22,26,27),"_inverse"),paste0("ielc_",c(1,2,4:9,14,16:21,23:25,28)))

# LOTR
lotr_to_recode <- paste0("lotr_",c(3,7,9))
predictive_df <- recode(predictive_df,lotr_to_recode,min=0,max=4)
lotr <- c(paste0("lotr_",c(3,7,9),"_inverse"),paste0("lotr_",c(1,4,10)))

# PAS
pas_torecode <- paste0("gpass_",c(1,8,10,13,16,18,20,26,27))
predictive_df <- recode(predictive_df,pas_torecode,min=1,max=4)
gpass <- c(paste0("gpass_",c(1,8,10,13,16,18,20,26,27),"_inverse"),paste0("gpass_",c(2,6,11,12,14)))

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
predictive_df[["bfi_sum"]] <- rowSums(predictive_df[,bfi],na.rm = TRUE)
predictive_df[["cdrisk_sum"]] <- rowSums(predictive_df[,cdrisk],na.rm = TRUE)
predictive_df[["cerq_sum"]] <- rowSums(predictive_df[,cerq],na.rm = TRUE)
predictive_df[["able"]] <- rowSums(predictive_df[,cope_able],na.rm = TRUE)
predictive_df[["verl"]] <- rowSums(predictive_df[,cope_verl],na.rm = TRUE)
predictive_df[["emu"]] <- rowSums(predictive_df[,cope_emu],na.rm = TRUE)
predictive_df[["ruck"]] <- rowSums(predictive_df[,cope_ruck],na.rm = TRUE)
predictive_df[["poum"]] <- rowSums(predictive_df[,cope_poum],na.rm = TRUE)
predictive_df[["hum"]] <- rowSums(predictive_df[,cope_hum],na.rm = TRUE)
predictive_df[["akbe"]] <- rowSums(predictive_df[,cope_akbe],na.rm = TRUE)
predictive_df[["aldro"]] <- rowSums(predictive_df[,cope_aldro],na.rm = TRUE)
predictive_df[["insun"]] <- rowSums(predictive_df[,cope_insun],na.rm = TRUE)
predictive_df[["ause"]] <- rowSums(predictive_df[,cope_ause],na.rm = TRUE)
predictive_df[["plan"]] <- rowSums(predictive_df[,cope_plan],na.rm = TRUE)
predictive_df[["akze"]] <- rowSums(predictive_df[,cope_akze],na.rm = TRUE)
predictive_df[["sebe"]] <- rowSums(predictive_df[,cope_sebe],na.rm = TRUE)
predictive_df[["reli"]] <- rowSums(predictive_df[,cope_reli],na.rm = TRUE)
predictive_df[["ctq_sum"]] <- rowSums(predictive_df[,ctq],na.rm = TRUE)
predictive_df[["fsozu_sum"]] <- rowSums(predictive_df[,fsozu],na.rm = TRUE)/14
predictive_df[["gpass_sum"]] <- rowSums(predictive_df[,gpass],na.rm = TRUE)
predictive_df[["pas_content_sum"]] <- rowSums(predictive_df[,pas_content],na.rm = TRUE)
predictive_df[["gse_sum"]] <- rowSums(predictive_df[,gse],na.rm = TRUE)
predictive_df[["ielc_sum"]] <- rowSums(predictive_df[,ielc],na.rm = TRUE)
predictive_df[["le_sum"]] <- rowSums(predictive_df[,le],na.rm = TRUE)
predictive_df[["lotr_sum"]] <- rowSums(predictive_df[,lotr],na.rm = TRUE)

# 2. Build new_d with inly the people who have enough predictive data (ie id in predictive_df) to get the pss, dh, ghq
# and also have people who participated all 12 first weeks
new_d <- d %>% filter(id %in% predictive_df$id) %>% dplyr::select(id,week,present,pss_sum,dh_sum,ghq_sum)

# 3. Redo the residualization since we have less people
new_d[["residuals_ghq_pss"]] <- NA
new_d[["residuals_ghq_dh"]] <- NA

for(adversity in c("pss_sum","dh_sum")){
  for(week_number in 1:25){
    # Get the people present that week
    index_present <- new_d$week==week_number&new_d$present
    # Check if there are enough people (ie > 2) to do the regression.
    if(sum(index_present)>2){
      data_regression <- new_d[index_present,c(adversity,"ghq_sum")]
      
      # Unadjusted linear model
      lm_unadjusted <- lm(as.formula(paste("ghq_sum", "~", adversity)), data = data_regression)
      
      # Identification of influencial points using Cook's D.
      used_data <- model.frame(lm_unadjusted)
      influencial_points <- which(cooks.distance(lm_unadjusted) > 4 / nrow(used_data))
      used_rows <- as.numeric(rownames(used_data))
      original_indices <- used_rows[influencial_points]
      
      # Cleaned data for the LM
      df_clean <- data_regression[-original_indices, ]
      
      # Verify that there are at least two points for the adjusted regression
      if(nrow(df_clean)>2){
        # Adjusted linear model
        lm_adjusted <- lm(as.formula(paste("ghq_sum", "~", adversity)), data = df_clean)
        
        # Residuals of the adjusted linear model
        predicted_all <- predict(lm_adjusted, newdata = data_regression)
        residuals_all <- data_regression[["ghq_sum"]] - predicted_all
        
        # Put the residual in the result df
        if(adversity=="pss_sum"){
          new_d[index_present,"residuals_ghq_pss"] <- residuals_all
        }
        else{
          new_d[index_present,"residuals_ghq_dh"] <- residuals_all
        }
      }
    }
  }
}

# 4. Redo the spike and slab
# GHQ~PSS
new_alpha_pss <- ss_ranef_alpha(y=new_d$residuals_ghq_pss, unit=new_d$id)
pct_PIP_new_alpha_pss<- pct_PIP(new_alpha_pss,seq(from=0,to=1,by=0.05))


# GHQ~DH
new_alpha_dh <- ss_ranef_alpha(y=new_d$residuals_ghq_dh, unit=new_d$id)
pct_PIP_new_alpha_dh<- pct_PIP (new_alpha_dh,seq(from=0,to=1,by=0.05))

## Evaluation of predictive power : Classification depending on the threshold ####

# Function to get the classification for each individual
get_classification <- function(d,ranef_summary,list_thresholds=seq(from = 0, to = 1, by = 0.05)){
  id_list <- unique(d$id)
  classifications <- tibble(id = id_list)
  
  for (PIP_threshold in list_thresholds) {
    
    col_name <- paste0("class_", PIP_threshold)
    
    ranef_df <- ranef_summary %>%
      mutate(id = id_list) %>%
      mutate(
        !!col_name := case_when(
          PIP < PIP_threshold ~ "average",
          PIP >= PIP_threshold & Post.mean < 0 ~ "resilient",
          PIP >= PIP_threshold & Post.mean > 0 ~ "non-resilient",
          TRUE ~ "average"
        )
      ) %>%
      dplyr::select(id, !!col_name)
    
    classifications <- classifications %>%
      left_join(ranef_df, by = "id")
  }
  return(classifications)
}

# Application to our data
#DH
ranef_summary_DH <- ranef_summary(new_alpha_dh, ci = 0.95, digits = 2)
classification_dh <- get_classification(new_d,ranef_summary_DH)
# PSS
ranef_summary_PSS <- ranef_summary(new_alpha_pss, ci = 0.95, digits = 2)
classification_pss <- get_classification(new_d,ranef_summary_PSS)

## Evaluation of predictive power : Functions for estimation and visualization of results ####

# Function to get the metrics (accuracy, precision, recall ...)
classification_metrics <- function(true_labels, predicted_labels,levels=c("resilient", "average", "non-resilient")) {
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

# Function to apply VIF methodology
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

# Function to dp the estimation of the predictive power of different thresholds
estimation_classification_cv <- function(predictive_df, classification_df, predictors,list_group_names, k = 5, seed = 123){
  set.seed(seed)
  res <- data.frame()
  df <- predictive_df
  
  # Remove high VIF predictors
  predictors <- remove_high_vif(predictive_df, predictors, threshold = 5)
  
  for(i in 1:length(list_group_names)){
    group_name <- list_group_names[[i]]
    print(paste("Running CV for group:", group_name))
    df[["groups"]] <- as.factor(classification_df[[group_name]])
    print(head(df[["groups"]]))
    
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
    
    support_total <- average_metric("support", "resilient") + average_metric("support", "average") + average_metric("support", "non-resilient")
    null_model <- max(average_metric("support", "resilient"), average_metric("support", "average"), average_metric("support", "non-resilient")) / support_total
    
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
      precision_non_resilient = average_metric("precision_per_class", "non-resilient"),
      recall_non_resilient = average_metric("recall_per_class", "non-resilient"),
      f1score_non_resilient = average_metric("f1_per_class", "non-resilient"),
      support_resilient = average_metric("support", "resilient"),
      support_average = average_metric("support", "average"),
      support_non_resilient = average_metric("support", "non-resilient")
    ))
  }
  
  return(res)
}

# Visualization of the accuracy of thre classifier and the null model as a function of the average group size
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

# Visualization of recall and precision for the resilient group
visualization_recall_precision <- function(df_perf,list_groups){
  df_perf_class_comparison <- df_perf %>% filter(group_name %in% list_groups)
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

## Evaluation of predictive power : Estimation and results ####
predictors_items <- c(bfi,cdrisk,cerq,cope,ctq,fsozu,gpass,pas_content,gse,ielc,le,lotr,"id","age","income","employment_status","gender")
predictors_sum <- c("bfi_sum","cdrisk_sum","cerq_sum","ctq_sum","fsozu_sum","gpass_sum","pas_content_sum","gse_sum","ielc_sum","le_sum","lotr_sum","id","age","income","employment_status","gender",
                    "able","verl","emu","ruck","poum","hum","akbe","aldro","insun","ause","plan","akze","sebe","reli")

to_test <- paste0("class_",seq(from=0,to=1,by=0.05))

# sums and subscales
df_perf_dh <- estimation_classification_cv(predictive_df,classification_dh,predictors_sum,list_group_names = to_test ,k=10)
df_perf_pss <- estimation_classification_cv(predictive_df,classification_pss,predictors_sum,list_group_names = to_test ,k=10)

# items
df_perf_dh_items <- estimation_classification_cv(predictive_df,classification_dh,predictors_items,list_group_names = to_test ,k=10)
df_perf_pss_items <- estimation_classification_cv(predictive_df,classification_pss,predictors_items,list_group_names = to_test ,k=10)


# Visualization of results / Figures
# sums and subscales
# DH
comparison_accuracy_null_model_classifier(df_perf_dh)
visualization_recall_precision(df_perf_dh,to_test)
# PSS
comparison_accuracy_null_model_classifier(df_perf_pss)
visualization_recall_precision(df_perf_pss,to_test)

# items
# DH
comparison_accuracy_null_model_classifier(df_perf_dh_items)
visualization_recall_precision(df_perf_dh_items,to_test)
# PSS
comparison_accuracy_null_model_classifier(df_perf_pss_items)
visualization_recall_precision(df_perf_pss_items,to_test)
