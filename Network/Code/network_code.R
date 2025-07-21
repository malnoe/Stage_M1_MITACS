# LORA Network Analysis

## Packages ####
library(dplyr) # Dataframe managment and selection
library(ggplot2) # Visualization
library(kml)
library(kml3d)
library(rpart) # Regression tree
library(tidyLPA) # LPA
library(tidyr) # Dataframe managment 

## Get back to work ####
load("~/Ecole/M1/Stage/Internship_repo/Network/Code/network_rdata.RData")

## Set-up ####
# Process and content variables
process_vars <- c("gpass_1", "gpass_2", "gpass_6", "gpass_8", "gpass_10", 
                  paste0("gpass_", 11:14), paste0("gpass_", 16:18), 
                  "gpass_20", "gpass_26", "gpass_27")
content_vars <- paste0("pass_content", 1:14)
# Loading data
original_LORA <- readRDS("C:/Users/garan/Documents/Ecole/M1/Stage/Internship_repo/Network/Data/ds_forJan.rds")
# Taking only the weeks of interest
LORA <- original_LORA[original_LORA$t%in%c(1,7,13,19,25),]
# Formating
LORA <- as.data.frame(lapply(LORA[,-2], function(x) as.numeric(as.character(x))))
# Replacing -991 and -993 by NAs
LORA[LORA < 0] <- NA
# Creating "present" variable indicating wether or not an individual is present that week
# All question must be NAs to be marked as not present
LORA$present <- !apply(is.na(LORA[ , -(1:3)]), 1, all)

dim(LORA) # 5955 = 5 x 1191

## Contents ####

# Sample size analysis -> look into the available data at each time point

# Week and data selection -> Depending on the results of the sample size analysis, get the weeks needed and also get the variables needed.

# Profiling -> Find latent profiles in the answers to the 
# Positive Appraisal Style Questionnaire in different ways :
## Cross-sectional LPA
## Longitudinal Joint KML3D
## Longitudinal Content
## Longitudinal Processes

# Network Analysis ->


## Sample size Analysis #####

# Function to get how many people are present in the asked week (possibilities: 1,7,13,19,25)
# People must have answered all the question of the PAS questionnaire
sample_size <- function(df, list_week) {
  required_vars <- c(process_vars, content_vars)
  
  # Remove individuals with NA in required variables
  df_filtered <- df %>%
    filter(across(all_of(required_vars), ~ !is.na(.x)))
  
  # Check presence across selected weeks
  df_filtered %>%
    filter(t %in% list_week) %>%
    group_by(id) %>%
    summarise(n_weeks_present = sum(present), .groups = "drop") %>%
    filter(n_weeks_present == length(list_week)) %>%
    nrow()
}

sample_size(LORA,list(1,7,13)) # 622  individuals participated in first 3 timepoints.
sample_size(LORA,list(7,13,19)) # 119  individuals partipated in the 3 timepoints in the middle
sample_size(LORA,list(1,7,13,19)) # 119 individuals participated in first 4 timepoints
sample_size(LORA,list(1,7,13,19,25)) # No one participated in all 5 timepoints
sample_size(LORA,list(25)) # Actually no one participated in timepoint 5
# -> There is a huge drop off at 4 weeks, only 119 people remaining and nobody participated in week 5
# Therefore : focus on the first 3 time points.

## Week and data selection ####
# Function to recode the variable that are negatively worded
recode <- function(df,variables_to_recode,min,max){
  for(variable in variables_to_recode){
    inverse_var <- paste0(variable, "_inverse")
    df[[inverse_var]] <- ifelse(df[[variable]] %in% min:max, max + min - df[[variable]], NA)
  }
  return(df)
}

# Function to select the data and build the sums
select_data <- function(df,list_week,variables_needed){
  required_vars <- c(process_vars, content_vars)
  # Filter to weeks of interest
  df_week <- df %>% filter(t %in% weeks_needed)
  # Keep only individuals with no missing values on required variables
  ids_with_complete_data <- df_week %>%
    group_by(id) %>%
    summarise(all_answered = all(!is.na(unlist(across(all_of(required_vars))))),
              .groups = "drop") %>%
    filter(all_answered) %>%
    pull(id)
  
  # Check who is present in all weeks
  ids_with_all_weeks_present <- df_week %>%
    filter(id %in% ids_with_complete_data) %>%
    group_by(id) %>%
    summarise(n_weeks_present = sum(present), .groups = "drop") %>%
    filter(n_weeks_present == length(weeks_needed)) %>%
    pull(id)
  
  # 4. Final data frame with selected individuals
  df_selected <- df %>% filter(id %in% ids_with_all_weeks_present & t%in%weeks_needed)
  
  # 5. Select the right variables
  df_selected <- df_selected[,variables_needed]
  
  # 6. Recode the variables
  pas_torecode <- paste0("gpass_",c(1,8,10,13,16,18,20,26,27))
  df_selected <- recode(df_selected,pas_torecode,min=1,max=4)
  process_vars <- c(paste0("gpass_",c(1,8,10,13,16,18,20,26,27),"_inverse"),paste0("gpass_",c(2,6,11,12,14)))
  
  # 7. Build the sums
  df_selected[["PAS_content"]] <- rowSums(df_selected[,content_vars],na.rm = TRUE)
  df_selected[["PAS_process"]] <- rowSums(df_selected[,process_vars],,na.rm = TRUE)
  
  return(df_selected)
}

# Define weeks and required variables
weeks_needed <- c(1, 7, 13)
variables_needed <- colnames(LORA)
LORAr <- select_data(LORA, weeks_needed,variables_needed)

## Profiling - Cross-sectional LPA ####
# LPA with baseline data
LORAr %>%
  filter(t==1) %>%
  dplyr::select(PAS_content, PAS_process) %>%
  single_imputation() %>%
  estimate_profiles(2:6, 
                    variances = c("equal","equal","varying"),
                    covariances = c("equal","zero","zero"),nrep = 5) %>%
  compare_solutions(statistics = c("AIC","BIC"))

# Model 3 with 2 classes -> AIC=7296, BIC=7331

# Plot best result for AIC, BIC and analytical process
LORAr %>%
  filter(t==1) %>%
  dplyr::select(PAS_content, PAS_process) %>%
  single_imputation() %>%
  estimate_profiles(2, 
                    variances = c("equal"),
                    covariances = c("equal"),nrep = 5) %>%
  plot_profiles()

LPA_LORA_profiles <- LORAr %>%
  filter(t==1) %>%
  dplyr::select(PAS_content, PAS_process) %>%
  single_imputation() %>%
  estimate_profiles(2, 
                    variances = c("equal"),
                    covariances = c("equal"),nrep = 5)

# Profile 1 (red): average-high on contents and average-higher on process
# Profile 2 (blue): average-low on contents and low on process


## Profiling - Cross-sectional Step function ####
LORAr[["ghq_sum"]] <- rowSums(LORAr[,paste0("ghq_",1:28)],na.rm = TRUE)
LORA_profiling <- LORAr %>%
    dplyr::select(id,t,PAS_content, PAS_process,ghq_sum)
LORA_profiling <- reshape(LORA_profiling,direction="wide",v.names=c("PAS_content","PAS_process","ghq_sum"),idvar="id",timevar="t")

handmade_cutoff <- function(k,type){
  if(type=="process"){
    df <- LORA_profiling %>%
      mutate(PAS_bin = cut(PAS_process.1, breaks = k))
    
    step_data <- df %>%
      group_by(PAS_bin) %>%
      summarize(mean_depression = mean(ghq_sum.1, na.rm = TRUE),
                PAS_mid = mean(PAS_process.1, na.rm = TRUE))
    
    plot <- ggplot(step_data, aes(x = PAS_mid, y = mean_depression)) +
      geom_step(direction = "hv") +
      labs(x = "PAS_process", y = "Mean Depression") +
      theme_minimal()
  }
  else{
    df <- LORA_profiling %>%
      mutate(PAS_bin = cut(PAS_content.1, breaks = k))
    
    step_data <- df %>%
      group_by(PAS_bin) %>%
      summarize(mean_depression = mean(ghq_sum.1, na.rm = TRUE),
                PAS_mid = mean(PAS_content.1, na.rm = TRUE))
    
    plot <- ggplot(step_data, aes(x = PAS_mid, y = mean_depression)) +
      geom_step(direction = "hv") +
      labs(x = "PAS_content", y = "Mean Depression") +
      theme_minimal()
  }
  return(plot)
}

rpart_cutoff <- function(df, type){
  if(type=="process"){
    # Fit tree
    tree <- rpart(ghq_sum.1 ~ PAS_process.1, data = df, control = rpart.control(cp = 0.01))
    
    # Assign each observation to a terminal node and get prediction
    df$node <- tree$where
    df$prediction <- predict(tree)
    
    # Total number of individuals
    n_total <- nrow(df)
    
    # Get bin info: prediction, range, count, and proportion
    step_data <- df %>%
      group_by(node, prediction) %>%
      summarize(x_min = min(PAS_process.1),
                x_max = max(PAS_process.1),
                count = n(),
                .groups = "drop") %>%
      mutate(percentage = round(100 * count / n_total, 1)) %>%
      arrange(x_min)
    
    # Plot
    plot <- ggplot() +
      # Raw data points
      geom_point(data = df, aes(x = PAS_process.1, y = ghq_sum.1), alpha = 0.3, color = "gray40", size = 1) +
      
      # Step segments
      geom_segment(data = step_data,
                   aes(x = x_min, xend = x_max,
                       y = prediction, yend = prediction),
                   size = 1.2, color = "steelblue") +
      
      # Percentage labels
      geom_text(data = step_data,
                aes(x = (x_min + x_max) / 2,
                    y = prediction + 0.1 * sd(df$ghq_sum.1, na.rm = TRUE),
                    label = paste0(percentage, "%")),
                size = 3.5, vjust = 0) +
      
      labs(
        x = "PAS_process score at T1",
        y = "GHQ Score at T1",
        title = "Step function with tree-based cutoffs and raw data",
      ) +
      theme_minimal()
    
  }
  else if(type=="content"){
    # Fit tree
    tree <- rpart(ghq_sum.1 ~ PAS_content.1, data = df, control = rpart.control(cp = 0.01))
    
    # Assign each observation to a terminal node and get prediction
    df$node <- tree$where
    df$prediction <- predict(tree)
    
    # Total number of individuals
    n_total <- nrow(df)
    
    # Get bin info: prediction, range, count, and proportion
    step_data <- df %>%
      group_by(node, prediction) %>%
      summarize(x_min = min(PAS_content.1),
                x_max = max(PAS_content.1),
                count = n(),
                .groups = "drop") %>%
      mutate(percentage = round(100 * count / n_total, 1)) %>%
      arrange(x_min)
    
    # Plot
    plot <- ggplot() +
      # Raw data points
      geom_point(data = df, aes(x = PAS_content.1, y = ghq_sum.1), alpha = 0.3, color = "gray40", size = 1) +
      
      # Step segments
      geom_segment(data = step_data,
                   aes(x = x_min, xend = x_max,
                       y = prediction, yend = prediction),
                   size = 1.2, color = "steelblue") +
      
      # Percentage labels
      geom_text(data = step_data,
                aes(x = (x_min + x_max) / 2,
                    y = prediction + 0.1 * sd(df$ghq_sum.1, na.rm = TRUE),
                    label = paste0(percentage, "%")),
                size = 3.5, vjust = 0) +
      
      labs(
        x = "PAS_content score at T1",
        y = "GHQ Score at T1",
        title = "Step function with tree-based cutoffs and raw data"
      ) +
      theme_minimal()
  }
  return(plot)
}

handmade_cutoff(6,type="process")
handmade_cutoff(6,type="content")
rpart_cutoff(df=LORA_profiling,type="content")
rpart_cutoff(df=LORA_profiling,type="process")


## Profiling - Longitudinal Joint KML3D ####
# Data Preparation
cld3dPregTemp <- cld3d(LORA_profiling,timeInData=list(PAS_content=c(2,5,8),PAS_process=c(3,6,9)))
# Building "optimal" clusteration
kml3d(cld3dPregTemp,2:7,nbRedrawing=50,toPlot="nothing")
choice(cld3dPregTemp)
# Visualizing in 3D
plotMeans3d(cld3dPregTemp,2)
plotMeans3d(cld3dPregTemp,3)
plotMeans3d(cld3dPregTemp,4)
plotMeans3d(cld3dPregTemp,5)
plotMeans3d(cld3dPregTemp,6)
plotMeans3d(cld3dPregTemp,7)


## Profiling - Longitudinal Content ####
cldLORA_content <- cld(LORA_profiling,timeInData =c(2,5,8))
kml(cldLORA_content,2:10,nbRedrawing = 15,toPlot="nothing")
plotAllCriterion(cldLORA_content)
plot(cldLORA_content,2)
plot(cldLORA_content,3)
plot(cldLORA_content,4)
plot(cldLORA_content,5)
plot(cldLORA_content,6)
plot(cldLORA_content,7)

## Profiling - Longitudinal Processes ####
cldLORA_process <- cld(LORA_profiling,timeInData =c(3,6,9))
kml(cldLORA_process,2:10,nbRedrawing = 15,toPlot="nothing")
plotAllCriterion(cldLORA_process)
plot(cldLORA_process,2)
plot(cldLORA_process,3)
plot(cldLORA_process,4)
plot(cldLORA_process,5)
plot(cldLORA_process,6)
plot(cldLORA_process,10)

## Network Analysis ####

