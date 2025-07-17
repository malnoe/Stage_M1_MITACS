# LORA Network Analysis

## Packages ####
library(dplyr) # Dataframe managment and selection
library(ggplot2) # Visualization
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
  
  # 6. Build the sums
  df_selected[["PAS_content"]] <- rowSums(df_selected[,content_vars],na.rm = TRUE)
  df_selected[["PAS_process"]] <- rowSums(df_selected[,process_vars],,na.rm = TRUE)
  
  return(df_selected)
}

# Define weeks and required variables
weeks_needed <- c(1, 7, 13)
variables_needed <- colnames(LORA)
LORAr <- select_data(LORA, weeks_needed,variables_needed)

## Profiling - Cross-sectional LPA ####

## Profiling - Longitudinal Joint KML3D ####

## Profiling - Longitudinal Content ####

## Profiling - Longitudinal Processes ####

## Network Analysis ####