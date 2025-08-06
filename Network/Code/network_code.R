# LORA Network Analysis

## Packages ####
library(dplyr) # Dataframe managment and selection
library(ggplot2) # Visualization
library(kml) # Profiling
library(kml3d) # Profiling
library(rpart) # Regression tree
library(tidyLPA) # Profiling LPA
library(tidyr) # Dataframe managment 
library(qgraph) # Plot the networks
library(glmnet) # Network
library(lavaan) # Network
library(networktools) # Item selection
library(tidyr)


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
  required_vars <- c(process_vars, content_vars,paste0("ghq_",1:28))
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
variables_needed <- c("t","id",process_vars,content_vars,paste0("ghq_",1:28))
LORAr <- select_data(LORA, weeks_needed,variables_needed)

# Scale
LORAr_scale1 <- scale(LORAr[LORAr$t==1,-c(1,2)])
LORAr_scale2 <- scale(LORAr[LORAr$t==7,-c(1,2)])
LORAr_scale3 <- scale(LORAr[LORAr$t==13,-c(1,2)])
LORAr_scale <- data.frame(rbind(LORAr_scale1,LORAr_scale2,LORAr_scale3))
LORAr_scale$t <- LORAr$t
LORAr_scale$id <- LORAr$id

# Get the right variables for the network analysis
LORA_network <- LORAr_scale[,c("id","t",paste0("ghq_",1:28))]


## Profiling - Cross-sectional LPA ####
LORA_mplus <- LORAr[LORAr$t==1,c("PAS_content","PAS_process")]
write.table(LORA_mplus, file = "LORAr.dat", row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE, na = ".")

# LPA with baseline data
LORAr_scale %>%
  filter(t==1) %>%
  dplyr::select(PAS_content, PAS_process) %>%
  single_imputation() %>%
  estimate_profiles(2:6, 
                    variances = c("equal","equal","varying"),
                    covariances = c("equal","zero","zero"),nrep = 5) %>%
  compare_solutions(statistics = c("AIC","BIC"))

# Model 3 with 2 classes best BIC, other lead to very small classes.

LORAr_scale %>%
  filter(t==1) %>%
  dplyr::select(PAS_content, PAS_process) %>%
  single_imputation() %>%
  estimate_profiles(2, 
                    variances = c("equal"),
                    covariances = c("equal"),nrep = 5) %>%
  plot_profiles()

LPA_LORA_profiles <- LORAr_scale %>%
  filter(t==1) %>%
  dplyr::select(PAS_content, PAS_process) %>%
  single_imputation() %>%
  estimate_profiles(2, 
                    variances = c("equal"),
                    covariances = c("equal"),nrep = 5)

# Profile 1 (red): average-high on contents and average-higher on process
# Profile 2 (blue): average-low on contents and low on process


## Profiling - Cross-sectional Step function ####
LORAr_scale[["ghq_sum"]] <- rowSums(LORAr_scale[,paste0("ghq_",1:28)],na.rm = TRUE)
LORA_profiling <- LORAr_scale %>%
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
kml3d(cld3dPregTemp,2:7,nbRedrawing=200,toPlot="nothing")
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
kml(cldLORA_content,2:10,nbRedrawing = 40,toPlot="nothing")
plotAllCriterion(cldLORA_content)
plot(cldLORA_content,2)
plot(cldLORA_content,3)
plot(cldLORA_content,4)
plot(cldLORA_content,5)

## Profiling - Longitudinal Processes ####
cldLORA_process <- cld(LORA_profiling,timeInData =c(3,6,9))
kml(cldLORA_process,2:10,nbRedrawing = 40,toPlot="nothing")
plotAllCriterion(cldLORA_process)
plot(cldLORA_process,2)
plot(cldLORA_process,3)
plot(cldLORA_process,4)
plot(cldLORA_process,5)

## Get clusters commands ####
getClusters(cldLORA_content,4,asInteger = 1)
getClusters(cldLORA_process,4,asInteger = 1)
getClusters(cld3dPregTemp,4,asInteger = 1)

## Functions - Load Functions for Network analysis ##########

getAdjMatList <- function(designMat, data){
  AdjMatList <- NULL
  lambdaList <- NULL
  k <- nrow(designMat)
  
  for (t in 1:(ncol(designMat) -1)){
    
    predictors <- as.matrix(data[, designMat[, t]])
    
    adjMat <- matrix(0, k, k)
    colnames(adjMat) <- designMat[, (t+1)]
    rownames(adjMat) <- colnames(predictors)
    
    lambdaVec <- rep(0,k)
    
    for (i in 1:k){
      
      set.seed(100)
      lassoreg <- cv.glmnet(x = predictors, 
                            y = data[, designMat[i , t+1 ]], 
                            family = "gaussian", alpha = 1, standardize=TRUE)
      
      lambdaVec[i] <- lassoreg$lambda.min
      
      adjMat[1:k,i] <- coef(lassoreg, s = lambdaVec[i], exact = FALSE)[2:(k+1)]
    }
    
    AdjMatList[[t]] <- adjMat
    lambdaList[[t]] <- lambdaVec
  }
  return(list(B = AdjMatList, lambdas = lambdaList))
}

getLavaanSyntax <- function(designMat, model = NULL){
  
  k <- nrow(designMat)
  
  regressions <- ""
  resVariances <- ""
  
  for (t in 1:(ncol(designMat) -1)){
    
    
    for (i in 1:k){
      
      predictors <- designMat[, t]
      
      if(!is.null(model)){
        
        predictors <- predictors[which(model[[t]][,i] != 0 )]
        
      }
      
      #regress variable on variables at the previous time point
      regressions <- paste(regressions, paste(designMat[i,(t+1)], "~", sep = ""), 
                           paste(predictors, collapse = "+"), "\n")
      
      resVariances <- paste(resVariances, paste(designMat[i,(t+1)], "~~", paste(designMat[(i : k), (t+1)], collapse = "+"), sep = ""), "\n")
    }
  }
  
  
  return(c(regressions, resVariances))
}

CreateSigB <- function(B, nonSigParam){
  
  sigB <- B
  
  if(nrow(nonSigParam) == 0){
    return(B)
  }else{
    for(i in 1: nrow(nonSigParam)){

      sigB[nonSigParam$rhs[i], nonSigParam$lhs[i]] <- 0
      
    }
    
    return(sigB)
  }
}

CreateSeparateB <- function(B, designMat){
  
  BList <- NULL
  
  for (t in 1:(ncol(designMat) -1)){
    
    k <- nrow(designMat)
    
    outcomes <- designMat[, (t + 1)]
    
    predictors <- designMat[, t]
    
    BList[[t]] <- B[predictors, outcomes]
    
  }
  
  return(BList)
}

# This function takes a lavaan object and returns lavaan syntax for two nested
# models (specifics about each model is described above).
modelComparisonSyntax <- function(fit, designMat){
  
  B <- lavInspect(fit, what = 'std.all')$beta
  
  k <- nrow(designMat)
  nwaves <- ncol(designMat)
  ncoef <- k^2*(nwaves-1)
  
  ZeroesMat <- matrix(0, k^2, (nwaves-1))
  
  
  # Create a matrix (separated by timepoint ) that tells us which where the 
  # 0s are in the beta matrix
  for(i in 1:(nwaves - 1)){
    ZeroesMat[, i] <- as.vector(t(B)[designMat[,i], designMat[,(i+1)] ]) == 0
  }
  
  #conservative list of zeroes (only if all are zero):
  Zeroes.Cons <- apply(ZeroesMat, 1, sum) >= (nwaves - 1) #set to zero only if all (nwaves-1) paths are zero
  
  # Create a vector with information about constraints 
  constraints <- rep(0, length(Zeroes.Cons))
  nlabels <- length(which(Zeroes.Cons == FALSE))
  paramLabels <- make.unique(rep(letters, length.out = nlabels), sep='')
  
  constraints[which(Zeroes.Cons == FALSE)] <- paramLabels
  constraints[which(Zeroes.Cons == TRUE)] <- 0
  
  # Now let's create the lavaan syntax 
  # Unconstrained Syntax
  regressions <- ""
  resVariances <- ""
  
  for (t in 1:(ncol(designMat) -1)){
    
    constraint.i <- 1
    for (i in 1:k){
      
      predictors <- designMat[, t]
      
      #regress variable on variables at the previous time point
      
      syn <- rep(0, length(predictors))
      
      for(l in 1:length(predictors)){
        
        if(constraints[constraint.i] == "0"){
          syn[l] <-  paste(constraints[constraint.i], "*", predictors[l])
        }else{
          
          syn[l] <-  predictors[l]
          
        }
        constraint.i <- constraint.i + 1
      }
      
      regressions <- paste(regressions, paste(designMat[i,(t+1)], "~", sep = ""), 
                           paste(syn, collapse = " + "), "\n")
      
      resVariances <- paste(resVariances, paste(designMat[i,(t+1)], "~~", paste(designMat[(i : k), (t+1)], collapse = "+"), sep = ""), "\n")
    }
  }
  
  unconsSyntax <- c(regressions, resVariances)
  
  # Constrained syntax
  regressions <- ""
  resVariances <- ""
  
  for (t in 1:(ncol(designMat) -1)){
    
    constraint.i <- 1
    for (i in 1:k){
      
      predictors <- designMat[, t]
      
      #regress variable on variables at the previous time point
      
      syn <- rep(0, length(predictors))
      
      for(l in 1:length(predictors)){
        syn[l] <-  paste(constraints[constraint.i], "*", predictors[l])
        constraint.i <- constraint.i + 1
      }
      
      regressions <- paste(regressions, paste(designMat[i,(t+1)], "~", sep = ""), 
                           paste(syn, collapse = " + "), "\n")
      
      resVariances <- paste(resVariances, paste(designMat[i,(t+1)], "~~", paste(designMat[(i : k), (t+1)], collapse = "+"), sep = ""), "\n")
    }
  }
  
  consSyntax <- c(regressions, resVariances)
  
  return(list(UnconstrainedSyntax = unconsSyntax, ConstrainedSyntax = consSyntax))
  
}


## Functions - Personal function to get the fit for the different groups ####
network_analysis_groups <- function(designMat,data,groups,subscale=FALSE){
  list_fit_cons <- list()
  list_fit_uncons <- list()
  
  # We verify that the number of individuals corresponds in the data and the given groups
  if(length(groups)==nrow(data)){
    n <- length(unique(groups))
    # For each group, we build the 
    for(i in 1:n){
      print(paste0("Starting to work on group ",i," out of ",n,"."))
      
      # 0. Select the lines corresponding to the group
      individuals <- groups==unique(groups)[i]
      sub_data <- data[individuals,]
      
      # 1. Regularized Regression Step -> LASSO.
      glmModel <- getAdjMatList(designMat = designMat, data = sub_data)
      
      # 2. Get the corresponding syntax
      hybridlavaanSyntax_base <- getLavaanSyntax(designMat = designMat, model = glmModel$B)
      model_lines <- unlist(strsplit(hybridlavaanSyntax_base, "\n"))
      model_clean <- model_lines[!grepl("~\\s*$", model_lines)]
      hybridlavaanSyntax <- paste(model_clean, collapse = "\n")
      
      # 3. Estimate SEM model
      print("Estimating the base model")
      hybridModel <- sem(model = hybridlavaanSyntax, data = sub_data, fixed.x = FALSE, missing = "FIML")
      
      # 4. Fix designMat
      missing_items <- setdiff(as.vector(designMat), rownames(t(inspect(hybridModel, "std.all")$beta)))
      designMat_clean <- designMat[
          apply(designMat, 1, function(row) all(!row %in% missing_items)),]
      # 5. Comparison between unconstrained and constrained
      print("Start comparison syntax")
      comparisonSyntax <- modelComparisonSyntax(fit = hybridModel, designMat = designMat_clean)
      
      # 6. Fit both models
      print("Estimating the constrained model")
      fit.cons <- sem(comparisonSyntax$ConstrainedSyntax, data = sub_data, 
                      missing = "FIML", fixed.x = F) 
      print("Estimating the unconstrained model")
      fit.uncons <- sem(comparisonSyntax$UnconstrainedSyntax, data = sub_data, 
                        missing = "FIML", fixed.x = F) 
      
      # 7. Compare the models with Chi square test and other indicators
      print(anova(fit.uncons, fit.cons)) 
      
      # 8. Save the fit in the list
      list_fit_uncons[[i]] <- fit.uncons
      list_fit_cons[[i]] <- fit.cons
      
    }
    return(list(cons=list_fit_cons,uncons=list_fit_uncons))
  }
  else{
    print("The length of the group doesn't correspond to the data.")
  }
}

get_final_matrix <- function(model,designMat,pvalue=0.05){
  allRegressions <- parameterestimates(model,  standardized = TRUE)[grep("~~", parameterestimates(model)$op, invert = TRUE),]
  allRegressions <- allRegressions[which(allRegressions$op != "~1"),]
  
  # Initialize to 0
  row_vars <- unique(designMat[, 1])
  col_vars <- unique(designMat[, 2])
  result_mat <- matrix(0, nrow = length(row_vars), ncol = length(col_vars),
                       dimnames = list(row_vars, col_vars))
  
  # Loop on designMatItems
  for(i in 1:nrow(designMat)) {
    rhs_var <- designMat[i, 1]
    lhs_var <- designMat[i, 2]
    
    # Look for corresponding line in allRegressions
    match_row <- allRegressions %>%
      filter(lhs == lhs_var, rhs == rhs_var, pvalue <= pvalue, est != 0)
    
    # And add it to the matrix
    if(nrow(match_row) == 1) {
      result_mat[rhs_var, lhs_var] <- match_row$est
    }
  }
  return(result_mat)
}

## Network analysis - ITEMS - Choose Nodes and build design matrix ####
# In this section we apply goldbricker to select which nodes to colides for
# the network analysis on the items-level

# 1. Apply Goldbricker on the whole dataset
goldbricker(LORA_network[,-c(1,2)])
# Res : ghq_16 & ghq_15, ghq_27 & ghq_25, ghq_28 & ghq_25, ghq_18 & ghq_17 

# 2. Per time point
goldbricker(LORA_network[LORA_network$t==1,-c(1,2)])
goldbricker(LORA_network[LORA_network$t==7,-c(1,2)])
goldbricker(LORA_network[LORA_network$t==13,-c(1,2)])
# Res : 14&11 and 16&15

# We choose to colide 14&11 and 16&15
LORA_network_items <- LORA_network[,c("id","t",paste0("ghq_",c(1:10,12,13,17:28)))]
LORA_network_items$ghq_1615 <- LORA_network$ghq_15+LORA_network$ghq_16/2
LORA_network_items$ghq_1411 <- LORA_network$ghq_14+LORA_network$ghq_11/2


# Build the wide dataframe
LORA_wide_items <- LORA_network_items %>%
  pivot_wider(
    id_cols = id,
    names_from = t,
    values_from = starts_with("ghq_"),
    names_sep = "."
  )
LORA_wide_items <- as.data.frame(LORA_wide_items[,-c(1)])

# Design Matrix
designMatItems <- matrix(colnames(LORA_wide_items),ncol=3,nrow=26,byrow=TRUE)



## Network analysis - ITEMS - Results ####
# KML3D
networks_kml3D_2 <- network_analysis_groups(designMatItems,LORA_wide_items,getClusters(cld3dPregTemp,2,asInteger = 1))

# KML Content
networks_content_2 <- network_analysis_groups(designMatItems,LORA_wide_items,getClusters(cldLORA_content,2,asInteger = 1))
networks_content_3 <-network_analysis_groups(designMatItems,LORA_wide_items,getClusters(cldLORA_content,3,asInteger = 1))

# KML Process
networks_process_2 <- network_analysis_groups(designMatItems,LORA_wide_items,getClusters(cldLORA_process,2,asInteger = 1))
networks_process_3 <- network_analysis_groups(designMatItems,LORA_wide_items,getClusters(cldLORA_process,3,asInteger = 1))

# Get the final matrix ?
model <- networks_kml3D_2$cons[[1]]
model <- networks_kml3D_2$cons[[2]]

model <- networks_content_2$cons[[1]]
model <- networks_content_2$cons[[2]]

model <- networks_process_2$cons[[1]]
model <- networks_process_2$cons[[2]]


sum(get_final_matrix(model,designMatItems)!=0)
qgraph(get_final_matrix(model,designMatItems))


## Network analysis - SUBSCALES - Build subscale and design matrix ####
LORA_network_subscale <- data.frame(id=LORA_network$id,t=LORA_network$t)
LORA_network_subscale[["somatic"]] <- rowSums(LORA_network[,paste0("ghq_",1:7)],na.rm = TRUE)
LORA_network_subscale[["anxiety"]] <- rowSums(LORA_network[,paste0("ghq_",8:14)],na.rm = TRUE)
LORA_network_subscale[["social_dysfunction"]] <- rowSums(LORA_network[,paste0("ghq_",15:21)],na.rm = TRUE)
LORA_network_subscale[["depression"]] <- rowSums(LORA_network[,paste0("ghq_",22:28)],na.rm = TRUE)

LORA_wide_subscale <- LORA_network_subscale %>%
  pivot_wider(
    id_cols = id,
    names_from = t,
    values_from = c("somatic","anxiety","social_dysfunction","depression"),
    names_sep = "."
  )
LORA_wide_subscale <- as.data.frame(LORA_wide_subscale[,-c(1)])

# Design Matrix
designMatSubscale <- matrix(colnames(LORA_wide_subscale),ncol=3,nrow=4,byrow=TRUE)

## Network analysis - SUBSCALES - Results ####
data <- LORA_wide_subscale
groups <- getClusters(cldLORA_process,2,asInteger = 1)
i <- 2
designMat <- designMatSubscale

individuals <- groups==unique(groups)[i]
sub_data <- data[individuals,]

# 1. Regularized Regression Step -> LASSO.
glmModel <- getAdjMatList(designMat = designMat, data = sub_data)

# 2. Get the corresponding syntax
hybridlavaanSyntax_base <- getLavaanSyntax(designMat = designMat, model = glmModel$B)
model_lines <- unlist(strsplit(hybridlavaanSyntax_base, "\n"))
model_clean <- model_lines[!grepl("~\\s*$", model_lines)]
hybridlavaanSyntax <- paste(model_clean, collapse = "\n")

# 3. Estimate SEM model
hybridModel <- sem(model = hybridlavaanSyntax, data = sub_data, fixed.x = FALSE, missing = "FIML")
qgraph(get_final_matrix(hybridModel,designMat))