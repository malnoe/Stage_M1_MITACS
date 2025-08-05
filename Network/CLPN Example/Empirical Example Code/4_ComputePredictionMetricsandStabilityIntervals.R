################################################################################
# Compute Node-wise Prediction Metrics and Stability Intervals 
# Code to accompany Wysocki, van Bork, Cramer, & Rhemtulla
# "Cross-Lagged Network Models"
################################################################################


### Conceptual Overview ########################################################

# A CLPN can be summarized into two measures of variable centrality: 
# In-prediction is the extent to which each variable is predicted by other 
# variables in the network, and out-prediction is the extent to which each 
# variable predicts other variables in the network. 

# Researchers may be particularly interested in two more specific in-prediction 
# measures: cross-lagged in-prediction is the proportion of variance accounted 
# for by all other variables at the previous measurement occasion, that is, 
# cross-lagged in-prediction excludes the auto-regressive path,
# within-construct in-prediction is the proportion of variance accounted for by 
# all the variables within the same construct, excluding the auto-regressive 
# effect, and cross-construct in-prediction is the proportion of variance 
# accounted for by all variables that belong to a different construct at the 
# previous measurement occasion. 

# Knowing the stability of an estimate (e.g., its fluctuation in response to 
# sampling variability) is important for calibrating interpretations. 
# We can use bootstrapping to obtain stability intervals around each of the 
# in- and out-prediction estimates. 

################################################################################


### Code Overview ##############################################################

# This code estimates prediction metrics for each node across each timepoint. 
# The prediction metrics capture how well a variable predicts other variables 
# and how well a variable is predicted by other variables. 

# The following code includes a function that will calculate the types of 
# in-prediction and out-prediction (total, cross-lagged, within-construct 
# and cross-construct). 

# We provide a function that uses bootstrapping to get stability intervals
# around each of the prediction metrics 

################################################################################


### Set-Up #####################################################################


# Install and load required package(s) ########

#install.packages("lavaan") #uncomment to install packages for the 1st time
#install.packages("dplyr")
#install.packages("ggplot2")
library(lavaan)
library(dplyr)
library(ggplot2)


### Load function(s) ############################

# This function takes a dataset, a Beta matrix, a designMat, and a vector 
# with group names and calculates the in and out prediction metrics 

InOutPredict <- function(data, B, designMat, groups){
  
  k <- nrow(designMat)
  predictionDf <- NULL
  
  for(t in 1:(ncol(designMat) - 1)){ # Loops over timepoints
    
    inpred0 <- rep(0, k) #set up empty vectors to hold predictability coefficients
    inpred1 <- rep(0, k) 
    inpred2 <- rep(0, k)
    inpred3 <- rep(0, k)
    
    # Calculate In Prediction
    for(i in 1:k){ 
      
      name.i <- designMat[i, (t + 1)] #what is the name of the DV
      group.i <- groups[i] #what construct does the DV belong to
      
      # Subset predictors for each of the in prediction metrics
      include.mod0 <- designMat[, t ] # include self at previous time point
      include.mod1 <- designMat[, t][-i] # exclude self at previous time point
      include.mod2 <- designMat[, t ][groups != group.i] #exclude predictors belonging to same construct
      include.mod3 <- designMat[, t ][groups == group.i] #exclude predictors belonging to different construct
      
      
      betas.mod0 <- B[name.i, include.mod0] #fix regression coefficients to the ones obtained by lavaan
      betas.mod1 <- B[name.i, include.mod1]
      betas.mod2 <- B[name.i, include.mod2] 
      betas.mod3 <- B[name.i, include.mod3]
      
      denom <- var(data[, name.i], na.rm = TRUE) # variance of outcome
      num0 <- var(as.matrix(data[,include.mod0]) %*% betas.mod0, na.rm = T)
      num1 <- var(as.matrix(data[,include.mod1]) %*% betas.mod1, na.rm = T) 
      num2 <- var(as.matrix(data[,include.mod2]) %*% betas.mod2, na.rm = T)
      num3 <- var(as.matrix(data[,include.mod3]) %*% betas.mod3, na.rm = T)
      
      inpred0[i] <- num0/denom 
      inpred1[i] <- num1/denom  
      inpred2[i] <- num2/denom  
      inpred3[i] <- num3/denom 
      
    } # End in-prediction loop
    
    
    # Calculate out prediction 
    outpred0 <- rep(0, k) 
    outpred1 <- rep(0, k) 
    outpred2 <- rep(0, k)
    outpred3 <- rep(0, k)
    
    for (i in 1:k){  
      
      name.i <- designMat[i, t ] #what is the name of the IV
      group.i <- groups[i] #what construct does the IV belong to
      
      include.mod0 <- designMat[, (t+1)] # include self at previous time point
      include.mod1 <- designMat[, (t+1)][-i] #exclude self at previous time point
      include.mod2 <- designMat[, (t+1)][groups != group.i] #exclude outcomes belonging to same construct
      include.mod3 <- designMat[, (t+1)][groups == group.i] 
      
      betas.mod0 <- B[include.mod0, name.i]
      betas.mod1 <- B[include.mod1, name.i] 
      betas.mod2 <- B[include.mod2, name.i]
      betas.mod3 <- B[include.mod3, name.i]
      
      
      yvars.mod0 <- apply(as.data.frame(data[, include.mod0]), 2, var, na.rm =TRUE)
      yvars.mod1 <- apply(as.data.frame(data[, include.mod1]), 2, var , na.rm =TRUE)
      yvars.mod2 <- apply(as.data.frame(data[, include.mod2]), 2, var, na.rm =TRUE)
      yvars.mod3 <- apply(as.data.frame(data[, include.mod3]), 2, var, na.rm =TRUE)
      
      xvar <- var(data[,name.i], na.rm = TRUE)
      
      rsq.mod0 <- (betas.mod0^2 * xvar)/yvars.mod0
      rsq.mod1 <- (betas.mod1^2 * xvar)/yvars.mod1
      rsq.mod2 <- (betas.mod2^2 * xvar)/yvars.mod2
      rsq.mod3 <- (betas.mod3^2 * xvar)/yvars.mod3
      
      outpred0[i] <- mean(rsq.mod0)
      outpred1[i] <- mean(rsq.mod1) 
      outpred2[i] <- mean(rsq.mod2)
      outpred3[i] <- mean(rsq.mod3)
      
    } #out prediction loop
    
    
    predictionDf <- rbind(predictionDf,
                          data.frame(t = t, 
                                     PredictionIndex = rep(c("InPred", "CrossLagInPred", "OtherConstructInPred", "SameConstructInPred",
                                                             "OutPred", "CrossLagOutPred", "OtherConstructOutPred", "SameConstructOutPred"), each = k), 
                                     PredictionValue = c(inpred0, inpred1, inpred2, inpred3,
                                                         outpred0, outpred1, outpred2, outpred3), 
                                     Variable = c(rep(designMat[, (t+1)], times = 4) , rep(designMat[, t], times = 4))))
    
  } # end t for loop
  
  
  return(predictionDf)
}


# This function bootstraps the data and calculates prediction metrics from each of the 
# bootstrapped datasets
bootstrapPredictionMetrics <- function(model, designMat, groups, samples = 10, parallel = "no"){
  
  # Bootstrap a lavaan model
  # Save the beta matrix from each sample
  bootRes <- bootstrapLavaan(model, R = samples, 
                             FUN = function(x){
                               require(lavaan)
                               as.vector(lavInspect(x, what = "est")$beta)}, 
                             keep.idx = TRUE, parallel = parallel)
  
  # the bootstrapLavaan() function returns a object where each row represents 
  # the beta values for a different sample in vector form
  # I need to re-create the B matrix
  
  BMatrices <- vector(mode = "list", length = nrow(bootRes)) 
  p <- nrow(inspect(model, "sampstat")$cov)
  
  for(i in 1:nrow(bootRes)){
    
    bootB <- matrix(bootRes[i, ],
                    nrow = p, ncol = p)
    
    colnames(bootB) <- colnames(lavInspect(model, what = "est")$beta)
    rownames(bootB) <- rownames(lavInspect(model, what = "est")$beta)
    
    BMatrices[[i]] <- bootB
    
  }
  
  k <- nrow(designMat)
  predictionDf <- NULL
  
  # Now I can calculate the in and out prediction metrics for each
  # timepoint and each bootstrapped sample
  
  for(b in 1:samples){ # Loops over the bootstrapped samples
    
    # subsetting the data rows that was used for this sample
    bootRows <- attr(bootRes, 'boot.idx')[[1]][b,]
    bootData <- lavInspect(model, what = "data")[bootRows,]
    
    # Subsetting the B matrix estimated from this sample
    bootB <- BMatrices[[b]]
    
    predictionDf <- rbind(predictionDf, 
                          cbind(b = b, InOutPredict(data = bootData, B = bootB, designMat = designMat, groups = groups))
    )
    
  } # end b for loop
  
  return(predictionDf)
}


### Load model fit and data ##################
# Load model/fit object you would like to calculate prediction metrics for
#load(file = "C:/Users/annawy/Box/Research/Davis/In Progress/CLPN/R/Applied Example Code/Fit Objects/modelComparisons.Rdata") # loads object modelComparisonFits

#load(file = "C:/Users/annawy/Box/Research/Davis/In Progress/CLPN/Data/EmpiricalExampleData.Rda") # loads object mydata and object mylegend




#### Calculate Prediction Indices ##############################################



### Overview of InOutPredict function ###############

# The function to calculate prediction metrics has 4 arguments: 
# InOutPredict(data, B, designMat, groups)

# data:      The data 
# B:         The estimated beta/regression matrix from the CLPN (i.e., the estimated weights)
#            The beta matrix should have column and row names
# designMat: This is a matrix with the names of each of the variables separated into 
#            different columns based on timepoints (see below for an example)
# groups:    A character vector with the group information for each variable. 
#            Used to calculate the within and between construct prediction metrics


### Creating the objects for each of these arguments ########

# 1. data argument
head(mydata) # I'm using the non-standardized data

# 2. B argument
consModel <- modelComparisonFits[[1]] # The final model in our applied example 
                                      # is the constrained model estimated in 
                                      # in lavaan. 
                                      # consModel is the lavaan object of that 
                                      # estimated model

B <- lavInspect(consModel, what = 'est')$beta # extract the estimated 
                                              # beta matrix

# 3. designMat argument  
designMat <- matrix(colnames(mydata),
                    nrow = 17, ncol = 4, 
                    byrow = FALSE)

designMat # Column 1 has the names of all the time 1 variables
          # Column 2 has the names of all the time 2 variables etc...
          # If we had a 5th timepoint, the designMat object would have a 
          # 5th column. 
          # If we had more variables in each time point, designMat would have 
          # more rows. 
          # If the variable names in the designMat don't match the variable names
          # in the data set and the beta matrix, the functions won't work!

# 4. groups argument
# A vector identifying which group each variable belongs to 
# Variables Needs to be in the same order as the beta matrix & data input
groups <- c(rep("Commitment to School", 7), rep("Self-Esteem", 10))



### Get prediction metrics ####################

predictionEstimates  <- InOutPredict(data = mydata, B = B, 
                                     designMat = designMat, groups = groups)

head(predictionEstimates)
# column t: tells you the time interval. For example t = 1 is looking at the 
#           in prediction of variables in timepoint 2 and the out prediction of 
#           time 1 variables. t= 2 is looking at the in prediction of time 3
#           variables and the out prediction of time 2 variables etc. 
# column PredictionIndex: tells you which prediction metric is depicted in that 
#                         row (options: In Prediction, CrossLag In Prediction, 
#                         Other Construct In Prediction, Same Construct In 
#                         Prediction, Out Prediction, CrossLag Out Prediciton, 
#                         Other Construct Out Prediction, Same Construct Out 
#                         Prediction)
# column: PredictionValue: The value for the prediction metric 
# Variable: The variable name 

# For example: 
predictionEstimates[1, ] # row 1 tells us that the in prediction value for SC01 at
                         # time 2 (time interval 1) is .36



### Compute Stability Intervals ################################################


### Overview of bootstrapPredictionMetrics function ##########

# The function to bootstrap prediction metrics has 5 arguments: 
# bootstrapPredictionMetrics(model, designMat, groups, samples = 10,
#                            parallel = "no")

# model:     A lavaan fit object (model estimated in lavaan)
# designMat: This is a matrix with the names of each of the variables separated into 
#            different columns based on timepoints (see below for an example)
# groups:    A character vector with the group information for each variable. 
#            Used to calculate the within and between construct prediction metrics
# samples:   the number of samples to bootstrap from the data. Default is 10. 
# parallel:  Do you want to run the bootstrap operation on multiple computer cores?
#            Default is "no". Other options are "multicore" and "snow"


## Creating the objects for each of these arguments
# 1. model argument
consModel # The final model in our applied example

# 2. designMat argument (we created this earlier)

designMat # Column 1 has the names of all the time 1 variables
# Column 2 has the names of all the time 2 variables etc...
# If we had a 5th timepoint, the designMat object would have a 
# 5th column. 
# If we had more variables in each time point, designMat would have 
# more rows. 
# If the variable names in the designMat don't match the variable names
# in the data set and the beta matrix, the functions won't work!

# 3. groups argument (we also created this earlier)
groups 

# 4. samples argument (I use 800 samples in the paper but 
#                       you can use a different number)

# 5. parallel arugment (bootstrapping is computationally intensive. I recommend 
                        # using parallel computing )




# Bootstrap Prediction Indices #####################################
predictionDF <- bootstrapPredictionMetrics(model = consModel, 
                           designMat = designMat, groups = groups, 
                           samples = 5, parallel = "multicore") # run 5 samples just to test


head(predictionDF)

# column b: tells you which bootstrapped sample the prediction metric in that 
#           row is calculated from. If you set the samples argument to 100
#           b should be from 1-100
# column t: tells you the time interval. For example t = 1 is looking at the 
#           in prediction of variables in timepoint 2 and the out prediction of 
#           time 1 variables. t= 2 is looking at the in prediction of time 3
#           variables and the out prediction of time 2 variables etc. 
# column PredictionIndex: tells you which prediction metric is depicted in that 
#                         row (options: In Prediction, CrossLag In Prediction, 
#                         Other Construct In Prediction, Same Construct In 
#                         Prediction, Out Prediction, CrossLag Out Prediciton, 
#                         Other Construct Out Prediction, Same Construct Out 
#                         Prediction)
# column: PredictionValue: The value for the prediction metric 
# Variable: The variable name 

# For example: 
predictionDF[1, ] # row 1 tells us that for the first boot strapped sample
                  # the in prediction value for SC01 at time 2 (time interval 1) is
                  # .33


# Calculate Stability Intervals #############################

# Bootstrapping means we have a set of values for each prediction metric 
# at each timepoint. We can create stability intervals around each prediction metric
# by calculating the standard deviation of the bootstrapped prediction metrics for each
# variable, and then +/- the standard deviation from the prediction value calculated
# from the original dataset

predictionDF %>%
  group_by(PredictionIndex, Variable) %>% 
  summarise(mu = mean(PredictionValue), s = sd(PredictionValue),
            max = max(PredictionValue), min = min(PredictionValue)) -> bootPredSummaryStat








