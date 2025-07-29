################################################################################
# Estimating a CLPN Network with an Hybrid Approach
# Code to accompany Wysocki, van Bork, Cramer & Rhemtulla 
# "Cross-Lagged Network Models"
################################################################################


### Conceptual Overview ########################################################

# We estimate the initial CLPN using lasso regression (a penalized regression
# approach that applies an l_1 penalty on the estimated regression coefficients)
# The penalty is selected using cross-validation. 
# This estimation technique has the effect of shrinking small 
# regression paths to exactly zero, thus achieving a sparse solution.
# Lasso regression is particularly useful in high-dimensional settings. 
# However lasso regression also biases the non-zero edges towards zero. 
# Additionally, the cross-validation approach to selecting the penalty
# has a high false positive rate. 

# Given the high false positive rate and the biased estimates, we include a 
# second estimation step that further prunes the model and provides
# non-regularized estimates for interpretation. 

# In the second step, the paths that are estimated to be zero in the regularized 
# regression results are fixed to zero. And the model is re-estimated as a 
# regular (non-regularized) regression model. We can look at the path 
# significance values and only interpret the significant paths
################################################################################

### Code Overview ##############################################################

# We use the glmnet R package to estimate the initial network structure (via 
# regularized regression). Then we estimate the model selected by glm within an
# SEM framework (non-regularized model) using the R package lavaan. Finally, 
# we can further prune the model (but do not re-estimate) 
# using the parameter significance from the SEM model 
################################################################################



### Set-Up #####################################################################


### Install and load required package(s) ########

#install.packages("glmnet")#uncomment to install packages for the 1st time
#install.packages("lavaan")

library(glmnet)
library(lavaan)


### Load Functions ###################

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


### Load Data ###########################
# If you're following along with the paper example, load the data file
# we cleaned and saved in the CleanData.R script

load(file = "mydata.Rda") 

#remove cases with missing data (glmnet doesn't allow them)
myfulldata <- mydata[complete.cases(mydata),]

### Create the inputs ########################


# Create a design matrix with the name of each variables separated by
# time point
# each column contains the set of variable names for a different timepoint

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

################################################################################


### Regularized Regression Step ################################################


# getAdjMatList is a function that uses glm regression to get a beta matrix for 
# t waves of data: This function takes the data, uses cv.glmnet to run 
# regularized regressions on each node at time t on all nodes at (t-1),
# saves (t-1) adjacency matrices (with coefficients and zeroes)

glmModel <- getAdjMatList(designMat = designMat , data = myfulldata)

glmModel[["B"]] # The edge weight matrices for the t-1 networks

################################################################################


### NonRegularized Regression Step #############################################


# The function (getLavaanSyntax()) takes the model selected by the glm
# regression and writes lavaan syntax that corresponds to that selected model

hybridlavaanSyntax <- getLavaanSyntax(designMat = designMat, model = glmModel$B)

# Estimate SEM model
hybridModel <- sem(model = hybridlavaanSyntax, data = mydata,
                   fixed.x = F, missing = 'FIML')

# Look at fit statistics
fitStats <- fitMeasures(hybridModel)

# Look at estimated regression parameters
allRegressions <- parameterestimates(hybridModel, 
                                     standardized = TRUE)[grep("~~", parameterestimates(hybridModel)$op,
                                                       invert = TRUE),]
allRegressions <- allRegressions[which(allRegressions$op != "~1"),]


hybridBStan <- t(inspect(hybridModel, "std.all")$beta)
hybridBUnstan <- t(inspect(hybridModel, "est")$beta)

################################################################################


# Further Prune Model Based on Significance ####################################

# The beta matrix we extracted from the lavaan object has the values for both 
# the significant and non-significant paths. If we only want the significant
# paths then we need to do some further pruning. 
# The function--CreateSigB()--does this. 
hybridSigBStan <- CreateSigB(hybridBStan, 
                             nonSigParam = allRegressions[allRegressions$pvalue > .05,])

# Also the beta matrix from the lavaan model has the parameter estimates for 
# all the timepoints in a single matrix. If we want to visualize separate networks
# for each timepoint then we need to split up the beta matrix based on the 
# timepoints
# The function--CreateSeparateB()--does this. 

hybridBList <- CreateSeparateB(hybridSigBStan, designMat)

################################################################################


### Save Objects ###############################################################

hybridModelList <- list(hybridBMatrixStandardized = hybridBStan,
                        hybridBMatrixUnstandardized = hybridBUnstan,
                        hybridBMatrixSep =hybridBList, 
                        fitStatistics = fitStats, 
                        lavaanObject = hybridModel, 
                        glmObject = glmModel)


save(hybridModelList, 
     file = "ModelList.Rda")
