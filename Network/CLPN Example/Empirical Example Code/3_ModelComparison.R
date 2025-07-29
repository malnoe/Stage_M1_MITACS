################################################################################
# Model Comparisons for Cross-time Constraints with a CLPN Network
# Code to accompany Wysocki, van Bork, Cramer, & Rhemtulla 
# "Cross-Lagged Network Models"
################################################################################


### Conceptual Overview ########################################################

# When more than two measurement occasions are available to model, it may be of 
# interest to investigate whether the predictive relations across successive 
# occasions are equal (e.g., if the paths from 7th to 8th grade are equivalent 
# to those from 8th to 9th grade). 

# This file contains code to fit two nested CLPNs within a SEM framework and 
# compare their fit
# Model 1: A CLPN model that is unconstrained (path weights are allowed to vary) 
#          across timepoints 
# Model 2: A CLPN with the same paths estimates as model 1, but with cross-time
#          constraints imposed (e.g,.CLxy_T1-T2 == CLxyT2-T3, etc.)



# Zeroes need to be imposed consistently across each measurement occasion for
# both of these models (i.e., the same set of effects need to be estimated 
# at T1-T2 and T2-T3 etc) to ensure that they are nested.  
# So we need some kind of decision rule to go from our hybrid model (where 
# there are likely different paths estimated at each time point) to our 
# unconstrained model (model 1). 

# In the paper and the following code, we use a conservative decision rule 
# where a path is constrained to 0 only if it is estimated to be zero across
# every neighboring pair of occasions in the hybrid model. 
################################################################################


### Code Overview ##############################################################

# This code uses the hybrid model estimated from the EstimationHybridCLPN code
# to create two nested models (the unconstrained [model 1] and constrained
# [model 2] models) that can then be compared to assess whether cross-time 
# constraints significantly worsens model fit

# The following code includes a function that will write the syntax for the 
# nested model based on another lavaan object (we use the hybrid model that
# was estimated in the previous step). 

# We then estimate the two models in lavaan and compare them using a chi
# square test. 
################################################################################ 



### Set-Up #####################################################################


# Install and load required package(s) ########
#install.packages("lavaan") #uncomment to install packages for the 1st time
library(lavaan)


### Load function(s) ############################

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


### Load hybrid model and data ##################

load(file = "ModelList.Rda") # Load model that you estimated with the EstimationHybridCLPN code
                                 # this code loads the object hybridModelList
                                 # which contains information about the CLPN 
                                 # we estimated in the previous step

load(file = "mydata.Rda") # Load your data


### Other Inputs #######################

# Create a design matrix with the name of each variable separated by the time
#           point i.e., in different columns

# For example if I only had 4 variables across 3 timepoints I would 
# specify the following design matrix
# designMat <- matrix(c('SC01.1','SC02.1', 'SC03.1', 'SE04.1',
#                      'SC01.2','SC02.2', 'SC03.2', 'SE04.2',
#                      'SC01.3','SC02.3', 'SC03.3', 'SE04.3'),
#                                 nrow = 4, ncol = 3, 
#                                 byrow = FALSE)

designMat <- matrix(colnames(mydata), 
                    nrow = 17, ncol = 4)

# I'm basing the nested model on the hybrid network I estimated with the 
# EstimationHybridCLPN code 
# So I need to extract that model
fit <- hybridModelList$lavaanObject


### Fitting and Comparing the Nested Models ####################################


### Create syntax for each model ###########
comparisonSyntax <- modelComparisonSyntax(fit = fit, designMat = designMat)

comparisonSyntax$UnconstrainedSyntax # Model 1
comparisonSyntax$ConstrainedSyntax # Model 2

### Estimate each model #####################
fit.cons <- sem(comparisonSyntax$ConstrainedSyntax, data = mydata, 
                missing = "FIML", fixed.x = F) 

fit.uncons <- sem(comparisonSyntax$UnconstrainedSyntax, data = mydata, 
                  missing = "FIML", fixed.x = F) 


### Save and compare the models ##################

summary(fit.cons)
summary(fit.uncons)

# Compare the models
anova(fit.uncons, fit.cons) # constrained model doesn't significantly worsen model fit 

unconsFitStats <- fitMeasures(fit.uncons)
consFitStats <- fitMeasures(fit.cons)

modelSave <- list("Unconstrained" = fit.uncons, "Constrained" = fit.cons,
                  "UnconstrainedFitStats = unconsFitStats",
                  "ConstrainedFitStats = consFitStats")
save(modelSave, file = "ModelComparisonFit.Rdat")
