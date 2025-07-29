GenerateBeta <- function(p, sparsity, weightBank){
  
  varNames <- LETTERS[1:p]
  
  # Number of paths to be set to a non-zero value
  numbNonzeroPaths <- round(p * p * sparsity, 0)
  
  # Create a list of all possible paths to randomly sample from this bank 
  # later to select with paths to set to a non-zero value
  pathBank <- expand.grid(varNames, varNames) # col 1 = predictor
  # col 2 = outcome
  
  B <- matrix(0, p, p)
  colnames(B) <- rownames(B) <- varNames
  
  
  # Populate B matrix with non-zero paths 
  pathIndex <- sample(1:nrow(pathBank), numbNonzeroPaths, replace = F)
  
  
  for(i in 1:numbNonzeroPaths){
    B[pathBank[pathIndex[i], 1], 
      pathBank[pathIndex[i], 2]] <- sample(weightBank, 1) }
  
  return(B)
  
}



GenerateNonstationaryCov <- function(BList, ResidualBank = NULL){
  
  p <- ncol(BList[[1]])
  timepoints <- length(BList) + 1
  
  saveMatrices <- vector(mode='list', length = timepoints) # we save the within time cov matrices here
  
  
  LongitudinalSigma <- matrix(0, p * timepoints, p * timepoints)
  rownames(LongitudinalSigma) <- colnames(LongitudinalSigma) <- paste0(LETTERS[1:p], 
                                                                       rep(1:timepoints, 
                                                                           each = p))
  # Create starting covariance matrix 
  start <- diag(p)
  
  BListPlus <- c(list(BList[[1]], BList[[1]]), BList)
  
  for (j in 1:(timepoints + 1)){
    
    # covariance matrix from j =1 and 2 are just a set-up covariance matrix
    # so that none of our covariance matrices have covariances of 0 
    # between the variables
    # we won't be saving this covariance matrix  
    
    psi <- diag(diag(1 - t(BListPlus[[j]]) %*% start %*% BListPlus[[j]]), p) # Creates a psi matrix
    # so the variance = 1
    
    if(is.null(ResidualBank) == FALSE){ 
      
      # Select residual covariances 
      residualCovariances <- sample(ResidualBank, (p * (p-1)/2), replace = TRUE)
      
      # Add residual covariances to psi matrix
      psi[upper.tri(psi)] <- residualCovariances
      psi[lower.tri(psi)] <- t(psi)[lower.tri(psi)]
      
    }
    
    if( sum(diag(psi) < 0) == 0){ # Only does the following steps if the psi values are
      # positive 
      
      save <- t(BListPlus[[j]]) %*% start %*% BListPlus[[j]] + psi
      start <- save
      
    }else{return("Negative Residual Variances")}
    
    if(j != 1){
      
      crossCovCode <- "save" # Equation for the covariances between time 
      # T variables and time T+1 variables
      
      if(matrixcalc::is.positive.definite(round(save, 6)) == TRUE){
        
        saveMatrices[[j - 1]] <- save
        
        
        LongitudinalSigma[grep((j - 1), colnames(LongitudinalSigma)), 
                          grep((j - 1), colnames(LongitudinalSigma))] <- save
        
        if(j != (timepoints + 1)){
          for(i in j:timepoints){
            
            # Calculate and save covariances between variables in different time points
            
            crossCovCode <- paste0(crossCovCode, "%*% BListPlus[[", (i + 1), "]]")
            
            LongitudinalSigma[grep((j - 1), colnames(LongitudinalSigma)), grep((i), colnames(LongitudinalSigma))]  <- eval(parse(text = crossCovCode))
            LongitudinalSigma[grep((i), colnames(LongitudinalSigma)), grep((j - 1), colnames(LongitudinalSigma))]  <- t(eval(parse(text = crossCovCode)))
            
          }
          
        }
        
        
      }else{
        saveMatrices[[j - 1]] <- "Not Positive Definite"
      }
    }
    
    
    
    
  } # end j loop
  
  
  
  return(list(SigmaList = saveMatrices,
              LongCov = LongitudinalSigma))
  
}

#function fixed 1/9/2025
CreateLavaanSyntax <- function(longVarNames, timepoints, model = NULL){
  regressions <- NULL
  covs <- NULL
  
  for(l in 1:(timepoints-1)){

    out <- paste(longVarNames[grep(as.character(l + 1), longVarNames)], "~")
    cov <- paste(longVarNames[grep(as.character(l + 1), longVarNames)], sep = "+",
                 collapse = "+ ")
    covs <- c(covs, paste0(out, "~", cov))
        
    if(is.null(model)){
      pred <- paste(longVarNames[grep(as.character(l), longVarNames)], sep = "+",
                    collapse = "+ ")
      regressions <- c(regressions, paste(out, pred))
      
      
    }else{
      
      p <- length(longVarNames) / timepoints
      
      for(i in 1:p){
        predictors <- longVarNames[grep(as.character(l), longVarNames)]
        pruned.predictors <- predictors[which(model[[l]][,i ] != 0)]
        
        if(length(pruned.predictors) != 0){
          pruned.predictors.syntax <- paste(pruned.predictors, sep = "+", collapse = "+")
          out <- paste(longVarNames[grep(as.character(l + 1), longVarNames)][i], "~")
          regressions <- c(regressions, paste(out, pruned.predictors.syntax)) 
        }
        
        
      }}
  }
  
  return(c(regressions, covs))
  
}


CreateSigB <- function(B, nonSigParam, alpha){
  
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

CalculatePerformanceCounts <- function(true, est){
  # arguments: 
  # true = True/Population Adjacency Matrix
  # est = Estimated matrix
  
  true <- as.matrix(true)
  est <- as.matrix(est)
  
  #Convert partial correlation matrix
  #to adjacency matrix 
  estSig <- as.vector(ifelse(est == 0, 0, 1))
  trueSig <- as.vector(ifelse(true == 0, 0, 1))
  
  
  #Count number of true/false positive and negatives
  #in estimated matrix 
  
  TP <- sum(trueSig == 1   & estSig == 1)
  TN <- sum(trueSig == 0   & estSig == 0)
  FP <- sum(trueSig == 0   & estSig == 1)
  FN <- sum(trueSig == 1   & estSig == 0)
  
  
  
  return(c(TP = TP, TN = TN, FP = FP, FN = FN))
}

CalculatePerformanceMeasures <- function(PerformanceCounts){
  
  Spec <- PerformanceCounts["TN"]/ (PerformanceCounts["TN"] + PerformanceCounts["FP"])
  
  FPR <- 1 - Spec
  
  Sens <- PerformanceCounts["TP"]/ (PerformanceCounts["TP"] + PerformanceCounts["FN"])
  
  TDR <- PerformanceCounts["TP"]/ (PerformanceCounts["TP"] + PerformanceCounts["FP"])
  
  
  return(c(Spec = Spec, FPR= FPR, Sens = Sens, TDR = TDR))
  
}
