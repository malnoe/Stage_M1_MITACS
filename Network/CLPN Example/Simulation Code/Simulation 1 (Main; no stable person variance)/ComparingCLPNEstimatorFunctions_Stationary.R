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

GenerateStationaryCov <- function(B, ResCov = NULL){
  
  saveMatrices <- vector(mode='list', length = 1000)
  
  p <- ncol(B)
  
  # Create starting covariance matrix 
  start <- diag(p)
  
  
  if(is.null(ResCov) == TRUE){
    
    for (j in 1:1000){
      psi <- diag(diag(1 - t(B) %*% start %*% B), p) # Creates a psi matrix
                                                     # so the variance = 1
      
      
      if( sum(psi < 0) == 0){ # Only does the following steps if the psi values are
                             # positive 
        
        save <- t(B) %*% start %*% B + psi
        saveMatrices[[j]] <- save
        start <- save
        
      }else{return("Negative Residual Variances")}
      
    } #end of j for loop
    
    
  }else{
    
    for (j in 1:1000){
      psi <- diag(diag(1 - t(B) %*% start %*% B), p) # Creates a psi matrix
                                                     # so the variance = 1
      
      # Add residual covariances to psi matrix
      psi[upper.tri(psi)] <- residualCovariances
      psi[lower.tri(psi)] <- t(psi)[lower.tri(psi)]
      
      if( sum(diag(psi) < 0) == 0){ # Only does the following steps if the psi values are
        # positive 
        
        save <- t(B) %*% start %*% B + psi
        saveMatrices[[j]] <- save
        start <- save
        
      }else{return("Negative Residual Variances")}
      
    } #end of j for loop
    
    
    
    } # end of CovRes/else command
  

    # Checks if the previous for loop returned 1000 matrices and if the last 
    # two matrices in the list are identical
  if(length(saveMatrices) == 1000 & 
     all.equal(saveMatrices[[999]], saveMatrices[[1000]])){ 
    
    Sigma1 <- round(save, 6)
    rownames(Sigma1) <- colnames(Sigma1) <- LETTERS[1:p]
    
  }else{return("Nonstationary Process")}
    
    # Check if matrix is positive definite
  if(matrixcalc::is.positive.definite(Sigma1) == TRUE){
      
      return(Sigma1)
    }else{
      return("Matrix is Positive Non-definite")
    }
  
  }

GenerateLongitudinalCov <- function(BList, Sigma, ResCov = NULL,
                                    timepoints){
  
  p <- ncol(Sigma)
  
  LongitudinalSigma <- matrix(0, p * timepoints, p * timepoints)
  rownames(LongitudinalSigma) <- colnames(LongitudinalSigma) <- paste0(LETTERS[1:p], 
                                                                       rep(1:timepoints, 
                                                                           each = p))
  
  if(length(BList) == 1){ # If process is stationary
    
    B <- BList[[1]]
    
    for(i in 1:timepoints){
      
      crossCovCode <- "Sigma2 %*% B" # Equation for the covariances between time 
                                     # T variables and time T+1 variables
      
      # Calculate and save covariances between variables in the same timepoint
      psi <- diag(diag(1 - t(B) %*% Sigma %*% B), p)
      
      if(is.null(ResCov) == FALSE){
        
        # Add residual covariances to psi matrix
        psi[upper.tri(psi)] <- residualCovariances
        psi[lower.tri(psi)] <- t(psi)[lower.tri(psi)]
        
      }
      
      Sigma2 <- t(B) %*% Sigma %*% B + psi # Equation for the covariances between time 
                                            # T variables
      
      LongitudinalSigma[grep(i, colnames(LongitudinalSigma)), 
                        grep(i, colnames(LongitudinalSigma))] <- Sigma2
      
    
   
       if(i != timepoints){
      
      for(j in (i+1):timepoints){
        
        # Calculate and save covariances between variables in different time points
        
        LongitudinalSigma[grep(i, colnames(LongitudinalSigma)), grep((j), colnames(LongitudinalSigma))]  <- eval(parse(text = crossCovCode))
        LongitudinalSigma[grep((j), colnames(LongitudinalSigma)), grep(i, colnames(LongitudinalSigma))]  <- t(eval(parse(text = crossCovCode)))
        
        crossCovCode <- paste0(crossCovCode, "%*% B")
        
      }}
      

    Sigma <- Sigma2
    }
    
  } # End stationary commands
  
  return(LongitudinalSigma)
}

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
