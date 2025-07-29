library(foreach)
library(doParallel)
library(doSNOW)
# install.packages("psychonetrics")
# library(psychonetrics)

## revert to previous version of psychonetrics
# PackageUrl <- "http://cran.r-project.org/src/contrib/Archive/psychonetrics/psychonetrics_0.10.tar.gz"
# install.packages(PackageUrl, repos=NULL, type="source")


## Load Functions
source("C:/Users/mijke/Box/CLPN/Run on Mijke Computer/ComparingCLPNEstimatorFunctions_Stationary.R")

#Create conditions/condition matrix 
pVec <- c(5, 9, 16)
nVec <- c(500, 2000)
sparsityVec <- c(.5, .8)

conditions <- expand.grid(p = pVec, n = nVec, sparsity = sparsityVec)

## which conditions * models have not already been done? ##########################
mlist <- NULL
for(c in 1:nrow(conditions)){
  mvec <- rep(0, 20)
  for (m in 1:20){
    p <- conditions[c, "p"]
    n <- conditions[c, "n"]
    sparsity <- conditions[c, "sparsity"]
    file_path <- paste('C:/Users/mijke/Box/CLPN/Run on Mijke Computer/BetweenPersonStableVariance/Results/Condition Results/Stationary/ResultsDf', n, p, sparsity, m, '.Rdata', sep = '.')
    mvec[m] <- file.exists(file_path)
  } #\end m loop
  mlist[[c]] <- which(mvec == 0)
} #\end c loop
mlist



## Set non-varying variables
timepoints <- 4 
alpha <- .01

## Create a weight bank that we can randomly sample from to populate the 
#  non-zero parameters in the beta and psi population matrices
weightBank <- c(seq(-.3, -.1, .005), seq(.1, .3, .005)) 
resCovBank <- rnorm(1000, mean =0, sd = .05)

#Parallel Set-up
cores <- detectCores()
cl <- makeCluster(8)
registerDoSNOW(cl)

exp_pkg <- c("glmnet", "psychonetrics", "lavaan", "matrixcalc", "dplyr", "clusterGeneration")

ResultList <- foreach(c = 1:nrow(conditions),  .errorhandling = "stop") %:% 
  foreach (m = mlist[[c]], .packages = exp_pkg) %dopar% {
    n <- conditions[c,"n"]
    p <- conditions[c,"p"]
    sparsity <- conditions[c, "sparsity"]
    model <- m
    
    ResultsDf <- NULL
    # Randomly generate a beta matrix with p variables in each timepoint 
    # and (1 - sparsity)% of the parameters set as non-zero in the 
    # population 
    # Then from the beta matrix generate the a stationary covariance matrix. 
    # Repeat these steps if we end up with a nonstationary covariance matrix
    # or a positive non-definite covariance matrix
    timesrepeat <- 0
    repeat{
      
      timesrepeat <- timesrepeat + 1
      
      # Select residual covariances
      residualCovariances <- sample(resCovBank, (p * (p-1)/2), replace = TRUE)

      B <- GenerateBeta(p = p, sparsity = sparsity, weightBank = weightBank)
      Sigma1 <- GenerateStationaryCov(B, ResCov = residualCovariances)
      
      #if the generatestationarycov function produces a positive definite matrix
      if (is.character(Sigma1) == FALSE) {
        
        #try to get the longitudinal sigma
        longitudinalSigma <- GenerateLongitudinalCov(list(B), Sigma1, ResCov = residualCovariances, 
                                                     timepoints)
        
        #since determinant of positive definite matrix is 0, then break the 
        #loop if it's greater than 0
        #if (round(det(longitudinalSigma), 20) > 0) {break}
        if(min(eigen(longitudinalSigma)$values) > 0) {break}
        
        #if the determinant is not positive, try again
        else {
          print("Non-positive definite, trying again"); flush.console() 
          next
        }
        
      }
      
      if(is.character(Sigma1) == FALSE){break}
      
      if(timesrepeat == 50){next}
      
    }#end repeat command
    
    # Generate a covariance matrix for all timepoints
    longitudinalSigmaWithin <- GenerateLongitudinalCov(list(B), Sigma1, ResCov = residualCovariances, 
                                                 timepoints)
      
    # set between-person stable effects
    betweenCovMat <- cov2cor(genPositiveDefMat(dim = p)$Sigma)

    #lambda matrix relating each observed var to random intercept
    lambda <- matrix(0, p, p*(timepoints))
    for(t in 1:(timepoints)){
      diag(lambda[1:p, ((t-1)*p+1):(t*p)]) <- 1
    }
    
    longitudinalSigmaBetween <- t(lambda) %*% betweenCovMat %*% lambda 

    longitudinalSigma <- longitudinalSigmaWithin + longitudinalSigmaBetween
    
    for(g in 1:100){ # Generate this many samples from each population covariance matrix
      
      longitudinalData <- MASS::mvrnorm(n = n, mu = rep(0, ncol(longitudinalSigma)), 
                                        Sigma = longitudinalSigma,
                                        empirical = FALSE)
      
      ########################################
      ##      Estimate Models               ##
      ########################################
      
      # We're comparing 4 estimation methods 
      # 1. lvar Model: Sacha approach which estimates a model and then prunes
      #                it using BIC
      # 2. GLM Model: Currently what we use in the CLPN paper. It uses regularized
      #               regression
      # 3. SEM Model: Unregularized approach where we use the significance levels 
      #               to prune the model 
      # 4. Hybrid Model: We estimate the network structure using the GLM approach
      #                  Then we feed that structure into lavaan to get the 
      #                  unregularized parameter estimates
      
      
      # lvar Model
      ##############################
      
      # The design matrix, with a row indicating a variable and a 
      # column a wave of measurements. 
      designMat <- matrix(paste0(colnames(B), rep(1:timepoints, each = p)),
                          nrow = p, ncol = timepoints, byrow = FALSE)
      rownames(designMat) <- colnames(B)
      
      
      lvarModel <- panelvar(
        data = longitudinalData, 
        vars = designMat) 
      
      lvarStart <- Sys.time()
      # Run model:
      pruneModel <- lvarModel %>% runmodel %>% 
        prune(alpha = alpha, adjust = "none", recursive = FALSE)  
      
      lvarTime <- Sys.time() - lvarStart
      
      # Get estimated B matrix
      lvarB <- t(getmatrix(pruneModel, "beta"))
      
      
      
      # GLM Model
      ######################################
      
      #the following code loops the regularized regression through each
      #variable at time 2 as a DV, regressing on each variable at time 1
      #as predictors. It assumes that the first p columns in your data
      #correspond to the variables at time 1 and the next p correspond
      #to the variables at time 2. These variables must be in the same 
      #order at each time point.
      
      longVarNames <- colnames(longitudinalData)
      adjList <- NULL
      
      for(l in 1:(timepoints-1)){
        
        adjMat <- matrix(0, p, p) #set up empty matrix of coefficients
        for (i in 1:p){
          
          glmStart <- Sys.time()
          
          lassoreg <- cv.glmnet(as.matrix(longitudinalData[,grep(as.character(l), longVarNames)]), 
                                longitudinalData[, grep(as.character(l + 1), longVarNames)[i]], 
                                family = "gaussian", alpha = 1, standardize=TRUE)
          
          glmTime <- Sys.time() - glmStart
          
          lambda <- lassoreg$lambda.min 
          
          adjMat[1:p, i] <- coef(lassoreg, s = lambda, exact = FALSE)[2:(p+1)]
          
        }
        
        adjList[[l]] <- adjMat 
      } 
      
      
      # Unregularized/SEM Model
      #########################################
      
      # Create lavaan syntax for regressing each T+1 variable on each T variable
      lavaanSyntax <- CreateLavaanSyntax(longVarNames, timepoints)
      
      semStart <- Sys.time()
      
      unregularizedModel <- sem(model = lavaanSyntax, data = longitudinalData,
                                fixed.x = F)
      
      semTime <- Sys.time() - semStart
      
      allRegressions <- parameterestimates(unregularizedModel)[grep("~~", parameterestimates(unregularizedModel)$op,
                                                                    invert = TRUE),]
      
      sigRegressions <- allRegressions[allRegressions$pvalue < alpha,]
      
      
      
      ## Hybrid Model
      ##########################################
      
      # Create lavaan syntax for a specific model
      hybridLavaanSyntax <- CreateLavaanSyntax(longVarNames, timepoints, model = adjList)
      
      hybridStart <- Sys.time()
      
      hybridModel <- sem(model = hybridLavaanSyntax, data = longitudinalData,
                         fixed.x = F)
      
      hybridTime <- Sys.time() - hybridStart
      
      allHybridRegressions <- parameterestimates(hybridModel)[grep("~~", parameterestimates(hybridModel)$op,
                                                                   invert = TRUE),]
      
      sigHybridRegressions <- allHybridRegressions[allHybridRegressions$pvalue < alpha,]
      
      
      
      ##################################################
      ##      Calculate Performance Measures          ##
      ##################################################
      
      hybridB <- t(inspect(hybridModel, "est")$beta)
      semB <- t(inspect(unregularizedModel, "est")$beta)
      
      # hybridB might be a different size than the other B matrices
      # if there's a variable that isn't involved in any of the 
      # significant effects it won't be included in hybridB
      # so we need to add columns and rows to make it the same size
      if(ncol(hybridB) != ncol(semB)){
        
        missingVars <- setdiff(colnames(semB), colnames(hybridB))
        
        for(v in 1:length(missingVars)){
          
          newNames <- c(colnames(hybridB), missingVars)
          
          hybridB <- cbind(hybridB, rep(0, nrow(hybridB)))
          hybridB <- rbind(hybridB, rep(0, ncol(hybridB)))
          
          colnames(hybridB) <- rownames(hybridB) <- newNames
          
        }}
      
      # The beta matrices for our sem and hybrid model have the weights for 
      # both the significant and nonsignificant parameters. Since we're using
      # significance levels to prune these models, I only want the significant
      # weights in my beta matrix when I calculate the performance metrics
      hybridSigB <- CreateSigB(B = hybridB, 
                               nonSigParam = allHybridRegressions[allHybridRegressions$pvalue > alpha,], 
                               alpha = alpha)
      
      semSigB <- CreateSigB(B = semB, 
                            nonSigParam = allRegressions[allRegressions$pvalue > alpha,], 
                            alpha = alpha)
      # Save population and estimates betas
      parameterResultsList <- list(Population = B,
                                   Sacha = lvarB,
                                   GlmReg = adjList,
                                   SEM = semSigB, 
                                   Hybrid = hybridSigB)
      
      save(parameterResultsList, file = paste('C:/Users/mijke/Box/CLPN/Run on Mijke Computer/BetweenPersonStableVariance/Results/Model Results/Stationary/ModelResults', n, p, sparsity, model, g, '.Rdata', sep = '.'))
      
      
      hybridPM <- NULL
      semPM <- NULL
      glmPM <- NULL
      
      # Caclculates true and false positives and true and false negatives. 
      # i can use these to calculate sensitivity, specifict, false negative rate, 
      # and true dectection rate
      
      for(i in 1:(timepoints-1)){
        
        glmPM <- cbind(glmPM, CalculatePerformanceCounts(true = B, est = adjList[[i]]))
        
        semBT <- round(semSigB[grep(i, rownames(semB)), grep(i + 1, colnames(semB))], 4)
        semPM <-  cbind(semPM, CalculatePerformanceCounts(true = B, est= semBT))
        
        
        hybridBT <- round(hybridSigB[grep(i, rownames(hybridSigB)), grep(i + 1, colnames(hybridSigB))], 4)
        hybridBT <- hybridBT[rownames(semBT), colnames(semBT)]
        hybridPM <- cbind(hybridPM, CalculatePerformanceCounts(true = B, est= hybridBT))
        
      }
      
      # The lvar model assumes stationarity so it only returns one beta matrix 
      # for all the time points. 
      # So I am multiplying the TP, TN, FP, and FN counts by the number of beta 
      # matrcies the other models have (The models that don't assume stationarity)
      lvarPM <- CalculatePerformanceCounts(true = B, est= lvarB) * (timepoints-1)
      
      # Sum TP, TN, FP, and FN across the timpoints
      hybridPM <- rowSums(hybridPM)
      semPM <- rowSums(semPM)
      glmPM <- rowSums(glmPM)
      
      
      # Now calculate sensitivity, FPR, TDR, and specificity
      Res <-rbind(c(glmPM, CalculatePerformanceMeasures(glmPM)),
                  c(lvarPM, CalculatePerformanceMeasures(lvarPM)),
                  c(semPM, CalculatePerformanceMeasures(semPM)),
                  c(hybridPM, CalculatePerformanceMeasures(hybridPM)))
      
      # Calculate correlation and distance between population and estimated true edge weights
      glmCor <- semCor <- hybridCor <- lvarCor <- rep(0, times= (timepoints-1))
      glmDist <- semDist <- hybridDist <- lvarDist <- rep(0, times= (timepoints-1))
      
      for(t in 1: (timepoints-1)){
        
        trueEdgeLoc <- which(B != 0)
        
        trueEdgeWeights <- as.vector(B)[trueEdgeLoc]
        
        #lvar
        lvarTP <- as.vector(lvarB)[trueEdgeLoc]
        lvarEstTP <- which(lvarTP !=0)
        lvarCor[t] <-cor(trueEdgeWeights[lvarEstTP], lvarTP[lvarEstTP])
        lvarDist[t] <- mean(abs(trueEdgeWeights[lvarEstTP] - lvarTP[lvarEstTP]))
      
        
        #glm
        
        glmTP <- as.vector(adjList[[t]])[trueEdgeLoc]
        glmEstTP <- which(glmTP !=0)
        glmCor[t] <-cor(trueEdgeWeights[glmEstTP], glmTP[glmEstTP])
        glmDist[t] <-  mean(abs(trueEdgeWeights[glmEstTP] - glmTP[glmEstTP]))
        
        # sem
        
        rnames <- paste0(colnames(B), t)
        cnames <- paste0(colnames(B), t+1)
        
        semTP <- as.vector(semB[rnames, cnames])[trueEdgeLoc]
        semEstTP <- which(semTP !=0)
        semCor[t] <-cor(trueEdgeWeights[semEstTP], semTP[semEstTP])
        semDist[t] <-mean(abs(trueEdgeWeights[semEstTP] - semTP[semEstTP]))
        
        # hybrid
        hybridTP <- as.vector(hybridB[rnames, cnames])[trueEdgeLoc]
        hybridEstTP <- which(hybridTP !=0)
        hybridCor[t] <-cor(trueEdgeWeights[hybridEstTP], hybridTP[hybridEstTP])
        hybridDist[t] <-mean(abs(trueEdgeWeights[hybridEstTP] - hybridTP[hybridEstTP]))
      }
      
      TrueEdgeCor <- c(mean(glmCor), mean(lvarCor), mean(semCor), mean(hybridCor))
      TrueEdgeDist <- c(mean(glmDist), mean(lvarDist), mean(semDist), mean(hybridDist))
      
      
      
      # Create results dataframe
      ResultsDf <-  rbind(ResultsDf, data.frame(conditions[c,], model = m, sample = g, 
                                    method = c("glm", "lvar", "sem", "hybrid"), Res, 
                                    TrueEdgeCor, TrueEdgeDist, 
                                    runtime = c(glmTime,
                                                lvarTime,
                                                semTime,
                                                hybridTime)))
      
    }#end g for loop 
    
    save(ResultsDf, file = paste('C:/Users/mijke/Box/CLPN/Run on Mijke Computer/BetweenPersonStableVariance/Results/Condition Results/Stationary/ResultsDf', n, p, sparsity, model, '.Rdata', sep = '.'))
    
    ResultsDf
  }
  

stopCluster(cl)

SimResults <- NULL
for(l in 1:length(ResultList)){
  SimResults <- rbind(SimResults, do.call(rbind, ResultList[[l]]))
}

date <- Sys.Date()

save(SimResults, file = paste('C:/Users/mijke/Box/CLPN/Run on Mijke Computer/BetweenPersonStableVariance/Results/StationRes', date, '.Rdata'))
