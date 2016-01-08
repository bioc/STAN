



#'  
#' Initialization of bidirectional hidden Markov models
#' 
#' @title Initialization of bidirectional hidden Markov models
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param nStates The number of states. 
#' @param method Emission distribution of the model. One out of c("NegativeBinomial", "PoissonLogNormal", "NegativeMultinomial", "ZINegativeBinomial", "Poisson", "Bernoulli", "Gaussian", "IndependentGaussian")
#' @param directedObs Integer vector defining the directionality (or strand-specificity) of the data tracks. Undirected (non-strand-specific) data tracks (e.g. ChIP) are indicated indicated by '0'. Directed (strand-specific) data tracks are indicated by increasing pairs of integers. For instance c(0,0,0,1,1,2,2): The first three data tracks are undirected, followed by two pairs of directed measurements.
#' @param sizeFactors Library size factors for Emissions PoissonLogNormal or NegativeBinomial as a length(obs) x ncol(obs[[1]]) matrix.
#' @param sharedCov If TRUE, (co-)variance of (Independent)Gaussian is shared over states. Only applicable to 'Gaussian' or 'IndependentGaussian' emissions. Default: FALSE. 
#' 
#' @return A HMM object.
#' @usage initBdHMM(obs, nStates, method, directedObs = rep(0, ncol(obs[[1]])), sizeFactors = matrix(1, nrow = length(obs), ncol = ncol(obs[[1]])), sharedCov = FALSE)
#'  
#' @examples
#' 
#' data(example)
#' hmm_ex = initHMM(observations, nStates=3, method="Gaussian") 
#'
#' @export initHMM
initBdHMM = function(obs, nStates, method, directedObs = rep(0, ncol(obs[[1]])), 
    sizeFactors = matrix(1, nrow = length(obs), ncol = ncol(obs[[1]])), 
    sharedCov = FALSE) {
    
    rev.operation = 1:ncol(obs[[1]])
    for (i in 1:length(directedObs)) {
        if (directedObs[i] == 0) {
            rev.operation[i] = i
        } else {
            myPair = which(directedObs == directedObs[i])
            myPair = myPair[myPair != i]
            rev.operation[i] = myPair
        }
    }
    
    myDirCols = tapply(1:length(directedObs), INDEX = directedObs, identity)
    myDirCols = myDirCols[-which(names(myDirCols) == "0")]
    
    if (length(myDirCols) > 0) {
        for (i in 1:length(myDirCols)) {
            for (k in 1:length(obs)) {
                obs[[k]] = obs[[k]][apply(obs[[k]], 1, function(x) all(!is.na(x))), 
                  ]
                obs[[k]][, myDirCols[[i]]] = t(apply(obs[[k]][, myDirCols[[i]]], 
                  1, sort, decreasing = TRUE))
            }
        }
    }
    
    myMat = do.call("rbind", obs)
    
    km = clusterMat(myMat, nStates, method)
    
    myEmissions = NULL
    if (method == "NegativeBinomial") {
        myEmissions = initNB(km, obs, sizeFactor = sizeFactors, stateLabels = c(paste("F", 
            1:(nStates), sep = ""), paste("R", 1:(nStates), sep = "")), 
            directedObs = directedObs)
    }
    if (method == "ZINegativeBinomial") {
        myEmissions = initNB(km, obs, sizeFactor = sizeFactors, zeroInflated = TRUE, 
            stateLabels = c(paste("F", 1:(nStates), sep = ""), paste("R", 
                1:(nStates), sep = "")), directedObs = directedObs)
    }
    if (method == "PoissonLogNormal") {
        myEmissions = initPoiLog(km, obs, sizeFactor = sizeFactors, stateLabels = c(paste("F", 
            1:(nStates), sep = ""), paste("R", 1:(nStates), sep = "")), 
            directedObs = directedObs)
    }
    if (method == "Gaussian") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = lapply(state2ind, function(x) apply(myMat[x,, drop=FALSE], 2, 
            mean))
        covs_Init = lapply(1:length(means_Init), function(x) matrix(cov(myMat), ncol=ncol(myMat)))
      
		myEmissions = list(HMMEmission(type = "Gaussian", parameters = list(mu = c(means_Init, 
            lapply(means_Init, function(x) x[rev.operation])), cov = c(covs_Init, 
            lapply(covs_Init, function(x) matrix(x[rev.operation, rev.operation], ncol=ncol(x)))), 
            sharedCov = sharedCov), nStates = nStates * 2))
    }
    if (method == "IndependentGaussian") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = lapply(state2ind, function(x) apply(myMat[x,, drop=FALSE], 2, 
            mean))
        covs_Init = lapply(1:length(means_Init), function(x) matrix(cov(myMat), ncol=ncol(myMat)))
        
        for (i in (nStates + 1):(nStates * 2)) {
            means_Init[[i]] = means_Init[[(i - nStates)]][rev.operation]
            covs_Init[[i]] = covs_Init[[(i - nStates)]][rev.operation, 
                rev.operation]
        }
        
        myEmissions = list()
        for (i in 1:ncol(obs[[1]])) {
            myEmissions[[i]] = HMMEmission(type = "Gaussian", parameters = list(mu = lapply(means_Init, 
                function(x) x[i]), cov = lapply(covs_Init, function(x) matrix(var(as.vector(obs[[1]])), 
                ncol = 1)), sharedCov = sharedCov), nStates = nStates * 
                2)
        }
    }
    if (method == "Bernoulli") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = lapply(state2ind, function(x) apply(myMat[x,, drop=FALSE], 2, 
            mean))
        
        myEmissions = list()
        
        for (i in 1:ncol(obs[[1]])) {
            myEmissions[[i]] = HMMEmission(type = "Bernoulli", parameters = list(p = lapply(means_Init, 
                function(x) x[i])), nStates = nStates)
        }
    }
    
    
    nStates = nStates * 2
    initProb = rep(1/nStates, nStates)
    transMat = matrix(1/nStates, nrow = nStates, ncol = nStates)
    bdhmm = bdHMM(initProb = initProb, transMat = transMat, emission = myEmissions, 
        nStates = nStates, stateNames = c(paste("F", 1:(nStates/2), sep = ""), 
            paste("R", 1:(nStates/2), sep = "")), status = "initial", transitionsOptim = "analytical", 
        directedObs = directedObs, dimNames = colnames(obs[[1]]))
    bdhmm
    
} 
