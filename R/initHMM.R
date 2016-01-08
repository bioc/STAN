


#'  
#' Initialization of hidden Markov models
#' 
#' @title Initialization of hidden Markov models
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param nStates The number of states. 
#' @param method Emission distribution of the model. One out of c("NegativeBinomial", "PoissonLogNormal", "NegativeMultinomial", "ZINegativeBinomial", "Poisson", "Bernoulli", "Gaussian", "IndependentGaussian")
#' @param sizeFactors Library size factors for Emissions PoissonLogNormal or NegativeBinomial as a length(obs) x ncol(obs[[1]]) matrix.
#' @param sharedCov If TRUE, (co-)variance of (Independent)Gaussian is shared over states. Only applicable to 'Gaussian' or 'IndependentGaussian' emissions. Default: FALSE. 
#' 
#' @return A HMM object.
#' @usage initHMM(obs, nStates, method, sizeFactors = matrix(1, nrow = length(obs), ncol = ncol(obs[[1]])), sharedCov = FALSE)
#' 
#' @examples
#' 
#' data(example)
#' hmm_ex = initHMM(observations, nStates=3, method="Gaussian") 
#'
#' @export initHMM
initHMM = function(obs, nStates, method, sizeFactors = matrix(1, nrow = length(obs), 
    ncol = ncol(obs[[1]])), sharedCov = FALSE) {

    ncols = sapply(obs, ncol)
    for (k in 1:length(obs)) {
        obs[[k]] = obs[[k]][apply(obs[[k]], 1, function(x) all(!is.na(x))), ]
		if(ncols[k] == 1) {
			obs[[k]] = matrix(obs[[k]], ncol=1)
		}
    }
    myMat = do.call("rbind", obs)
	km = clusterMat(myMat, nStates, method)
    
    myEmissions = NULL
    if (method == "NegativeBinomial") {
        myEmissions = initNB(km, obs, sizeFactor = sizeFactors)
    }
    if (method == "ZINegativeBinomial") {
        myEmissions = initNB(km, obs, sizeFactor = sizeFactors, zeroInflated = TRUE)
    }
    if (method == "PoissonLogNormal") {
        myEmissions = initPoiLog(km, obs, sizeFactor = sizeFactors)
    }
    if (method == "IndependentGaussian") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = lapply(state2ind, function(x) apply(myMat[x, 1:ncol(myMat), drop=FALSE], 2, 
            mean))
        covs_Init = lapply(state2ind, function(x) apply(myMat[x, 1:ncol(myMat), drop=FALSE], 2, 
            var))
        
        myEmissions = list()
        
        for (i in 1:ncol(obs[[1]])) {
            myEmissions[[i]] = HMMEmission(type = "Gaussian", parameters = list(mu = lapply(means_Init, 
                function(x) x[i]), cov = lapply(covs_Init, function(x) matrix(var(as.vector(obs[[1]])), 
                ncol = 1)), sharedCov = sharedCov), nStates = nStates)
        }
    }
    if (method == "Gaussian") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = lapply(state2ind, function(x) apply(myMat[x, 1:ncol(myMat), drop=FALSE], 2, 
            mean))
        covs_Init = lapply(1:length(means_Init), function(x) cov(myMat))
        
        myEmissions = list(HMMEmission(type = "Gaussian", parameters = list(mu = means_Init, 
            cov = covs_Init, sharedCov = sharedCov), nStates = nStates))
    }
    if (method == "NegativeMultinomial") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = do.call("rbind", lapply(state2ind, function(x) apply(myMat[x, 1:ncol(myMat), drop=FALSE], 2, mean)))
        
        multiNomInit = lapply(1:nrow(means_Init), function(x) means_Init[x, 
            2:ncol(means_Init)]/sum(means_Init[x, 2:ncol(means_Init)]))
        multiNomEmission = HMMEmission(type = "Multinomial", parameters = list(p = multiNomInit), 
            nStates = nStates)
        myEmissionsNB = initNB(km, obs)[1]
        
        myEmissions = c(myEmissionsNB, multiNomEmission)
    }
    if (method == "Bernoulli") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = lapply(state2ind, function(x) apply(myMat[x, 1:ncol(myMat), drop=FALSE], 2, 
            mean))
        
        myEmissions = list()
        
        for (i in 1:ncol(obs[[1]])) {
            myEmissions[[i]] = HMMEmission(type = "Bernoulli", parameters = list(p = lapply(means_Init, 
                function(x) x[i])), nStates = nStates)
        }
    }
    if (method == "Poisson") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = lapply(state2ind, function(x) apply(myMat[x, 1:ncol(myMat), drop=FALSE], 2, 
            mean))
        
        myEmissions = list()
        
        for (i in 1:ncol(obs[[1]])) {
            myEmissions[[i]] = HMMEmission(type = "Poisson", parameters = list(lambda = lapply(means_Init, 
                function(x) x[i])), nStates = nStates)
        }
    }
    
    
    initProb = rep(1/nStates, nStates)
    transMat = matrix(1/nStates, nrow = nStates, ncol = nStates)
    hmm = HMM(initProb = initProb, transMat = transMat, emission = myEmissions, 
        nStates = nStates)
    hmm
    
}



#' Internal function that does k-means clustering for initialization. The function is called by initHMM and initBdHMM.
#' @param myMat Data matrix.
#' @param nStates Number of states.
#' @param method Emission distribution of the model. One out of c("NegativeBinomial", "PoissonLogNormal", "NegativeMultinomial", "ZINegativeBinomial", "Poisson", "Bernoulli", "Gaussian", "IndependentGaussian").
#' 
#' @keywords internal
#' @noRd
clusterMat = function(myMat, nStates, method) {
    km = NULL
    if (method %in% c("NegativeBinomial", "PoissonLogNormal", "NegativeMultinomial", 
        "ZINegativeBinomial", "Poisson")) {
        myLambdas = c()
        for (i in 1:ncol(myMat)) {
            myLambdas[i] = mean(myMat[, i])
        }
        mycutoffs = sapply(myLambdas, function(x) min(which(1 - ppois(1:50000, 
            x) < 1e-04)))
		for(i in 1:ncol(myMat)) {
			if(mycutoffs[i] > max(myMat[,i])) {
				mycutoffs[i] = quantile(myMat[,i], 0.25)
			}
		}			

        for (i in 1:ncol(myMat)) {
			if(mycutoffs[i] > 0) {
				myMat[(myMat[, i] <= mycutoffs[i]), i] = 0
			}
        }
		myBckg = which(apply(myMat, 1, function(x) all(x == 0)))
        mySignal = setdiff(1:nrow(myMat), myBckg)
		
        suppressWarnings(km <- kmeans(log(myMat[mySignal, drop=FALSE] + 1), centers = nStates - 
            1, iter.max = 1000, nstart = 10))
        
        km$centers = matrix(NA, nrow = nStates, ncol = ncol(myMat))
        signalClusters = km$cluster
        km$cluster = rep(NA, nrow(myMat))
        km$cluster[myBckg] = nStates
        km$cluster[mySignal] = signalClusters
        if (any(is.na(km$cluster))) 
            stop("NA in clustering!\n")
    } else if (method %in% c("Bernoulli", "Gaussian", "IndependentGaussian")) {
        suppressWarnings(km <- kmeans(myMat, centers = nStates, iter.max = 1000, 
            nstart = 10))
    } else {
        stop("Method ", method, " is not supported!\n", sep = "")
    }
    km
}
 
