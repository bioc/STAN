



#'  
#' Initialization of bidirectional Clustering
#' 
#' @title Initialization of bidirectional Clustering
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param dStates The number of directed states.
#' @param uStates The number of undirected states. 
#' @param method Emission distribution of the model. One out of c("NegativeBinomial", "PoissonLogNormal", "NegativeMultinomial", "ZINegativeBinomial", "Poisson", "Bernoulli", "Gaussian", "IndependentGaussian")
#' @param dirFlags The flag sequence is needed when a bdHMM is fitted on undirected data (e.g.) ChIP only. It is a \code{list} of \code{character} vectors indication for each position its knwon directionality. U allows all states. F allows undirected states and states in forward direction. R allows undirected states and states in reverse direction. 
#' @param directedObs Integer vector defining the directionality (or strand-specificity) of the data tracks. Undirected (non-strand-specific) data tracks (e.g. ChIP) are indicated indicated by '0'. Directed (strand-specific) data tracks are indicated by increasing pairs of integers. For instance c(0,0,0,1,1,2,2): The first three data tracks are undirected, followed by two pairs of directed measurements.
#' @param sizeFactors Library size factors for Emissions PoissonLogNormal or NegativeBinomial as a length(obs) x ncol(obs[[1]]) matrix.
#' @param sharedCov If TRUE, (co-)variance of (Independent)Gaussian is shared over states. Only applicable to 'Gaussian' or 'IndependentGaussian' emissions. Default: FALSE. 
#' 
#' @return A HMM object.
#' @usage initBdClust(obs, dStates = 0, uStates = 0, method, directedObs = rep(0, ncol(obs[[1]])), sizeFactors = matrix(1, nrow = length(obs), ncol = ncol(obs[[1]])), sharedCov = FALSE, dirFlags = NULL)
#'  
#' @export initbdClust


initBdClust = function(obs, dStates = 0 , uStates = 0, method, 
                       directedObs = rep(0, ncol(obs[[1]])), 
                       sizeFactors = matrix(1, nrow = length(obs), ncol = ncol(obs[[1]])), 
                       sharedCov = FALSE, dirFlags = NULL ){
    initBdHMM(obs, dStates = dStates, uStates = uStates, method = method, directedObs = directedObs, sizeFactor = sizeFactors, sharedCov = sharedCov, dirFlags = dirFlags)
}

   

