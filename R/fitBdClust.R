#'  
#' The function is used to fit (bidirectional) Clusters, given one or more observation sequence. 
#' 
#' @title Fit a bidirectional Clustering
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param BdClust The initial Bidirectional Cluster.
#' @param convergence Convergence cutoff for EM-algorithm (default: 1e-6).
#' @param maxIters Maximum number of iterations.
#' @param dirFlags The flag sequence is needed when a bdHMM is fitted on undirected data (e.g.) ChIP only. It is a \code{list} of \code{character} vectors indication for each position its knwon directionality. U allows all states. F allows undirected states and states in forward direction. R allows undirected states and states in reverse direction. 
#' @param emissionProbs List of precalculated emission probabilities of emission function is of type 'null'.
#' @param effectiveZero Transitions below this cutoff are analytically set to 0 to speed up comptuations.
#' @param verbose \code{logical} for printing algorithm status or not.
#' @param nCores Number of cores to use for computations.
#' @param incrementalEM When TRUE, the incremental EM is used to fit the model, where parameters are updated after each iteration over a single observation sequence.
#' @param updateTransMat Wether transitions should be updated during model learning, default: TRUE.
#' @param sizeFactors Library size factors for Emissions PoissonLogNormal or NegativeBinomial as a length(obs) x ncol(obs[[1]]) matrix.
#' 
#' @return A list containing the trace of the log-likelihood during EM learning and the fitted HMM model.
#' @usage fitBdClust(obs=list(), BdClust , convergence=1e-6, maxIters=1000, dirFlags=list(), emissionProbs=list(), effectiveZero=0, verbose=FALSE, nCores=1, incrementalEM=FALSE, updateTransMat=TRUE, sizeFactors=matrix(1, nrow=length(obs), ncol=ncol(obs[[1]])))
#' 
#' @examples 
#' 
#' data(example)
#' bdclust_ex = initBdClust(observations, dStates=3, method="Gaussian")
#' bdclust_fitted = fitBdClust(observations, bdclust_ex)
#'
#' @export fitBdClust


fitBdClust <- function(obs=list(), BdClust = BdClust, convergence=1e-6, maxIters=1000, dirFlags=list(), emissionProbs=list(), effectiveZero=0, verbose=FALSE, nCores=1, incrementalEM=FALSE, updateTransMat=TRUE, sizeFactors=matrix(1, nrow=length(obs), ncol=ncol(obs[[1]]))){
    
    fitHMM(obs = obs,hmm =  BdClust, convergence, maxIters, dirFlags, emissionProbs, effectiveZero, verbose, nCores, incrementalEM, updateTransMat, sizeFactors, clustering = TRUE)
}