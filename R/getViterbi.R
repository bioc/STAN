

#'  
#' Given a Hidden Markov Model, the function calculates the most likely state path (viterbi) for one or more observation sequence.
#' 
#' @title Calculate the most likely state path
#' 
#' @param hmm The initial Hidden Markov Model.
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param NAtol Successive positions having NAs longer than this threshold are masked in the viterbi path. 
#' @param emissionProbs List of precalculated emission probabilities of emission function is of type 'null'.
#' @param verbose \code{logical} for printing algorithm status or not.
#' @param sizeFactors Library size factors for Emissions PoissonLogNormal or NegativeBinomial as a length(obs) x ncol(obs[[1]]) matrix.
#' 
#' @return A list containint the vterbi paths.
#' @usage getViterbi(hmm, obs=list(), NAtol=5, emissionProbs=list(), verbose=FALSE, sizeFactors=matrix(1, nrow=length(obs), ncol=ncol(obs[[1]])))
#' 
#' @examples
#' 
#' data(example)
#' hmm_ex = initHMM(observations, nStates=3, method="Gaussian") 
#' hmm_fitted = fitHMM(observations, hmm_ex)
#' viterbi = getViterbi(hmm_fitted, observations)
#'
#' @export getViterbi
getViterbi = function(hmm, obs = list(), NAtol = 5, emissionProbs = list(), 
    verbose = FALSE, sizeFactors = matrix(1, nrow = length(obs), ncol = ncol(obs[[1]]))) {
    if (class(obs) != "list") 
        stop("Observations must be submitted as a list of matrices!")
    if (!all(sapply(obs, class) == "matrix")) 
        stop("Observations must be submitted as a list of matrices!")
    if (!class(hmm) %in% c("HMM", "bdHMM")) 
        stop("Parameter hmm is not of class HMM or bdHMM!")
    if (length(emissionProbs) > 0) {
        if (length(emissionProbs) != length(obs)) 
            stop("Pre-computed emission probabilities and observations must have the same length!")
        if (sum(sapply(emissionProbs, nrow)) != sum(sapply(obs, nrow))) 
            stop("Pre-computed emission probabilities and observations must have the same length!")
        if (!all(unique(sapply(emissionProbs, ncol)) == hmm@nStates)) 
            stop("Number of columns in emissionProbs object must equal number of HMM states.")
    }
    if (length(unique(sapply(obs, ncol))) != 1) 
        stop("All entries of observation sequences must have the same number of columns!")
    if (ncol(obs[[1]]) != sum(sapply(hmm@emission, function(y) y@dim))) 
        stop("Number of columns in observation sequence must match dimensions of observations sequence!")
    
    
    myParameters = list(emissions = hmm@emission)
    jiEmission = HMMEmission(type = "JointlyIndependent", parameters = myParameters, 
        nStates = as.integer(hmm@nStates))
    hmm@emission = list(jiEmission)
    emissionTypes = unlist(sapply(hmm@emission[[1]]@parameters$emissions, 
        function(x) rep(x@type, x@dim)))
    matchEmission2Data(obs, emissionTypes, verbose)
    hmm@emission[[1]]@parameters = prepareEmission(hmm@emission[[1]]@parameters, 
        hmm@emission[[1]]@type)
    
    observationEmissionType = unlist(sapply(hmm@emission[[1]]@parameters$emissions, 
        function(x) rep(x@type, x@dim)))
    if (any(length(observationEmissionType) != sapply(obs, ncol))) 
        stop("Number of data tracks (columns in observation matrix) do not match dimension of (bd)HMM.")
    
    for (currType in c("NegativeBinomial", "PoissonLogNormal")) {
        myD = which(observationEmissionType == currType)
        if (length(myD) > 0) {
            for (i in 1:length(hmm@emission[[1]]@parameters$emissions)) {
                hmm@emission[[1]]@parameters$emissions[[i]]@parameters$sizeFactor = lapply(1:hmm@nStates, 
                  function(x) sizeFactors[, i])
                # print(emissionParams$emissions[[i]]@parameters)
            }
        }
        
    }
    
    if (!all(unlist(lapply(obs, function(x) apply(x, 2, class))) == "numeric")) {
        stop("All observations must be of class numeric!")
    }
    
    
    if (is.null(obs)) 
        stop("Zero observations\n")
    
    D = NULL
    if (length(obs) > 0) {
        D = as.integer(dim(obs[[1]])[2])
        D = hmm@emission[[1]]@dim
    } else if (length(emissionProbs) > 0) {
        D = as.integer(dim(emissionProbs[[1]])[2])
    } else {
        stop("Must either submit observations of pre-calculated emissions!\n")
    }
    
    mySplit = list()
    if (hmm@emission[[1]]@type == "JointlyIndependent") {
        emissionTypes = unlist(sapply(hmm@emission[[1]]@parameters$emissions, 
            function(x) rep(x@type, x@dim)))
        for (currType in c("NegativeBinomial", "PoissonLogNormal")) {
            mySplit[[currType]] = list()
        }
        hmm@emission[[1]]@parameters$mySplit = mySplit
    }
    
    
    hmm@emission[[1]]@parameters$updateCov = FALSE
    viterbi = .Call("RHMMVITERBI", SEXPobs = obs, SEXPpi = as.numeric(hmm@initProb), 
        SEXPtransProb = as.numeric(hmm@transMat), SEXPemission = hmm@emission[[1]]@parameters, 
        SEXPtype = as.character(hmm@emission[[1]]@type), SEXPdim = D, SEXPnStates = as.integer(hmm@nStates), 
        SEXPverbose = as.integer(verbose), SEXPfixedEmission = emissionProbs, 
        PACKAGE = "STAN")
    if (length(emissionProbs) == 0) {
        for (n in 1:length(viterbi)) {
            wherena = which(apply(obs[[n]], 1, function(x) any(is.na(x))))
            if (length(wherena) > 0 & (length(wherena) > 1)) {
                nablocks = diff(wherena)
                currblock = 1
                for (i in 1:length(nablocks)) {
                  if (nablocks[i] == 1) {
                    nablocks[i] = currblock
                  } else {
                    currblock = currblock + 1
                    nablocks[i] = currblock
                  }
                }
                nablocks = c(1, nablocks)
                nablocks = tapply(wherena, INDEX = nablocks, identity)
                
                
                for (i in 1:length(nablocks)) {
                  if (length(nablocks[[i]]) > NAtol) {
                    viterbi[[n]][nablocks[[i]]] = NA
                  }
                }
            }
        }
        names(viterbi) = names(obs)
    } else if (length(emissionProbs) > 0) {
        names(viterbi) = names(emissionProbs)
    }
    
    if (class(hmm) == "bdHMM") {
        myLevels = hmm@stateNames
        viterbi = lapply(viterbi, function(x) factor(myLevels[x], levels = myLevels))
    } else {
        myLevels = 1:hmm@nStates
        viterbi = lapply(viterbi, function(x) factor(myLevels[x], levels = myLevels))
    }
    
    
    return(viterbi)
}
 
