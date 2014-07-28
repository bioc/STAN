

#'  
#' Given a Hidden Markov Model, the function calculates the most likely state path (viterbi) for one or more observation sequence.
#' 
#' @title Calculate the most likely state path
#' 
#' @param HMM The initial Hidden Markov Model. This is a \code{\linkS4class{HMM}}.
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param NAtol Successive positions having NAs longer than this threshold are masked in the viterbi path. 
#' @param emissionProbs List of precalculated emission probabilities of emission function is of type 'null'.
#' @param verbose \code{logical} for printing algorithm status or not.
#' 
#' @return A list containint the vterbi paths.
#' @usage getViterbi(HMM, obs=list(), NAtol=5, emissionProbs=list(), verbose=FALSE)
#' 
#' @seealso \code{\linkS4class{HMM}}
#' @examples
#' data(example)
#' hmm_fitted = fitHMM(observations, hmm_ex)
#' viterbi_hmm = getViterbi(hmm_fitted$hmm, observations)
#'
#' @export getViterbi
getViterbi = function(HMM, obs=list(), NAtol=5, emissionProbs=list(), verbose=FALSE) {
    if (class(obs) != "list") 
        stop("submit obs as list!\n")
    
    if (is.null(obs)) 
        stop("Zero observations\n")
    
    D = NULL
    if (length(obs) > 0) {
        D = as.integer(dim(obs[[1]])[2])
        D = HMM@emission@dim
    } else if (length(emissionProbs) > 0) {
        D = as.integer(dim(emissionProbs[[1]])[2])
    } else {
        stop("Must either submit observations of pre-calculated emissions!\n")
    }
    
    viterbi = .Call("RHMMVITERBI", SEXPobs = obs, SEXPpi = as.numeric(HMM@initProb), 
        SEXPtransProb = as.numeric(HMM@transMat), SEXPemission = HMM@emission@parameters, 
        SEXPtype = as.character(HMM@emission@type), SEXPdim = D, SEXPnStates = as.integer(HMM@nStates), 
        SEXPverbose = as.integer(verbose), SEXPfixedEmission = emissionProbs, PACKAGE = "STAN")
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
    
    if (class(HMM) == "bdHMM") {
        myLevels = HMM@stateLabel
        viterbi = lapply(viterbi, function(x) factor(myLevels[x], levels = myLevels))
    } else {
        myLevels = 1:HMM@nStates
        viterbi = lapply(viterbi, function(x) factor(myLevels[x], levels = myLevels))
    }
    
    
    return(viterbi)
} 
