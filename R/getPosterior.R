

#'  
#' The function calculates posterior state probabilities for one or more observation sequence. 
#' 
#' @title Calculate posterior state distribution.
#' 
#' @param hmm The initial Hidden Markov Model. This is a \code{\linkS4class{HMM}}.
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param emissionProbs List of precalculated emission probabilities of emission function is of type 'null'.
#' @param dirFlags The flag sequence is needed when a bdHMM is fitted on undirected data (e.g.) ChIP only. It is a \code{list} of \code{character} vectors indication for each position its knwon directionality. U allows all states. F allows undirected states and states in forward direction. R allows undirected states and states in reverse direction. 
#' @param verbose \code{logical} for printing algorithm status or not.
#' @param nCores Number of cores to use for computations.
#' 
#' @return A list containing for the observation sequences the posterior state (col) distribution at each position (row).
#' @usage getPosterior(hmm, obs=list(), emissionProbs=list(), dirFlags=list(), verbose=FALSE, nCores=1)
#' 
#' @seealso \code{\linkS4class{HMM}}
#' @examples 
#' data(example)
#' hmm_fitted = fitHMM(observations, hmm_ex)
#' posterior_hmm = getPosterior(hmm_fitted$hmm, observations)
#' 
#' @export getPosterior
getPosterior = function(hmm, obs=list(), emissionProbs=list(), dirFlags=list(), 
    verbose=FALSE, nCores=1) {
    emission = hmm@emission
    emissionParams = emission@parameters
    emissionParams = prepareEmission(emissionParams, hmm@emission@type)
    D = hmm@emission@dim
    initProb = hmm@initProb
    transMat = hmm@transMat
    nStates = hmm@nStates
    
    state2flag = NULL
    if (class(hmm) == "bdHMM") {
        stateLabel = hmm@stateLabel
        state2flag[grep("F", stateLabel)] = 1
        state2flag[grep("R", stateLabel)] = -1
        state2flag[which(state2flag == 0)] = 100
        state2flag = -state2flag
    }
    if (length(dirFlags) > 0) {
        flag2int = c(1, -1, 0)
        names(flag2int) = c("F", "R", "U")
        dirFlags = lapply(dirFlags, function(x) as.integer(flag2int[x]))
        
    }
    
    
    
    myGamma = .Call("RGETPOSTERIOR", SEXPobs = obs, SEXPpi = initProb, SEXPA = transMat, 
        SEXPemission = emissionParams, SEXPtype = as.character(hmm@emission@type), 
        SEXPdim = as.integer(D), SEXPk = as.integer(nStates), SEXPverbose = as.integer(verbose), 
        SEXPfixedEmission = emissionProbs, SEXPnCores = as.integer(nCores), SEXPflags = lapply(dirFlags, 
            as.integer), SEXPstate2flag = as.integer(state2flag))
    names(myGamma) = names(obs)
    
    lapply(myGamma, matrix, ncol = nStates)
} 
