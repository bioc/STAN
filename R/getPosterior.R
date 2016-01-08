

#'  
#' The function calculates posterior state probabilities for one or more observation sequence. 
#' 
#' @title Calculate posterior state distribution.
#' 
#' @param hmm The Hidden Markov Model.
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param emissionProbs List of precalculated emission probabilities of emission function is of type 'null'.
#' @param dirFlags The flag sequence is needed when a bdHMM is fitted on undirected data (e.g.) ChIP only. It is a \code{list} of \code{character} vectors indication for each position its knwon directionality. U allows all states. F allows undirected states and states in forward direction. R allows undirected states and states in reverse direction. 
#' @param verbose \code{logical} for printing algorithm status or not.
#' @param nCores Number of cores to use for computations.
#' @param sizeFactors Library size factors for Emissions PoissonLogNormal or NegativeBinomial as a length(obs) x ncol(obs[[1]]) matrix.
#' 
#' @return A list containing for the observation sequences the posterior state (col) distribution at each position (row).
#' @usage getPosterior(hmm, obs=list(), emissionProbs=list(), dirFlags=list(), verbose=FALSE, nCores=1, sizeFactors=matrix(1, nrow=length(obs), ncol=ncol(obs[[1]])))
#' 
#' @examples 
#' 
#' data(example)
#' hmm_ex = initHMM(observations, nStates=3, method="Gaussian") 
#' hmm_fitted = fitHMM(observations, hmm_ex)
#' posterior = getPosterior(hmm_fitted, observations)
#' 
#' @export getPosterior
getPosterior = function(hmm, obs = list(), emissionProbs = list(), dirFlags = list(), 
    verbose = FALSE, nCores = 1, sizeFactors=matrix(1, nrow=length(obs), ncol=ncol(obs[[1]]))) {
    mySplit = list()
    if (class(obs) != "list") 
        stop("Observations must be submitted as a list of matrices!")
    if (!all(sapply(obs, class) == "matrix")) 
        stop("Observations must be submitted as a list of matrices!")
    if (!class(hmm) %in% c("HMM", "bdHMM")) 
        stop("Parameter hmm is not of class HMM or bdHMM!")
    if (length(emissionProbs) > 0) {
        if (length(emissionProbs) != length(obs)) 
            stop("Pre-computed emission probabilities and observations must have the same length!")
        if (sum(sapply(emissionProbs, length)) != sum(sapply(obs, nrow))) 
            stop("Pre-computed emission probabilities and observations must have the same length!")
        if (!all(unique(sapply(emissionProbs, ncol)) == hmm@nStates)) 
            stop("Number of columns in emissionProbs object must equal number of HMM states.")
    }
    if (length(dirFlags) > 0) {
        if (length(dirFlags) != length(obs)) 
            stop("Flag sequence and observations must have the same length!")
        if (sum(sapply(dirFlags, length)) != sum(sapply(obs, nrow))) 
            stop("Flag sequence and observations must have the same length!")
    }
    if (length(unique(sapply(obs, ncol))) != 1) 
        stop("All entries of observation sequences must have the same number of columns!")
    if (ncol(obs[[1]]) != sum(sapply(hmm@emission, function(y) y@dim))) 
        stop("Number of columns in observation sequence must match dimensions of observations sequence!")
    
	for(i in 1:length(hmm@emission)) {
		if(!all(sizeFactors == 1) & (! hmm@emission[[i]]@type %in% c("NegativeBinomial", "PoissonLogNormal"))) stop("Size factor correction is only implemented for NegativeBinomial and PoissonLogNormal distributions!\n")
	}
	
    myParameters = list(emissions = hmm@emission)
    jiEmission = HMMEmission(type = "JointlyIndependent", parameters = myParameters, 
        nStates = as.integer(hmm@nStates))
    hmm@emission = list(jiEmission)
    
    emissionTypes = unlist(sapply(hmm@emission[[1]]@parameters$emissions, 
        function(x) rep(x@type, x@dim)))
    matchEmission2Data(obs, emissionTypes, verbose)
    
	observationEmissionType = unlist(sapply(hmm@emission[[1]]@parameters$emissions, function(x) rep(x@type, x@dim)))
	if(hmm@emission[[1]]@type == "JointlyIndependent") {
		if(length(observationEmissionType) == 0) stop("Please specify observationEmissionType to ensure that JointlyIndependent emission match observation data types!")
		if(length(observationEmissionType) != dim(obs[[1]])[2]) stop("Length of observationDataType does not match number of observation cols (data tracks)!")
		emissionTypes = unlist(sapply(hmm@emission[[1]]@parameters$emissions, function(x) rep(x@type, x@dim)))
		matchEmission2Data(obs, emissionTypes, verbose)
		
		if(! all(emissionTypes == observationEmissionType)) stop("Types of emission functions does not match emission types of observations!")
		## TODO: check what happens when dimensionality does not match
		for(currType in c("NegativeBinomial", "PoissonLogNormal")) {
			
			myD = which(observationEmissionType == currType)
			if(length(myD) > 0) {
				for(i in 1:length(hmm@emission[[1]]@parameters$emissions)) {
					hmm@emission[[1]]@parameters$emissions[[i]]@parameters$sizeFactor = lapply(1:hmm@nStates, function(x) sizeFactors[,i])
				}
			}
		}
	}
	
	for (currType in c("NegativeBinomial", "PoissonLogNormal")) {
        mySplit[[currType]] = list()
    }
    hmm@emission[[1]]@parameters$mySplit = mySplit
    
    
    emission = hmm@emission[[1]]
    emissionParams = emission@parameters
    emissionParams = prepareEmission(emissionParams, emission@type)
    D = emission@dim
    initProb = hmm@initProb
    transMat = hmm@transMat
    nStates = hmm@nStates
    
    state2flag = NULL
    if (class(hmm) == "bdHMM") {
        stateLabel = hmm@stateNames
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
    
    
    
    
    
    myGamma = .Call("RGETPOSTERIOR", SEXPobs = obs, SEXPpi = initProb, 
        SEXPA = transMat, SEXPemission = emissionParams, SEXPtype = as.character(emission@type), 
        SEXPdim = as.integer(D), SEXPk = as.integer(nStates), SEXPverbose = as.integer(verbose), 
        SEXPfixedEmission = emissionProbs, SEXPnCores = as.integer(nCores), 
        SEXPflags = lapply(dirFlags, as.integer), SEXPstate2flag = as.integer(state2flag), 
        PACKAGE = "STAN")
    names(myGamma) = names(obs)
    
    lapply(myGamma, matrix, ncol = nStates)
} 
