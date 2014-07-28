
#'  
#' The function is used to fit (bidirectional) Hidden Markov Models, given one or more observation sequence. 
#' 
#' @title Fit a Hidden Markov Model
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param hmm The initial Hidden Markov Model. This is a \code{\linkS4class{HMM}}.
#' @param convergence Convergence cutoff for EM-algorithm (default: 1e-6).
#' @param maxIters Maximum number of iterations.
#' @param dirFlags The flag sequence is needed when a bdHMM is fitted on undirected data (e.g.) ChIP only. It is a \code{list} of \code{character} vectors indication for each position its knwon directionality. U allows all states. F allows undirected states and states in forward direction. R allows undirected states and states in reverse direction. 
#' @param emissionProbs List of precalculated emission probabilities of emission function is of type 'null'.
#' @param effectiveZero Transitions below this cutoff are analytically set to 0 to speed up comptuations.
#' @param verbose \code{logical} for printing algorithm status or not.
#' @param nCores Number of cores to use for computations.
#' @param incrementalEM When TRUE, the incremental EM is used to fit the model, where parameters are updated after each iteration over a single observation sequence.
#' @param observationEmissionType Only needed when HMMEmission is 'JointlyIndependent'. Defines for each dimension (columns in obs) of the data the type of emission to be used.
#' 
#' @return A list containing the trace of the log-likelihood during EM learning and the fitted HMM model.
#' @usage fitHMM(obs=list(), hmm, convergence=1e-06, maxIters=1000, dirFlags=list(), emissionProbs=list(), effectiveZero=0, verbose=FALSE, nCores=1, incrementalEM=FALSE, observationEmissionType=c())
#' 
#' @seealso \code{\linkS4class{HMM}}
#' @examples 
#' data(example)
#' hmm_fitted = fitHMM(observations, hmm_ex)
#'
#' @export fitHMM
fitHMM = function(obs=list(), hmm, convergence=1e-06, maxIters=1000, dirFlags=list(), 
    emissionProbs=list(), effectiveZero=0, verbose=FALSE, nCores=1, incrementalEM=FALSE, 
    observationEmissionType=c()) {
    emission = hmm@emission
    emissionParams = emission@parameters
    emissionParams = prepareEmission(emissionParams, hmm@emission@type)
    D = hmm@emission@dim
    initProb = hmm@initProb
    transMat = hmm@transMat
    nStates = hmm@nStates
    
    if (hmm@emission@type == "null" & length(emissionProbs) == 0) 
        stop("Must supply pre-calculated emission probabilities when using emission function of type 'null'!")
    if (hmm@emission@type == "null" & (length(emissionProbs) != length(obs))) 
        stop("Observation sequence must have same length as pre-calculated emissions!")
    
    if (emission@type == "JointlyIndependent") {
        if (length(observationEmissionType) == 0) 
            stop("Please specify observationEmissionType to ensure that JointlyIndependent emission match observation data types!")
        if (length(observationEmissionType) != dim(obs[[1]])[2]) 
            stop("Length of observationDataType does not match number of observation cols (data tracks)!")
        emissionTypes = sapply(emission@parameters$emissions, function(x) rep(x@type, 
            x@dim))
        if (!all(emissionTypes == observationEmissionType)) 
            stop("Types of emission functions does not match emission types of observations!")
    }
    
    bdHMM.settings = list(dirFlags = dirFlags, stateLabel = character(), bidirectional.mc = TRUE, 
        optim.method = "", rev.operation = NULL, bidirOptimParams = list())
    
    if (class(hmm) == "bdHMM") {
        bdHMM.settings$optim.method = hmm@transitionsOptim
        if (length(bdHMM.settings$dirFlags) > 0 & (!all(unlist(bdHMM.settings$dirFlags) %in% 
            c("F", "R", "U")))) {
            stop("dirFlags must be one of c(\"F\", \"R\", \"U\")!\n")
        } else {
            flag2int = c(1, -1, 0)
            names(flag2int) = c("F", "R", "U")
            bdHMM.settings$dirFlags = lapply(bdHMM.settings$dirFlags, function(x) as.integer(flag2int[x]))
            # print(dirFlags)
        }
        
        bdHMM.settings$stateLabel = hmm@stateLabel
        bdHMM.settings$bidirectional.mc = hmm@transitionsOptim %in% c("rsolnp", "analytical")
        bdHMM.settings$rev.operation = 1:sum(hmm@emission@dim)  ##changed Julia
        if (any(table(hmm@directedObs[hmm@directedObs != 0]) > 2)) 
            stop("Matching of directed observations is not possible!")
        for (i in 1:length(hmm@directedObs)) {
            if (hmm@directedObs[i] == 0) {
                bdHMM.settings$rev.operation[i] = i
            } else {
                myPair = which(hmm@directedObs == hmm@directedObs[i])
                myPair = myPair[myPair != i]
                bdHMM.settings$rev.operation[i] = myPair
            }
        }
    }
    
    
    couples = NULL
    
    if (class(hmm) == "bdHMM") {
        bdHMM.settings$state2flag = rep(0, nStates)
        bdHMM.settings$couples = NULL
        bdHMM.settings = bdHMM_get_info(bdHMM.settings)
        if (any(bdHMM.settings$rev.operation < 1)) 
            stop("Reverse oparations contain dimensions < 0!\n")
        if (!is.null(bdHMM.settings$rev.operation)) {
            emissionParams = emissionRevOp(emissionParams, hmm@emission@type, bdHMM.settings, 
                nStates, hmm@directedObs)
            bdHMM.settings$rev.operation = bdHMM.settings$rev.operation - 1
        }
        if (bdHMM.settings$optim.method == "analytical") {
            bdHMM.settings$bidirOptimParams = list()
        }
        else {
            bdHMM.settings$bidirOptimParams$c2optimize = c2optimize
        }
    }
   
    
    
    emissionPrior = emissionParams$emissionPrior
    if (is.null(emissionPrior)) {
        emissionPrior = list()
    }
    else {
        emissionPrior$calldiwish = calldiwish
    }
    
    if(nCores < 1) {
        nCores = 1
    }
    hmm_out = .Call("RHMMFit", SEXPobs = obs, SEXPpi = initProb, SEXPA = transMat, 
        SEXPemission = emissionParams, SEXPtype = as.character(hmm@emission@type), 
        SEXPdim = D, SEXPregularize = as.numeric(0), SEXPk = as.integer(nStates), 
        SEXPmaxIters = as.integer(maxIters), SEXPparallel = as.integer(nCores), SEXPflags = lapply(bdHMM.settings$dirFlags, 
            as.integer), SEXPstate2flag = as.integer(bdHMM.settings$state2flag), 
        SEXPcouples = as.integer(bdHMM.settings$couples), SEXPrevop = as.integer(bdHMM.settings$rev.operation), 
        SEXPverbose = as.integer(verbose), SEXPupdateTransMat = as.integer(1), SEXPfixedEmission = emissionProbs, 
        SEXPbidiroptim = bdHMM.settings$bidirOptimParams, SEXPemissionPrior = emissionPrior, 
        SEXPeffectivezero = as.numeric(effectiveZero), SEXPconvergence = as.numeric(convergence), 
        SEXPincrementalEM = as.integer(incrementalEM), PACKAGE = "STAN")
    
    
    hmm@transMat = matrix(hmm_out$transMat, nrow = nStates, ncol = nStates, byrow = TRUE)
    hmm@initProb = hmm_out$initProb
    # print(hmm_out$emission)
    fixed.emission = list()
    if (length(fixed.emission) == 0) {
        hmm@emission = finishEmission(hmm_out$emission, hmm@emission@type, D, nStates)  #hmm_out$emission#
        
        # hmm@emission$invsigma = lapply(hmm@emission$invsigma, function(x) matrix(x,
        # nrow=D, ncol=D, byrow=TRUE))
    }
    if (!is.null(bdHMM.settings$rev.operation)) {
        bdHMM.settings$rev.operation = bdHMM.settings$rev.operation + 1
        hmm@emission@parameters = emissionRevOp(hmm@emission@parameters, hmm@emission@type, 
            bdHMM.settings, nStates, hmm@directedObs)
    }
    hmm@status = "EM"
    # hmm@converged = (fit$loglik[length(fit$loglik)] -
    # fit$loglik[length(fit$loglik)-1]) < 1e-6
    
    return(list(loglik = hmm_out$loglik, hmm = hmm))
}




#' Prepares and reformats emission for model fitting. 
#' 
#' @keywords internal
#' @noRd
prepareEmission = function(emissionParams, type) {
    if (type == "JointlyIndependent") 
        warning("emissionRevOp not implemented for JointlyIndependent!")
    if (type == "Gaussian" & (all(c("mean", "cov") %in% names(emissionParams)))) {
        emissionParams$cov = lapply(emissionParams$cov, function(x) as.matrix(x, 
            mode = "numeric"))
        emissionParams$mean = lapply(emissionParams$mean, as.numeric)
    }
    emissionParams
}

#' Processing and reformatting of emission after model fitting. 
#' 
#' @keywords internal
#' @noRd
finishEmission = function(emissionParams, type, D, nStates) {
    myEmission = list()
    if (type == "Gaussian" & (all(c("mean", "cov") %in% names(emissionParams)))) {
        emissionParams = emissionParams[c("mean", "cov")]
        emissionParams$cov = lapply(emissionParams$cov, function(x) matrix(x, nrow = D, 
            ncol = D, byrow = TRUE))
        emissionParams$mean = emissionParams$mean
        myEmission = HMMEmission(type = "Gaussian", parameters = emissionParams, 
            nStates = nStates)
    }
    if (type == "JointlyIndependent") {
        types = emissionParams$types
        
        # for(i in 1:l)
        for (i in 1:length(emissionParams$emissions)) {
            if (types[i] == "Bernoulli") {
                emissionParams$emissions[[i]] = unlist(emissionParams$emissions[[i]], 
                  recursive = FALSE)
                emissionParams$emissions[[i]] = list(p = emissionParams$emissions[[i]])
                # print(emissionParams$emissions[[i]])
            } else if (types[i] == "Gaussian") {
                # TODO names(emissionParams$emissions[[i]]) = c('mean', 'cov')
            }
            
        }
        
        myEmissionList = lapply(1:D, function(x) HMMEmission(type = types[x], parameters = emissionParams$emissions[[x]], 
            nStates = nStates))
        myEmission = HMMEmission(type = "JointlyIndependent", parameters = list(emissions = myEmissionList), 
            nStates = as.integer(nStates))
        
    }
    myEmission
}


#' Reverse operation on observations before and after model fitting.
#' 
#' @keywords internal
#' @noRd
emissionRevOp = function(emissionParams, type, bdHMM.settings, nStates, directedObs) {
    if (type == "JointlyIndependent") 
        warning("emissionRevOp not implemented for JointlyIndependent!")
    if (type == "Gaussian" & (all(c("mean", "cov") %in% names(emissionParams)))) {
        for (i in 1:nStates) {
            if (bdHMM.settings$state2flag[i] == 1) {
                emissionParams$mean[[i]] = emissionParams$mean[[i]][bdHMM.settings$rev.operation]
                emissionParams$cov[[i]] = emissionParams$cov[[i]][bdHMM.settings$rev.operation, 
                  bdHMM.settings$rev.operation]
            }
        }
    }
    emissionParams
}

#' Prepares settings for bdHMM transition optimization. 
#' 
#' @keywords internal
#' @noRd
bdHMM_get_info = function(bdHMM.settings) {
    
    
    stateLabel = bdHMM.settings$stateLabel
    myF = sapply(strsplit(stateLabel[grep("F", stateLabel)], "F"), function(x) x[2])
    R = sapply(strsplit(stateLabel[grep("R", stateLabel)], "R"), function(x) x[2])
    U = sapply(strsplit(stateLabel[grep("U", stateLabel)], "U"), function(x) x[2])
    nStates = length(stateLabel)
    dirStates = (nStates - length(U))/2
    undirStates = length(U)
    myU = ""
    if (length(U) > 0) {
        myU = paste("U", 1:undirStates, sep = "")
    }
    if (length(myF) != length(R)) {
        stop("Unequal numbers of states with forward and reverse orientation!")
    }
    
    bdHMM.settings$couples = (1:nStates) - 1
    if (dirStates > 0) {
        bdHMM.settings$couples[1:(2 * dirStates)] = c((1:dirStates) + dirStates, 
            1:dirStates) - 1
    }
    
    
    bdHMM.settings$state2flag[grep("F", stateLabel)] = 1
    bdHMM.settings$state2flag[grep("R", stateLabel)] = -1
    bdHMM.settings$state2flag[which(bdHMM.settings$state2flag == 0)] = 100
    bdHMM.settings$state2flag = -bdHMM.settings$state2flag
    
    
    
    if (bdHMM.settings$bidirectional.mc) {
        
        eqB = c(rep(0, nStates), 0)
        statD = rep(1/nStates, nStates)
        xsi_sum = matrix(as.numeric(NA), nrow = nStates, ncol = nStates)
        initGamma = as.numeric(rep(0, nStates))
        constraints = getConstraints(stateLabel)
        x0 = rep(1/nStates, length(constraints) + ((nStates - length(U))/2) + length(U))
        LB = rep(0, length(x0))
        UB = rep(1, length(x0))
        control = list(trace = 0, tol = 1e-06)  ## => Likelihood may decrease if estimate is not accurate enough ( solution: choose relative convergence criterion for optimizer )
        transMat = matrix(1/nStates, nrow = nStates, ncol = nStates)
        bdHMM.settings$bidirOptimParams = list(x0 = x0, eqB = eqB, xsi_sum = xsi_sum, 
            initGamma = initGamma, statD = statD, constraints = constraints, nStates = nStates, 
            LB = LB, UB = UB, control = control, transMat = transMat, stateLabel = stateLabel, 
            method = bdHMM.settings$optim.method, objective = as.numeric(Inf), nrm = as.integer(0), 
            couples = bdHMM.settings$couples, update = as.integer(0))
        
    }
    
    bdHMM.settings
} 
