



#'
#' This class is a generic container for different emission functions of Hidden Markov Models.
#' 
#' @slot type The type of emission function c('Gaussian').
#' @slot parameters A list containing the the parameters for each state.
#' @slot dim Number of dimensions.
#' @slot nStates The number of states.
#' @examples
#' nStates = 5
#' means = list(4,11,4,11,-1)
#' Sigma = lapply(list(4,4,4,4,4), as.matrix)
#' transMat = matrix(1/nStates, nrow=nStates, ncol=nStates)
#' initProb = rep(1/nStates, nStates)
# 
#' HMMEmission(type='Gaussian', parameters=list(mean=means, cov=Sigma), nStates=length(means))
#' @exportClass HMMEmission
.HMMEmission <- setClass("HMMEmission", representation(type = "character", parameters = "list", 
    dim = "integer", nStates = "integer"))

#' Internal function used to inititialize an HMMEmission object.
#' 
#' @keywords internal
#' @noRd
checkParameters = function(type, parameters, nStates) {
    val = FALSE
    if (class(parameters) != "list") 
        stop("Emission parameters must be submitted as a list!")
    
    if (type == "Gaussian") {
        if (!(all(c("mean", "cov") %in% names(parameters)))) 
            stop("Gaussian emissions must contain parameters \"mean\" and \"cov\"!")
        means = parameters$mean
        covs = parameters$cov
        if (!(class(means) == "list" & class(covs) == "list")) 
            stop("Means and covariances need to be of class list!")
        if (!(nStates == length(means) & nStates == length(covs))) 
            stop("Number of states does not match dimensions of means and covariances!")
        for (i in 1:length(means)) {
            stopifnot(dim(covs[[i]])[1] == dim(covs[[i]])[2] & length(means[[i]]) == 
                dim(covs[[i]])[1])
            stopifnot((is.numeric(covs[[i]])) & (is.numeric(means[[i]])))
        }
        val = TRUE
    } else if (type == "Bernoulli") {
        if (class(parameters$p) != "list") 
            stop("Bernoulli emission parameters must be submitted as a list (one entry per state)!")
        if (!"p" %in% names(parameters)) 
            stop("Bernoulli emission parameters must contain entry p storing emission probabilities!")
        myLens = sapply(parameters$p, function(x) length(x) == 1)
        if (!all(myLens)) 
            stop(paste("Dimensionality of Bernoulli probabilities unequal to one for states: ", 
                which(!myLens), sep = ""))
        
        val = TRUE
    } else if (type == "JointlyIndependent") {
        if (class(parameters$emissions) != "list") 
            stop("Jointly Independent emission parameters must be submitted as a list (one entry per emission function)!")
        myClasses = sapply(parameters$p, function(x) class(x) == "HMMEmission")
        if (!all(myClasses)) 
            stop(paste("Parameters of Jointly Independent emissions are not of class HMMEmission: ", 
                which(!myLens), sep = ""))
        
        val = TRUE
    } else if (type == "null") {
        # cat('You are not using an Emission Function! Must supply fixed Emission
        # probabilities!')
        val = TRUE
    } else {
        stop("Unknown emission function specified: ", type, "\n", sep = "")
    }
    val
}

#' Internal function used to convert Bernoulli representation in Independent Model
#' 
#' @keywords internal
#' @noRd
changeParams = function(type, parameters) {
    if (type == "JointlyIndependent") {
        myDs = sapply(parameters$emissions, function(x) x@dim)
        if (!all(myDs > 0)) 
            stop("Dimensions must be greater than 0!")
        myDs = cumsum(myDs)
        emissionDims = list()
        for (i in 1:length(myDs)) {
            if (i == 1) {
                emissionDims[[i]] = 1:myDs[i]
            } else {
                emissionDims[[i]] = (myDs[i - 1] + 1):myDs[i]
            }
            
        }
        parameters$emissionDim = emissionDims
        parameters$types = sapply(parameters$emissions, function(x) x@type)
    }
    return(parameters)
    
}



#' Get dimensions of the emission function
#' 
#' @keywords internal
#' @noRd
getDim = function(type, parameters) {
    dim = NULL
    if (type == "Gaussian") {
        dim = as.integer(length(parameters$mean[[1]]))
    }
    if (type == "Bernoulli") {
        dim = as.integer(length(parameters$p[[1]]))
    }
    if (type == "JointlyIndependent") {
        dim = as.integer(sum(sapply(parameters$emissions, function(x) x@dim)))
    }
    if (type == "null") {
        if (is.null(parameters$dim)) 
            stop("Must supply number of dimensions in parameter list if emission type is 'null'!")
        dim = parameters$dim
    }
    dim
}

#'
#' This function creates a HMMEmission object.
#' 
#' @title Create a HMMEmission object
#' 
#' @param type The type of emission function c('Gaussian').
#' @param parameters A list containing the the parameters for each state.
#' @param nStates The number of states.
#' @examples
#' nStates = 5
#' means = list(4,11,4,11,-1)
#' Sigma = lapply(list(4,4,4,4,4), as.matrix)
#' transMat = matrix(1/nStates, nrow=nStates, ncol=nStates)
#' initProb = rep(1/nStates, nStates)
# 
#' HMMEmission(type='Gaussian', parameters=list(mean=means, cov=Sigma), nStates=length(means))
#' @export HMMEmission
HMMEmission = function(type=character(), parameters=list(), nStates=integer()) {
    if (length(nStates) == 0) 
        stop("Must specify number of states for emission!")
    if (!checkParameters(type, parameters, nStates)) 
        stop("Parameters not correctly specified for distribution ", type)
    param = changeParams(type, parameters)
    dim = getDim(type, parameters)
    .HMMEmission(type = type, parameters = param, dim = dim, nStates = as.integer(nStates))
}


#' Prints description of HMMEmission object
#' 
#' @keywords internal
#' @noRd
setMethod(f="show", signature=c("HMMEmission"), function(object) {
    cat("An object of class \"", class(object), "\"\n", " type: ", object@type, "\n dim: ", 
        paste(rbind(names(object@parameters$order), if (length(object@parameters$order) > 
            0) {
            ": "
        }, object@dim, if (length(object@parameters$order) > 0) {
            "; "
        })), "\n nStates: ", object@nStates, "\n", sep = "")
})


#'
#' This class is a generic container for Hidden Markov Models.
#' 
#' @slot initProb Initial state probabilities.
#' @slot transMat Transition probabilities
#' @slot emission Emission parameters as an HMMEmission object.
#' @slot nStates Number of states.
#' @slot status of the HMM. On of c('initial', 'EM').
#' @examples 
#' nStates = 5
#' means = list(4,11,4,11,-1)
#' Sigma = lapply(list(4,4,4,4,4), as.matrix)
#' transMat = matrix(1/nStates, nrow=nStates, ncol=nStates)
#' initProb = rep(1/nStates, nStates)
# 
#' HMM(initProb=initProb, transMat=transMat, emission=HMMEmission(type='Gaussian', parameters=list(mean=means, cov=Sigma), nStates=length(means)), nStates=nStates, status='initial')
# 
#' @seealso \code{\linkS4class{HMMEmission}}
#' @exportClass HMM
.HMM <- setClass("HMM", representation(initProb = "numeric", transMat = "matrix", 
    emission = "HMMEmission", nStates = "numeric", status = "character"))

#'
#' This function creates a HMM object.
#' @title Create a HMM object
#' 
#' @param initProb Initial state probabilities.
#' @param transMat Transition probabilities
#' @param emission Emission parameters as an HMMEmission object.
#' @param nStates Number of states.
#' @param status of the HMM. On of c('initial', 'EM').
#' 
#' 
#' @examples 
#' nStates = 5
#' means = list(4,11,4,11,-1)
#' Sigma = lapply(list(4,4,4,4,4), as.matrix)
#' transMat = matrix(1/nStates, nrow=nStates, ncol=nStates)
#' initProb = rep(1/nStates, nStates)
# 
#' HMM(initProb=initProb, transMat=transMat, emission=HMMEmission(type='Gaussian', parameters=list(mean=means, cov=Sigma), nStates=length(means)), nStates=nStates, status='initial')
# 
#' @seealso \code{\linkS4class{HMMEmission}}
#' @export HMM
HMM <- function(initProb=numeric(), transMat=matrix(numeric(), ncol=1, nrow=1), 
    emission, nStates = numeric(), status = character()) {
    .HMM(initProb = initProb, transMat = transMat, emission = emission, nStates = nStates, 
        status = status)
}




#' Prints description of HMM object
#' 
#' @keywords internal
#' @noRd
setMethod(f="show", signature=c("HMM"), function(object) {
    cat("An object of class \"", class(object), "\" \n", sep = "")
    cat(" nStates: ", object@nStates, "\n", sep = "")
    cat(" emission function: ", object@emission@type, "\n", sep = "")
    cat(" dim: HMM", paste(rbind(names(object@emission@parameters$order), if (length(object@emission@parameters$order) > 
        0) {
        ": "
    }, object@emission@dim, if (length(object@emission@parameters$order) > 0) {
        "; "
    })), "\n", sep = "")
    cat(" status: ", object@status, "\n", sep = "")
})


#'
#' This class is a generic container for bidirectional Hidden Markov Models.
#' 
#' @slot initProb Initial state probabilities.
#' @slot transMat Transition probabilities
#' @slot emission Emission parameters as an HMMEmission object.
#' @slot nStates Number of states.
#' @slot stateLabel Indicates directinality of states. States can be forward (F1, F2, ..., Fn), reverse (R1, R2, ..., Rn) or undirectional (U1, U2, ..., Um). Number of F and R states must be equal and twin states are indicated by integers in id (e.g. F1 and R1 and twins).
#' @slot transitionsOptim There are three methods to choose from for fitting the transitions. Bidirectional transition matrices (invariant under reversal of time and direction) can be fitted using c('rsolnp', 'ipopt'). 'None' uses standard update formulas and the resulting matrix is not constrained to be bidirectional.
#' @slot directedObs An integer indicating which dimensions are directed. Undirected dimensions are 0. Directed observations must be marked as unique integer pairs. For instance c(0,0,0,0,0,1,1,2,2,3,3) contains 5 undirected observations, and thre pairs (one for each direction) of directed observations.
#' @slot status Status of the bdHMM. 'Initial' means that the model was not fitted yet. 'EM' means that the model was optimized using Expectation maximization.
#' @examples 
#' nStates = 5
#' stateLabel = c('F1', 'F2', 'R1', 'R2', 'U1')
#' means = list(4,11,4,11,-1)
#' Sigma = lapply(list(4,4,4,4,4), as.matrix)
#' transMat = matrix(1/nStates, nrow=nStates, ncol=nStates)
#' initProb = rep(1/nStates, nStates)
#' 
#' bdHMM(initProb=initProb, transMat=transMat, emission=HMMEmission(type='Gaussian', parameters=list(mean=means, cov=Sigma), nStates=length(means)), nStates=nStates, status='initial', stateLabel=stateLabel, transitionsOptim='none', directedObs=as.integer(0))
#' @seealso \code{\linkS4class{HMMEmission}}
#' @exportClass bdHMM
.bdHMM <- setClass("bdHMM", contains = "HMM", representation(stateLabel = "character", 
    transitionsOptim = "character", directedObs = "integer"))


#' Checks equality of emission of twin states
#' 
#' @keywords internal
#' @noRd
checkEmissionTwins = function(twin1, twin2, emission, rev.operation, stateLabel) {
    val = FALSE
    if (emission@type == "Gaussian" & (all(c("mean", "cov") %in% names(emission@parameters)))) {
        means = emission@parameters$mean
        covs = emission@parameters$cov
        if (!all(means[[twin1]] - means[[twin2]][rev.operation] < 1e-12)) 
            stop("Means of twins (", stateLabel[twin1], ", ", stateLabel[twin2], 
                ") are not equal!")
        # if(!all(covs[[twin1]] == covs[[twin2]][rev.operation, rev.operation]))
        # stop('Covariances of twins (', stateLabel[twin1], ', ', stateLabel[twin2], ')
        # are not equal!')
        
        val = TRUE
    }
    if (emission@type == "Independent") {
        order = emission@parameters$order
        
        if (all(c("mean", "cov") %in% names(emission@parameters$Gaussian))) {
            means = emission@parameters$Gaussian$mean
            covs = emission@parameters$Gaussian$cov
            ## recalculate rev.operation like function emissionRevOp in fitHMM
            ## if(!all(means[[twin1]] == means[[twin2]][rev.operation])) stop('Means of twins
            ## (', stateLabel[twin1], ', ', stateLabel[twin2], ') are not equal!')
            ## if(!all(covs[[twin1]] == covs[[twin2]][rev.operation, rev.operation]))
            ## stop('Covariances of twins (', stateLabel[twin1], ', ', stateLabel[twin2], ')
            ## are not equal!')
            val = TRUE
        }
        if (c("p") %in% names(emission@parameters$Bernoulli)) {
            ## recalculate changing partner from rev.operation like function emissionRevOp in
            ## fitHMM
            p = emission@parameters$Bernoulli$p
            val = TRUE
        }
    }
    val
}


#' Checks validity of a bdHMM object
#' 
#' @keywords internal
#' @noRd
valid.bdHMM = function(mybdHMM) {
    myF = grep("F", mybdHMM@stateLabel)
    myR = grep("R", mybdHMM@stateLabel)
    F2R = gsub("F", "R", mybdHMM@stateLabel[myF])
    if (!all(F2R %in% mybdHMM@stateLabel[myR])) 
        stop("No matching between twin states possible!")
    couples = rep(-1, length(mybdHMM@stateLabel))
    if (length(myF) > 0) {
        for (i in 1:length(myF)) {
            myTwin = which(mybdHMM@stateLabel == F2R[i])
            couples[myF[i]] = myTwin
            couples[myTwin] = myF[i]
        }
    }
    
    dim = sum(getDim(mybdHMM@emission@type, mybdHMM@emission@parameters))
    rev.operation = 1:dim
    mapping.rev.op = c()
    changed <- c()
    
    # rev.operation = 1:mybdHMM@emission@dim
    
    if (any(table(mybdHMM@directedObs[mybdHMM@directedObs != 0]) > 2)) 
        stop("Matching of directed observations is not possible!")
    for (i in 1:length(mybdHMM@directedObs)) {
        if (mybdHMM@directedObs[i] == 0) {
            rev.operation[i] = i
            
            if (mybdHMM@emission@type == "Independent") {
                mapping.rev.op <- append(mapping.rev.op, which(mybdHMM@emission@parameters$order$Gaussian == 
                  i))
            }
        } else {
            myPair = which(mybdHMM@directedObs == mybdHMM@directedObs[i])
            myPair = myPair[myPair != i]
            rev.operation[i] = myPair
            if (mybdHMM@emission@type == "Independent") {
                if (!(myPair %in% changed) && !(i %in% changed)) {
                  mapping.rev.op <- append(mapping.rev.op, which(mybdHMM@emission@parameters$order$Gaussian == 
                    myPair))
                  mapping.rev.op <- append(mapping.rev.op, which(mybdHMM@emission@parameters$order$Gaussian == 
                    i))
                  changed <- append(changed, c(myPair, i))
                }
            }
        }
    }
    if (mybdHMM@emission@type == "Independent") {
        rev.operation = mapping.rev.op
    }
    for (i in 1:length(couples)) {
        if (couples[i] != -1) {
            checkEmissionTwins(i, couples[i], mybdHMM@emission, rev.operation, mybdHMM@stateLabel)
        }
    }
}

#'
#' This function creates a bdHMM function.
#' 
#' @title Create a bdHMM object
#' 
#' @param initProb Initial state probabilities.
#' @param transMat Transition probabilities
#' @param emission Emission parameters as an HMMEmission object.
#' @param nStates Number of states.
#' @param stateLabel Indicates directinality of states. States can be forward (F1, F2, ..., Fn), reverse (R1, R2, ..., Rn) or undirectional (U1, U2, ..., Um). Number of F and R states must be equal and twin states are indicated by integers in id (e.g. F1 and R1 and twins).
#' @param transitionsOptim There are three methods to choose from for fitting the transitions. Bidirectional transition matrices (invariant under reversal of time and direction) can be fitted using c('rsolnp', 'ipopt'). 'None' uses standard update formulas and the resulting matrix is not constrained to be bidirectional.
#' @param directedObs An integer indicating which dimensions are directed. Undirected dimensions are 0. Directed observations must be marked as unique integer pairs. For instance c(0,0,0,0,0,1,1,2,2,3,3) contains 5 undirected observations, and thre pairs (one for each direction) of directed observations.
#' @param status Status of the bdHMM. 'Initial' means that the model was not fitted yet. 'EM' means that the model was optimized using Expectation maximization.
#' @examples 
#' nStates = 5
#' stateLabel = c('F1', 'F2', 'R1', 'R2', 'U1')
#' means = list(4,11,4,11,-1)
#' Sigma = lapply(list(4,4,4,4,4), as.matrix)
#' transMat = matrix(1/nStates, nrow=nStates, ncol=nStates)
#' initProb = rep(1/nStates, nStates)
#' 
#' bdHMM(initProb=initProb, transMat=transMat, emission=HMMEmission(type='Gaussian', parameters=list(mean=means, cov=Sigma), nStates=length(means)), nStates=nStates, status='initial', stateLabel=stateLabel, transitionsOptim='none', directedObs=as.integer(0))
#' @seealso \code{\linkS4class{HMMEmission}}
#' @export bdHMM
bdHMM = function(initProb=numeric(), transMat=matrix(numeric(), ncol=1, nrow=1), 
    emission, nStates = numeric(), status = character(), stateLabel = character(), 
    transitionsOptim = character(), directedObs = integer()) {
    mybdHMM = .bdHMM(HMM(initProb = initProb, transMat = transMat, emission = emission, 
        nStates = nStates, status = status), stateLabel = stateLabel, transitionsOptim = transitionsOptim, 
        directedObs = directedObs)
    valid.bdHMM(mybdHMM)
    mybdHMM
}

#' Prints description of a bdHMM
#' 
#' @keywords internal
#' @noRd
setMethod(f="show", signature=c("bdHMM"), function(object) {
    cat("An object of class \"", class(object), "\" \n", sep = "")
    cat(" nStates: ", object@nStates, "\n", sep = "")
    cat(" emission function: ", object@emission@type, "\n", sep = "")
    cat(" dim.obs: ", paste(object@emission@dim, collapse = "; "), "\n", sep = "")
    cat(" directedObs: ", paste(object@directedObs, collapse = " "), "\n", sep = "")
    cat(" status: ", object@status, "\n", sep = "")
    cat(" transitionsOptim: ", object@transitionsOptim, "\n", sep = "")
    cat(" stateLabel: ", paste(object@stateLabel, collapse = " "), "\n", sep = "")
}) 
