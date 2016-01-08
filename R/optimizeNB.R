
#' Initializes NegativeBinomial emission functions (HMM or bdHMM) given a k-means clustering of the data.
#' 
#' @keywords internal
#' @noRd
initNB = function(km, signalbychrom, celltypes = NULL, stateLabels = NULL, 
    directedObs = NULL, sizeFactor = NULL, zeroInflated = FALSE) {
    celltypes = list(CD4T = 1)
    myAvg = apply(t(sapply(celltypes, function(x) apply(do.call("rbind", 
        signalbychrom[x]), 2, sum))), 2, mean)
    
    if (is.null(celltypes)) {
        celltypes = list(1:length(signalbychrom))
    }
    
    nStates = nrow(km$centers)
    myD = ncol(km$centers)
    myMat = do.call("rbind", signalbychrom)
    
    if (is.null(sizeFactor)) {
        sizeFactor = matrix(1, ncol = myD, nrow = length(signalbychrom))
    }
    # print(sizeFactor)
    
    cl2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
    cl2pos = list()
    cl2val = list()
    cl2init = list()
    initMeans_humanCD4T = list()
    initCovs_humanCD4T = list()
    for (i in 1:length(cl2ind)) {
        # print(i)
        cl2pos[[i]] = cl2ind[[i]]
        cl2val[[i]] = myMat[cl2pos[[i]], ]
        if (ncol(myMat) == 1) {
            cl2val[[i]] = matrix(cl2val[[i]], ncol = 1)
        }
        mySplit = list(lapply(1:ncol(cl2val[[i]]), function(d) as.list(tapply(1:nrow(cl2val[[i]]), 
            INDEX = cl2val[[i]][, d], identity))))
        # mySplit = list(lapply(1:ncol(cl2val[[i]]), function(d)
        # as.list(tapply(1:nrow(cl2val[[i]]), INDEX=cl2val[[i]][,d],
        # identity)))) print(str(mySplit)) cl2init[[i]] =
        # sapply(1:ncol(cl2val[[i]]), function(x)
        # optimizePoiLog(list(mu=mean(log(c(1e-3, cl2val[[i]][(cl2val[[i]][,x]
        # != 0),x]))), sigma=1, gamma=rep(1,nrow(myMat)), d=x,
        # sizeFactor=mySF[,x], uniqueCountSplit=mySplit)))
        if (zeroInflated) {
            cl2init[[i]] = sapply(1:ncol(cl2val[[i]]), function(x) optimizeNBInit(list(mu = mean(cl2val[[i]][, 
                x]), size = mean(cl2val[[i]][, x]), pi = max(1e-06, length(which(cl2val[[i]][, 
                x] == 0))/length(cl2val[[i]][, x])), gamma = rep(1, nrow(myMat)), 
                sizeFactor = sizeFactor, uniqueCountSplit = mySplit, d = x)))
        } else {
            cl2init[[i]] = sapply(1:ncol(cl2val[[i]]), function(x) optimizeNBInit(list(mu = mean(cl2val[[i]][, 
                x]), size = mean(cl2val[[i]][, x]), pi = 0, gamma = rep(1, 
                nrow(myMat)), sizeFactor = sizeFactor, uniqueCountSplit = mySplit, 
                d = x)))
        }
        # print(cl2init[[i]])
    }
    
    
    myEmissions = list()
    
    if (is.null(stateLabels)) {
        for (currDim in 1:myD) {
            currParameters = list(mu = lapply(1:nStates, function(state) cl2init[[state]][1, 
                currDim]), size = lapply(1:nStates, function(state) cl2init[[state]][2, 
                currDim]), sizeFactor = lapply(1:nStates, function(x) sizeFactor[, 
                currDim]), pi = lapply(1:nStates, function(state) cl2init[[state]][3, 
                currDim]))
            myEmissions[[currDim]] = HMMEmission(type = "NegativeBinomial", 
                parameters = currParameters, nStates = nStates)
        }
    } else {
        if (!length(stateLabels) == 2 * nStates) 
            stop("Length if object stateLabels must be 2*nStates")
        
        nStates = length(stateLabels)
        rev.operation = c()
        for (i in 1:length(directedObs)) {
            if (directedObs[i] == 0) {
                rev.operation[i] = i
            } else {
                myPair = which(directedObs == directedObs[i])
                myPair = myPair[myPair != i]
                rev.operation[i] = myPair
            }
        }
        cl2init = cl2init[c(rep(grep("F", stateLabels), 2), grep("U", stateLabels))]
        myR = grep("R", stateLabels)
        cl2init[myR] = lapply(cl2init[myR], function(x) x[, rev.operation])
        for (currDim in 1:myD) {
            currParameters = list(mu = lapply(1:nStates, function(state) cl2init[[state]][1, 
                currDim]), size = lapply(1:nStates, function(state) cl2init[[state]][2, 
                currDim]), sizeFactor = lapply(1:nStates, function(x) sizeFactor[, 
                currDim]), pi = lapply(1:nStates, function(state) cl2init[[state]][3, 
                currDim]))
            myEmissions[[currDim]] = HMMEmission(type = "NegativeBinomial", 
                parameters = currParameters, nStates = nStates)
        }
    }
    
    myEmissions
}


#' Optimizes parameters of NegativeBinomial emission functions during EM-learning.
#' 
#' @keywords internal
#' @noRd
optimizeNB = function(myPars) {
    # print(myPars$sizeFactor)
    mySplit = myPars$uniqueCountSplit
    out = NULL
    myFile = file.path(tempdir(), paste("STAN.temp.params.", Sys.getpid(), 
        ".rda", sep = ""))
    if (myPars$d != 1) {
        load(myFile)
    } else {
        load(myFile)
        out = mclapply(1:length(mySplit[[1]]), function(currD) {
            gammaSums = list()
            myUniques = list()
            sizeFactors = list()
            for (i in 1:length(mySplit)) {
                gammaSums[[i]] = sapply(mySplit[[i]][[currD]], function(x) sum(myPars$gamma[x]))
                myUniques[[i]] = as.integer(names(mySplit[[i]][[currD]]))
                gammaSums[[i]] = unlist(gammaSums[[i]])
                # print(i) print(currD) bdHMM
                if (i > nrow(myPars$sizeFactor)) {
                  sizeFactors[[i]] = rep(myPars$sizeFactor[i/2, currD], 
                    length(gammaSums[[i]]))
                } else {
                  sizeFactors[[i]] = rep(myPars$sizeFactor[i, currD], length(gammaSums[[i]]))
                }
                
            }
            gammaSums = unlist(gammaSums)
            myUniques = unlist(myUniques)
            sizeFactors = unlist(sizeFactors)
            
            pos2fac = tapply(1:length(sizeFactors), INDEX = sizeFactors, 
                identity)
            fac2unique = lapply(pos2fac, function(x) myUniques[x])
            fac2gamma = lapply(pos2fac, function(x) gammaSums[x])
            pos2fac = tapply(1:length(sizeFactors), INDEX = sizeFactors, 
                identity)
            sf_collapsed = list()
            gamma_collapsed = list()
            unique_collapsed = list()
            for (i in 1:length(pos2fac)) {
                gamma_collapsed[[i]] = tapply(fac2gamma[[i]], INDEX = fac2unique[[i]], 
                  sum)
                unique_collapsed[[i]] = as.numeric(names(gamma_collapsed[[i]]))
                sf_collapsed[[i]] = rep(as.numeric(names(pos2fac)[i]), 
                  length(gamma_collapsed[[i]]))
            }
            gammaSums = unlist(gamma_collapsed)
            myUniques = unlist(unique_collapsed)
            sizeFactors = unlist(sf_collapsed)
            # print(sd(log(myUniques+1)))
            
            currPi = myparams$pi[myPars$currstate, currD]
            if (currPi == 0) {
                currPi = Inf
            }
            
            initPars = c(myparams$mu[myPars$currstate, currD], myparams$size[myPars$currstate, 
                currD], currPi)
            # initPars = c(sum(gammaSums*(myUniques))/sum(gammaSums),
            # sum(gammaSums*(myUniques))/sum(gammaSums),
            # currPi)#max(c(log(sum(gammaSums*myUniques)/sum(gammaSums))))
            # print(initPars) before = myQPoiLog(par=initPars,myGammas=gammaSums,
            # sizeFactor=sizeFactors, x=myUniques) print(initPars)
            out = optim(par = initPars, fn = myQNBinom, zeroInflation = ifelse(initPars[3] == 
                Inf, FALSE, TRUE), myGammas = gammaSums, sizeFactor = sizeFactors, 
                x = myUniques, lower = c(1e-06, 1e-06, 1e-06), upper = c(1e+06, 
                  1e+06, 1), method = "L-BFGS-B")$par
            if (currPi == Inf) {
                out[3] = 0
            }
            out
        }, mc.cores = min(c(myPars$ncores, length(mySplit[[1]]))))
        # print(unlist(sapply(out, function(x) x[1])))
        myparams$mu[myPars$currstate, ] = unlist(sapply(out, function(x) x[1]))
        myparams$size[myPars$currstate, ] = unlist(sapply(out, function(x) x[2]))
        myparams$pi[myPars$currstate, ] = unlist(sapply(out, function(x) x[3]))
        # print(myparams)
        save(myparams, file = myFile)
    }
    # print(myparams$mu)
    out = c(myparams$mu[myPars$currstate, myPars$d], myparams$size[myPars$currstate, 
        myPars$d], myparams$pi[myPars$currstate, myPars$d])
    # print(out) out = optim(par=c(myPars$mu, myPars$sigma), fn=myQPoiLog,
    # myGammas=gammaSums, sizeFactor=sizeFactors, x=myUniques,
    # method='BFGS')$par#, control=list(reltol=1e-8, maxit=10000))$par
    # cat('mean: ', log(sum(gammaSums*myUniques)/sum(gammaSums)), '
    # par[2]=', out[1], '\n', sep='') print(out) out[2] = abs(out[2])
    # after = myQPoiLog(par=out,myGammas=gammaSums, sizeFactor=sizeFactors,
    # x=myUniques) print(before-after) if(c(before-after) < 0.1) {
    # print(sum(gammaSums)) out = c(myPars$mu, myPars$sigma) } if(myPars$d
    # == length(mySplit[[1]])) { system(paste('rm ', myFile, sep='')) }
    
    out
}


#' Optimizes parameters of NegativeBinomial emission functions during initialization. Called by initNB.
#' 
#' @keywords internal
#' @noRd
optimizeNBInit = function(myPars) {
    gammaSums = list()
    myUniques = list()
    sizeFactors = list()
    mySplit = myPars$uniqueCountSplit
    for (i in 1:length(mySplit)) {
        gammaSums[[i]] = sapply(mySplit[[i]][[myPars$d]], function(x) sum(myPars$gamma[x]))
        myUniques[[i]] = as.integer(names(mySplit[[i]][[myPars$d]]))
        gammaSums[[i]] = unlist(gammaSums[[i]])
        # print(i)
        sizeFactors[[i]] = rep(myPars$sizeFactor[i, myPars$d], length(gammaSums[[i]]))
    }
    # print(myPars$gamma)
    gammaSums = unlist(gammaSums)
    myUniques = unlist(myUniques)
    # print(gammaSums) print(myUniques)
    sizeFactors = unlist(sizeFactors)
    
    
    pos2fac = tapply(1:length(sizeFactors), INDEX = sizeFactors, identity)
    fac2unique = lapply(pos2fac, function(x) myUniques[x])
    fac2gamma = lapply(pos2fac, function(x) gammaSums[x])
    pos2fac = tapply(1:length(sizeFactors), INDEX = sizeFactors, identity)
    sf_collapsed = list()
    gamma_collapsed = list()
    unique_collapsed = list()
    for (i in 1:length(pos2fac)) {
        gamma_collapsed[[i]] = tapply(fac2gamma[[i]], INDEX = fac2unique[[i]], 
            sum)
        unique_collapsed[[i]] = as.numeric(names(gamma_collapsed[[i]]))
        sf_collapsed[[i]] = rep(as.numeric(names(pos2fac)[i]), length(gamma_collapsed[[i]]))
    }
    gammaSums = unlist(gamma_collapsed)
    myUniques = unlist(unique_collapsed)
    sizeFactors = unlist(sf_collapsed)
    
    if (myPars$pi == 0) {
        myPars$pi = Inf
    }
    
    out = optim(par = c(myPars$mu, myPars$size, myPars$pi), fn = myQNBinom, 
        zeroInflation = ifelse(myPars$pi == Inf, FALSE, TRUE), myGammas = gammaSums, 
        sizeFactor = sizeFactors, x = myUniques, lower = c(1e-06, 1e-06, 
            1e-06), upper = c(1e+06, 1e+06, 1), method = "L-BFGS-B")$par
    if (myPars$pi == Inf) {
        out[3] = 0
    }
    
    out
}


#' Optimizes parameters of NegativeBinomial emission functions during EM-learning or initialization. Called by optimizeNB and optimizeNBInit.
#' 
#' @keywords internal
#' @noRd
myQNBinom <- function(x, par, myGammas = rep(1, length(x)), sizeFactor = 1, 
    zeroInflation = FALSE) {
    mynegbinom = dnbinom(x, mu = par[1]/sizeFactor, size = par[2])
    if (zeroInflation) {
        myZeros = which(x == 0)
        notZeros = which(x != 0)
        mynegbinom[myZeros] = log(par[3] + ((1 - par[3]) * mynegbinom[myZeros]) + 
            1e-300)
        mynegbinom[notZeros] = log(((1 - par[3]) * mynegbinom[notZeros]) + 
            1e-300)
    } else {
        mynegbinom = log(mynegbinom + 1e-300)
    }
    
    myval = -sum(myGammas * mynegbinom)
    myval
}
 
