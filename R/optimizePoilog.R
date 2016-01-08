
#' Initializes PoissonLogNormal emission functions (HMM or bdHMM) given a k-means clustering of the data.
#' 
#' @keywords internal
#' @noRd
initPoiLog = function(km, signalbychrom, celltypes = NULL, stateLabels = NULL, 
    directedObs = NULL, sizeFactor = NULL) {
    celltypes = list(CD4T = 1)
    myAvg = apply(t(sapply(celltypes, function(x) apply(do.call("rbind", 
        signalbychrom[x]), 2, sum))), 2, mean)
    
    if (is.null(celltypes)) {
        celltypes = list(1:length(signalbychrom))
    }
    # mySF = matrix(NA, nrow=length(signalbychrom),
    # ncol=ncol(signalbychrom[[1]])) for(cell in 1:length(celltypes)) {
    # for(currDim in 1:ncol(mySF)) { mySF[celltypes[[cell]],currDim] =
    # 1#myAvg[currDim]/sum(do.call('rbind',
    # signalbychrom[celltypes[[cell]]])[,currDim]) } }
    
    myD = ncol(km$centers)
    if (is.null(sizeFactor)) {
        sizeFactor = matrix(1, ncol = myD, nrow = length(signalbychrom))
    }
    
    nStates = nrow(km$centers)
    myMat = do.call("rbind", signalbychrom)
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
        # print(str(mySplit))
        cl2init[[i]] = sapply(1:ncol(cl2val[[i]]), function(x) optimizePoiLogInit(list(mu = mean(log(c(0.001, 
            cl2val[[i]][(cl2val[[i]][, x] != 0), x]))), sigma = 1, gamma = rep(1, 
            nrow(myMat)), d = x, sizeFactor = sizeFactor, uniqueCountSplit = mySplit)))
    }
    
    
    myEmissions = list()
    
    if (is.null(stateLabels)) {
        for (currDim in 1:myD) {
            currParameters = list(mu = lapply(1:nStates, function(state) cl2init[[state]][1, 
                currDim]), sigma = lapply(1:nStates, function(state) cl2init[[state]][2, 
                currDim]))  #, sizeFactor=lapply(1:nStates, function(x) sizeFactor[,currDim]))
            myEmissions[[currDim]] = HMMEmission(type = "PoissonLogNormal", 
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
                currDim]), sigma = lapply(1:nStates, function(state) cl2init[[state]][2, 
                currDim]))  #, sizeFactor=lapply(1:nStates, function(x) sizeFactor[,currDim]))
            myEmissions[[currDim]] = HMMEmission(type = "PoissonLogNormal", 
                parameters = currParameters, nStates = nStates)
        }
    }
    
    myEmissions
}



#' Optimizes parameters of PoissonLogNormal emission functions during initialization. Called by initNB.
#' 
#' @keywords internal
#' @noRd
optimizePoiLogInit = function(myPars) {
    gammaSums = list()
    myUniques = list()
    sizeFactors = list()
    mySplit = myPars$uniqueCountSplit
    for (i in 1:length(mySplit)) {
        gammaSums[[i]] = sapply(mySplit[[i]][[myPars$d]], function(x) sum(myPars$gamma[x]))
        myUniques[[i]] = as.integer(names(mySplit[[i]][[myPars$d]]))
        gammaSums[[i]] = unlist(gammaSums[[i]])
        sizeFactors[[i]] = rep(myPars$sizeFactor[i, myPars$d], length(gammaSums[[i]]))
    }
    gammaSums = unlist(gammaSums)
    myUniques = unlist(myUniques)
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
    # print(sd(log(myUniques+1)))
    initPars = c(myPars$mu, myPars$sigma)  #max(c(log(sum(gammaSums*myUniques)/sum(gammaSums))))
    # before = myQPoiLog(par=initPars,myGammas=gammaSums,
    # sizeFactor=sizeFactors, x=myUniques) save(initPars, myQPoiLog,
    # gammaSums, sizeFactors, myUniques, file=paste('optim.data.',
    # myPars$d, '.rda',sep=''))
    out = optim(par = initPars, fn = myQPoiLog, myGammas = gammaSums, sizeFactor = sizeFactors, 
        x = myUniques, lower = c(log(1e-06), 1e-06), upper = c(log(1e+06), 
            10000), method = "L-BFGS-B")$par
    # out = optim(par=c(myPars$mu, myPars$sigma), fn=myQPoiLog,
    # myGammas=gammaSums, sizeFactor=sizeFactors, x=myUniques,
    # method='BFGS')$par#, control=list(reltol=1e-8, maxit=10000))$par
    # cat('mean: ', log(sum(gammaSums*myUniques)/sum(gammaSums)), '
    # par[2]=', out[1], '\n', sep='') print(out) out[2] = abs(out[2])
    # after = myQPoiLog(par=out,myGammas=gammaSums, sizeFactor=sizeFactors,
    # x=myUniques) print(before-after) if(c(before-after) < 0.1) {
    # print(sum(gammaSums)) out = c(myPars$mu, myPars$sigma) }
    
    out
}




#' Optimizes parameters of PoissonLogNormal emission functions during EM-learning.
#' 
#' @keywords internal
#' @noRd
optimizePoiLog = function(myPars) {
    
    mySplit = myPars$uniqueCountSplit
    out = NULL
    myFile = file.path(tempdir(), paste("STAN.temp.params.", Sys.getpid(), 
        ".rda", sep = ""))
    if (myPars$d != 1) {
        load(myFile)
    } else {
        # print(myPars$d) print(myPars$currstate)
        
        # if(file.exists(myFile)) { system(paste('rm ', myFile, sep='')) }
        load(myFile)
        # print(myparams) print(myFile) print(sapply(mySplit, length))
        out = mclapply(1:length(mySplit[[1]]), function(currD) {
            gammaSums = list()
            myUniques = list()
            sizeFactors = list()
            for (i in 1:length(mySplit)) {
                gammaSums[[i]] = sapply(mySplit[[i]][[currD]], function(x) sum(myPars$gamma[x]))
                myUniques[[i]] = as.integer(names(mySplit[[i]][[currD]]))
                gammaSums[[i]] = unlist(gammaSums[[i]])
                if (i > nrow(myPars$sizeFactor)) {
                  # bdHMM
                  sizeFactors[[i]] = rep(myPars$sizeFactor[i/2, currD], 
                    length(gammaSums[[i]]))
                } else {
                  sizeFactors[[i]] = rep(myPars$sizeFactor[i, currD], length(gammaSums[[i]]))
                }
                # sizeFactors[[i]] = rep(myPars$sizeFactor[i,currD],
                # length(gammaSums[[i]]))
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
            initPars = c(myparams$mu[myPars$currstate, currD], myparams$sigma[myPars$currstate, 
                currD])
            # print(initPars) before = myQPoiLog(par=initPars,myGammas=gammaSums,
            # sizeFactor=sizeFactors, x=myUniques) print(initPars)
            optim(par = initPars, fn = myQPoiLog, myGammas = gammaSums, 
                sizeFactor = sizeFactors, x = myUniques, lower = c(log(1e-06), 
                  1e-06), upper = c(log(1e+06), 10000), method = "L-BFGS-B")$par
            
        }, mc.cores = min(c(myPars$ncores, length(mySplit[[1]]))))
        # print(unlist(sapply(out, function(x) x[1])))
        myparams$mu[myPars$currstate, ] = unlist(sapply(out, function(x) x[1]))
        myparams$sigma[myPars$currstate, ] = unlist(sapply(out, function(x) x[2]))
        # print(myparams)
        save(myparams, file = myFile)
    }
    # print(myparams$mu)
    out = c(myparams$mu[myPars$currstate, myPars$d], myparams$sigma[myPars$currstate, 
        myPars$d])
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



#' Optimizes parameters of PoissonLogNormal emission functions during EM-learning or initialization. Called by optimizePoiLog and optimizePoiLogInit.
#' 
#' @keywords internal
#' @noRd
myQPoiLog <- function(x, par, myGammas = rep(1, length(x)), sizeFactor = 1) {
    # cat(par[1], '\n', sep=' ') myMeans = log((exp(par[1])/sizeFactor))
    # par[2] = abs(par[2])
    myMeans = par[1] - log(sizeFactor)
    # print(x)
    mypl = NULL
    # if(all(sizeFactor == 1)) { mypln =log(dpoilog(x, mu=myMeans,
    # sig=par[2])+1e-300) } else {
    mypln = sapply(1:length(x), function(i) log(dpoilog(x[i], mu = myMeans[i], 
        sig = par[2]) + 1e-300))
    # mypln = ifelse(mypln==-Inf,-1e12,mypln) }
    
    # print(range(mypln))
    
    myval = -sum(myGammas * mypln)
    myval
}




#'  
#' Calculate density of the Poisson-Log-Normal distribution.
#' 
#' @title Calculate density of the Poisson-Log-Normal distribution.
#' 
#' @param x A vector c(n, mu, sigma), where n is the number of observed counts, mu the mean of the Log-Normal distribution and sigma its variance.
#' 
#' @return Density of the Poisson-Log-Normal distribution.
#' @usage call_dpoilog(x)
#' 
#' @examples
#' call_dpoilog(c(5, 2, 1))
#' 
#' @export call_dpoilog
call_dpoilog = function(x) {
    dpoilog(x[1], x[2], x[3])
} 
