



#'  
#' Initialization of bidirectional hidden Markov models
#' 
#' @title Initialization of bidirectional hidden Markov models
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param dStates The number of directed states.
#' @param uStates The number of undirected states. 
#' @param method Emission distribution of the model. One out of c("NegativeBinomial", "PoissonLogNormal", "NegativeMultinomial", "ZINegativeBinomial", "Poisson", "Bernoulli", "Gaussian", "IndependentGaussian")
#' @param dirFlags The flag sequence is needed when a bdHMM is fitted on undirected data (e.g.) ChIP only. It is a \code{list} of \code{character} vectors indication for each position its knwon directionality. U allows all states. F allows undirected states and states in forward direction. R allows undirected states and states in reverse direction. 
#' @param directedObs Integer vector defining the directionality (or strand-specificity) of the data tracks. Undirected (non-strand-specific) data tracks (e.g. ChIP) are indicated by '0'. Directed (strand-specific) data tracks are indicated by increasing pairs of integers. For instance c(0,0,0,1,1,2,2): The first three data tracks are undirected, followed by two pairs of directed measurements.
#' @param sizeFactors Library size factors for Emissions PoissonLogNormal or NegativeBinomial as a length(obs) x ncol(obs[[1]]) matrix.
#' @param sharedCov If TRUE, (co-)variance of (Independent)Gaussian is shared over states. Only applicable to 'Gaussian' or 'IndependentGaussian' emissions. Default: FALSE. 
#' 
#' @return A HMM object.
#' @usage initBdHMM(obs, dStates = 0, uStates = 0, method, dirFlags = NULL, directedObs = rep(0, ncol(obs[[1]])), sizeFactors = matrix(1, nrow = length(obs), ncol = ncol(obs[[1]])), sharedCov = FALSE)
#'  
#' @examples
#' 
#' data(example)
#' bdHMM_ex = initBdHMM(observations, dStates=3, method="Gaussian") 
#'
#' @export initBdHMM
initBdHMM = function(obs, dStates = 0, uStates = 0, method, dirFlags = NULL, directedObs = rep(0, ncol(obs[[1]])), 
    sizeFactors = matrix(1, nrow = length(obs), ncol = ncol(obs[[1]])), 
    sharedCov = FALSE) {
    if(is.null(dirFlags) & uStates != 0 & !sum(directedObs != 0 ) & dStates != 0)  stop("It is impossible to estimate appropiate initial parameters without any directinal information")
  if(is.null(directedObs)){
    if(length(dirFlags) != length(obs) ) stop("length of Flags and observation must match!")
  }  
  nStates = dStates + uStates
    rev.operation = 1:ncol(obs[[1]])
    for (i in 1:length(directedObs)) {
        if (directedObs[i] == 0) {
            rev.operation[i] = i
        } else {
            myPair = which(directedObs == directedObs[i])
            myPair = myPair[myPair != i]
            rev.operation[i] = myPair
        }
    }
    
    myDirCols = tapply(1:length(directedObs), INDEX = directedObs, identity)
    myDirCols = myDirCols[-which(names(myDirCols) == "0")]
    
    if (length(myDirCols) > 0) {
        for (i in 1:length(myDirCols)) {
            for (k in 1:length(obs)) {
                obs[[k]] = obs[[k]][apply(obs[[k]], 1, function(x) all(!is.na(x))), 
                  ]
                # obs[[k]][, myDirCols[[i]]] = t(apply(obs[[k]][, myDirCols[[i]]], 
                  # 1, sort, decreasing = TRUE))
            }
        }
    } ## Makes always one column of the observation matrices to have the higher number when strand-specific data is included
#     #print(paste("rev.operation ", rev.operation))
    myMat = do.call("rbind", obs)
    stateLabels = c()
    indexStates = c()
    if(uStates != 0){
        stateLabels = c(stateLabels, paste("U", 1:uStates, sep = ""))
        indexStates = c(indexStates, paste("U", 1:uStates, sep = "") )
    }
    if(dStates != 0){
        stateLabels = c( paste("F", 1:(dStates), sep = ""), paste("R", 1:(dStates), sep = ""), stateLabels)
        indexStates = c( paste("D", 1:(dStates), sep = ""), indexStates)
    } 
    #print(indexStates)
    if(!is.null(dirFlags) & dStates != 0 & uStates != 0){
        
        myFlags <- unlist(dirFlags)
        UMat <- myMat[myFlags == "U",]
        kmU <- clusterMat(UMat, uStates, method)
#         return(kmU)
        for(s in uStates:1){
            #print(paste("in loop s cluster = ", grep("U", indexStates)[s]))
            #print(paste("length cluster ", grep("U", indexStates)[s], sum(kmU$cluster == s)))
            kmU$cluster[kmU$cluster == s] = grep("U", indexStates)[s]
        }
        #print(paste("line 78, levels kmU$cluster", levels(as.factor(kmU$cluster))))
        DMat <- myMat[myFlags != "U",]
        kmD <- clusterMat(DMat, dStates, method)
        km <- clusterMat(myMat, nStates, method)
        km$cluster[myFlags == "U"] = kmU$cluster
        km$cluster[myFlags != "U"] = kmD$cluster
        km$centers[as.numeric(levels(as.factor(kmD$cluster))),] = kmD$centers
        km$centers[as.numeric(levels(as.factor(kmU$cluster))),] = kmU$centers
        
    }
    else{
        km = clusterMat(myMat, nStates, method)
    }
#     return(km)
    
    myEmissions = NULL
    if (method == "NegativeBinomial") {
        #print(paste("inside method ==  NB", indexStates))
        #print("He entrado en NegativeBinomial")
        
        myEmissions = initNB(km, obs, sizeFactor = sizeFactors, stateLabels = stateLabels, 
            directedObs = directedObs, indexStates = indexStates)
    }
    if (method == "ZINegativeBinomial") {
        myEmissions = initNB(km, obs, sizeFactor = sizeFactors, zeroInflated = TRUE, 
                             stateLabels = stateLabels, directedObs = directedObs, indexStates = indexStates)
    }
    if (method == "PoissonLogNormal") {
#         print(paste("I enter in PoissonLogNormal"))
        myEmissions = initPoiLog(km, obs, sizeFactor = sizeFactors, stateLabels = stateLabels, 
            directedObs = directedObs, indexStates = indexStates)
    }
    stateLabels = c()
    indexStates = c()
    if(uStates != 0){
        stateLabels = c(stateLabels, paste("U", 1:uStates, sep = ""))
        indexStates = c(indexStates, paste("U", 1:uStates, sep = "") )
    }
    if(dStates != 0){
        stateLabels = c( paste("F", 1:(dStates), sep = ""), paste("R", 1:(dStates), sep = ""), stateLabels)
        indexStates = c( paste("D", 1:(dStates), sep = ""), indexStates)
    } 
    
    if (method == "Gaussian") {
        # print(km)
        # print(myDirCols)
        if(length(myDirCols) > 0){
            means_Init = as.list(as.data.frame(t(km$centers[order(km$centers[,myDirCols[[1]][1]], decreasing = TRUE ),])))
            print(means_Init)
        }else{
            means_Init = as.list(as.data.frame(t(km$centers)))
        }
        
        
        covs_Init = lapply(1:length(means_Init), function(x) matrix(cov(myMat), ncol=ncol(myMat)))
      
		myEmissions = list(HMMEmission(type = "Gaussian", parameters = list(mu = c(means_Init[grep("D", indexStates)], 
            lapply(means_Init[grep("D", indexStates)], function(x) x[rev.operation]), means_Init[grep("U", indexStates)]), 
            cov = c(covs_Init[grep("D", indexStates)], 
            lapply(covs_Init[grep("D", indexStates)], function(x) matrix(x[rev.operation, rev.operation], ncol=ncol(x))), 
            covs_Init[grep("U", indexStates)]), 
            sharedCov = sharedCov), nStates = length(stateLabels)))
    }
    if (method == "IndependentGaussian") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, mean)
        means_Init = lapply(state2ind, function(x) apply(myMat[x,, drop=FALSE], 2,
            mean))
        # print(means_Init)
        means_Init <- do.call(rbind, means_Init)
        # print(means_Init)
        
        if(length(myDirCols) > 0){
            means_Init =means_Init[order(means_Init[,myDirCols[[1]][1]], decreasing = TRUE ),]
        }
        # print(means_Init)
        
        means_Init = split(means_Init, rep(1:nrow(means_Init), each = ncol(means_Init)))
        
        # print(means_Init)
        
        covs_Init = lapply(1:length(means_Init), function(x) matrix(cov(myMat), ncol=ncol(myMat)))
        
        
        Rmeans_Init <- list()
        Rcovs_Init <- list ()
        for (i in 1:dStates) {
            Rmeans_Init[[i]] = means_Init[[i]][rev.operation]
            Rcovs_Init[[i]] = covs_Init[[i]][rev.operation, 
                rev.operation]
        }
        means_Init <- c(means_Init[grep("D", indexStates)], Rmeans_Init, means_Init[grep("U", indexStates)])
        covs_Init <- c(covs_Init[grep("D", indexStates)], Rcovs_Init, covs_Init[grep("U", indexStates)])
        myEmissions = list()
        for (i in 1:ncol(obs[[1]])) {
            myEmissions[[i]] = HMMEmission(type = "Gaussian", parameters = list(mu = lapply(means_Init, 
                function(x) x[i]), cov = lapply(covs_Init, function(x) matrix(var(as.vector(obs[[1]])), 
                ncol = 1)), sharedCov = sharedCov), nStates = length(stateLabels) )
        }
    }
    if (method == "Bernoulli") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = lapply(state2ind, function(x) apply(myMat[x,, drop=FALSE], 2, 
            mean))
        
        means_Init = c(means_Init[grep("D", indexStates)], means_Init[grep("D", indexStates)], means_Init[grep("U", indexStates)])
        names(means_Init) = stateLabels
        
        for( i in grep("R", stateLabels)){
            means_Init[[i]] <- means_Init[[i]][rev.operation]
        }
        
        myEmissions = list()
        
        for (i in 1:ncol(obs[[1]])) {
            myEmissions[[i]] = HMMEmission(type = "Bernoulli", parameters = list(p = lapply(means_Init, 
                function(x) x[i])), nStates = length(stateLabels))
        }
    }
    
    if (method == "Poisson") {
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = lapply(state2ind, function(x) apply(myMat[x, 1:ncol(myMat), drop=FALSE], 2, 
                                                         mean))
        
        for(s in 1:length(means_Init)){
            for(v in 1:length(means_Init[[s]])){
                if(means_Init[[s]][v] == 0 ) means_Init[[s]][v] = 1e-6
            }
        }
        means_Init = c(means_Init[grep("D", indexStates)], means_Init[grep("D", indexStates)], means_Init[grep("U", indexStates)])
        
        names(means_Init) = stateLabels
        
        for( i in grep("R", stateLabels)){
            means_Init[[i]] <- means_Init[[i]][rev.operation]
        }
        
        myEmissions = list()
        
        for (i in 1:ncol(obs[[1]])) {
            myEmissions[[i]] = HMMEmission(type = "Poisson", parameters = list(lambda = lapply(means_Init, 
                                                                                               function(x) x[i])), nStates = length(stateLabels))
        }
    }
    
    if (method == "NegativeMultinomial") {
        
        state2ind = tapply(1:length(km$cluster), INDEX = km$cluster, identity)
        means_Init = do.call("rbind", lapply(state2ind, function(x) apply(myMat[x, 1:ncol(myMat), drop=FALSE], 2, mean)))
        
        means_Init = rbind(means_Init[grep("D", indexStates),], means_Init[grep("D", indexStates),rev.operation], means_Init[grep("U", indexStates),])
        
        multiNomInit = lapply(1:nrow(means_Init), function(x) means_Init[x, 
                                                                         2:ncol(means_Init)]/sum(means_Init[x, 2:ncol(means_Init)]))
        multiNomEmission = HMMEmission(type = "Multinomial", parameters = list(p = multiNomInit), 
                                       nStates = length(stateLabels))
        myEmissionsNB = initNB(km, obs, sizeFactor = sizeFactors, stateLabels = stateLabels, 
                             directedObs = directedObs, indexStates = indexStates)[1]
        
        # print(myEmissionsNB)
        # print("-----------------------")
        # print( multiNomEmission)
        # print("-----------------------")
        myEmissions = c(myEmissionsNB, multiNomEmission)
        # print(myEmissions)
    }
    
    
    nStates = length(stateLabels)
    
    initProb = rep(1/nStates, nStates)
    transMat = matrix(1/nStates, nrow = nStates, ncol = nStates)
    # cat("TransMat:", "\n", transMat, "\n", "My Emissions:", "\n", myEmissions)
    bdhmm = bdHMM(initProb = initProb, transMat = transMat, emission = myEmissions, 
        nStates = nStates, stateNames = stateLabels, status = "initial", transitionsOptim = "analytical", 
        directedObs = directedObs, dimNames = colnames(obs[[1]]))
    bdhmm
    
} 
