

#'  
#' Compute size factors
#' 
#' @title Compute size factors
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param celltypes Indicates the cell type/tissue for each entry in obs. 
#' 
#' @return A celltype/tissue x data tracks matrix containing the size factors.
#' @usage getSizeFactors(obs, celltypes)
#' 
#' @examples
#' 
#' data(trainRegions)
#' celltypes = list("E123"=grep("E123", names(trainRegions)), 
#'         "E116"=grep("E116", names(trainRegions)))
#' sizeFactors = getSizeFactors(trainRegions, celltypes)
#' sizeFactors
#'
#' @export getSizeFactors
getSizeFactors = function(obs, celltypes) {
    myAvg = apply(t(sapply(celltypes, function(x) apply(do.call("rbind", 
        obs[x]), 2, sum))), 2, mean)
    sizeFactors = matrix(NA, nrow = length(obs), ncol = ncol(obs[[1]]))
    rownames(sizeFactors) = sapply(names(celltypes), function(x) rep(x, 
        length(celltypes[[x]])))
    colnames(sizeFactors) = colnames(obs[[1]])
    for (cell in 1:length(celltypes)) {
        for (currDim in 1:ncol(sizeFactors)) {
            sizeFactors[celltypes[[cell]], currDim] = myAvg[currDim]/sum(do.call("rbind", 
                obs[celltypes[[cell]]])[, currDim])
        }
    }
    # list('sizeFactors'=sizeFactors, 'mean'=myAvg)
    sizeFactors
}


#'  
#' Smooth data with running mean
#' 
#' @title Smooth data with running mean
#' 
#' @param x A vector with the data.
#' @param winHalfSize The smoothing window half size.
#' 
#' @return A vector containing the smoothed data.
#' @usage runningMean(x, winHalfSize = 2)
#' 
#' @examples
#' 
#' data(trainRegions)
#' celltypes = list("E123"=grep("E123", names(trainRegions)), 
#'         "E116"=grep("E116", names(trainRegions)))
#' sizeFactors = getSizeFactors(trainRegions, celltypes)
#' sizeFactors
#'
#' @export runningMean
runningMean = function(x, winHalfSize = 2) {
    out = x
    for (i in 1:length(x)) {
        start = (max(c(1, i - winHalfSize)))
        end = (min(c(length(x), i + winHalfSize)))
        out[i] = mean(x[start:end])
    }
    out
}


#'  
#' Binarize Sequencing data with the default ChromHMM binarization
#' 
#' @title Binarize Sequencing data with the default ChromHMM binarization
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param thresh Upper tail probability to find a value equal or higher than Y (P(Y >= y))
#' 
#' @return Binarized observation sequences as a list.
#' @usage binarizeData(obs, thresh = 1e-04)
#' 
#' @examples
#' 
#' data(trainRegions)
#' binData = binarizeData(obs = trainRegions, thresh = 1e-4)
#'
#' @export binarizeData
binarizeData = function(obs, thresh = 1e-4){
  myMat = do.call("rbind", obs)
  
  myLambdas = c()
  for (i in 1:ncol(myMat)) {
    myLambdas[i] = mean(myMat[,i])
  }
  myMax <- apply(myMat, 2, max)
  mycutoffs = sapply(1:length(myLambdas), function(x) min(which(1 - ppois(1:myMax[x], 
                                                                          myLambdas[x]) < thresh)))
  for (i in 1:length(obs)) {
    for (j in 1:ncol(obs[[i]])) {
      obs[[i]][, j] = ifelse(obs[[i]][, j] <= mycutoffs[j], 
                             0, 1)
    }
  }
  obs
}



#'  
#' Convert the viterbi path to a GRanges object
#' 
#' @title Convert the viterbi path to a GRanges object
#' 
#' @param viterbi A list containing the viterbi paths as factors. The output from getViterbi.
#' @param regions GRanges object of the regions (e.g. chromosomes) stored in the viterbi path.
#' @param binSize The bin size of the viterbi path. 
#' 
#' @return The viterbi path as GRanges object.
#' @usage viterbi2GRanges(viterbi, regions, binSize)
#'
#' @examples
#' library(GenomicRanges) 
#' data(yeastTF_databychrom_ex)
#' dStates = 6
#' dirobs = as.integer(c(rep(0,10), 1, 1))
#' bdhmm_gauss = initBdHMM(yeastTF_databychrom_ex, dStates = dStates, method ="Gaussian", directedObs=dirobs)
#' bdhmm_fitted_gauss = fitHMM(yeastTF_databychrom_ex, bdhmm_gauss)
#' viterbi_bdhmm_gauss = getViterbi(bdhmm_fitted_gauss, yeastTF_databychrom_ex)
#' yeastGRanges = GRanges(IRanges(start=1214616, end=1225008), seqnames="chrIV")
#' names(viterbi_bdhmm_gauss) = "chrIV"
#' viterbi_bdhmm_gauss_gr = viterbi2GRanges(viterbi_bdhmm_gauss, yeastGRanges, 8)
#' 
#' @export viterbi2GRanges
viterbi2GRanges = function(viterbi, regions, binSize) {
    names(viterbi) = as.character(seqnames(regions))
    vit = lapply(viterbi, as.character)
    myRle = lapply(vit, Rle)
    mySegm = list()
    for (i in 1:length(myRle)) {
        myTake = 1:length(myRle[[i]]@values)
        mySeqNames = rep(names(myRle)[i], length(myTake))
        
        mySegm[[i]] = GRanges(seqnames = mySeqNames, ranges = IRanges(start = (start(myRle[[i]])[myTake] - 
            1) * binSize + 1, end = end(myRle[[i]])[myTake] * binSize), 
            name = as.character(myRle[[i]]@values))
        end(mySegm[[i]]) = end(mySegm[[i]]) + start(regions[i])
        start(mySegm[[i]]) = start(mySegm[[i]]) + start(regions[i]) - 1
    }
    
    suppressWarnings(mySegm <- do.call(c, mySegm))
    mySegm
}


#'  
#' Compute average signal in state segmentation
#' 
#' @title Compute average signal in state segmentation
#' 
#' @param viterbi A list containing the viterbi paths as factors. The output from getViterbi.
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param fct The averaging function, default: mean.
#' 
#' @return A state x data track matrix containing the average signal.
#' @usage getAvgSignal(viterbi, obs, fct=mean)
#' 
#' @examples
#' data(yeastTF_databychrom_ex)
#' dStates = 6
#' dirobs = as.integer(c(rep(0,10), 1, 1))
#' bdhmm_gauss = initBdHMM(yeastTF_databychrom_ex, dStates = dStates, method = "Gaussian", directedObs=dirobs)
#' bdhmm_fitted_gauss = fitHMM(yeastTF_databychrom_ex, bdhmm_gauss)
#' viterbi_bdhmm_gauss = getViterbi(bdhmm_fitted_gauss, yeastTF_databychrom_ex)
#' avg_signal = getAvgSignal(viterbi_bdhmm_gauss, yeastTF_databychrom_ex)
#' 
#' @export getAvgSignal
getAvgSignal = function(viterbi, obs, fct=mean) {
	viterbi = unlist(viterbi)
    mySignals = do.call("rbind", obs)
    
    state2pos = tapply(1:nrow(mySignals), INDEX = viterbi, function(x) mySignals[x, 
        ])
    state2val = lapply(state2pos, function(x) apply(x, 2, fct, na.rm=TRUE))
    vit_means = do.call("rbind", state2val)
 
	vit_means
}



#'  
#' Convert data for plotting with Gviz
#' 
#' @title Convert data for plotting with Gviz
#' 
#' @param obs The observations. A list of one or more entries containing the observation matrix (\code{numeric}) for the samples (e.g. chromosomes).
#' @param regions GRanges object of the regions (e.g. chromosomes) stored in the viterbi path.
#' @param binSize The bin size of the viterbi path. 
#' @param gen The geome id, e.g. hg19, hg38 for human.
#' @param col The color of the data tracks.
#' @param type Type of plot (See Gviz DataTrack documentation).
#' @param chrom Chromosome in chich to create the object.
#' 
#' @return A list containing the data tracks converted to Gviz objects for plotting.
#' @usage data2Gviz(obs, regions, binSize, gen, col = "black", type = "h", chrom)
#' 
#' @export data2Gviz
data2Gviz = function (obs, regions, binSize, gen, col = "black", type = "h", chrom) 
{
  dlist = list()
  mySignals = obs
  regions = regions[seqnames(regions) == chrom]
  for (n in colnames(mySignals)) {
    dlist[[n]] = DataTrack(data = mySignals[, n], start = seq(start(regions), 
                                                              end(regions), by = binSize)[1:nrow(mySignals)], end = seq(start(regions) + binSize - 1, end(regions) + binSize, by = binSize)[1:nrow(mySignals)]
                           , chromosome = as.character(seqnames(regions)), 
                           genome = gen, name = n, type = type, col = col)
  }
  dlist
}


#'  
#' Convert state segmentation for plotting with Gviz
#' 
#' @title Convert state segmentation for plotting with Gviz
#' 
#' @param viterbi A list containing the viterbi paths as factors. The output from getViterbi.
#' @param chrom The chromosome/sequence if to convert.
#' @param gen The geome id, e.g. hg19, hg38 for human.
#' @param from Genomic start poistion.
#' @param to Genomic end poistion.
#' @param statecols Named vector with state colors.
#' @param col Background color.
#' 
#' @return A list containing the viterbi path converted to Gviz objects for plotting.
#' @usage viterbi2Gviz(viterbi, chrom, gen, from, to, statecols, col = NULL)
#' 
#' @export viterbi2Gviz
viterbi2Gviz = function(viterbi, chrom, gen, from, to, statecols, col = NULL) {
    if(is.null(col)) col = "transparent"
    myViterbiPanels = list()
    viterbi = viterbi[as.character(seqnames(viterbi)) == chrom]
    for (dir in "*") {
        currGRange = viterbi[end(viterbi) >= from & start(viterbi) <= to]
        
        myViterbiPanels[[dir]] = AnnotationTrack(range = currGRange, genome = gen, 
            chromosome = chrom, name = "", id = currGRange$name, shape = "box", 
            fill = statecols[currGRange$name], col = "transparent", stacking = "dense", background.panel = col)
    }
    myViterbiPanels
}
 
