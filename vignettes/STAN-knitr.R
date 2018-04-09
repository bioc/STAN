## ----style, eval=TRUE, echo=FALSE, results="asis"------------------------
BiocStyle::latex()

## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE,dev="png",fig.show="hide",
               fig.width=4,fig.height=4.5,
               message=FALSE, highlight=FALSE)

## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----Quick_start, results="hide", eval=FALSE-----------------------------
## ## Loading library and data
## library(STAN)
## data(trainRegions)
## data(pilot.hg19)
## 
## ## Model initialization
## hmm_nb = initHMM(trainRegions[1:3], nStates=10, "NegativeBinomial")
## 
## ## Model fitting
## hmm_fitted_nb = fitHMM(trainRegions[1:3], hmm_nb, maxIters=10)
## 
## ## Calculate state path
## viterbi_nb = getViterbi(hmm_fitted_nb, trainRegions[1:3])
## 
## ## Convert state path to GRanges object
## viterbi_nb_gm12878 = viterbi2GRanges(viterbi_nb, pilot.hg19, 200)

## ----Loading_library, results="hide"-------------------------------------
library(STAN)

## ----Loading_example_data------------------------------------------------
data(trainRegions)
names(trainRegions)
str(trainRegions[c(1,4)])

## ----Loading_example_data_regions----------------------------------------
data(pilot.hg19)
pilot.hg19

## ----Calculate size factors----------------------------------------------
celltypes = list("E123"=grep("E123", names(trainRegions)), 
        "E116"=grep("E116", names(trainRegions)))
sizeFactors = getSizeFactors(trainRegions, celltypes)
sizeFactors

## ----STAN-PoiLog---------------------------------------------------------
nStates = 10
hmm_poilog = initHMM(trainRegions, nStates, "PoissonLogNormal", sizeFactors)
hmm_fitted_poilog = fitHMM(trainRegions, hmm_poilog, sizeFactors=sizeFactors, maxIters=10)
viterbi_poilog = getViterbi(hmm_fitted_poilog, trainRegions, sizeFactors)
str(viterbi_poilog)

## ----STAN-PoiLog viterbi-------------------------------------------------
viterbi_poilog_gm12878 = viterbi2GRanges(viterbi_poilog[1:3], regions=pilot.hg19, binSize=200)
viterbi_poilog_gm12878

## ----STAN-NB, results="hide"---------------------------------------------
hmm_nb = initHMM(trainRegions, nStates, "NegativeBinomial", sizeFactors)
hmm_fitted_nb = fitHMM(trainRegions, hmm_nb, sizeFactors=sizeFactors, maxIters=10)
viterbi_nb = getViterbi(hmm_fitted_nb, trainRegions, sizeFactors=sizeFactors)
viterbi_nb_gm12878 = viterbi2GRanges(viterbi_nb[1:3], pilot.hg19, 200)

## ----STAN coverage, results="hide"---------------------------------------
avg_cov_nb = getAvgSignal(viterbi_nb, trainRegions)
avg_cov_poilog = getAvgSignal(viterbi_poilog, trainRegions)

## ----STAN_coverage_plotting_pdf, echo=FALSE, eval=TRUE, results="hide"----

library(gplots)
heat = c("dark blue", "dodgerblue4", "darkred", "red", "orange", "gold", "yellow")
colfct = colorRampPalette(heat)
colpal_statemeans = colfct(200)
ord_nb = order(apply(avg_cov_nb,1,max), decreasing=TRUE)
statecols_nb = rainbow(nStates)
names(statecols_nb) = ord_nb

pdf("nb_avg_cov.pdf")
heatmap.2(log(avg_cov_nb+1)[as.character(ord_nb),], margins=c(8,7),srtCol=45, RowSideColors=statecols_nb[as.character(ord_nb)], dendrogram="none", Rowv=FALSE, Colv=FALSE, col=colpal_statemeans, trace="none", cellnote=round(avg_cov_nb,1)[as.character(ord_nb),], notecol="black")
dev.off()

ord_poilog = order(apply(avg_cov_poilog,1,max), decreasing=TRUE)
statecols_poilog = rainbow(nStates)
names(statecols_poilog) = ord_poilog
pdf("poilog_avg_cov.pdf")
heatmap.2(log(avg_cov_poilog+1)[ord_poilog,], RowSideColors=statecols_poilog[as.character(ord_poilog)], margins=c(8,7),srtCol=45, dendrogram="none", Rowv=FALSE, Colv=FALSE, col=colpal_statemeans, trace="none", cellnote=round(avg_cov_poilog,1)[ord_poilog,], notecol="black")
dev.off()


## ----STAN_coverage_plotting, results="hide"------------------------------
## specify color palette
library(gplots)
heat = c("dark blue", "dodgerblue4", "darkred", "red", "orange", "gold", "yellow")
colfct = colorRampPalette(heat)
colpal_statemeans = colfct(200)

## define state order and colors
ord_nb = order(apply(avg_cov_nb,1,max), decreasing=TRUE)
statecols_nb = rainbow(nStates)
names(statecols_nb) = ord_nb
heatmap.2(log(avg_cov_nb+1)[as.character(ord_nb),], margins=c(8,7), srtCol=45, 
        RowSideColors=statecols_nb[as.character(ord_nb)], dendrogram="none", 
        Rowv=FALSE, Colv=FALSE, col=colpal_statemeans, trace="none", 
        cellnote=round(avg_cov_nb,1)[as.character(ord_nb),], notecol="black")


## define state order and colors
ord_poilog = order(apply(avg_cov_poilog,1,max), decreasing=TRUE)
statecols_poilog = rainbow(nStates)
names(statecols_poilog) = ord_poilog
heatmap.2(log(avg_cov_poilog+1)[as.character(ord_poilog),], margins=c(8,7), srtCol=45, 
        RowSideColors=statecols_poilog[as.character(ord_poilog)], dendrogram="none", 
        Rowv=FALSE, Colv=FALSE, col=colpal_statemeans, trace="none", 
        cellnote=round(avg_cov_poilog,1)[as.character(ord_poilog),], notecol="black")

## ----STAN_convert_to_Gviz, results="hide"--------------------------------
library(Gviz)
from = start(pilot.hg19)[3]
to = from+300000
gvizViterbi_nb = viterbi2Gviz(viterbi_nb_gm12878, "chr11", "hg19", from, to, statecols_nb)
gvizViterbi_poilog = viterbi2Gviz(viterbi_poilog_gm12878, "chr11", "hg19", from, to, 
        statecols_poilog)
gvizData = data2Gviz(trainRegions[[3]], pilot.hg19[3], binSize = 200, gen = "hg19", col="black", chrom = "chr11")

## ----STAN_plot_with_Gviz, eval=FALSE, results="hide"---------------------
## gaxis = GenomeAxisTrack()
## data(ucscGenes)
## mySize = c(1,rep(1.2,9), 0.5,0.5,3)
## plotTracks(c(list(gaxis), gvizData,gvizViterbi_nb,gvizViterbi_poilog,ucscGenes["chr11"]),
##         from=from, to=to, showFeatureId=FALSE, featureAnnotation="id", fontcolor.feature="black",
##         cex.feature=0.7, background.title="darkgrey", lwd=2, sizes=mySize)

## ----STAN_plot_with_Gviz_pdf, echo=FALSE, eval=TRUE, results="hide"------

gaxis = GenomeAxisTrack()
data(ucscGenes)
mySize = c(1,rep(1.2,9), 0.5,0.5,3)
pdf("gviz_example.pdf", width=7*1.5)
plotTracks(c(list(gaxis), gvizData,gvizViterbi_nb,gvizViterbi_poilog,ucscGenes["chr11"]),
        from=from, to=to, showFeatureId=FALSE, featureAnnotation="id", fontcolor.feature="black",
        cex.feature=0.7, background.title="darkgrey", lwd=2, sizes=mySize)#, stacking="dense")#, ylim=c(0,100))
dev.off()

## ----STAN_poisson_emissions, results="hide"------------------------------
hmm_pois = initHMM(trainRegions, nStates, "Poisson")
hmm_fitted_pois = fitHMM(trainRegions, hmm_pois, maxIters=10)
viterbi_pois = getViterbi(hmm_fitted_pois, trainRegions)

## ----STAN_nmn_emissions, results="hide"----------------------------------
simData_nmn = lapply(trainRegions, function(x) cbind(apply(x,1,sum), x))
hmm_nmn = initHMM(simData_nmn, nStates, "NegativeMultinomial")
hmm_fitted_nmn = fitHMM(simData_nmn, hmm_nmn, maxIters=10)
viterbi_nmn = getViterbi(hmm_fitted_nmn, simData_nmn)

## ----STAN_gaussian_emissions, results="hide"-----------------------------
trainRegions_smooth = lapply(trainRegions, function(x) 
            apply(log(x+sqrt(x^2+1)), 2, runningMean, 2))
hmm_gauss = initHMM(trainRegions_smooth, nStates, "IndependentGaussian", sharedCov=TRUE)
hmm_fitted_gauss = fitHMM(trainRegions_smooth, hmm_gauss, maxIters=10)
viterbi_gauss = getViterbi(hmm_fitted_gauss, trainRegions_smooth)

## ----STAN_bernoulli_emissions, results="hide"----------------------------
trainRegions_binary = binarizeData(trainRegions)
hmm_ber = initHMM(trainRegions_binary, nStates, "Bernoulli")
hmm_fitted_ber = fitHMM(trainRegions_binary, hmm_ber, maxIters=10)
viterbi_ber = getViterbi(hmm_fitted_ber, trainRegions_binary)

## ----STAN_other_emissions_avg_cov, results="hide"------------------------
avg_cov_gauss = getAvgSignal(viterbi_gauss, trainRegions)
avg_cov_nmn = getAvgSignal(viterbi_nmn, trainRegions)
avg_cov_ber = getAvgSignal(viterbi_ber, trainRegions)
avg_cov_pois = getAvgSignal(viterbi_pois, trainRegions)

## ----STAN_other_emissions_avg_cov_plot, eval=FALSE, results="hide"-------
## heatmap.2(log(avg_cov_gauss+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE,
##         Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1,
##         cellnote=round(avg_cov_gauss,1), notecol="black")
## heatmap.2(log(avg_cov_nmn+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE,
##         Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1,
##         cellnote=round(avg_cov_nmn,1), notecol="black")
## heatmap.2(log(avg_cov_ber+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE,
##         Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1,
##         cellnote=round(avg_cov_ber,1), notecol="black")
## heatmap.2(log(avg_cov_pois+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE,
##         Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1,
##         cellnote=round(avg_cov_pois,1), notecol="black")

## ----STAN_other_emissions_avg_cov_plot_pdf, echo=FALSE, eval=TRUE, results="hide"----
pdf("avg_cov_gauss.pdf")
heatmap.2(log(avg_cov_gauss+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE, Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1, cellnote=round(avg_cov_gauss,1), notecol="black")
dev.off()
pdf("avg_cov_nmn.pdf")
heatmap.2(log(avg_cov_nmn+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE, Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1, cellnote=round(avg_cov_nmn,1), notecol="black")
dev.off()
pdf("avg_cov_ber.pdf")
heatmap.2(log(avg_cov_ber+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE, Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1, cellnote=round(avg_cov_ber,1), notecol="black")
dev.off()
pdf("avg_cov_pois.pdf")
heatmap.2(log(avg_cov_pois+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE, Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1, cellnote=round(avg_cov_pois,1), notecol="black")
dev.off()

## ----bdHMM_yeast_fit, results="hide"-------------------------------------
data(yeastTF_databychrom_ex)
dStates = 6
dirobs = as.integer(c(rep(0,10), 1, 1))
bdhmm_gauss = initBdHMM(yeastTF_databychrom_ex, dStates = dStates, method = "Gaussian", directedObs=dirobs)
bdhmm_fitted_gauss = fitHMM(yeastTF_databychrom_ex, bdhmm_gauss)
viterbi_bdhmm_gauss = getViterbi(bdhmm_fitted_gauss, yeastTF_databychrom_ex)

## ----bdHMM_yeast_params_plot, eval=FALSE, results="hide"-----------------
## statecols_yeast = rep(rainbow(nStates), 2)
## names(statecols_yeast) = StateNames(bdhmm_fitted_gauss)
## means_fitted = EmissionParams(bdhmm_fitted_gauss)$mu
## heatmap.2(means_fitted, col=colpal_statemeans,
##         RowSideColors=statecols_yeast[rownames(means_fitted)],
##         trace="none", cexCol=0.9, cexRow=0.9,
##         cellnote=round(means_fitted,1), notecol="black", dendrogram="row",
##         Rowv=TRUE, Colv=FALSE, notecex=0.9)

## ----bdHMM_yeast_params_plot_pdf, echo=FALSE, eval=TRUE, results="hide"----
pdf("yeast_means_gauss.pdf")
statecols_yeast = rep(rainbow(nStates), 2)
names(statecols_yeast) = StateNames(bdhmm_fitted_gauss)
means_fitted = EmissionParams(bdhmm_fitted_gauss)$mu
heatmap.2(means_fitted, col=colpal_statemeans, RowSideColors=statecols_yeast[rownames(means_fitted)], trace="none", cexCol=0.9, cexRow=0.9,
        cellnote=round(means_fitted,1), notecol="black", dendrogram="row",
        Rowv=TRUE, Colv=FALSE, notecex=0.9)
dev.off()

## ----bdHMM_convert_GRanges-----------------------------------------------
yeastGRanges = GRanges(IRanges(start=1214616, end=1225008), seqnames="chrIV")
names(viterbi_bdhmm_gauss) = "chrIV"
viterbi_bdhmm_gauss_gr = viterbi2GRanges(viterbi_bdhmm_gauss, yeastGRanges, 8)
viterbi_bdhmm_gauss_gr

## ----bdHMM_gviz_plot, eval=FALSE, results="hide"-------------------------
## chr = "chrIV"
## gen = "sacCer3"
## gtrack <- GenomeAxisTrack()
## 
## from=1217060
## to=1225000
## forward_segments = grep("F", viterbi_bdhmm_gauss_gr$name)
## reverse_segments = grep("R", viterbi_bdhmm_gauss_gr$name)
## gvizViterbi_yeast = viterbi2Gviz(viterbi_bdhmm_gauss_gr[forward_segments],
##         "chrIV", "sacCer3", from, to, statecols_yeast)
## gvizViterbi_yeast2 = viterbi2Gviz(viterbi_bdhmm_gauss_gr[reverse_segments],
##         "chrIV", "sacCer3", from, to, statecols_yeast)
## 
## gvizData_yeast = data2Gviz(obs = yeastTF_databychrom_ex[[1]], regions = yeastGRanges, binSize = 8, gen = "sacCer3", col="black", chrom = chr)
## gaxis = GenomeAxisTrack()
## data(yeastTF_SGDGenes)
## mySize = c(1,rep(1,12), 0.5,0.5,3)
## 
## plotTracks(c(list(gaxis), gvizData_yeast,gvizViterbi_yeast,gvizViterbi_yeast2,
##         list(yeastTF_SGDGenes)), cex.feature=0.7, background.title="darkgrey", lwd=2,
##         sizes=mySize, from=from, to=to, showFeatureId=FALSE, featureAnnotation="id",
##         fontcolor.feature="black", cex.feature=0.7, background.title="darkgrey",
##         showId=TRUE)

## ----bdHMM_gviz_plot_pdf, echo=FALSE, eval=TRUE, results="hide"----------
yeastGRanges = GRanges(IRanges(start=1214616, end=1225008), seqnames="chrIV")
names(viterbi_bdhmm_gauss) = "chrIV"
viterbi_bdhmm_gauss_gr = viterbi2GRanges(viterbi_bdhmm_gauss, yeastGRanges, 8)
chr = "chrIV"
gen = "sacCer3"
gtrack <- GenomeAxisTrack()

from=1217060
to=1225000
forward_segments = grep("F", viterbi_bdhmm_gauss_gr$name)
reverse_segments = grep("R", viterbi_bdhmm_gauss_gr$name)
gvizViterbi_yeast = viterbi2Gviz(viterbi_bdhmm_gauss_gr[forward_segments], "chrIV", "sacCer3", from, to, statecols_yeast)
gvizViterbi_yeast2 = viterbi2Gviz(viterbi_bdhmm_gauss_gr[reverse_segments], "chrIV", "sacCer3", from, to, statecols_yeast)

gvizData_yeast = data2Gviz(obs = yeastTF_databychrom_ex[[1]], regions = yeastGRanges, binSize = 8, gen = "sacCer3", col="black", chrom = chr)
gaxis = GenomeAxisTrack()
data(yeastTF_SGDGenes)
mySize = c(1,rep(1,12), 0.5,0.5,3)
pdf("gviz_example_yeast.pdf", width=7*1.5)
plotTracks(c(list(gaxis), gvizData_yeast,gvizViterbi_yeast,gvizViterbi_yeast2,list(yeastTF_SGDGenes)),
        cex.feature=0.7, background.title="darkgrey", lwd=2, sizes=mySize, from=from, to=to,
        showFeatureId=FALSE, featureAnnotation="id",
        fontcolor.feature="black", cex.feature=0.7, background.title="darkgrey", showId=TRUE)#, stacking="dense")#, ylim=c(0,100))
plot(1:10, (1:10)^3, type = "l")
dev.off()

## ----sessInfo, results="asis", echo=FALSE--------------------------------
toLatex(sessionInfo())

