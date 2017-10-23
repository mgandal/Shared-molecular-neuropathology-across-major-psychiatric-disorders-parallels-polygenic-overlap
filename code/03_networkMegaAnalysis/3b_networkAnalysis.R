#3b_networkAnalysis.R

##CrossDisorder: step 4_WGCNA on mega analysis
##--------------------------------------------
##This script takes a mega-analysis of gene expression studies, datExpr,
##and runs WGCNA to find modules of co-expressed genes


rm(list=ls())
library(WGCNA)
rootdir = "~/Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/"
setwd(rootdir)
load("./working_data/NetworkAnalysis/All_datasets_combined_092017_11245x625.RData")


## NETWORK ANALYSIS
## ----------------
multiExpr = vector(mode="list", length=1)
multiExpr[[1]] = list(data=as.data.frame(t(datExpr)))
bsize = 5000
nSets = 1
powers = c(seq(1,9,by=1),seq(10,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()


## Compute soft-threshold for mega-analysis

if(FALSE) {   ## only need to do this once
  pdf("./results/figures/WGCNA/WGCNA-softthresh.pdf")
  par(mfrow=c(1,2))
  n = 1
  multiExpr[[n]]$softThresh = pickSoftThreshold(data= multiExpr[[n]]$data, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = bsize)
  
  sft = multiExpr[[n]]$softThresh
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2")
  abline(h=0.8, col="black")
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
  dev.off()
}


wgcna_parameters = list(powers =  9)
wgcna_parameters$minModSize = 40  
wgcna_parameters$minHeight = 0.1
wgcna_parameters$bsize = 26000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 2  ##deep split parameter contorls number of modules
wgcna_parameters$networkType = "signed"    ## using signed networks
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = TRUE

n = 1

## ROBUST WGCNA goes here
## --> we want to generate a "consensus TOM" based on resampling datExpr so as to make 
## --> module definitions robust to outliers
if(FALSE) multiExpr[[n]]$netData = blockwiseModules(datExpr=multiExpr[[n]]$data, maxBlockSize=wgcna_parameters$bsize, 
                                                   networkType=wgcna_parameters$networkType, corType = wgcna_parameters$corFnc ,  
                                                   power = wgcna_parameters$powers[n], mergeCutHeight= wgcna_parameters$minHeight, 
                                                   nThreads=4, saveTOMFileBase=paste("./processed_data/WGCNA/network_signed", dataset_params, "_exprSet", as.character(n),sep=""), 
                                                   saveTOMs=TRUE, minModuleSize= wgcna_parameters$minModSize, pamStage=wgcna_parameters$pamStage, 
                                                   reassignThreshold=1e-6, verbose = 3, deepSplit=wgcna_parameters$ds)




