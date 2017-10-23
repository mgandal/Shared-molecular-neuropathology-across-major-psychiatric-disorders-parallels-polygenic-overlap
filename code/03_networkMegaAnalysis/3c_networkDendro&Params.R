#3c_networkDendroAndParams.R

rm(list=ls()); options(stringsAsFactors = F)

library(WGCNA)
rootdir = "~/Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap//"
setwd(rootdir)

#Load combined datasets
load("./working_data/NetworkAnalysis/All_datasets_combined_092017_11245x625.RData")

#load TOM comptued by rWGCNA
load("./working_data/NetworkAnalysis/consensusTOM_final.rda")
geneTree = hclust(1-as.dist(consensusTOM_final), method="average")


# Iterate WGCNA parameters for robustness -- this takes a while
if(FALSE) {
   colors = vector(mode="list")
  labels = vector(mode="list")
  for (pam in c(FALSE,TRUE)) {
    for (minModSize in c(50,100, 200)) {
      for (dthresh in c(0.1, 0.2)) {
        for(ds in c(0:4)) { 
            print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
            
            tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pam, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-as.dist(consensusTOM_final)))
            merged = mergeCloseModules(exprData= t(datExpr), colors = tree$labels, cutHeight=dthresh)
            colors = cbind(colors, labels2colors(merged$colors))
            
            labels = c(labels, paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
        }
      }
    }
  }
  
  plotDendroAndColors(geneTree, colors, groupLabels=labels, addGuide= TRUE, dendroLabels=FALSE, main="Dendrogram", cex.colorLabels=0.5)
  save(file="./working_data/NetworkAnalysis//WGCNA_diffParams.rda", geneTree, colors, labels)
  
  colors2 = colors
  colors2[,seq(1:30)] = colors[,seq(1,60,by=2)]
  colors2[,seq(31:60)] = colors[,seq(2,60,by=2)]
  
  plotDendroAndColors(geneTree,colors,addGuide=T,dendroLabels=F)
  
  for (i in 1:60){
    ci = as.character(colors[,i])
    c_new = matchLabels(ci, c_ref)
    colors[,i] = c_new
  }
  
  colors = cbind(colors[,15], colors)
  labels = c("Final Modules", labels)
  
  pdf("./figures/WGCNA/WGCNA_diffParams.pdf",width=6,height=8)
  plotDendroAndColors(geneTree,colors,groupLabels = labels,addGuide=T,dendroLabels=F,cex.colorLabels=0.3)
  dev.off()
  
}



# Finalized Parameters
# --------------------
# Parameters to Use: "DS=4,MMS=100,DCOR=0.1,PAM=FALSE"
wgcna_parameters = list(powers =  9)
wgcna_parameters$minModSize = 100
wgcna_parameters$minHeight = 0.1
wgcna_parameters$bsize = 26000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 4  ##deep split parameter contorls number of modules
wgcna_parameters$networkType = "signed"    ## using signed networks
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = FALSE

tree = cutreeHybrid(dendro = geneTree, minClusterSize= wgcna_parameters$minModSize, pamStage=wgcna_parameters$pamStage, cutHeight = 0.999, 
                    deepSplit=wgcna_parameters$ds, distM=as.matrix(1-as.dist(consensusTOM_final)))
merged = mergeCloseModules(exprData= t(datExpr), colors = tree$labels, cutHeight=wgcna_parameters$minHeight)
colors = labels2colors(merged$colors)
table(colors)
length(table(colors))

plotDendroAndColors(geneTree,colors,groupLabels = "mod",cex.colorLabels = 0.5,addGuide=T,dendroLabels=F)


MEs = moduleEigengenes(expr = t(datExpr), colors, softPower = wgcna_parameters$powers)
kMEtable = signedKME(t(datExpr),MEs$eigengenes)
tableS1 = data.frame(kMEtable[,paste0("kME", labels2colors(1:13))])
colnames(tableS1) = paste0("kME.CD", 1:13, ".", labels2colors(1:13))

tableS1 = cbind(datProbes, data.frame(Module.Color=colors, Module.name = paste0("CD",merged$colors)), tableS1)

write.csv(file="./results/tables/Manuscript/TableS1 - kME table.csv", tableS1)
save(file="./working_data/NetworkAnalysis/finalizedNetwork_092017.RData", datExpr, datMeta, datProbes, multiExpr, geneTree, colors, wgcna_parameters, colors,MEs, kMEtable)
