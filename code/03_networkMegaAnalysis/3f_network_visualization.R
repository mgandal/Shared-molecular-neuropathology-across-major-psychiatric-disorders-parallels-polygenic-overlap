#3f_network_visualization.R

rm(list=ls())
options(stringsAsFactors = F)
#source("http://bioconductor.org/biocLite.R"); biocLite("igraph")
library(WGCNA);library(ggplot2); library(reshape); library(igraph); library(RColorBrewer); library(WGCNA)
rootdir = "~//Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap//"; setwd(rootdir)

load("./working_data/NetworkAnalysis/finalizedNetwork_092017.RData")

## Make module eigengene-MDS plot
eigmat = MEs$eigengenes
eigmat = eigmat[,paste0("ME", labels2colors(1:13))]
adj = bicor(eigmat)
mds = cmdscale(dist(t(eigmat)), eig = T);   
mds$points[,1] = -1* mds$points[,1]
g1 <- graph.adjacency(as.matrix(adj),mode="undirected",weighted=T,diag=FALSE)
layoutFR <- mds$points
c = paste0("CD", 1:13, ".",labels2colors(1:13))
  
pdf("./results/figures/Manuscript/Fig3B.pdf",height=4,width=4,useDingbats=FALSE)
edgecolors = numbers2colors(E(g1)$weight, colors = redWhiteGreen(10, gamma=3), signed=T, centered=T, lim=c(-1,1))
plot.igraph(g1, vertex.label = c,
            vertex.label.dist=0.3, 
            vertex.size=6,
            vertex.label.color="black",
            vertex.label.family = "sans",
            vertex.color = labels2colors(1:13),
            vertex.label.cex=0.6,
            layout=layoutFR,
            edge.color=edgecolors,
            edge.width=2,asp=1, main="Module Eigengene MDS")
labs = seq(1,-1,by=-.25)
str = paste(labs,sep="\n")
#text(-1.25,-1, labels = paste(labs,collapse='\n'),pos = 4,cex = .5)
p = matrix(NA,nrow=9,ncol=4)
p[,1]=-1.35; p[,2]=-1.3
p[,3]=p[,4] = -.6-.6*seq(0,1,by=.12)
for(i in 1:9) {
  lines(x=p[i,1:2],y=p[i,3:4], lwd = 2, col=numbers2colors(as.numeric(labs[i]), colors = redWhiteGreen(10, gamma=3), signed=T, centered=T, lim=c(-1,1)))
  text(x=-1.3,y=p[i,3],labs[i],cex=.4,pos=4)
}
dev.off()






## Part 2) 
## Make Module MDS plot
## --------------------
moduleColors = colors
modColors = data.frame(color= moduleColors,row.names = rownames(datExpr))

cons_kme = signedKME(t(datExpr), MEs$eigengenes)
cons_kme = cons_kme[,-which(colnames(cons_kme)=="kMEgrey")]

maxsize = 20  #plot top 10 hub genes for each module

# [1] "kMEblack"       "kMEblue"        "kMEbrown"       "kMEgreen"       "kMEgreenyellow"
# [6] "kMEmagenta"     "kMEpink"        "kMEpurple"      "kMEred"         "kMEsalmon"     
# [11] "kMEtan"         "kMEturquoise"   "kMEyellow"         

gene_idx = order(cons_kme[which(colors=="green"),1], decreasing = T)[1:maxsize]
for(i in c(2, 5,8,10,11,12,13)){
  gene_idx = c(gene_idx, order(cons_kme[,i], decreasing = T)[1:maxsize])
}

hubGenes = character()
hubGenes.kme = numeric()
for(col in c("blue","green","greenyellow","purple","salmon","tan","turquoise","yellow")) {
  #  col = labels2colors(i)
  
  modgenes = rownames(datExpr)[which((colors) == col)]
  kmes = cons_kme[modgenes, paste("kME", col,sep="")]
  top_hubs = modgenes[order(kmes, decreasing=T)[1:maxsize]]
  top_hubs.kme = kmes[order(kmes, decreasing=T)[1:maxsize]]
  hubGenes = c(hubGenes,top_hubs)
  hubGenes.kme = c(hubGenes.kme, top_hubs.kme)
}
gene_idx = match(hubGenes,rownames(datExpr))

adjMat = bicor(t(datExpr))
keepgenes = rownames(cons_kme)[gene_idx]
adjMat = adjMat[gene_idx,gene_idx]
topcors=0.65^9
adjMat[adjMat< topcors]=0
  
  geneSymbols = datProbes$external_gene_id[match(keepgenes, datProbes$ensembl_gene_id)]
  g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=FALSE)
  #
  
  
  # mds = cmdscale(dist(t(adjMat)), eig = T)
  # layoutFR = mds$points
  edgecolors = numbers2colors(E(g1)$weight, colors = redWhiteGreen(100, gamma=4), signed=T, centered=T, lim=c(-1,1))
  
  
  pdf("./results/figures/Manuscript/Fig3D.pdf",useDingbats = F, width=12,height=12)
  plot.igraph(g1, vertex.label = geneSymbols,
              vertex.label.dist=0, edge.width=0.25,
              vertex.size=1, vertex.frame.color="black",
              vertex.label.color="black",
              vertex.color = colors[gene_idx],
              vertex.label.cex=0.3,
              layout=layout.fruchterman.reingold(g1),
              edge.color="green")
  dev.off()
