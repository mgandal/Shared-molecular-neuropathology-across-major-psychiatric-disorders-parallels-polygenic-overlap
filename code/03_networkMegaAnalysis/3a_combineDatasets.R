##3a_combineDatasets
## This script combines microarray studies together and performs ComBat to normalize them
## QC plots are made


rm(list=ls()); options(stringsAsFactors=F)

#source("http://bioconductor.org/biocLite.R")
library(ggplot2);library(sva); library(WGCNA)
rootdir = "~/Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap//"
setwd(rootdir)

par(mfrow=c(2,2))

files = dir("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned//", pattern="_CR_"); files = files[!grepl("IBD",files)]
multiExpr = vector(mode="list",length = length(files))
for( i in 1:length(files)) {
  load(paste("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/",files[[i]],sep=""))
  multiExpr[[i]]$datExpr= datExpr
  multiExpr[[i]]$datMeta= datMeta
  rm(datExpr); rm(datMeta); 
}

genes = rownames(multiExpr[[1]]$datExpr)
for(i in 2:length(multiExpr)) genes = intersect(genes, rownames(multiExpr[[i]]$datExpr))
idx = match(genes, rownames(datProbes))
datProbes = datProbes[idx, c("ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")]

all_datExpr = data.frame(row.names = genes)
all_datMeta = data.frame(matrix(NA, nrow=0, ncol=10));

for(i in 1:length(multiExpr)) {
  all_datExpr = cbind(all_datExpr, multiExpr[[i]]$datExpr[match(genes, rownames(multiExpr[[i]]$datExpr)),])
  datMetaNA = as.data.frame(matrix(NA, ncol=10, nrow=nrow(multiExpr[[i]]$datMeta)))
  idx = pmatch(c("Study", "Subject", "Group", "Region", "Age", "Sex", "PMI", "pH", "RIN", "RNAdeg"),colnames(multiExpr[[i]]$datMeta))
  datMetaNA[,which(!is.na(idx))] = multiExpr[[i]]$datMeta[,na.omit(idx)]
  all_datMeta=rbind(all_datMeta, datMetaNA)
}

colnames(all_datMeta) = c("Study", "Subject", "Group", "Region", "Age", "Sex", "PMI", "pH", "RIN", "RNAdeg")
all_datMeta$RNAdeg[all_datMeta$Study=="SCZ.BD.Chen"]= NA #3' bias is not compatible from this study to the others (array platforM)


##QC Pre-Combat
all_datMeta$Study = as.factor(all_datMeta$Study)
sex_col = rep("blue", times = nrow(all_datMeta))
sex_col[all_datMeta$Sex=="F"] = "pink"
age_col = numbers2colors(all_datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(all_datMeta$Age, na.rm=T),max(all_datMeta$Age, na.rm=T)))
ph_col = numbers2colors(all_datMeta$pH, blueWhiteRed(100), signed=F, centered=T, lim=c(min(all_datMeta$pH, na.rm=T),max(all_datMeta$pH, na.rm=T)))
pmi_col = numbers2colors(all_datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(all_datMeta$PMI, na.rm=T),max(all_datMeta$PMI, na.rm=T)))
rin_col = numbers2colors(all_datMeta$RIN, blueWhiteRed(100), signed=F, centered=T, lim=c(min(all_datMeta$RIN, na.rm=T),max(all_datMeta$RIN, na.rm=T)))
rna_col = numbers2colors(all_datMeta$RNAdeg, blueWhiteRed(100), signed=F, centered=T, lim=c(min(all_datMeta$RNAdeg, na.rm=T),max(all_datMeta$RNAdeg, na.rm=T)))

plot(density(all_datExpr[,1]), xlim=c(-5,20), ylim=c(0, 0.5), col = as.numeric(all_datMeta$Study[1]), xlab="Intensity (log2)", ylab="Density", main="Mega-Analysis: Pre-Combat")
for(i in 2:dim(all_datExpr)[[2]])
  lines(density(all_datExpr[,i]), xlim=c(0,20), col = as.numeric(all_datMeta$Study[i]))  
legend("topleft", (levels(all_datMeta$Study)), col=c(1:8), pch=16, cex=0.5)

mds = cmdscale(dist(t(all_datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(all_datMeta$Study)), pch=20, main="Multidimensional Scaling Plot\nPre-ComBat", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
legend("topleft", (levels(all_datMeta$Study)), col=c(1:8), pch=16, cex=0.5)

tree = hclust(dist(t(all_datExpr)), method = "average")
par(mfrow=c(1,1))
plotDendroAndColors(tree, cbind(as.numeric(all_datMeta$Group), as.numeric(all_datMeta$Study), sex_col, age_col, ph_col, pmi_col, rin_col, rna_col), 
                    groupLabels = c("Group", "Study", "Sex", "Age", "pH", "PMI", "RIN", "RNA"), cex.colorLabels=0.8, cex.dendroLabels=0.15,
                    main="Dendrogram\nPre-Combat")



#Normalize by Study
mod = model.matrix(~Group+Age+Sex, data=all_datMeta)
batch = as.factor(all_datMeta$Study)
datExpr = ComBat(all_datExpr, batch=batch, mod=mod, prior.plots = F)


##QC - PostCombat
par(mfrow=c(2,2))
plot(density(datExpr[,1]),xlim=c(0,17), ylim=c(0, 0.4), col = as.numeric(all_datMeta$Study[1]), xlab="", ylab="", main="")
for(i in 2:dim(datExpr)[[2]])
  lines(density(datExpr[,i]), xlim=c(0,16), ylim=c(0,0.3), col = as.numeric(all_datMeta$Study[i]))  
legend("topright", levels(all_datMeta$Study), col=c(1:8), pch=16,cex=0.7)

#MDS Plot
mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(all_datMeta$Study)), pch=20, main="MDS: Study", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
plot(mds$points, col=as.numeric(as.factor(all_datMeta$Group)), pch=16, main="MDS: Group",  xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
legend("bottomleft", levels(all_datMeta$Group), col=c(1:length(levels(all_datMeta$Group))), pch=16, cex=0.8)
plot(mds$points, col=sex_col, pch=16, main="MDS - Sex", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
plot(mds$points, col=age_col, pch=16, main="MDS - Age", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
#plot(mds$points, col=pmi_col, pch=16, main="MDS - PMI", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
#plot(mds$points, col=ph_col, pch=16, main="MDS - pH", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
#plot(mds$points, col=rin_col, pch=16, main="MDS - RIN", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
#plot(mds$points, col=rna_col, pch=16, main="MDS - RNA", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 

#Dendrogram
par(mfrow=c(1,1))
tree = hclust(dist(t(datExpr)), method = "average")
plotDendroAndColors(tree, cbind(as.numeric(all_datMeta$Group), as.numeric(all_datMeta$Study), sex_col, age_col, ph_col, pmi_col, rin_col, rna_col), 
                    groupLabels = c("Group", "Study", "Sex", "Age", "pH", "PMI", "RIN", "RNA"), cex.colorLabels=0.6, cex.dendroLabels=0.2,
                    main="Dendrogram\nPost-Combat")

dev.off()


datMeta = all_datMeta
save(file = "./working_data//NetworkAnalysis///All_datasets_combined_092017_11245x625.RData", datExpr,datMeta,multiExpr,datProbes)
