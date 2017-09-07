##1b_Microarray_ASD_chow.R

rm(list=ls()); options(stringsAsFactors = F)

library(limma); library(WGCNA); library(biomaRt); library(sva)

rootdir = "~/Dropbox/GeschwindLab/Projects/CrossDisorder3/"
setwd(rootdir)


if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_ASD_Chow_normalized.RData")) {
  
  #Load & clean metaData
  datMeta = read.csv("./raw_data/Microarray/Chow_GSE28475/Chow_GSE28475_datMeta.csv")
  rownames(datMeta) = datMeta$GSM
  datMeta$Group = factor(datMeta$Group,levels =  c("CTL", "ASD"))
  datMeta$Sex = as.factor(datMeta$SEX)
  datMeta$Batch = as.factor(datMeta$Batch)
  datMeta$Age = as.numeric(datMeta$AGE)
  datMeta$ChipID = as.factor(datMeta$ChipID)
  datMeta$ChipPosition = as.factor(datMeta$ChipPosition)
  datMeta$ID = gsub(" ", ".", datMeta$Sample2)
  datMeta$Study="ASD.chow"
  
  #Load quantile normaliezed micorarray data
  #Downlaod from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE28475&format=file&file=GSE28475%5Fquantile%5Fnormalized%2Etxt%2Egz
  datExpr = read.delim(file="./raw_data/Microarray/Chow_GSE28475/Chow_GSE28475_quantile_normalized.txt")
  rownames(datExpr) = datExpr[,1]
  datExpr=datExpr[,-1]
  datExpr= datExpr[,seq(1,65,by=2)]
  colnames(datExpr) = gsub("DASL_Frozen_","", colnames(datExpr))
  
  idx = match(colnames(datExpr), datMeta$ID)
  datMeta = datMeta[idx,]
  
  #Initial QC Plots
  par(mar=c(5,4,2,2), mfrow=c(2,2))
  boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot", ylab = "Intensity");  legend("bottomleft", legend=levels(datMeta$Group), col = as.numeric(datMeta$Group), pch=19)
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Batch[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))
  mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="Batch", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""))
  plot(mds$points, col=as.numeric(as.factor(datMeta$ChipID)), pch=16, main="ChipID", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""))
  par(mfrow=c(1,1))
  plot(hclust(dist(t(datExpr)),method="average"))
  
  #Assess initial levels of counfounding
  par(mfrow=c(3,3))  
  plot(datMeta$Group, ylim=c(0,20), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ as.factor(datMeta$Batch), main=paste("Batch, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
  A = anova(lm(as.numeric(datMeta$Group) ~ as.factor(datMeta$ChipID))); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ as.factor(datMeta$ChipID), main=paste("Chip, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
  A = anova(lm(as.numeric(datMeta$Group) ~ as.factor(datMeta$ChipPosition))); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ as.factor(datMeta$ChipPosition), main=paste("Chip Position, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
  #A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$Group, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="");

  ##Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "illumina_humanref_8_v3"
  getinfo <- c("illumina_humanref_8_v3", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  
  idx = match(rownames(datExpr), geneDat$illumina_humanref_8_v3)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  
  colnames(datExpr) = rownames(datMeta)
  save(file="./working_data/Microarray/01_Normalized/Microarray_ASD_Chow_normalized.RData", datMeta, datExpr, datProbes)
  
}

load("./working_data/Microarray/01_Normalized/Microarray_ASD_Chow_normalized.RData")


#pdf("./results/figures/MicroarrayQC/ASD_Chow_QC.pdf", width=11,height=8.5)
par(mfrow=c(3,4))

sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
legend("bottomleft",pch=16, legend = levels(datMeta$Group), col = 1:2)
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

datMeta$RIN[is.na(datMeta$RIN)]= mean(datMeta$RIN,na.rm=T)
model = model.matrix(~Group+Age+PMI+RIN+Batch, data=datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_ASD_Chow_normalized_balanced.RData", datMeta, datExpr, datProbes,model)


#Batch Correction
table(datMeta$Batch)
mod = model.matrix(~Group, data=datMeta)
batch = factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat


#Post-Normalization QC
boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot Post-Normalization", ylab = "Intensity");   
#legend("topright", legend=levels(datMeta$Group), col = 1:2, pch=19)
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp", xlim=c(3, 16), ylim=c(0,0.2));   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:2)

mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=as.numeric(datMeta$Group), pch=16, cex=0.8)
plot(mds$points, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="MDS: Batch", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
plot(mds$points, col=as.numeric(as.factor(datMeta$ChipID)), pch=16, main="ChipID", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))

# Plot potential Colinearities
plot(datMeta$Group, ylim=c(0,20), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Group) ~ as.factor(datMeta$ChipID))); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ as.factor(datMeta$ChipID), main=paste("Chip, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
A = anova(lm(as.numeric(datMeta$Group) ~ as.factor(datMeta$ChipPosition))); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ as.factor(datMeta$ChipPosition), main=paste("Chip Position, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=1"), ylab="", xlab="")
A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$Group, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="");

# Cluster Dendrogram
par(mfrow=c(1,1))
tree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI),max(datMeta$PMI)))
rna_qual = numbers2colors(datMeta$RIN, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RIN,na.rm=T), max(datMeta$RIN,na.rm=T)))
plotDendroAndColors(tree, colors = cbind(as.numeric(datMeta$Group), as.numeric(datMeta$Batch), as.numeric(datMeta$ChipID),  sex_col, age_col, pmi_col, rna_qual), groupLabels = c("Group", "Batch", "Chip", "Sex", "Age", "PMI",  "RIN"), cex.colorLabels=0.6, cex.dendroLabels=0.5)
dev.off()


# CollapseRows Probes --> Genes
realGenes = !is.na(datProbes$ensembl_gene_id)  #
table(realGenes)
datExpr = datExpr[realGenes,]; datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$illumina_humanref_8_v3) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$illumina_humanref_8_v3)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

## Regress all Technical and Non-DX biological covariates
datMeta$RIN[is.na(datMeta$RIN)]= mean(datMeta$RIN,na.rm=T)
X = model.matrix(~Group+Age+PMI+RIN, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,3:5]) %*% (as.matrix(beta[3:5,]))) 
datExpr = datExpr - t(to_regress)

save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_ASD_chow_normalized_CR_cleaned.RData", datExpr, datMeta, datProbes)