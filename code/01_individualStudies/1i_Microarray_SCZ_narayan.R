##1i) Microarray_SCZ_narayan
#--------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

home= ## insert your GitHub home directory, ie. "C://Users/me/GitHub/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)

#1) Load Data
#------------
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_SCZ_Narayan_normalized.RData")) {
  
  #Download and ProcesPhenotypic Data
  datMeta = read.csv("./raw_data/Microarray/Narayan_GSE21138/Narayan_GSE21138_datMeta.csv")
  rownames(datMeta) = datMeta$GE
  
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Group = factor(datMeta$Group, levels = c("CTL", "SCZ"))
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$PMI = as.numeric(datMeta$PMI)
  datMeta$pH = as.numeric(datMeta$pH)
  datMeta$Duration = as.numeric(datMeta$Duration)
  datMeta$Study = "Narayan_SCZ"
  datMeta$Subject = paste("VictorianBrainBank", datMeta$Group, datMeta$Age, datMeta$Sex,datMeta$PMI, sep="_")
  
  #Read in Raw Expression Data
  data.affy = ReadAffy(celfile.path="./raw_data/Microarray/Narayan_GSE21138/Narayan_GSE21138_RAW/")
  
  RNAdeg = AffyRNAdeg(data.affy)
  plotAffyRNAdeg(RNAdeg)
  datMeta$RNAdeg = RNAdeg$slope
  sd = protocolData(data.affy)$ScanDate; sd = substring(sd,1,nchar(sd)-9)
  datMeta$Batch = as.factor(sd)
  
  ## RMA Normalize
  datExpr = rma(data.affy, background =T, normalize=T, verbose=T)
  datExpr = exprs(datExpr)
  colnames(datExpr) = substr(colnames(datExpr),1,9)
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]
  
  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "affy_hg_u133_plus_2"
  getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene", "external_gene_id","chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  
  idx = match(rownames(datExpr), geneDat$affy_hg_u133_plus_2)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_SCZ_Narayan_normalized.RData", datExpr, datMeta, datProbes)
  
  
  #QC Plots -- Prenormalization
  datExpr = log2(exprs(data.affy))
  
  #Boxplot
  par(mfrow=c(2,2))    
  boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot", ylab = "Intensity");   legend("topright", legend=c("CTL", "SCZ"), col = 1:length(levels(datMeta$Group)), pch=19)
  
  # Histogram
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))
  
  #MDS Plots
  mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="Batch", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
  
  #Dendrogram
  tree = hclust(as.dist(1-bicor(datExpr)), method = "average");   
  plot(tree, xlab="", cex=0.3)
  
  #Covariate Plots
  par(mfrow=c(3,3))
  plot(datMeta$Group, ylim=c(0,40), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="");
  A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")
  
}

load("./working_data/Microarray/01_Normalized/Microarray_SCZ_Narayan_normalized.RData")

pdf("./results/figures/MicroarrayQC/SCZ_Narayan_QC.pdf",width=11,height=8.5)
par(mfrow=c(3,4))


#Balance Groups by covariates, remove singular batches
table(datMeta$Batch)
to_keep = (datMeta$Batch != "06/14/06") & (datMeta$RNAdeg < 5.5)
datMeta = datMeta[to_keep, ]; datExpr = datExpr[,to_keep]
datMeta$Batch = factor(datMeta$Batch)

##Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

model = model.matrix(~Group+Sex+Age+PMI+pH+RNAdeg+Batch, data=datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_SCZ_Narayan_normalized_balanced.RData", datExpr, datMeta, datProbes,model)


## Batch Correction
mod = model.matrix(~Group, data=datMeta)
batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat

#QC Plots post-normalization: Boxplot
boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot Post-Normalization", ylab = "Intensity");   legend("topright", legend=levels(datMeta$Group), col = 1:length(levels(datMeta$Group)), pch=19)

#Histogram
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))

#MDS
mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)
plot(mds$points, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="MDS: Batch", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))

#Covariates
plot(datMeta$Group, ylim=c(0,40), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="");
A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")

#Cluster Dendrogram
par(mfrow=c(1,1))
sampleTree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI),max(datMeta$PMI)))
ph_col = numbers2colors(datMeta$pH, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$pH),max(datMeta$pH)))
rna_qual = numbers2colors(datMeta$RNAdeg, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RNAdeg), max(datMeta$RNAdeg)))
duration = numbers2colors(datMeta$Duration, blueWhiteRed(100), signed=F, centered=T, lim=c(1,3))
plotDendroAndColors(sampleTree, colors = cbind(as.numeric(datMeta$Group), as.numeric(datMeta$Batch), duration, sex_col, age_col, pmi_col, ph_col, rna_qual), groupLabels = c("Group", "Batch", "Duration", "Sex", "Age", "PMI", "pH",  "RNAqual"), cex.colorLabels=0.6, cex.dendroLabels=0.5)
dev.off()


# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)  #
table(realGenes)
datExpr = datExpr[realGenes,]; datProbes = datProbes[realGenes,]
CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$affy_hg_u133_plus_2) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hg_u133_plus_2)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

#Regress covariates
X = model.matrix(~Group+Sex+Age+PMI+pH+RNAdeg, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,3:7]) %*% (as.matrix(beta[3:7,])))  # Technical Covariates
datExpr = datExpr - t(to_regress)

save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_SCZ_Narayan_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)

