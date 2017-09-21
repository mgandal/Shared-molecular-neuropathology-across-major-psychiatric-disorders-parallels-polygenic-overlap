##1j) Microarray_MDD_Sibille
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
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_MDD_Sibille1_normalized.RData")) {
  
  #Download and ProcesPhenotypic Data
  ## GSE54565, GSE54567, GSE54568, GSE54571, GSE54572 -- GPL570
  datMeta = read.csv("./raw_data/Microarray/Sibille/Sibille_datMeta.csv")
  rownames(datMeta) = datMeta$GSM
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Group = factor(gsub("Control", "CTL", datMeta$GROUP), levels = c("CTL", "MDD"))
  datMeta$Region = as.factor(datMeta$REGION)
  datMeta$StudyGSE = as.factor(datMeta$GSE)
  datMeta$Study = "MDD_Sibille"
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$PMI = as.numeric(datMeta$PMI)
  datMeta$pH = as.numeric(datMeta$pH)
  datMeta$Race = as.factor(datMeta$Race)
  
  #Read in Raw Expression Data
  data.affy = ReadAffy(celfile.path="./raw_data/Microarray/Sibille/Sibille/Sibille_GSE54565_GSE54567_GSE54568_GSE54571_GSE54572_RAW/")
  idx = match(substr(colnames(exprs(data.affy)),1,10), rownames(datMeta))
  datMeta = datMeta[idx,]
  
  RNAdeg = AffyRNAdeg(data.affy)
  plotAffyRNAdeg(RNAdeg)
  datMeta$RNAdeg = RNAdeg$slope
  sd = protocolData(data.affy)$ScanDate; sd = substring(sd,1,nchar(sd)-9)
  datMeta$Batch = as.factor(sd)
  
  ## RMA Normalize
  datExpr = rma(data.affy, background =T, normalize=T, verbose=T)
  datExpr = exprs(datExpr)
  colnames(datExpr) = substr(colnames(datExpr),1,10)
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]
  
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  identifier <- "affy_hg_u133_plus_2"
  getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene", "external_gene_id","chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  
  idx = match(rownames(datExpr), geneDat$affy_hg_u133_plus_2)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_MDD_Sibille1_normalized.RData", datExpr, datMeta, datProbes)
  
  #QC Plots -- Prenormalization
  datExpr = log2(exprs(data.affy))
  
  par(mfrow=c(3,3))    
  boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot", ylab = "Intensity");   legend("topright", legend=c("CTL", "SCZ"), col = 1:length(levels(datMeta$Group)), pch=19)
  
  # Histogram
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))
  
  #MDS Plots
  mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="Batch", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
  sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
  plot(mds$points, col=sex_col, pch=16, main="Sex", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
  plot(mds$points, col=as.numeric(datMeta$Region), pch=16, main="Region", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
  plot(mds$points, col=as.numeric(datMeta$Study), pch=16, main="Study", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))
  plotAffyRNAdeg(RNAdeg)
  
  #Covariate Plots
  plot(datMeta$Group, ylim=c(0,80), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="");
  A = anova(lm(as.numeric(datMeta$Group) ~ datMeta$Region)); p = A$"Pr(>F)"[1];   plot(datMeta$Region ~ datMeta$Group, main=paste("Region, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="");
  A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$Group, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="");
  A = anova(lm(as.numeric(datMeta$Group) ~ datMeta$Race));p  = A$"Pr(>F)"[1];   plot(datMeta$Group ~ datMeta$Race, main=paste("Race, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Group) ~ datMeta$Study)); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ datMeta$Study, main=paste("Study, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")
  
  #Dendrogram
  par(mfrow=c(1,1));   tree = hclust(as.dist(1-bicor(datExpr)), method = "average");   plot(tree, xlab="", cex=0.3)
  

} 

load("./working_data/Microarray/01_Normalized/Microarray_MDD_Sibille1_normalized.RData")

pdf("./results/figures/MicroarrayQC/MDD_Sibille1_QC.pdf",width=11,height=8.5)
par(mfrow=c(3,4))


## Remove singular batch
table(datMeta$Batch)
to_keep = (datMeta$Batch != "12/03/09")
datMeta = datMeta[to_keep, ]; datExpr = datExpr[,to_keep]
datMeta$Batch = factor(datMeta$Batch)

datMeta$Region = factor(datMeta$Region)

## Remove Group x covariate confound --> none
model = model.matrix(~Group+Age+Region+Race+Sex+PMI+pH+RIN+RNAdeg+Batch, data=datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_MDD_Sibille1_normalized_balanced.RData", datExpr, datMeta, datProbes,model)



## Batch Correction
mod = model.matrix(~Group, data=datMeta)
batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat




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




#QC Plots post normalization
#Boxplot
boxplot(datExpr, range=0,col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot Post-Normalization", ylab = "Intensity");   legend("topright", legend=levels(datMeta$Group), col = 1:length(levels(datMeta$Group)), pch=19)
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))

mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)
plot(mds$points, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="MDS: Batch", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))



## Plot potential Colinearities
plot(datMeta$Group, ylim=c(0,80), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="");
A = anova(lm(as.numeric(datMeta$Group) ~ datMeta$Region)); p = A$"Pr(>F)"[1];   plot(datMeta$Region ~ datMeta$Group, main=paste("Region, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="");
A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$Group, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="");
A = anova(lm(as.numeric(datMeta$Group) ~ datMeta$Race));p  = A$"Pr(>F)"[1];   plot(datMeta$Group ~ datMeta$Race, main=paste("Race, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")


#Cluster Dendrogram
par(mfrow=c(1,1))
sampleTree = hclust(as.dist(1-bicor(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI),max(datMeta$PMI)))
ph_col = numbers2colors(datMeta$pH, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$pH),max(datMeta$pH)))
rin_col = numbers2colors(datMeta$RIN, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RIN), max(datMeta$RIN)))
rna_qual = numbers2colors(datMeta$RNAdeg, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RNAdeg), max(datMeta$RNAdeg)))
plotDendroAndColors(sampleTree, colors = cbind(as.numeric(datMeta$Group), as.numeric(datMeta$Batch), as.numeric(datMeta$Region), as.numeric(datMeta$Race), sex_col, age_col, pmi_col, ph_col, rin_col, rna_qual), groupLabels = c("Group", "Batch", "Region", "Race", "Sex", "Age", "PMI", "pH", "RIN", "RNAqual"), cex.colorLabels=0.6, cex.dendroLabels=0.5)
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

datMeta$Region = factor(datMeta$Region)
X = model.matrix(~Group+Age+Region+Race+Sex+PMI+pH+RIN+RNAdeg, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,3:10]) %*% (as.matrix(beta[3:10,])))  # Technical Covariates
datExpr = datExpr - t(to_regress)

save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_MDD_Sibille1_Normalized_CR_regressed.RData", datExpr, datMeta, datProbes)


