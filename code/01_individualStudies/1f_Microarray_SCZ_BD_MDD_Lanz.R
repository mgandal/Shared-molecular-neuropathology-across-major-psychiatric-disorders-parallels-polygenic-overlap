##1f) Microarray_SCZ_BD_MDD_lanz
#--------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

home= "~/Github/" ## insert your GitHub home directory, ie. "C://Users/me/GitHub/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)

set.seed(100)

if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_MDD_Lanz_normalized.RData")) {
  
  pfc_samples = list.files("./raw_data/Microarray/Lanz_GSE53987/Lanz_GSE53987_RAW/Lanz_GSE53987_RAW/")
  pfc_samples = pfc_samples[grep("B46", pfc_samples)]
  data.affy = ReadAffy(celfile.path="./raw_data/Microarray/Lanz_GSE53987/Lanz_GSE53987_RAW/Lanz_GSE53987_RAW/", filenames=pfc_samples)
  datExpr=exprs(data.affy)
  
  datMeta = read.csv("./raw_data/Microarray/Lanz_GSE53987/Lanz_GSE53987_datMeta.csv")
  rownames(datMeta) = datMeta$Chip
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Race = as.factor(datMeta$Race)
  datMeta$Region[datMeta$BrainRegion=="hippocampus"] = "HC";   datMeta$Region[datMeta$BrainRegion=="Pre-frontal cortex (BA46)"] = "PFC";   datMeta$Region[datMeta$BrainRegion=="Associative striatum"] = "STR"
  datMeta$Region = as.factor(datMeta$Region)
  datMeta$Group[datMeta$Disorder=="bipolar disorder"]="BD"
  datMeta$Group[datMeta$Disorder=="control"]="CTL"
  datMeta$Group[datMeta$Disorder=="major depressive disorder"]="MDD"
  datMeta$Group[datMeta$Disorder=="schizophrenia"]="SCZ"
  datMeta$Group = factor(datMeta$Group, levels = c("CTL", "BD", "SCZ", "MDD"))
  datMeta$Study = "SCZ.BD.MDD.Lanz"
  datMeta$Subject = paste(datMeta$BrainBank, datMeta$Group, datMeta$Sex, datMeta$Age, datMeta$PMI, datMeta$Race, sep=".")
  
  samples = substr(colnames(datExpr),1,10)
  datMeta = datMeta[samples,]
  sd = protocolData(data.affy)$ScanDate; sd = substring(sd,1,nchar(sd)-9)
  datMeta$Batch = sd
  
  RNAdeg = AffyRNAdeg(data.affy)
  plotAffyRNAdeg(RNAdeg)
  datMeta$RNAdeg = RNAdeg$slope
  
  # Normalize
  datExpr = rma(data.affy, background =T, normalize=T, )
  datExpr = exprs(datExpr)
  colnames(datExpr) = substr(colnames(datExpr),1,10)
  
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "affy_hg_u133_plus_2"
  getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  
  idx = match(rownames(datExpr), geneDat$affy_hg_u133_plus_2)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_MDD_Lanz_normalized.RData", datExpr, datMeta, datProbes)
  
  
  ##QC Plots -- Prenormalization
  ##----------------------------
  datExpr = log2(exprs(data.affy))
  par(mfrow=c(2,2))
  boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot", ylab = "Intensity");   legend("topright", legend=levels(datMeta$Group), col = 1:length(levels(datMeta$Group)), pch=19)

  # Histogram
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))
  
  #MDS Plots
  mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="Batch", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))  
  
  #Covariate Plots
  plot(datMeta$Group, ylim=c(0,20), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$Group, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")
  
  #Dendrogram
  par(mfrow=c(1,1));   tree = hclust(as.dist(1-bicor(datExpr)), method = "average");   
  plot(tree, xlab="", cex=0.5)
  
} 

load("./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_MDD_Lanz_normalized.RData")

pdf("./results/figures/MicroarrayQC/SCZ_BD_MDD_Lanz_QC.PDF",width=11,height=8.5)
par(mfrow=c(3,4))

#Remove singular batches
table(datMeta$Batch)

#Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*cor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(z.ku, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Euclidean distance (z score)")
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

#Balance groups
to_remove = (datMeta$Group=="SCZ") & (datMeta$RNAdeg > 5)
datExpr = datExpr[,!to_remove]
datMeta = datMeta[!to_remove,]

model = model.matrix(~Group+Sex+Age+PMI+pH+Race+RIN+RNAdeg+Batch, data=datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection//Microarray_SCZ_BD_MDD_Lanz_normalized_balanced.RData", datExpr, datMeta, datProbes,model)


#Batch Correction
library(sva)
mod = model.matrix(~datMeta$Group)
batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch, mod)
datExpr = datExpr.combat


#QC-postNormalization

#Boxplot
boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot Post-Normalization", ylab = "Intensity");   legend("topright", legend=levels(datMeta$Group), col = 1:length(levels(datMeta$Group)), pch=19)
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))

mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)
plot(mds$points, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="MDS: Batch", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""))

## Plot potential Colinearities
plot(datMeta$Group, ylim=c(0,20), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Race) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Race ~ datMeta$Group, main=paste("Race, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(as.factor(datMeta$Batch)) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(as.factor(datMeta$Batch) ~ datMeta$Group, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$RIN ~ datMeta$Group, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")


#Cluster Dendrogram
par(mfrow=c(1,1))
sampleTree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI),max(datMeta$PMI)))
ph_col = numbers2colors(datMeta$pH, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$pH),max(datMeta$pH)))
rin_col = numbers2colors(datMeta$RIN, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RIN),max(datMeta$RIN)))
rna_qual = numbers2colors(datMeta$RNAdeg, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RNAdeg), max(datMeta$RNAdeg)))

plotDendroAndColors(sampleTree, colors = cbind(as.numeric(datMeta$Group), sex_col, age_col, as.numeric(datMeta$Race), pmi_col, ph_col, rin_col, rna_qual), groupLabels = c("Group", "Sex", "Age", "Race", "PMI", "pH", "RIN", "RNAqual"), cex.colorLabels=0.6, cex.dendroLabels=0.5)
dev.off()



# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)  #
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$affy_hg_u133_plus_2) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hg_u133_plus_2)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id



X = model.matrix(~Group+Sex+Age+PMI+pH+Race+RIN+RNAdeg, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,5:11]) %*% (as.matrix(beta[5:11,])))  
datExpr = datExpr - t(to_regress)




idx=  which(datMeta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
datMeta$Group.BD = datMeta$Group.SCZ = datMeta$Group
datMeta$Group.SCZ[c(which(datMeta$Group=="BD"), idx[!subset.scz])] = NA
datMeta$Group.BD[c(which(datMeta$Group=="SCZ"), idx[subset.scz])] = NA

save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned//Microarray_SCZ_BD_MDD_Lanz_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)




  
  
  