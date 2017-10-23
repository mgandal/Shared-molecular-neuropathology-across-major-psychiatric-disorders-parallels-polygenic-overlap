##1h) Microarray_AAD_Mayfield
#---------------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva); library(lumi)

home= "~/Github/" ## insert your GitHub home directory, ie. "C://Users/me/GitHub/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)


#1) Load Data
#------------
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_AAD_mayfield_normalized.RData")) {
  
  ## Format MetaData
  ## ---------------
  datMeta = read.csv("./raw_data/Microarray/Mayfield_GSE29555/Mayfield_GSE29555_datMeta.csv")
  rownames(datMeta)=datMeta$ID
  datMeta$Group = as.factor(datMeta$Group)
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$Study = "AAD.mayfield"
  
  ## Format ExpressionData
  ## ---------------------
  data.lumi = lumiR("./raw_data/Microarray/Mayfield_GSE29555/Mayfield_GSE29555_non-normalized_region1.txt")
  datExpr = lumiN(data.lumi, method="quantile")
  datExpr = log2(exprs(datExpr))
  
  idx = match(colnames(datExpr),rownames(datMeta));
  datMeta = datMeta[idx,]
  
  ## Annotate Probes
  library(biomaRt)
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  identifier <- "illumina_humanht_12_v3"
  getinfo <- c("illumina_humanht_12_v3", "ensembl_gene_id","external_gene_id", "entrezgene","chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr),geneDat[,1])
  datProbes = geneDat[idx,]
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_AAD_mayfield_normalized.RData",datExpr,datProbes,datMeta)
  
  ## QC-PreNormalization
  par(mfrow=c(2,2))
  datExpr = log2(exprs(data.lumi))
  
  boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot", ylab = "Intensity");   legend("topright", legend=levels(datMeta$Group), col = 1:length(levels(datMeta$Group)), pch=19, cex=0.7)
  i = 1; plot(density((datExpr[,i]), na.rm=T), xlim = c(5,15), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   
  for(i in 2:dim(datExpr)[2]) {       lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   
  legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))
  
  mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)
  
  tree = hclust(as.dist(1-bicor(datExpr)), method = "average");   plot(tree, xlab="", cex=0.5)
  
  #Covariate Plots
  plot(datMeta$Group, ylim=c(0,20), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$RIN) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$Group, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="")
  
} 

load("./working_data/Microarray/01_Normalized/Microarray_AAD_mayfield_normalized.RData")
pdf("./results/figures/MicroarrayQC/AAD_Mayfield_QC.pdf",width=11,height=8.5)
par(mfrow=c(3,4))


# Balance groups
#--> No singular batches or group x covariate confounds

# Outlier Removal
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

model =  model.matrix(~Group+Age+Sex+PMI+pH+RIN, data = datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_AAD_mayfield_normalized_balanced.RData",datExpr,datProbes,datMeta,model)


#Boxplot
boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot Post-Normalization", ylab = "Intensity");   legend("topright", legend=levels(datMeta$Group), col = 1:length(levels(datMeta$Group)), pch=19)
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp"); 
for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))

mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)

#Covariate Plots
plot(datMeta$Group, ylim=c(0,20), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$RIN ~ datMeta$Group, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="")

# Cluster Dendrogram for outlier removal
par(mfrow=c(1,1))
sampleTree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI),max(datMeta$PMI)))
ph_col = numbers2colors(datMeta$pH, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$pH),max(datMeta$pH)))
rna_qual = numbers2colors(datMeta$RIN, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RIN), max(datMeta$RIN)))
plotDendroAndColors(sampleTree, colors = cbind(as.numeric(datMeta$Group), sex_col, age_col, pmi_col, ph_col, rna_qual), groupLabels = c("Group", "Sex", "Age", "PMI", "pH",  "RIN"), cex.colorLabels=0.6, cex.dendroLabels=0.5)
dev.off()


# Collapse Rows
realGenes = !is.na(datProbes[,2])  #
table(realGenes)
datExpr = datExpr[realGenes,]; datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes[,2], rowID = datProbes[,1]) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes[,1])
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes[,2]
dim(datExpr)



## Regress Covariates
X =  model.matrix(~Group+Age+Sex+PMI+pH+RIN, data = datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,3:7]) %*% (as.matrix(beta[3:7,])))  # PMI + Sex +  Age + RNAdeg
datExpr = datExpr - t(to_regress)




save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_AAD_mayfield_normalized_CR_regressed.RData", datExpr,datMeta, datProbes)




