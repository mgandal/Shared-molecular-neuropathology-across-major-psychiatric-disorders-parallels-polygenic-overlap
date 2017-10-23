##1g) Microarray_SCZ_Maycox
#--------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

home= "~/Github/" ## insert your GitHub home directory, ie. "C://Users/me/GitHub/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)

#1) Load Data
#------------
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_SCZ_Maycox_normalized.RData")) {
  
  data.affy = ReadAffy(celfile.path="./raw_data/Microarray/Maycox_GSE17612/Maycox_GSE17612_RAW/")
  datExpr = rma(data.affy, normalize=T, background=T, verbose=T)
  datExpr = exprs(datExpr)
  
  batch = as.factor(substr(protocolData(data.affy)$ScanDate,1,8))
  
  datMeta = read.csv("./raw_data/Microarray/Maycox_GSE17612/Maycox_GSE17612_datMeta.csv")
  rownames(datMeta) = datMeta$Chip
  datMeta$Group = factor(datMeta$Group, levels=c("CTL", "SCZ"))
  datMeta$Sex = as.factor(datMeta$Sex)
  RNAdeg = AffyRNAdeg(data.affy)
  plotAffyRNAdeg(RNAdeg)
  datMeta$RNAdegBias = RNAdeg$slope
  datMeta$Batch = batch
  datMeta$Study = "SCZ.maycox"
  datMeta$Subject = paste(datMeta$BrainBank, datMeta$ID,sep="_")
  
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  identifier <- "affy_hg_u133_plus_2"
  getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene", "external_gene_id",  "chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  
  idx = match(rownames(datExpr), geneDat$affy_hg_u133_plus_2)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  colnames(datExpr) = substring(colnames(datExpr),1,9)
  save(file="./working_data/Microarray/01_Normalized/Microarray_SCZ_Maycox_normalized.RData", datExpr, datMeta, datProbes)
  
  #QC Plots Pre-Norm
  datExpr = log2(exprs(data.affy))
  par(mfrow=c(2,2))
  
  #Boxplot
  boxplot(datExpr, col = as.numeric(datMeta$Group),range=0)
  
  #Histogram
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group)[i], main="Hist of Log2 Exp", xlab = "log2 exp",ylim=c(0,0.5))
  for(i in 2:dim(datExpr)[2])
    lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group)[i])
  
  #MDS Plot
  mds = cmdscale(dist(t(datExpr)))
  plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
  legend("topright", levels(datMeta$Group), col=c(1:2), pch=16, cex=0.8)
  
  #Dendrogram
  tree = hclust(dist(t(datExpr)),method="average")
  plot(tree, cex=0.5)
  
  #Covariates
  par(mfrow=c(3,3))
  plot(datMeta$Group, ylim=c(0,50), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")

  
} 

load("./working_data/Microarray/01_Normalized/Microarray_SCZ_Maycox_normalized.RData")

pdf("./results/figures/MicroarrayQC/SCZ_Maycox_QC.pdf",width=11,height=8.5)
par(mfrow=c(3,4))


## Remove singular batches
table(datMeta$Batch)

##Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(z.ku, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Connectivity (Z)")
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

#Balance Groups, Remove Singular Batches
#--> No singular batches or unbalanced covariates

model = model.matrix(~Group+Sex+Age+PMI+RNAdegBias+Batch, data = datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_SCZ_Maycox_normalized_balanced.RData", datExpr, datMeta, datProbes,model)


## Batch Correction
#plot(datMeta$Group  ~ datMeta$Batch, main="Batch x Group", xlab = "Batch", ylab="") # No group by batch confound
mod = model.matrix(~Group, data=datMeta)
batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat


#QC Plots
boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), main ="Boxplot", ylab = "Intensity")
legend("topright", legend=c("CTL", "SCZ"), col = 1:2, pch=19)

# Histogram
i = 1; plot(density((datExpr[,i]), na.rm=T), col = i, main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(datExpr)[2])
  lines(density((datExpr[,i]), na.rm=T), col = datMeta$Group[i])
legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:2)

#MDS Plots
mds = cmdscale(dist(t(datExpr)))
plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
legend("topleft", c("CTL", "SCZ"), col=c(1:2), pch=16, cex=0.8)
plot(mds, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="MDS Plot by Batch")

#Covariate Plots
plot(datMeta$Group, ylim=c(0,40), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")


#Cluster Dendrogram
par(mfrow=c(1,1))
sampleTree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI),max(datMeta$PMI)))
ph_col = numbers2colors(datMeta$pH, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$pH),max(datMeta$pH)))
rna_qual = numbers2colors(datMeta$RNAdegBias, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RNAdegBias), max(datMeta$RNAdegBias)))
batch_col = (as.numeric(as.factor(datMeta$Batch)))
plotDendroAndColors(sampleTree, colors = cbind(as.numeric(datMeta$Group), batch_col, sex_col, age_col, pmi_col, ph_col, rna_qual), groupLabels = c("Group", "Batch", "Sex", "Age", "PMI", "pH", "RNAqual"), cex.colorLabels=0.6, cex.dendroLabels=0.5)
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

#Regress Covariates
X = model.matrix(~Group+Sex+Age+PMI+RNAdegBias, data = datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
b = as.data.frame(t(beta))
to_regress = (as.matrix(X[,3:6]) %*% (as.matrix(beta[3:6,])))
datExpr = datExpr - t(to_regress)

save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_SCZ_Maycox_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)






