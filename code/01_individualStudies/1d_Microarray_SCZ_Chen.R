##1D) Microarray_SCZ_BD_chen
#--------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

home= '~/Github/' ## insert your GitHub home directory, ie. "C://Users/me/GitHub/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)

set.seed(100)

#1) Load Data
#------------
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_chen_normalized.RData")) {
  
  datMeta=read.csv("./raw_data/Microarray/Chen_GSE35978/GSE35978_datMeta.csv")
  rownames(datMeta)= datMeta$Chip
  datMeta$Region = as.factor(datMeta$Region)
  datMeta$Group[datMeta$Group=="Bipolar"] = "BD"; datMeta$Group[datMeta$Group=="Control"] = "CTL";  datMeta$Group[datMeta$Group=="Schizophrenia"] = "SCZ";   datMeta$Group[datMeta$Group=="Depression"] = "MDD";
  datMeta$Group = factor(datMeta$Group, levels = c("CTL", "MDD", "BD", "SCZ"))
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$PMI = as.numeric(datMeta$PMI)
  datMeta$pH = as.numeric(datMeta$pH)
  datMeta$Sex=as.factor(datMeta$Sex)
  datMeta$Study="SCZ.BD.Chen"

  #Remove NA group
  to_keep = !is.na(datMeta$Group)
  datMeta = datMeta[to_keep,]
  
  ## Get parietal cortex samples only
  ctx_only = rownames(datMeta)[datMeta$Region=="PCTX"]
  all_samples = list.files("./raw_data/Microarray/Chen_GSE35978/GSE35978_CEL/")
  ctx_sample_files = all_samples[pmatch(ctx_only, all_samples)]
  
  ## Read in Expression Data
  data.affy = ReadAffy(celfile.path="./raw_data/Microarray/Chen_GSE35978/GSE35978_CEL", filenames=ctx_sample_files)
  datExpr = affy::rma(data.affy,verbose=T,normalize=T,background=T)
  datExpr = exprs(datExpr)
  colnames(datExpr) = substring(colnames(datExpr), 1, 9)
  
  ## Get RNA degradation
  RNAdeg = AffyRNAdeg(data.affy)
  plotAffyRNAdeg(RNAdeg)
  RNAdeg$sample.names = substring(RNAdeg$sample.names,1,9)
  idx=  match(rownames(datMeta), RNAdeg$sample.names)
  datMeta$RNAdeg = RNAdeg$slope[idx]
  
  ## Get batch information
  batch = as.factor(substring(protocolData(data.affy)$ScanDate,1,10))
  idx = match(rownames(datMeta), colnames(datExpr))
  datMeta$Batch = batch[idx]
  
  ## Align expression and phenoData matricies
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]
  
  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  a = listAttributes(ensembl)
  identifier <- "affy_hugene_1_0_st_v1"
  getinfo <- c("affy_hugene_1_0_st_v1", "ensembl_gene_id","hgnc_symbol", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat$affy_hugene_1_0_st_v1)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_chen_normalized.RData", datExpr, datMeta, datProbes)
  
  
  #QC Plots before normalization
  datExpr = log2(exprs(data.affy))
  par(mfrow=c(2,2))
  boxplot(datExpr, range=0, col = as.numeric(datMeta$Group), main = "Array Boxplot")
  legend("topright", legend=c("CTL", "Depression", "Bipolar", "Schizophrenia"), fill = c(1:4), cex=0.7)
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp")
  for(i in 2:dim(datExpr)[2])
    lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]))
  legend("topright", legend=c("CTL", "Depression", "Bipolar", "Schizophrenia"), fill = c(1:4), cex=0.7)
  mds = cmdscale(dist(t(datExpr)))
  plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
  plot(mds, col=as.numeric(datMeta$Batch), pch=16, main = "MDS by batch")
  
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
  A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")
  
  par(mfrow=c(1,1))
  tree = hclust(dist(t(datExpr)),method="average")
  plot(tree, cex=0.1)
}


load("./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_chen_normalized.RData")

pdf("./results/figures/MicroarrayQC/SCZ_BD_Chen_QC.pdf",width=11,height=8.5)
par(mfrow=c(3,4))


## Remove singular batch
table(datMeta$Batch)
to_keep = datMeta$Batch != "2010-04-14"
datExpr = datExpr[,to_keep]
datMeta = datMeta[to_keep,]


## Remove MDD group (confounded by pH)
to_keep = datMeta$Group != "MDD"
datExpr = datExpr[,to_keep]
datMeta = datMeta[to_keep,]

## Remove samples to match PMI, pH
to_keep = datMeta$Group!="MDD" & datMeta$PMI < 70
to_keep = to_keep & ((datMeta$pH < 6.9) | !(datMeta$Group=="CTL")) & (datMeta$pH > 5.8)
table(to_keep)
datExpr = datExpr[,to_keep]
datMeta = datMeta[to_keep,]
datMeta$Group = factor(datMeta$Group)

## Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(z.ku, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Network Connectivity (z score)")
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

datMeta$Batch = factor(datMeta$Batch)
model = model.matrix(~Group +Sex + Age + pH + PMI + RNAdeg + Batch, datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_SCZ_BD_chen_normalized_balanced.RData", datExpr, datMeta, datProbes,model)


## Batch Correction
#plot(datMeta$Group ~ datMeta$Batch, main = "Batch Balance", ylab ="", xlab = "Batch")
mod = model.matrix(~Group, data=datMeta)
batch = factor(datMeta$Batch)
datExpr = ComBat(datExpr, batch=batch, mod=mod)

#QC Plots
boxplot(datExpr, range=2, col = as.numeric(datMeta$Group), main = "Array Boxplot")
legend("topright", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)

# Histogram
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(datExpr)[2])
  lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]))
legend("topright", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)

#MDS Plots
mds = cmdscale(dist(t(datExpr)))
plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
legend("topright", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)
plot(mds, col=as.numeric(datMeta$Batch), pch=16, main = "MDS by batch")

#Covariate Plots
plot(datMeta$Group, ylim=c(0,50), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")

#Cluster Dendrogram
par(mfrow=c(1,1))
sampleTree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI),max(datMeta$PMI)))
ph_col = numbers2colors(datMeta$pH, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$pH),max(datMeta$pH)))
batch_col = (as.numeric(as.factor(datMeta$Batch)))
rna_col = numbers2colors(datMeta$RNAdeg, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RNAdeg),max(datMeta$RNAdeg)))
plotDendroAndColors(sampleTree, colors = cbind(as.numeric(datMeta$Group), batch_col, sex_col, age_col, pmi_col, ph_col, rna_col), groupLabels = c("Group", "Batch", "sex", "age", "pmi", "pH", "RNAdeg"), cex.colorLabels=0.6, cex.dendroLabels=0.2)
#legend("topleft", c("CTL", "Depression", "Bipolar", "Schizophrenia"), col=c(1:4), pch=16, cex=0.8)

dev.off()

# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)  #
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$affy_hugene_1_0_st_v1) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hugene_1_0_st_v1)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id


# Regress Covariates
X = model.matrix(~Group +Sex + Age + pH + PMI + RNAdeg, datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,4:8]) %*% (as.matrix(beta[4:8,])))  # PMI + Sex +  Age + RNAdeg
datExpr = datExpr - t(to_regress)



idx=  which(datMeta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
datMeta$Group.BD = datMeta$Group.SCZ = datMeta$Group
datMeta$Group.SCZ[c(which(datMeta$Group=="BD"), idx[!subset.scz])] = NA
datMeta$Group.BD[c(which(datMeta$Group=="SCZ"), idx[subset.scz])] = NA


save(file="./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_SCZ_BD_Chen_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)

