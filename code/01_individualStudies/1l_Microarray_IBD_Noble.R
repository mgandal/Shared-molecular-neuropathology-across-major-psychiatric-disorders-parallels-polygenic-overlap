##1l) Microarray_IBD_Noble
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
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_IBD_noble_normalized.RData")) {
  
  datMeta=read.csv("./raw_data/Microarray/Noble_GSE11223/Noble_GSE11223_datMeta2.csv")
  rownames(datMeta)= datMeta$GSM
  datMeta$Sex = as.factor(datMeta$Gender)
  datMeta$Ethnicity = as.factor(datMeta$Ethnicity)
  datMeta$Group = as.factor(gsub("Normal", "CTL", gsub("UC", "IBD", datMeta$Disease)))
  datMeta$Tissue = as.factor(datMeta$Anatomic_Location)
  datMeta$Inflammation_State = as.factor(datMeta$Inflammation_State)
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$Study= "Noble.IBD"
  datMeta$Array = "Agilent_G4112A"
  datMeta$Batch = as.factor(datMeta$Run_Date)
  
  ## Read in Expression Data
  filenames=list.files("./raw_data/Microarray/Noble_GSE11223/Noble_GSE11223_RAW/")
  RG = read.maimages(files=filenames,source="agilent",path="./raw_data/Microarray/Noble_GSE11223/Noble_GSE11223_RAW")
  plotMD(RG)
  RGb = backgroundCorrect(RG,method= "normexp",offset=50)
  plotDensities(RGb)
  MA <- normalizeWithinArrays(RGb,method="loess")
  plotDensities(MA)
  MA.q <- normalizeBetweenArrays(MA,method="quantile")
  plotDensities(MA.q)
  datExpr = getEAWP(MA.q)$exprs
  rownames(datExpr) = getEAWP(MA.q)$probes$ProbeName
  colnames(datExpr) = gsub(".txt", "", colnames(datExpr))
  
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]
  
  to_remove = which(datMeta$Gender=="unknown")
  datMeta = datMeta[-to_remove,]
  datExpr = datExpr[,-to_remove]
  
  
  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  a = listAttributes(ensembl); f=listFilters(ensembl)
  identifier <- "efg_agilent_wholegenome_4x44k_v1"
  getinfo <- c("efg_agilent_wholegenome_4x44k_v1", "ensembl_gene_id","hgnc_symbol", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat$efg_agilent_wholegenome_4x44k_v1)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_IBD_noble_normalized.RData", datExpr, datMeta, datProbes)
  
  
  #QC Plots before normalization
  datExpr = log2(getEAWP(MA)$exprs)
  par(mfrow=c(2,2))
  boxplot(datExpr, range=0, col = as.numeric(datMeta$Group), main = "Array Boxplot")
  legend("topright", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp")
  for(i in 2:dim(datExpr)[2])
    lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]))
  legend("topright", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)
  mds = cmdscale(dist(t(datExpr)))
  plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
  plot(mds, col=as.numeric(datMeta$Batch), pch=16, main = "MDS by batch")
  
  par(mfrow=c(3,3))
  plot(datMeta$Group, ylim=c(0,150), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Ethnicity) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Group ~ datMeta$Ethnicity, main=paste("Ethnicity, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Tissue) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Group ~ datMeta$Tissue, main=paste("Tissue, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
  
  par(mfrow=c(1,1))
  tree = hclust(dist(t(datExpr)),method="average")
  plot(tree, cex=0.1)
}

load("./working_data/Microarray/01_Normalized/Microarray_IBD_noble_normalized.RData")

pdf("./results/figures/MicroarrayQC/IBD_noble_QC.pdf",width=11,height=8.5)
par(mfrow=c(3,4))


## Remove singular batch --> no singular batches
table(datMeta$Batch)

## Remove group x covariate confounds
to_keep = (datMeta$Ethnicity != "ASIAN") & (datMeta$Ethnicity != "JEWISH")

to_keep = to_keep & !((datMeta$Group=="IBD") & (datMeta$Sex=="M") & (datMeta$Age> 50))
to_keep = to_keep & !((datMeta$Group=="CTL") & (datMeta$Sex=="F") & (datMeta$Age< 30))
table(to_keep)
datExpr = datExpr[,to_keep]
datMeta = datMeta[to_keep,]
datMeta$Ethnicity = factor(datMeta$Ethnicity)
datMeta$Sex = factor(datMeta$Sex)

## Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(z.ku, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Network Connectivity (z score)")
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

model = model.matrix(~Group + Sex + Age + Tissue + Batch, datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_IBD_noble_normalized_balanced.RData", datExpr, datMeta, datProbes,model)


## Batch Correction
mod = model.matrix(~Group, data=datMeta)
batch = factor(datMeta$Batch)
datExpr = ComBat(datExpr, batch=batch, mod=mod)

#QC Plots
boxplot(datExpr, range=0, col = as.numeric(datMeta$Group), main = "Array Boxplot")
legend("topright", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)

# Histogram
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(datExpr)[2])
  lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]))
legend("topright", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)

#MDS Plots
mds = cmdscale(dist(t(datExpr)))
plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
legend("bottomleft", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)
plot(mds, col=as.numeric(datMeta$Batch), pch=16, main = "MDS by batch")

#Covariate Plots
plot(datMeta$Group, ylim=c(0,150), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Tissue) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Group ~ datMeta$Tissue, main=paste("Tissue, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")

#Cluster Dendrogram
par(mfrow=c(1,1))
sampleTree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
batch_col = (as.numeric(as.factor(datMeta$Batch)))
tissue_col = (as.numeric(as.factor(datMeta$Tissue)))
plotDendroAndColors(sampleTree, colors = cbind(as.numeric(datMeta$Group), batch_col, sex_col, age_col, tissue_col, datMeta$Ethnicity), groupLabels = c("Group", "Batch", "sex", "age", "tissue", "ethnicity"), cex.colorLabels=0.6, cex.dendroLabels=0.2)

dev.off()




# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)  & !duplicated(rownames(datExpr))
table(realGenes)
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$efg_agilent_wholegenome_4x44k_v1) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$efg_agilent_wholegenome_4x44k_v1)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id


# Regress Covariates
X = model.matrix(~Group + Sex + Age + Tissue, datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,3:7]) %*% (as.matrix(beta[3:7,])))
datExpr = datExpr - t(to_regress)

save(file="./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_IBD_noble_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)

