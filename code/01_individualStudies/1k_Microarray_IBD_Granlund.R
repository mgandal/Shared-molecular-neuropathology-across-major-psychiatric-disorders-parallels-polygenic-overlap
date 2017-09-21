##1k_Microarray_IBD_Granlund.R
### Granlund_Crohns: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0056818

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(lumi); library(limma); library(biomaRt); library(plyr); library(sva)

home= ## insert your GitHub home directory, ie. "C://Users/me/GitHub/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)

#1) Load Data
#------------
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_IBD_Granlund_normalized.RData")) {

  ## Load MetaData
  ## ---------------
  datMeta = read.delim("./raw_data/Microarray/Granlund/E-MTAB-184.sdrf.txt")
  datMeta$Group = gsub("unaffected", "CTL", datMeta$Characteristics.DiseaseState.)
  datMeta$Group = gsub("diseased", "IBD", datMeta$Group)
  datMeta$Group = factor(datMeta$Group, levels=c("CTL", "IBD"))
  datMeta$Sex = gsub("female", "F", datMeta$Characteristics.Sex.)
  datMeta$Sex = gsub("male", "M", datMeta$Sex)
  datMeta$Sex = factor(datMeta$Sex, levels=c("M", "F"))
  rownames(datMeta) = datMeta$Sample.Name
  datMeta$Study="Granlund.IBD"
  
  ## Load Expression Data
  ## ---------------------
  ## Split dataset into two to make it compatible with GitHub syncing
  #body = read.table("./raw_data/Microarray/Granlund/atle_mage_ml_atle_mage_ml_nonorm_nobkgd_ArrayExpress_DataFile.txt",skip=8)
  #body1 = body[1:24000,]
  #body2 = body[24001:dim(body)[1],]
  #head = read.table("./raw_data/Microarray/Granlund/atle_mage_ml_atle_mage_ml_nonorm_nobkgd_ArrayExpress_DataFile.txt",nrow=6,sep="\n")
  #colnames = read.table("./raw_data/Microarray/Granlund/atle_mage_ml_atle_mage_ml_nonorm_nobkgd_ArrayExpress_DataFile.txt",skip=7,nrow=1,sep="\t")
  #space = rep("",769)
  #head = cbind(head,matrix("",6,768))
  #names(space) <- names(head) <- names(body2) <- names(colnames) <- names(body1)
  #save(head,colnames,body1,file="Granlund_rawDat1.RData")
  #save(body2,file="Granlund_rawDat2.RData")
  
  load("Granlund_rawDat1.RData")
  load("Granlund_rawDat2.RData")
  rawDat = rbind(head,space,colnames,body1,body2)
  
  write.table(rawDat,file="raw_data/Microarray/Granlund/Granlund_rawDat.txt",sep="\t",col.names = FALSE,row.names = FALSE)
  data.lumi = lumiR("./raw_data/Microarray/Granlund/Granlund_rawDat.txt")
  datExpr = lumiN(data.lumi, method="quantile")
  datExpr = log2(exprs(datExpr))
  idx = match(colnames(datExpr), datMeta$Sample.Name)
  datMeta = datMeta[idx,]
  
  ## Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",  host="feb2014.archive.ensembl.org") 
  #f = listFilters(ensembl);   a = listAttributes(ensembl)
  identifier <- "illumina_humanht_12_v3"
  getinfo <- c("illumina_humanht_12_v3", "ensembl_gene_id","external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr),geneDat[,1])
  datProbes = geneDat[idx,]
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_IBD_Granlund_normalized.RData",datExpr,datProbes,datMeta)
  
  ## QC-PreNormalization
  datExpr = log2(exprs(data.lumi))
  par(mfrow=c(1,1))
  
  #Boxplot
  boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot", ylab = "Intensity")
  legend("topright", legend=levels(datMeta$Group), col = 1:length(levels(datMeta$Group)), pch=19, cex=0.7)
  
  #Histogram
  i = 1; 
  plot(density((datExpr[,i]), na.rm=T), xlim = c(5,15), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   
  for(i in 2:dim(datExpr)[2]) {     
    lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   
  legend("topright", levels(datMeta$Group), cex=0.7, col = 1:length(levels(datMeta$Group)), pch=19)
  
  #MDS Plot
  mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""))
  legend("bottomright", levels(datMeta$Group), col=1:length(levels(datMeta$Group)), pch=16, cex=0.8)
  
  #Dendrogram
  tree = hclust(as.dist(1-bicor(datExpr)), method = "average");   plot(tree, xlab="", cex=0.5)
  
  #Covariate Plot - we only have Gender information
  plot(datMeta$Group, ylim=c(0,150), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  
} 

load("./working_data/Microarray/01_Normalized/Microarray_IBD_Granlund_normalized.RData")

pdf("./results/figures/MicroarrayQC/IBD_Granlund_QC.pdf",width=11,height=8.5)
par(mfrow=c(3,4))


# Remove singular batches, remove confounders --> none 

# Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-2, lty=2)
datExpr_wOutliers = datExpr
datMeta_wOutliers = datMeta
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

## Batch Correction --> no batch information

to_keep = !is.na(datMeta$Sex)
datMeta = datMeta[to_keep,]
datExpr =datExpr[,to_keep]
model = model.matrix(~Group + Sex, datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_IBD_Granlund_normalized_balanced.RData",datExpr,datProbes,datMeta,model)


#QC Plots
boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), xaxt='n', xlab = "Array", main ="Boxplot Post-Normalization", ylab = "Intensity");   legend("topright", legend=levels(datMeta$Group), col = 1:length(levels(datMeta$Group)), pch=19)
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp"); 
for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),) } ;   legend("topright", levels(datMeta$Group), cex=0.7, text.col = 1:length(levels(datMeta$Group)))

mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)

#covariates
plot(datMeta$Group, ylim=c(0,150), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")

# Cluster Dendrogram
par(mfrow=c(1,1))
sampleTree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
plotDendroAndColors(sampleTree, colors = cbind(as.numeric(datMeta$Group), sex_col), groupLabels = c("DiseaseState", "Sex"), cex.colorLabels=0.6, cex.dendroLabels=0.5)

dev.off()

# Collapse Rows
realGenes = !is.na(datProbes[,2])  
table(realGenes)
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes[,2], rowID = datProbes[,1]) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes[,1])
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes[,2]
dim(datExpr)



#Regress covariate
X = model.matrix(~Group + Sex, datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = X[,3] %*% t(beta[3,])
datExpr = datExpr - t(to_regress)

save(file="./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_IBD_Granlund_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)
