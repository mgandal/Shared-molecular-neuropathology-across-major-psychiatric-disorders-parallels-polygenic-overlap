#1e_Microrray_SCZ_BD_Iwamoto.R

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

home= "~/Github/" ## insert your GitHub home directory, ie. "C://Users/me/GitHub/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)

set.seed(100)

if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_Iwa_normalized.RData")) {
  
  #Load Meta Data
  datMeta = read.csv(file="./raw_data/Microarray/Iwa_GSE12649/Iwa_GSE12649_datMeta.csv", head=T)
  rownames(datMeta) = paste(datMeta$Filename,".CEL.gz",sep="")
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Race = as.factor(datMeta$Race)
  datMeta$Group = datMeta$Profile
  datMeta$Group[datMeta$Group=="Bipolar"]="BD"
  datMeta$Group[datMeta$Group=="Schizophrenia"]="SCZ"
  datMeta$Group[datMeta$Group=="Control"]="CTL"
  datMeta$Group[datMeta$Group=="Depression"]="MDD"
  datMeta$Group = factor(datMeta$Group, levels = c("CTL", "BD", "SCZ"))

  datMeta$RD = as.factor(datMeta$RD)
  datMeta = datMeta[!is.na(datMeta$Group),]
  datMeta$Study="SCZ.BD.Iwa"
  
  #Load Expression Data
  data.affy = ReadAffy(celfile.path="./raw_data/Microarray/Iwa_GSE12649/Iwa_GSE12649_raw/Iwa_GSE12649_raw/", filenames = paste(datMeta$Filename,".CEL.gz",sep=""))
  datExpr = affy::rma(data.affy, background =T, normalize=T)
  datExpr = exprs(datExpr)
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]
  
  
  RNAdeg = AffyRNAdeg(data.affy)
  plotAffyRNAdeg(RNAdeg)
  datMeta$RNAdeg = RNAdeg$slope
  
  sd = protocolData(data.affy)$ScanDate; sd = substring(sd,1,nchar(sd)-9)
  datMeta$Batch = as.factor(sd)
  table(sd)
  
  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "affy_hg_u133a"
  getinfo <- c("affy_hg_u133a", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat$affy_hg_u133a)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_Iwa_normalized.RData", datExpr, datMeta, datProbes)
  
  
  #QC-PreNorm
  datExpr =log2(exprs(data.affy))
  
  #Boxplot
  boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), main ="Boxplot", ylab = "Intensity")
  legend("topright", legend=c("CTL", "BAD", "SCZ"), col = 1:3, pch=19)
  
  # Histogram
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp")
  for(i in 2:dim(datExpr)[2])
    lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),)
  legend("topright", colnames(datExpr), cex=0.7, text.col = 1:dim(datExpr)[2])
  
  #MDS Plots
  mds = cmdscale(dist(t(datExpr)))
  
  plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
  legend("bottomleft", c("CTL", "BAD", "SCZ"), col=c(1:4), pch=16, cex=0.8)
  plot(mds, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="MDS Plot")
  
  #Covariate plots
  plot(datMeta$Group, ylim=c(0,40), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Race) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Race ~ datMeta$Group, main=paste("Race, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")
  
}  

load("./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_Iwa_normalized.RData")


pdf("./results/figures/MicroarrayQC/SCZ_BD_Iwa_QC.pdf",width=11,height=8.5)
par(mfrow=c(3,4))


#Balance Groups by batch, RNAdeg
table(datMeta$Batch)
to_keep = !(datMeta$Batch == "07/09/03") & (datMeta$RNAdeg < 3.8) & (datMeta$RNAdeg > 2) & !is.na(datMeta$Group) & (datMeta$Filename != "77-HGU133A")
table(to_keep)
datMeta = datMeta[to_keep,]; datExpr = datExpr[,to_keep]


##Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(z.ku, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Euclidean distance (z score)")
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

datMeta$Batch=factor(datMeta$Batch)
model = model.matrix(~Group+Sex+Age+pH+PMI+RNAdeg+Batch,data=datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_SCZ_BD_Iwa_normalized_balanced.RData", datExpr, datMeta, datProbes,model)


## Batch Correction
mod = model.matrix(~Group, data=datMeta)
batch = factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat


#QC Plots
boxplot(datExpr, range = 0, col= as.numeric(datMeta$Group), main ="Boxplot", ylab = "Intensity")
legend("topright", legend=c("CTL", "BD", "SCZ"), col = 1:3, pch=19)

# Histogram
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(datExpr)[2])
  lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),)

#MDS Plots
mds = cmdscale(dist(t(datExpr)))
plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
legend("bottomleft", c("CTL", "BD", "SCZ"), col=c(1:4), pch=16, cex=0.8)
plot(mds, col=as.numeric(as.factor(datMeta$Batch)), pch=16, main="MDS Plot by Batch")



## Plot potential Colinearities
plot(datMeta$Group, ylim=c(0,40), ylab="Number", main="Subjects")
A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Race) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Race ~ datMeta$Group, main=paste("Race, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")

par(mfrow=c(1,1))
sampleTree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI),max(datMeta$PMI)))
ph_col = numbers2colors(datMeta$pH, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$pH),max(datMeta$pH)))
rna_qual = numbers2colors(datMeta$RNAdeg, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RNAdeg), max(datMeta$RNAdeg)))
plotDendroAndColors(sampleTree, colors = cbind(as.numeric(datMeta$Group), as.numeric(datMeta$Batch), sex_col, age_col, as.numeric(datMeta$Race), pmi_col, ph_col, rna_qual), groupLabels = c("Group", "Batch", "Sex", "Age", "Race", "PMI", "pH", "RNAqual"), cex.colorLabels=0.6, cex.dendroLabels=0.5)
dev.off()



# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)  #
table(realGenes)
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$affy_hg_u133a) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hg_u133a)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)


X = model.matrix(~Group+Sex+Age+pH+PMI+RNAdeg,data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,4:8]) %*% (as.matrix(beta[4:8,]))) 
datExpr = datExpr - t(to_regress)

idx=  which(datMeta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
datMeta$Group.BD = datMeta$Group.SCZ = datMeta$Group
datMeta$Group.SCZ[c(which(datMeta$Group=="BD"), idx[!subset.scz])] = NA
datMeta$Group.BD[c(which(datMeta$Group=="SCZ"), idx[subset.scz])] = NA

save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_SCZ_BD_Iwa_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)
