##1p_RNAseq_SCZ_BD_GVEX.R

library(biomaRt); library(WGCNA); library(cqn); library(corrplot); library(ggplot2); library(sva); library(limma)
options(stringsAsFactors = FALSE); theme_update(plot.title = element_text(hjust = 0.5))
plot_pdf=T

setwd("~/Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/")

##-----Load meta Data
if(!file.exists("./working_data/RNAseq/RNAseq_SCZ_BD_GVEX_normalized.Rdata")) {
  datMeta = read.csv("./raw_data/RNAseq_GVEX/RNAseq_SCZ_BD_GVEX_datMeta.csv")
  rownames(datMeta) = datMeta$Individual_ID..RNAseq.library.BID.
  datSeq = datMeta[,grep("Picard",colnames(datMeta))]
  rownames(datSeq) = rownames(datMeta)
  datMeta$Group = factor(datMeta$Group, levels=c("CTL", "BD", "SCZ"))
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$LibraryBatch = as.factor(datMeta$LibraryBatch); 
  datMeta$Ethnicity = as.factor(datMeta$Ethnicity)
  datMeta$BrainBank=as.factor(datMeta$BrainBank)
  
  #-----Create sequencing statistics
  datSeq[,c(1:3,5:6,11:12)] = log10(datSeq[,c(1:3,5:6,11:12)])
  PC <- prcomp(na.omit(t(scale((datSeq),scale=T))), center=T)
  topPC <- PC$rotation[,1:8];  colnames(topPC) = paste0("seq", colnames(topPC))
  varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
  topvar= varexp[1:8]; 
  for(i in 1:8) datMeta[,paste0("seqPC",i)] = topPC[,i]
  
  #corrplot(cor(cbind(topPC, datSeq),method="spearman"))
  
  ##----Load Expression Data
  datExpr = read.csv("./raw_data/RNAseq_GVEX//RNAseq_SCZ_BD_GVEX_datExpr.csv",row.names=1)
  colnames(datExpr)= gsub("[.]","-", gsub("X","",colnames(datExpr)))
  datExpr = datExpr[grep("ENSG",rownames(datExpr)),]
  datExpr = datExpr[,match(rownames(datMeta),colnames(datExpr))]
  
  ##------Annotate Probes
  getinfo <- c("ensembl_gene_id","external_gene_id","chromosome_name","start_position",
               "end_position","strand","band","gene_biotype","percentage_gc_content")
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="feb2014.archive.ensembl.org") ## Gencode v19,
  datProbes <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values= rownames(datExpr),mart=mart)
  datProbes <- datProbes[match(rownames(datExpr),datProbes[,1]),]
  datProbes = datProbes[match(rownames(datExpr), datProbes$ensembl_gene_id),]
  rownames(datProbes) = datProbes$ensembl_gene_id
  datProbes$length = datProbes$end_position - datProbes$start_position
  to_keep = !is.na(datProbes$length)
  datExpr = datExpr[to_keep,]; datProbes = datProbes[to_keep,]
  
  

  #-----CQN normalize
  cqn.dat <- cqn(datExpr,lengths = as.numeric(datProbes$length), x = as.numeric(datProbes$percentage_gc_content),
                 lengthMethod=c("smooth"),sqn=FALSE) ## Run cqn with specified depths and with no quantile normalization
  cqn.dat <- cqn.dat$y + cqn.dat$offset ## Get the log2(Normalized FPKM) values
  datExpr.preCQN = datExpr
  datExpr= cqn.dat
  
  ##-----Filter out genes with low counts: threshold of 1 FPKM in at least 50% of samples
  pres = apply(datExpr>1,1,sum) 
  to_keep = (pres > 0.5*ncol(datExpr)) 
  table(to_keep)
  datExpr = datExpr[to_keep,]
  datProbes = datProbes[to_keep,]
  
  
  save(file="./working_data/RNAseq/RNAseq_SCZ_BD_GVEX_normalized.Rdata",datExpr, datMeta, datProbes,datSeq)
}

load("./working_data/RNAseq/RNAseq_SCZ_BD_GVEX_normalized.Rdata")

#---Check Consistency
all(rownames(datProbes)==rownames(datExpr))
all(colnames(datExpr)==rownames(datMeta))
all(rownames(datMeta)==rownames(datSeq))

#---Balance Groups by covariates, remove singular batches 
##Remove singular batches
table(datMeta$LibraryBatch)
to_remove = (datMeta$LibraryBatch=="11/6/2014") #Remove singular batch
datExpr= datExpr[,!to_remove]; datMeta = datMeta[!to_remove,]; datSeq = datSeq[!to_remove,]
datMeta$LibraryBatch = factor(datMeta$LibraryBatch)

##Balance Groups
to_remove = FALSE
to_remove = to_remove | (datMeta$pH > 6.75 & datMeta$Group=="CTL" & datMeta$Sex=="M")
to_remove = to_remove | (datMeta$Group=="BD" & datMeta$Sex=="F" & datMeta$RIN<7) | is.na(datMeta$RIN)
table(to_remove)
datMeta =  datMeta[!to_remove,]
datExpr = datExpr[,!to_remove]
datSeq = datSeq[!to_remove,]

if(plot_pdf) pdf(file="./results/figures/RNAseqQC//SCZ_BD_RNAseq-GVEX_QC.pdf",width=15,height=8)
par(mfrow=c(3,8))

## Remove outliers based on network connectivity z-scores
normadj <- (0.5+0.5*bicor(datExpr))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,col=datMeta$Group,pch=19)
legend("bottomleft",legend = levels(datMeta$Group), col = 1:3,pch=19,cex=.7)
abline(h=-2, lty=2)
outliers = (z.ku < -2)
table(outliers)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]
datSeq = datSeq[!outliers,]
 
# Batch Correction 
batch = factor(datMeta$LibraryBatch)
mod = model.matrix(~datMeta$Group)
datExpr.combat = ComBat(dat = as.matrix(datExpr),batch = batch,mod = mod)
datExpr.cqn.preCombat = datExpr
datExpr = as.data.frame(datExpr.combat)


##----------------Post-normalization QC
boxplot(datExpr,col=as.numeric(datMeta$Group), range=0)

#Histogram
i = 1; plot(density((datExpr[,i]), na.rm=T),ylim=c(0,0.3), col = as.numeric(datMeta$Group[i]), main="Hist of Norm Exp By Dx", xlab = "norm exp")   
for(i in 2:dim(datExpr)[2]) {  lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]),ylim=c(0,0.21))}
legend("topright", levels(datMeta$Group), cex=0.7, col = 1:3, pch=19)

#MDS
mds = cmdscale(dist(t(datExpr)),eig=TRUE)
plot(mds$points,col=as.numeric(datMeta$Group),pch=19,main = "MDS Plot HTSeq Counts Gene Dx")

#Covariate plots
plot(datMeta$Group, ylim=c(0,60), ylab="Number", main="Subjects", col=1:3)
A = anova(lm((datMeta$AgeDeath) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$AgeDeath ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$Group, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="")
A = chisq.test(datMeta$LibraryBatch, datMeta$Group); p = A$p.value;   plot(datMeta$Group ~ datMeta$LibraryBatch, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
A = chisq.test(datMeta$Ethnicity, datMeta$Group); p = A$p.value;   plot(datMeta$Group ~ datMeta$Ethnicity, main=paste("Ethnicity, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$BrainBank) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$BrainBank ~ datMeta$Group, main=paste("BrainBank, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$ReadDepth) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$ReadDepth ~ datMeta$Group, main=paste("ReadDepth, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1];   plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNA 3' Bias, p=", signif(p,2)), ylab="", xlab="")
for(s in paste("seqPC",1:2,sep=""))  {
  A = anova(lm(as.numeric((datMeta[,s])) ~ datMeta$Group));   p = A$"Pr(>F)"[1];   plot((datMeta[,s]) ~ datMeta$Group, main=paste(s, ", p=", signif(p,2)), ylab="", xlab="")
}

#Dendrogram
par(mfrow=c(1,1))
tree= hclust(dist(t(datExpr)),method="average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
datMeta$Age = datMeta$AgeDeath; age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=F, lim=c(min(datMeta$Age),max(datMeta$Age)))
rin_col = numbers2colors(datMeta$RIN, blueWhiteRed(100), signed=F, centered=F, lim=c(min(datMeta$RIN,na.rm=T),max(datMeta$RIN,na.rm=T)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=F, lim=c(min(datMeta$PMI,na.rm=T),max(datMeta$PMI,na.rm=T)))
seqdepth_col = numbers2colors(datSeq[,"PicardQC_TOTAL_READS"], blueWhiteRed(100), signed=F, centered=F, lim=c(min(datSeq[,"PicardQC_TOTAL_READS"]),max(datSeq[,"PicardQC_TOTAL_READS"])))
rnadeg_col = numbers2colors(datSeq[,"PicardQC_MEDIAN_5PRIME_TO_3PRIME_BIAS"], blueWhiteRed(100), signed=F, centered=F, lim=c(min(datSeq[,"PicardQC_MEDIAN_5PRIME_TO_3PRIME_BIAS"]),max(datSeq[,"PicardQC_MEDIAN_5PRIME_TO_3PRIME_BIAS"])))
plotDendroAndColors(tree, colors = cbind(as.factor(datMeta$Group), as.numeric(datMeta$BrainBank), as.numeric(datMeta$LibraryBatch), as.numeric(datMeta$Ethnicity), sex_col, age_col, pmi_col, rin_col, seqdepth_col,rnadeg_col), 
                    groupLabels = c("Group","BrainBank", "Batch", "Ethnicity", "Sex", "Age", "PMI", "RIN", "seqDepth", "RNA 3' bias"), cex.colorLabels=0.6, cex.dendroLabels=0.5)

if(plot_pdf) dev.off()

#-------Calculate differential expression:  linear model
sumstats = vector(mode="list", length=2); names(sumstats)= c("SCZ", "BD")
datMeta$RIN2 = datMeta$RIN^2
mod = model.matrix(~Group+AgeDeath+Ethnicity+Sex+RIN+PMI+pH+RIN2+seqPC1+seqPC2, data=datMeta)
fit = eBayes(lmFit(datExpr, mod),trend = T, robust=T)
sumstats$BD = topTable(fit, coef=2, number = Inf, sort.by = "none",confint = T)
sumstats$SCZ = topTable(fit, coef=3, number = Inf, sort.by = "none", confint = T)
rnaseq.gvex = do.call("cbind", sumstats)
write.csv(file="./results/tables/RNAseq_SCZ_BD_GVEX.csv", rnaseq.gvex)


#qSVA assessment
degMatrix=read.csv("./raw_data/RNAseq_GVEX/RNAseq_GVEX_riboZero_degRegion_featureCounts.csv",row.names=1)
colnames(degMatrix) = gsub("\\.", "-", substr(colnames(degMatrix),2,10))
degMatrix = degMatrix[,match(colnames(datExpr), colnames(degMatrix))]
degMatrix=degMatrix[rowSums(degMatrix)>0,]
intervalLength=unlist(lapply(strsplit(rownames(degMatrix), "_"), function(x) { as.numeric(x[3]) - as.numeric(x[2]) } ))/1000
intervalLength = matrix(rep(intervalLength), ncol=ncol(degMatrix), nrow=nrow(degMatrix), byrow = FALSE)
normFactor=  matrix(rep(datMeta$PicardQC_TOTAL_READS/8e07),  ncol=ncol(degMatrix), nrow=nrow(degMatrix), byrow = TRUE)
degMatrix_fpk = degMatrix / normFactor
degMatrix_fpk80m = log2(1 + degMatrix_fpk / intervalLength)
boxplot(degMatrix_fpk80m)

degPca = prcomp(t(degMatrix_fpk80m))
k = num.sv(degMatrix_fpk80m, mod=rep(1,ncol(datExpr)), B = 100)   #using 11 qSVs

fit.qsva = eBayes(lmFit(datExpr, cbind(mod, degPca$x[,1:k])),trend = T,robust=T) 
sumstats.qsva = vector(mode="list", length=2); names(sumstats.qsva)= c("SCZ", "BD")
sumstats.qsva$BD = topTable(fit.qsva, coef=2, number = Inf, sort.by = "none", confint = T)
sumstats.qsva$SCZ = topTable(fit.qsva, coef=3, number = Inf, sort.by = "none", confint = T)

write.csv(file="./results/tables/RNAseq_SCZ_BD_GVEX.qsva.csv", do.call("cbind", sumstats.qsva))
