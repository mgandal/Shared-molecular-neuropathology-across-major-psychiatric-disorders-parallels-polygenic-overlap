#1o_RNAseq_ASD_4region.R
library(WGCNA); library(nlme); library(reshape); 
library(ggplot2); library(corrplot); library(biomaRt); library(cqn); library(limma)

options(stringsAsFactors = F)
theme_update(plot.title = element_text(hjust = 0.5))
plot_pdf=T
setwd("~/Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/")

##-----Load meta Data
if(!file.exists("./working_data/RNAseq/RNAseq_ASD_4region_normalized.Rdata")) {
  datMeta = read.csv("./raw_data/RNAseq_ASD/RNAseq_ASD_datMeta.csv")
  rownames(datMeta) = datMeta$Dissected_Sample_ID
  datMeta$Dx = factor(datMeta$Diagnosis_, levels=c("CTL", "ASD"))
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Brain_Bank = as.factor(datMeta$Brain_Bank)
  datMeta$Brain_Region = as.factor(datMeta$Region)
  datMeta$Brain_lobe = "Occipital"
  datMeta$Brain_lobe[datMeta$Brain_Region %in% c("BA4_6", "BA9", "BA24", "BA44_45")] = "Frontal"
  datMeta$Brain_lobe[datMeta$Brain_Region %in% c("BA3_1_2_5", "BA7")] = "Parietal"
  datMeta$Brain_lobe[datMeta$Brain_Region %in% c("BA38", "BA39_40", "BA20_37", "BA41_42_22")] = "Temporal"
  datMeta$Brain_lobe[datMeta$Brain_Region == "CBL"] = "CBL"
  datMeta$Brain_lobe=factor(datMeta$Brain_lobe, levels=c("Frontal", "Temporal", "Parietal", "Occipital", "CBL"))
  datMeta$RIN[is.na(datMeta$RIN)] = mean(datMeta$RIN, na.rm=T)
  
  ##----Create Sequencing Statistics
  datSeq = datMeta[,grep("Picard", colnames(datMeta))]
  rownames(datSeq) = rownames(datMeta)
  datSeq[,c(1:3,5:11, 17)] = log10(datSeq[,c(1:3,5:11, 17)])
  PC <- prcomp(na.omit(t(scale((datSeq),scale=T))), center=T)
  topPC <- PC$rotation[,1:8];
  colnames(topPC) = paste0("seq", colnames(topPC)) 
  varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
  topvar= varexp[1:8]
  for(i in 1:8) datMeta[,paste0("seqPC",i)] = topPC[,i]
  
  corrplot(cor(cbind(topPC, datSeq),method="spearman"))
  
  ##----Load Expression Data
  datExpr = read.csv("./raw_data/RNAseq_ASD/RNAseq_ASD_datExpr.csv", row.names = 1)
  colnames(datExpr)=gsub("[.]","-", gsub("X","",colnames(datExpr)))
  idx = match(rownames(datMeta),colnames(datExpr))
  datExpr = datExpr[,idx]
  
  ##------Annotate Probes
  getinfo <- c("ensembl_gene_id","external_gene_id","chromosome_name","start_position",
               "end_position","strand","band","gene_biotype","percentage_gc_content")
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="feb2014.archive.ensembl.org") ## Gencode v19
  datProbes <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values= rownames(datExpr),mart=mart)
  datProbes <- datProbes[match(rownames(datExpr),datProbes$ensembl_gene_id),]
  datProbes$length = datProbes$end_position - datProbes$start_position
  to_keep = !is.na(datProbes$length)
  datProbes = datProbes[to_keep,]
  datExpr = datExpr[to_keep,]
  rownames(datProbes) = datProbes$ensembl_gene_id
  
  #-----CQN normalize
  cqn.dat <- cqn(datExpr,lengths = as.numeric(datProbes$length), x = as.numeric(datProbes$percentage_gc_content),
                 lengthMethod=c("smooth"),sqn=FALSE) ## Run cqn with specified depths and with no quantile normalization
  cqn.dat <- cqn.dat$y + cqn.dat$offset ## Get the log2(Normalized FPKM) values
  datExpr.preCQN = datExpr
  datExpr= cqn.dat
  
  ##-----Filter out genes with low counts
  pres = apply(datExpr>1,1,sum) 
  to_keep = pres > 0.5 * ncol(datExpr)
  table(to_keep)
  datExpr = datExpr[to_keep,]
  datProbes = datProbes[to_keep,] 
  
  save(file="./working_data/RNAseq/RNAseq_ASD_4region_normalized.Rdata",datMeta,datExpr,datProbes,datSeq)
}

load("./working_data/RNAseq/RNAseq_ASD_4region_normalized.Rdata")

#---Check Consistency
all(rownames(datProbes)==rownames(datExpr))
all(colnames(datExpr)==rownames(datMeta))
all(rownames(datMeta)==rownames(datSeq))


#---Balance Groups by covariates, remove singular batches (none)
to_keep = (datMeta$Subject_ID != "AN03345") & !is.na(datMeta$Dx)
table(to_keep)
datMeta = datMeta[to_keep,]
datExpr = datExpr[,to_keep]
datSeq = datSeq[to_keep,]

if(plot_pdf) pdf(file="./results/figures//RNAseqQC/ASD_RNAseq-4Region_QC.pdf",width=15,height=8)
par(mfrow=c(3,8))

## Remove outliers based on network connectivity z-scores
normadj <- (0.5+0.5*bicor(datExpr))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,col=datMeta$Dx,pch=19)
legend("bottomleft",legend = levels(datMeta$Dx), col = 1:3,pch=19,cex=.7)
abline(h=-2, lty=2)
outliers = (z.ku < -2)
table(outliers)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]
datSeq = datSeq[!outliers,]

##----------------Post-normalization QC
boxplot(datExpr,col=as.numeric(datMeta$Dx), range=0)

#Histogram
i = 1; plot(density((datExpr[,i]), na.rm=T),ylim=c(0,0.3), col = as.numeric(datMeta$Dx[i]), main="Hist of Norm Exp By Dx", xlab = "norm exp")   
for(i in 2:dim(datExpr)[2]) {  lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Dx[i]),ylim=c(0,0.21))}
legend("topright", levels(datMeta$Dx), cex=0.7, col = 1:3, pch=19)

#MDS
mds = cmdscale(dist(t(datExpr)),eig=TRUE)
pc1=signif(100*mds$eig[1]^2/sum(mds$eig^2),2);pc2=signif(100*mds$eig[2]^2/sum(mds$eig^2), 2)
plot(mds$points,col=as.numeric(datMeta$Dx),pch=19,main = "MDS Plot HTSeq Counts Gene Dx", xlab=paste("PC1: ", pc1, "% variance", sep=""),ylab=paste("PC2: ", pc2, "% variance", sep=""))

## Covariate
plot(datMeta$Dx, col=1:2)
p = chisq.test(datMeta$Dx,datMeta$Brain_Region)$p.value; plot(datMeta$Dx ~ datMeta$Brain_Region, main=paste("Region, p=", signif(p,2)), ylab="", xlab="")

plot(datMeta$Dx ~ datMeta$Brain_Region, col=c(1,2))
A = anova(lm(as.numeric(datMeta$Age) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Dx, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Dx, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Dx, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   plot(as.numeric(datMeta$RIN) ~ datMeta$Dx, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datSeq$PicardQC_TOTAL_READS) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   plot(as.numeric(datSeq$PicardQC_TOTAL_READS) ~ datMeta$Dx, main=paste("ReadDepth, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datSeq$PicardQC_MEDIAN_5PRIME_TO_3PRIME_BIAS) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   plot(as.numeric(datSeq$PicardQC_MEDIAN_5PRIME_TO_3PRIME_BIAS) ~ datMeta$Dx, main=paste("RNA 3' Bias, p=", signif(p,2)), ylab="", xlab="")
for(s in paste("seqPC",1:2,sep=""))  {
  A = anova(lm(as.numeric((datMeta[,s])) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta[,s]) ~ datMeta$Dx, main=paste(s, ", p=", signif(p,2)), ylab="", xlab="")
}



##----Dendrogram
par(mfrow=c(1,1))
tree= hclust(as.dist(1-bicor(datExpr)),method="average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=F, lim=c(min(datMeta$Age),max(datMeta$Age)))
rin_col = numbers2colors(datMeta$RIN, blueWhiteRed(100), signed=F, centered=F, lim=c(min(datMeta$RIN,na.rm=T),max(datMeta$RIN,na.rm=T)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=F, lim=c(min(datMeta$PMI,na.rm=T),max(datMeta$PMI,na.rm=T)))
seqdepth_col = numbers2colors(datSeq[,"PicardQC_TOTAL_READS"], blueWhiteRed(100), signed=F, centered=F, lim=c(min(datSeq[,"PicardQC_TOTAL_READS"]),max(datSeq[,"PicardQC_TOTAL_READS"])))
plotDendroAndColors(tree, colors = cbind(as.factor(datMeta$Diagnosis_), as.numeric(datMeta$Brain_Region), sex_col, age_col, pmi_col, rin_col, seqdepth_col), 
                    groupLabels = c("Group","Region", "Sex", "Age", "PMI", "RIN", "seqDepth"), cex.colorLabels=0.6, cex.dendroLabels=0.5)

dev.off()

#-----Calculate differential gene expression (DGE) summary statistics
sumstats = vector(mode="list", length=5)
names(sumstats) = c("All", "Frontal", "Temporal", "Parietal", "Occipital")

datMeta$RIN2 = datMeta$RIN^2
mod=model.matrix(~Dx + Brain_Region + Age + Sex + PMI + RIN +RIN2+ seqPC1 + seqPC2, data=datMeta)
corfit <- duplicateCorrelation(datExpr,mod,block=datMeta$Subject_ID)
fit = eBayes(lmFit(datExpr, mod,block=datMeta$Subject_ID, correlation = corfit$consensus),trend=T, robust=T)
sumstats$All = topTable(fit,coef=2,number=Inf, sort.by="none")

for(region in c("Frontal", "Temporal", "Parietal", "Occipital")) {
  idx = which(datMeta$Brain_lobe==region)
  mod.region=model.matrix(~Dx + Age + Sex + PMI + RIN + RIN2+ seqPC1 + seqPC2, data=datMeta[idx,])
  fit.region = eBayes(lmFit(datExpr[,idx], mod.region), trend=T, robust=T)
  sumstats[[region]] = topTable(fit.region,coef=2,number=Inf, sort.by="none")
}

asd.rnaseq = do.call("cbind", sumstats)
write.csv(file="./results/tables/RNAseq_ASD_4region_sumstats.csv", asd.rnaseq)


#Assess qSVA
library(sva)
degMatrix=read.csv("./raw_data/RNAseq_ASD/RNAseq_ASD_riboZero_degRegion_featureCounts.csv",row.names=1)
colnames(degMatrix) = gsub("X", "", substr(colnames(degMatrix),2,10))
degMatrix = degMatrix[,match(colnames(datExpr), colnames(degMatrix))]
degMatrix=degMatrix[rowSums(degMatrix)>0,]
intervalLength=unlist(lapply(strsplit(rownames(degMatrix), "_"), function(x) { as.numeric(x[3]) - as.numeric(x[2]) } ))/1000
intervalLength = matrix(rep(intervalLength), ncol=ncol(degMatrix), nrow=nrow(degMatrix), byrow = FALSE)
normFactor=  matrix(rep(datMeta$PicardQC_TOTAL_READS/8e07),  ncol=ncol(degMatrix), nrow=nrow(degMatrix), byrow = TRUE)
degMatrix_fpk = degMatrix / normFactor
degMatrix_fpk80m = log2(1 + degMatrix_fpk / intervalLength)

degPca = prcomp(t(degMatrix_fpk80m))
k = num.sv(degMatrix_fpk80m, mod=rep(1,ncol(datExpr)),B = 100) #Using 5 qSVs


fit = eBayes(lmFit(datExpr, cbind(mod, degPca$x[, seq_len(k)]), block=datMeta$Subject_ID, correlation = corfit$consensus),trend=T,robust=T)
sumstats.qsva= topTable(fit, coef=2, number = Inf, sort.by = "none")
#corrplot(cor(cbind(degPca$x[,1:6], datSeq, datMeta$RIN, datMeta$seqPC1, datMeta$seqPC2, as.numeric(datMeta$Dx))))
write.csv(file="./results/tables/RNAseq_ASD_4region_sumstats_qSVA.csv", sumstats.qsva)



# Compare with microarray sumstats
asd = read.csv("~/CrossDisorder3//results/tables/Microarray_ASD_metaanalysis_030117.csv")
asd=asd[match(rownames(asd.rnaseq), asd$X),]
rho = cor.test(asd$beta, asd.rnaseq$All.logFC,method="spearman")
text = paste("rho=", signif(rho$estimate,2), "\n", "p=", signif(rho$p.value,2),sep="")

d= data.frame(ASD.array=asd$beta, Frontal=sumstats$Frontal$logFC, Temporal=sumstats$Temporal$logFC, Parietal=sumstats$Parietal$logFC, Occipital=sumstats$Occipital$logFC)
idx.sig = which(asd$fdr<.05)

d=d[idx.sig,]
colnames(d)[2] = paste("Frontal, rho=",signif(cor(d$Frontal, d$ASD.array,method="spearman",use='pairwise.complete.obs'),2))
colnames(d)[3] = paste("Temporal, rho=",signif(cor(d$Temporal, d$ASD.array,method="spearman",use='pairwise.complete.obs'),2))
colnames(d)[4] = paste("Parietal, rho=",signif(cor(d$Parietal, d$ASD.array,method="spearman",use='pairwise.complete.obs'),2))
colnames(d)[5] = paste("Occipital, rho=",signif(cor(d$Occipital, d$ASD.array,method="spearman",use='pairwise.complete.obs'),2))


d2 = melt(d,id.vars = 1)
colnames(d2)[2] = "Region"
plot.microarrayVsRNASeq4region=
  ggplot(d2, aes(x=ASD.array, y=value, color=Region)) + geom_point(alpha=0.5) + geom_smooth(method="lm") + 
  ggtitle("Region Specific RNAseq\nvs Cortical Array") + xlab("Microarray log2FC") + ylab("RNAseq log2FC") +
  geom_abline(slope=1,lty=2) + coord_fixed(ratio=1) + xlim(-2,2) + ylim(-2,2)
plot.microarrayVsRNASeq4region
if(plot_pdf) ggsave(filename = "./results/figures/Manuscript//Fig2D.pdf",plot.microarrayVsRNASeq4region, width=5,height=5)
