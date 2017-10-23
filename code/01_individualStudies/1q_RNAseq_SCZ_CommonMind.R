#1q_RNAseq_SCZ_CommonMind.R

rm(list=ls()); options(stringsAsFactors = FALSE)
library(biomaRt); library(WGCNA); library(cqn); library(corrplot); library(ggplot2); library(sva); library(limma)
plot_pdf=T

setwd("~/Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/")


##-----Load expression and meta data
if(!file.exists("./working_data/RNAseq/RNAseq_SCZ_BD_CMC_normalized.RData")) {
  datMeta1 = read.csv("./raw_data/RNAseq_CMC/CommonMind-release1//CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv")
  datMeta2 = read.csv("./raw_data/RNAseq_CMC//CommonMind-release1/CMC_MSSM-Penn-Pitt_Clinical.csv")
  datMeta3 = read.delim("./raw_data/RNAseq_CMC//CommonMind-release1/CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_GemToolsAncestry.tsv")
  datMeta = merge(datMeta1, datMeta2, by="Individual_ID")
  datMeta = merge(datMeta, datMeta3, by="Genotyping_Sample_ID")
  
  rownames(datMeta) = datMeta$DLPFC_RNA_Sequencing_Sample_ID.x
  datMeta$Dx=factor(datMeta$Dx, levels=c("Control", "SCZ", "BP"))
  datMeta$Batch = as.factor(datMeta$DLPFC_RNA_Sequencing_Library_Batch)
  datMeta$RIN2 = datMeta$DLPFC_RNA_isolation_RIN^2
  
  Age = datMeta$Age_of_Death
  Age[Age=="90+"] = 90
  Age = as.numeric(Age)
  datMeta$Age = Age
  
  ##----Calculate sequencing statistics
  datSeq = datMeta[,c(20:29)]; rownames(datSeq) = datMeta$DLPFC_RNA_Sequencing_Sample_ID.x
  picard= read.csv("./raw_data/RNAseq_CMC//CommonMind-release1/PicardToolsQC.csv",row.names = 1)
  idx = match(rownames(datSeq), rownames(picard))
  datSeq = cbind(datSeq, picard[idx,])
  a = apply(is.na(datSeq),1,any) | (rownames(datSeq)=="PITT_RNA_BP_PFC_686") #Duplicate sample
  datSeq = datSeq[!a,]; datMeta=datMeta[!a,]
  idx = which(datSeq[1,]>100)
  datSeq[,idx] = log10(datSeq[,idx])
  PC <- prcomp(t(scale(datSeq, scale=T)), center=T)
  topPC <- PC$rotation[,1:10];
  varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
  topvar= varexp[1:10]
  for(i in 1:10) datMeta[,paste0("seqPC",i)] = topPC[,i]
  
  #corrplot(cor(cbind(topPC, datSeq),method="spearman"),tl.cex = .5)
  
  
  ##----Load Expression Data: 56632 genes x 613 samples
  datExpr = read.table("./raw_data/RNAseq_CMC//CommonMind-release1/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw.tsv", skip =1, row.names = 1)
  id = read.table("./raw_data/RNAseq_CMC//CommonMind-release1/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw.tsv", nrow=1)
  colnames(datExpr) = id
  idx =match(rownames(datMeta),colnames(datExpr))
  datExpr=datExpr[,idx]
  
  
  ##------Annotate Probes
  getinfo <- c("ensembl_gene_id","external_gene_id","chromosome_name","start_position",
               "end_position","strand","band","gene_biotype","percentage_gc_content")
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="feb2014.archive.ensembl.org")
  datProbes <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values= rownames(datExpr),mart=mart)
  datProbes = datProbes[match(rownames(datExpr), datProbes$ensembl_gene_id),]
  datProbes$length = datProbes$end_position - datProbes$start_position
  to_keep = !is.na(datProbes$length)
  datProbes = datProbes[to_keep,]
  datExpr = datExpr[to_keep,]
  rownames(datProbes) = datProbes$ensembl_gene_id
  
  #Check consistency
  all(colnames(datExpr) == rownames(datMeta), (rownames(datSeq) == rownames(datMeta)))
  
  #-----CQN normalize
  cqn.dat <- cqn(datExpr,lengths = as.numeric(datProbes$length), x = as.numeric(datProbes$percentage_gc_content),
                 lengthMethod=c("smooth"),sqn=FALSE) ## Run cqn with specified depths and with no quantile normalization
  cqn.dat <- cqn.dat$y + cqn.dat$offset ## Get the log2(Normalized FPKM) values
  datExpr= as.data.frame(cqn.dat)
  
  ##-----Filter out genes with low counts: threshold of 1 FPKM in at least 50% of samples
  pres = apply(datExpr>1,1,sum) 
  to_keep = (pres > 0.5*ncol(datExpr)) ## 16780 genes
  table(to_keep)
  datExpr = datExpr[to_keep,]
  datProbes = datProbes[to_keep,]

  
  
  to_keep = !is.na(datMeta$Dx)
  datMeta = datMeta[to_keep,]; datExpr = datExpr[,to_keep]; datSeq = datSeq[to_keep,]
  
  save(file="./working_data/RNAseq/RNAseq_SCZ_BD_CMC_normalized.RData",datExpr,datMeta,datProbes,datSeq)
}

load("./working_data/RNAseq/RNAseq_SCZ_BD_CMC_normalized.RData")

if(plot_pdf) pdf("./results/figures/RNAseqQC///SCZ_BD_RNAseq-CommonMind_QC.pdf",width=15,height=8)
par(mfrow=c(3,4), mar=c(2,5,2,2), oma=c(0,0,0,0))

##----------------QC Post-Normalization, Outlier Removal ----------------
## Remove outliers based on network connectivity z-scores
normadj <- (0.5+0.5*bicor(datExpr))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,col=datMeta$Dx,pch=19, main="Outlier Detection", xlab="", ylab="Standardized Network\nConnectivity (Z score)")
legend("bottomright",legend = levels(datMeta$Dx), col = 1:3,pch=19,cex=.7)
abline(h=-2, lty=2)
outliers = (z.ku < -2)
table(outliers)
datExpr = datExpr[,!outliers]; datMeta = datMeta[!outliers,]; datSeq = datSeq[!outliers,]

mod = model.matrix(~datMeta$Dx)
combat = ComBat(datExpr, batch=factor(datMeta$Batch), mod=mod, prior.plots = F)
datExpr = combat


##Balance SCZ by Covariates
to_keep = !(datMeta$Dx=="BP")
to_keep = to_keep & !(datMeta$Dx=="SCZ" & datMeta$PMI_hrs > 35) & !(datMeta$Dx=="Control" & datMeta$PMI_hrs < 5); 
to_keep =to_keep & !(datMeta$Dx=="SCZ" & datMeta$Age > 89) & !(is.na(datMeta$pH)) & !(datMeta$Dx=="SCZ" & datMeta$pH < 6)
to_keep =to_keep & !(datMeta$Dx=="SCZ" & datMeta$DLPFC_RNA_isolation_RIN < 6) & !(datMeta$Dx=="Control" & datMeta$DLPFC_RNA_isolation_RIN > 9); 
to_keep =to_keep & !(datMeta$Dx=="Control" & datMeta$DLPFC_RNA_isolation_28S_18S_ratio > 2)
to_keep =to_keep & !(datMeta$Dx=="SCZ" & datMeta$DLPFC_RNA_Sequencing_Expression_Profiling_Efficiency < 30) & !(datMeta$Dx=="Control" & datMeta$DLPFC_RNA_Sequencing_Expression_Profiling_Efficiency > 68)
to_keep = to_keep & !(datMeta$Dx=="SCZ" & datMeta$seqPC1 < -.05) & !(datMeta$Dx=="Control" & datMeta$seqPC1 > .08)
to_keep = to_keep & !(datMeta$Dx=="SCZ" & datMeta$Institution=="Penn" & datMeta$DLPFC_RNA_isolation_RIN < 7.4)
to_keep = to_keep
table(to_keep) 
datMeta_scz = datMeta[to_keep,]; datExpr_scz = datExpr[,to_keep]
datMeta_scz$Dx = factor(datMeta_scz$Dx)


##Balance BD by Covariates
to_keep = (datMeta$Institution=="Pitt") & !(datMeta$Dx=="SCZ")
to_keep = to_keep & !(datMeta$Dx=="Control" & datMeta$DLPFC_RNA_isolation_RIN >8.7) 
table(to_keep)
datMeta_bd = datMeta[to_keep,]; datExpr_bd = datExpr[,to_keep]
datMeta_bd$Dx = factor(datMeta_bd$Dx)


boxplot(datExpr,range=0, col=as.numeric(datMeta$Dx), main="Expression Boxplot)",xaxt = "n")
legend("topright", levels(datMeta$Dx), col=c(1:length(levels(datMeta$Dx))), pch=16, cex=0.8)

plot(density(datExpr[,1]), col=as.numeric(datMeta$Dx)[1], main="Density", ylim=c(0,0.2))
for(i in 2:dim(datExpr)[2]) {
  lines(density(datExpr[,i]), col=as.numeric(datMeta$Dx)[i])  
}
legend("topright", levels(datMeta$Dx), col=c(1:length(levels(datMeta$Dx))), pch=16, cex=0.8)

par(mfrow=c(3,8), mar=c(2,5,2,2), oma=c(0,0,0,0))
#Covariate Plots
for(dz in c("SCZ", "BD")) {
  if(dz=="SCZ") {
    datMeta = datMeta_scz; datExpr = datExpr_scz
    plot_cols = c("black", "red")
  } else {
    datMeta = datMeta_bd; datExpr = datExpr_bd
    plot_cols = c("black", "green")
  }
  mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  plot(mds$points, col=plot_cols, pch=16, main="MDS: Dx", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""));  
  legend("bottomright", levels(datMeta$Dx), col=plot_cols, pch=16, cex=0.8)
  
  plot(datMeta$Dx, ylim=c(0,200),col = plot_cols, main="Subjects")
  A = anova(lm(as.numeric(as.factor(datMeta$Gender)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot(as.factor(datMeta$Gender) ~ datMeta$Dx, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$Age)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$Age) ~ datMeta$Dx, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$PMI_hrs)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$PMI_hrs) ~ datMeta$Dx, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$pH)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$pH) ~ datMeta$Dx, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$DLPFC_RNA_isolation_RIN)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$DLPFC_RNA_isolation_RIN) ~ datMeta$Dx, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="")
  A = chisq.test(datMeta$Batch,datMeta$Dx); p = A$p.value; plot(datMeta$Dx ~ as.factor(datMeta$Batch), col=plot_cols, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
  A = chisq.test(datMeta$Cluster,datMeta$Dx); p = A$p.value; plot(datMeta$Dx ~ as.factor(datMeta$Cluster), col=plot_cols, main=paste("Ancestry Cluster,\np=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$DLPFC_RNA_Sequencing_Mapped_Reads)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$DLPFC_RNA_Sequencing_Mapped_Reads) ~ datMeta$Dx, main=paste("Mapped Reads, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$DLPFC_RNA_Sequencing_Genes_Detected)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$DLPFC_RNA_Sequencing_Genes_Detected) ~ datMeta$Dx, main=paste("Genes Detected, p=", signif(p,2)), ylab="", xlab="")
  for(s in paste("seqPC",1:2,sep=""))   {
    A = anova(lm(as.numeric((datMeta[,s])) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta[,s]) ~ datMeta$Dx, main=paste(s, ", p=", signif(p,2)), ylab="", xlab="")
  }
  if(length(levels(factor(datMeta$Institution))) > 1) A = chisq.test(datMeta$Institution,datMeta$Dx); p = A$p.value; plot(datMeta$Dx ~ as.factor(datMeta$Institution), col=plot_cols, main=paste("Institution, p=", signif(p,2)), ylab="", xlab="")
}
if(plot_pdf) dev.off()




#-------Calculate differential expression:  linear model
sumstats = vector(mode="list", length=2); names(sumstats)= c("SCZ", "BD")
mod.scz=model.matrix(~Dx+Age+Gender+Institution+DLPFC_RNA_isolation_RIN+RIN2+PMI_hrs+seqPC1+seqPC2+EV.1+EV.2+EV.3+EV.4+EV.5,data=datMeta_scz)
mod.bd=model.matrix(~Dx+Age+Gender+DLPFC_RNA_isolation_RIN+RIN2+PMI_hrs+seqPC1+seqPC2+EV.1+EV.2+EV.3+EV.4+EV.5,data=datMeta_bd)
fit.scz = eBayes(lmFit(datExpr_scz, mod.scz), trend=T,robust=T)
fit.bd = eBayes(lmFit(datExpr_bd, mod.bd), trend=T, robust=T)
sumstats$SCZ = topTable(fit.scz, coef=2, number = Inf, sort.by = "none", confint = T)
sumstats$BD = topTable(fit.bd, coef=2, number = Inf, sort.by = "none", confint = T)

rnaseq.cmc = do.call("cbind", sumstats)
write.csv(file="./results/tables/RNAseq_SCZ_BD_CMC.csv", rnaseq.cmc)

cor.test(sumstats$SCZ$logFC, sumstats$BD$logFC,method="spearman")




#---qSVA------
degMatrix=read.csv("./raw_data/RNAseq_CMC/RNAseq_CMC_riboZero_degRegion_featureCounts.csv",row.names=1)
degMatrix=degMatrix[rowSums(degMatrix)>0,]
degMatrix.scz = degMatrix[,match(colnames(datExpr_scz), colnames(degMatrix))]
degMatrix.bd = degMatrix[,match(colnames(datExpr_bd), colnames(degMatrix))]

intervalLength=unlist(lapply(strsplit(rownames(degMatrix), "_"), function(x) { as.numeric(x[3]) - as.numeric(x[2]) } ))/1000
intervalLength = matrix(rep(intervalLength), ncol=ncol(degMatrix), nrow=nrow(degMatrix), byrow = FALSE)
 
normFactor.scz =  matrix(rep(datMeta_scz$DLPFC_RNA_Sequencing_Total_Reads/8e07),  ncol=ncol(degMatrix.scz), nrow=nrow(degMatrix.scz), byrow = TRUE)
normFactor.bd =  matrix(rep(datMeta_bd$DLPFC_RNA_Sequencing_Total_Reads/8e07),  ncol=ncol(degMatrix.bd), nrow=nrow(degMatrix.bd), byrow = TRUE)

degMatrix_fpk80m.scz = log2(1 + (degMatrix.scz /  normFactor.scz) / intervalLength)
degMatrix_fpk80m.bd = log2(1 + (degMatrix.bd /  normFactor.bd) / intervalLength)

degPca.scz = prcomp(t(degMatrix_fpk80m.scz))
k.scz = num.sv(degMatrix_fpk80m.scz, mod=rep(1,ncol(datExpr_scz)),B=100)  #Using 6 qSVs
degPca.bd = prcomp(t(degMatrix_fpk80m.bd))
k.bd = num.sv(degMatrix_fpk80m.bd, mod=rep(1,ncol(datExpr_bd)),B=100)  #Using 4 qSVs


sumstats.qsva = vector(mode="list", length=2); names(sumstats.qsva)= c("SCZ", "BD")
fit.scz = eBayes(lmFit(datExpr_scz, cbind(mod.scz, degPca.scz$x[,1:k.scz])))
sumstats.qsva$SCZ = topTable(fit.scz, coef=2, number=Inf,sort.by = "none", confint = T)
fit.bd = eBayes(lmFit(datExpr_bd, cbind(mod.bd, degPca.bd$x[,1:k.bd])))
sumstats.qsva$BD = topTable(fit.bd, coef=2, number=Inf,sort.by = "none", confint = T)
write.csv(file="./results/tables/RNAseq_SCZ_BD_CMC.qsva.csv", sumstats.qsva)
