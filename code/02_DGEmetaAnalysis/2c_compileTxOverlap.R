#2c_compileTxOverlap
rm(list=ls()); options(stringsAsFactors=F)
source("http://bioconductor.org/biocLite.R")
library(ggplot2); library(mada); library(reshape)
library(NMF); library(WGCNA)
rootdir = "~/Github//Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/"
setwd(rootdir)

PLOT_PDF=T

##L
asd_meta = read.csv("./results/tables//Microarray_ASD_metaanalysis_092017.csv", row.names=1)
scz_meta = read.csv("./results/tables/Microarray_SCZ_metaanalysis_092017.csv", row.names=1)
bd_meta = read.csv("./results/tables/Microarray_BD_metaanalysis_092017.csv", row.names=1)
mdd_meta = read.csv("./results/tables/Microarray_MDD_metaanalysis_092017.csv", row.names=1)
aad_meta = read.csv("./results/tables/Microarray_AAD_metaanalysis_092017.csv", row.names=1)
ibd_meta = read.csv("./results/tables/Microarray_IBD_metaanalysis_092017.csv", row.names=1)

all_genes = intersect(intersect(intersect(intersect(intersect(rownames(asd_meta), rownames(scz_meta)), rownames(bd_meta)), rownames(aad_meta)), rownames(mdd_meta)), rownames(aad_meta))
allmeta = matrix(NA,nrow=length(all_genes), 6)
allmeta[,1] = asd_meta$beta[match(all_genes, rownames(asd_meta))]
allmeta[,2] = scz_meta$beta[match(all_genes, rownames(scz_meta))]
allmeta[,3] = bd_meta$beta[match(all_genes, rownames(bd_meta))]
allmeta[,4] = mdd_meta$beta[match(all_genes, rownames(mdd_meta))]
allmeta[,5] = aad_meta$beta[match(all_genes, rownames(aad_meta))]
allmeta[,6] = ibd_meta$beta[match(all_genes, rownames(ibd_meta))]

colnames(allmeta) = c("ASD", "SCZ", "BD", "MDD", "AAD", "IBD")
rownames(allmeta) = all_genes
allmeta=  as.data.frame(allmeta)
cor(allmeta,use="pairwise.complete.obs",method="spearman")




#Fig2A_inset: ASD vs SCZ, ASD vs AAD
dat= data.frame(x=allmeta$ASD, y=allmeta$SCZ, Disease="SCZ")
dat = rbind(dat, data.frame(x=allmeta$ASD, y=allmeta$AAD, Disease="AAD"))

Fig2a_inset = ggplot(dat,aes(x,y,group=Disease,color=Disease)) +   
  xlim(c(-1.5,1.5)) + ylim(c(-1.5,1.5)) +
  geom_point(alpha=.5,size=2) +   
  geom_smooth(method="lm",fullrange=T) + 
  scale_colour_discrete(name="Disease") +
  labs(x="ASD log2FC", y="log2FC") +  
  geom_abline(intercept=0,slope=1,colour=scales::alpha("black",.5),lty=2)
Fig2a_inset

if(PLOT_PDF) ggsave(filename="./results/figures/Manuscript/Fig2A_inset.pdf",Fig2a_inset,width=4,height=2)

#Load null
null = data.frame(read.delim("./working_data//NullDistribution/null.txt",head=F))
null = sort(null$V1)
null = cbind(null, data.frame(prob=1-abs(2*seq(1:length(null))/length(null)-1)))
hist(null$null, 50, main="Null Distribution\n40,000 permutations", xlab="Spearman's Rho")

#Make Bargraph
comparisons = t(combn(seq(1,ncol(allmeta)),2))
barplot = data.frame(Mean = NA, SEM=NA, p.fdr=NA)
for (i in 1:dim(comparisons)[1]) {
  x = comparisons[i,1]
  y = comparisons[i,2]
  R = cor.test(allmeta[,x], allmeta[,y])
  rho =cor(allmeta[,x], allmeta[,y], method="spearman", use="pairwise.complete.obs")
  sem = (tanh(atanh(rho + 1.96/sqrt(nrow(allmeta)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(allmeta)-3))))/3.92
  
  barplot[i,] = c(rho, sem, R$p.value)
  rownames(barplot)[i] = paste(colnames(allmeta)[x],colnames(allmeta)[y],sep="-")
}
barplot$p.fdr = p.adjust(barplot$p.fdr,method="fdr")
barplot$p.bootstrap = null$prob[findInterval(barplot$Mean, null$null)]
barplot$p.symbol = ""
barplot$p.symbol[barplot$p.bootstrap<0.05] = "*"
barplot$p.symbol[barplot$p.bootstrap<0.01] = "**"
barplot$p.symbol[barplot$p.bootstrap<0.001] = "***"
barplot$Comparison = rownames(barplot)
barplot$modality="microarray"

Fig2_main = ggplot(barplot,aes(x = reorder(Comparison, -Mean), y=Mean, label=p.symbol)) + ylim(-0.17,1) +  
  geom_bar(stat="identity",fill="royalblue3",width=0.75) +
  geom_errorbar(aes(ymin=(Mean - SEM), ymax=(Mean + SEM)), position=position_dodge(width=0.8), width=0.25,size=0.25) +   
  #  ggtitle(title) +   
  theme(plot.title = element_text(size=20, face="bold", vjust=2)) +   
  labs(x="", y=expression(paste("Transcriptome correlation (", rho, ")", sep=""))) +     	
  theme(
    axis.text.x=element_text(angle=50, size=12, hjust=1),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.title = element_text(size=12), 
    plot.title = element_text(size=12, face="bold"),
    axis.title.x = element_text(size=10, vjust=-0.35, face="bold"),
    axis.title.y = element_text(size=12, vjust=0.5),
    plot.margin=unit(c(2,2,1,2),"mm")
  ) + geom_text(color="red",size=4,aes(y=Mean+ sign(Mean)*SEM + sign(Mean)*.02))
Fig2_main

if(PLOT_PDF) ggsave("./results/figures/Manuscript/Fig2A.pdf", Fig2_main,height=4, width=5)


## From Cross-Disorder Group of the Psychiatric Genomics Consortium et al. 
## Genetic relationship between five psychiatric disorders estimated from genome-wide SNPs. 
## Nat Genet 45, 984–994 (2013).
barplot["SCZ-BD", "genetic.correlation"] = 0.68; barplot["SCZ-BD", "genetic.correlation.SEM"] = 0.04
barplot["SCZ-MDD", "genetic.correlation"] = 0.43; barplot["SCZ-MDD", "genetic.correlation.SEM"] = 0.06
barplot["ASD-SCZ", "genetic.correlation"] = 0.16; barplot["ASD-SCZ", "genetic.correlation.SEM"] = 0.06
barplot["BD-MDD", "genetic.correlation"] = 0.47; barplot["BD-MDD", "genetic.correlation.SEM"] = 0.06
barplot["ASD-BD", "genetic.correlation"] = 0.04; barplot["ASD-BD", "genetic.correlation.SEM"] = 0.06
barplot["ASD-MDD", "genetic.correlation"] = 0.05; barplot["ASD-MDD", "genetic.correlation.SEM"] = 0.09
barplot["ASD-IBD", "genetic.correlation"] = -.07; barplot["ASD-IBD", "genetic.correlation.SEM"] = 0.06
barplot["SCZ-IBD", "genetic.correlation"] = -.01; barplot["SCZ-IBD", "genetic.correlation.SEM"] = 0.03
barplot["BD-IBD", "genetic.correlation"] = -.05; barplot["BD-IBD", "genetic.correlation.SEM"] = 0.04
barplot["MDD-IBD", "genetic.correlation"] = 0.02; barplot["MDD-IBD", "genetic.correlation.SEM"] = 0.05


barplot2 = barplot
barplot2$modality = "RNAseq"
barplot2$Mean=NA; barplot2$SEM=NA

rnaseq.asd = read.csv("./results/tables//RNAseq_ASD_4region_sumstats.csv")
rnaseq.asd$log2FC = rnaseq.asd$All.logFC
gvex = read.csv("./results/tables/RNAseq_SCZ_BD_GVEX.csv")
all_genes = intersect(rnaseq.asd$X,gvex$X)
rnaseq.scz = data.frame(X= gvex$X, log2FC=gvex$SCZ.logFC)
rnaseq.bd = data.frame(X= gvex$X, log2FC=gvex$BD.logFC)
rnaseq.asd = rnaseq.asd[match(all_genes,rnaseq.asd$X),]
rnaseq.scz = rnaseq.scz[match(all_genes,rnaseq.scz$X),]
rnaseq.bd = rnaseq.bd[match(all_genes,rnaseq.bd$X),]

barplot2 = barplot
barplot2$modality = "RNAseq"
barplot2$Mean=NA; barplot2$SEM=NA

c = cor.test(rnaseq.asd$All.logFC, rnaseq.scz$log2FC,use="pairwise.complete.obs", method="spearman")
n = sum(!is.na(rnaseq.scz$log2FC[match(rnaseq.asd$X,rnaseq.scz$X)]))
barplot2["ASD-SCZ", "Mean"] = c$estimate; 
barplot2["ASD-SCZ", "SEM"] = (tanh(atanh(c$estimate + 1.93/sqrt(n-3))) - tanh(atanh(c$estimate - 1.93/sqrt(n-3))))/3.92

c = cor.test(rnaseq.asd$log2FC, rnaseq.bd$log2FC,use="pairwise.complete.obs", method="spearman")
barplot2["ASD-BD", "Mean"] = c$estimate; 
barplot2["ASD-BD", "SEM"] = (tanh(atanh(c$estimate + 1.93/sqrt(n-3))) - tanh(atanh(c$estimate - 1.93/sqrt(n-3))))/3.92

c = cor.test(rnaseq.bd$log2FC, rnaseq.scz$log2FC, method="spearman")
n = nrow(rnaseq.bd)
barplot2["SCZ-BD", "Mean"] = c$estimate; 
barplot2["SCZ-BD", "SEM"] =(tanh(atanh(c$estimate + 1.93/sqrt(n-3))) - tanh(atanh(c$estimate - 1.93/sqrt(n-3))))/3.92
barplot = rbind(barplot,barplot2)



c = cor.test(barplot$Mean, barplot$genetic.correlation, method="spearman", use="pairwise.complete.obs")
CI.95 = CIrho(c$estimate, N=13)[2:3]
e = expression(paste(rho, "=0.79**"))
linearMod =  lm(formula = barplot$Mean ~ barplot$genetic.correlation)


plot2 = ggplot(barplot, aes(x=genetic.correlation,y=Mean,label=Comparison, color=modality)) + 
  xlim(c(-0.19,1)) + ylim(c(-.19,1)) +
  geom_point(size=3) +
  labs(x="Genetic (SNP-based) correlation", y=expression(paste("Transcriptome correlation (", rho, ")", sep=""))) +       
  geom_errorbar(aes(ymin=(Mean - SEM), ymax=(Mean + SEM),colour=modality), width=0,size=0.8) +
  geom_errorbarh(aes(xmin=(genetic.correlation - genetic.correlation.SEM), xmax=(genetic.correlation + genetic.correlation.SEM)), height=0,size=0.8) + 
  theme(
    axis.text.x=element_text(size=10, hjust=1),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.title = element_text(size=12), 
    plot.title = element_text(size=12, face="bold"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12, vjust=0.5),
    plot.margin=unit(c(2,2,16,2),"mm")
  ) + geom_text(size=4,aes(colour=modality),position = position_jitter(width = 0, height = .05)) + 
  stat_smooth(method="lm",alpha=0, fullrange=T) +
  geom_abline(intercept=linearMod$coefficients[1], slope=linearMod$coefficients[2],lty=2) + theme(legend.position="none") + annotate("text", label=as.character(e), x=Inf,y=Inf, vjust=1,hjust=1, parse=T)

plot2




if(PLOT_PDF) ggsave(filename="./results/figures/Manuscript/Fig2C.pdf",plot2,width=5,height=4)







##Calculate the slope of transcriptome overlap using principle components regression
pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  return(list(slope,intercept))
}

dat= data.frame(SCZ=allmeta$SCZ, ASD=allmeta$ASD, BD=allmeta$BD, MDD=allmeta$MDD)
dat2= melt(dat,id=1)
dat2$value = as.numeric(dat2$value)

fit.asd = pcreg(allmeta$SCZ, allmeta$ASD)
fit.bd = pcreg(allmeta$SCZ, allmeta$BD)
fit.mdd=pcreg(allmeta$SCZ, allmeta$MDD)

dat2$variable = as.character(dat2$variable)
dat2$variable = gsub("ASD", paste("ASD, slope=", signif(fit.asd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("BD", paste("BD, slope=", signif(fit.bd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("MDD", paste("MDD, slope=", signif(fit.mdd[[1]],2), sep=""), dat2$variable)


TxSlope.array=ggplot(dat2,aes(x=SCZ,y=value,color=variable)) + geom_point(alpha=.5) + 
  geom_abline(slope=1, lty=2) + xlim(-1,1) + ylim(-1,1) + 
  geom_abline(slope=fit.asd[[1]], intercept = fit.asd[[2]], color="#F8766D") + 
  geom_abline(slope=fit.bd[[1]], intercept = fit.bd[[2]], color="#00BA38") + 
  geom_abline(slope=fit.mdd[[1]], intercept = fit.mdd[[2]], color="#619CFF") +
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) 

if(PLOT_PDF) ggsave(TxSlope.array, filename="./results/figures/Manuscript/Fig2B.pdf",width=5,height=5)

TxSlope.array


dat= data.frame(SCZ=rnaseq.scz$log2FC, ASD=rnaseq.asd$log2FC, BD=rnaseq.bd$log2FC)
dat2= melt(dat,id=1)
dat2$value = as.numeric(dat2$value)

fit.asd = pcreg(rnaseq.scz$log2FC, rnaseq.asd$log2FC)
fit.bd = pcreg(rnaseq.scz$log2FC, rnaseq.bd$log2FC)


dat2$variable = as.character(dat2$variable)
dat2$variable = gsub("ASD", paste("ASD, slope=", signif(fit.asd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("BD", paste("BD, slope=", signif(fit.bd[[1]],2), sep=""), dat2$variable)

TxSlope.rnaseq=ggplot(dat2,aes(x=SCZ,y=value,color=variable)) + geom_point(alpha=.5) + 
  geom_abline(slope=1, lty=2) + xlim(-1,1) + ylim(-1,1) + 
  geom_abline(slope=fit.asd[[1]], intercept = fit.asd[[2]], color="#F8766D") + 
  geom_abline(slope=fit.bd[[1]], intercept = fit.bd[[2]], color="#00BA38") +
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) 






##Plot dendrogram fo the top Genes
all_genes = intersect(intersect(intersect(intersect(rownames(asd_meta), rownames(scz_meta)), rownames(bd_meta)), rownames(mdd_meta)), rownames(aad_meta))
gene.symbols= asd_meta$symbol[match(all_genes, rownames(asd_meta))]

all_beta = matrix(NA,nrow=length(all_genes), 5)
all_beta[,1] = asd_meta$beta[match(all_genes, rownames(asd_meta))]
all_beta[,2] = scz_meta$beta[match(all_genes, rownames(scz_meta))]
all_beta[,3] = bd_meta$beta[match(all_genes, rownames(bd_meta))]
all_beta[,4] = mdd_meta$beta[match(all_genes, rownames(mdd_meta))]
all_beta[,5] = aad_meta$beta[match(all_genes, rownames(aad_meta))]
all_beta=as.data.frame(all_beta)

all_pvals = matrix(NA,nrow=length(all_genes), 5)
all_pvals[,1] = asd_meta$fdr[match(all_genes, rownames(asd_meta))]
all_pvals[,2] = scz_meta$fdr[match(all_genes, rownames(scz_meta))]
all_pvals[,3] = bd_meta$fdr[match(all_genes, rownames(bd_meta))]
all_pvals[,4] = mdd_meta$fdr[match(all_genes, rownames(mdd_meta))]
all_pvals[,5] = aad_meta$fdr[match(all_genes, rownames(aad_meta))]
all_pvals=  as.data.frame(all_pvals)

colnames(all_beta) = colnames(all_pvals) = c("ASD", "SCZ", "BD", "MDD", "AAD")


rowsums=rowSums(all_beta)
idx = order(rowsums,decreasing = T)[1:25]
idx = c(idx, order(rowsums)[1:25])
gene.symbols[idx]
mat.plot = as.matrix(-sign(all_beta[idx,]) * log10(all_pvals[idx,]))

rownames(mat.plot) = gene.symbols[idx]
colnames(mat.plot) = c("ASD", "SCZ", "BD", "MDD", "AAD")

textMat = signif(all_beta,2)
#textMat[all_pvals.fdr>0.05] = ''
textMat = textMat[idx,]







pdf("./results/figures/Manuscript/FigS2.pdf", width=5,height=8)
aheatmap(mat.plot,color=blueWhiteRed(100)[90:0],cexRow=.7, cexCol=1,fontsize=8,border_color="grey",scale="none",treeheight = 10, txt=textMat,Colv=F)
dev.off()


## Compile Table S1 - Transcriptome Signatures for Each Disease
union_genes = union(union(union(union(union(rownames(asd_meta), rownames(scz_meta)), rownames(bd_meta)), rownames(ibd_meta)), rownames(mdd_meta)), rownames(aad_meta))
all_disorders = vector(mode="list",length=6)
all_disorders[[1]] = asd_meta
all_disorders[[2]] = scz_meta
all_disorders[[3]] = bd_meta
all_disorders[[4]] = mdd_meta
all_disorders[[5]] = aad_meta
all_disorders[[6]] = ibd_meta

tableS1 = as.data.frame(matrix(NA,nrow=length(union_genes), 18))
for(i in 1:6) {
  all_disorders[[i]] = all_disorders[[i]][match(union_genes, rownames(all_disorders[[i]])),]
  tableS1[,3*i-2] = all_disorders[[i]]$beta
  tableS1[,3*i-1] = all_disorders[[i]]$p
  tableS1[,3*i] = all_disorders[[i]]$fdr
}

i = 1
for(dx in c("ASD", "SCZ", "BD", "MDD", "AAD", "IBD")) {
  for(var in c("beta_log2FC", "P.value", "FDR")) {
    colnames(tableS1)[i] = paste(dx,var,sep=".")
    i=i+1
  }
}


library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",  host="feb2014.archive.ensembl.org") 
#f=listFilters(ensembl); a =listAttributes(ensembl)
bm = getBM(attributes = c("ensembl_gene_id", "external_gene_id", "hgnc_symbol", "entrezgene", "chromosome_name", "start_position", "end_position"),
          filters="ensembl_gene_id", values=union_genes, mart=ensembl)
bm = bm[match(union_genes, bm$ensembl_gene_id),]
tableS1 = cbind(bm, tableS1)
write.csv(file="./results/tables/Manuscript/TableS1 - Microarray Meta-Analyses.csv", tableS1)



#Plot SCZ-BD and SCZ-ASD Tx Overlap Colored by Cell-Type
load("./working_data//NetworkAnalysis/finalizedNetwork_092017.RData")
dat = data.frame(allmeta); dat$Colors = factor(colors[match(rownames(allmeta),datProbes$ensembl_gene_id)]); 
dat$Cell[dat$Colors == "yellow"]= "Astrocyte"
dat$Cell[dat$Colors == "greenyellow"]= "Microglia"
dat$Cell[dat$Colors == "brown"]= "Oligodendrocyte"
dat$Cell[dat$Colors == "tan"]= "Endothelial"
dat$Cell[dat$Colors %in% c("green", "purple", "salmon", "turquoise")]= "Neuron"

dat2 = melt(dat[,c(1:3)],id=2)
dat2$CellType = c(dat$Cell,dat$Cell)

FigS15=ggplot(dat2,aes(x=SCZ,y=value,color=CellType)) +   
  xlim(c(-.5,.5))  + 
  geom_point(alpha=.5,size=2) +   
  scale_colour_discrete(name="CellType") +
  labs(x="SCZ log2FC", y="Diseae2 log2FC") +  
  facet_wrap(~variable,scales="free")+
  geom_hline(yintercept = 0, lty=2) + geom_vline(xintercept = 0, lty=2) + 
  theme(strip.text.x = element_text(size = 12))
  
ggsave(file="./results/figures/Manuscript/FigS15.pdf", FigS15, width=8, height=4)
