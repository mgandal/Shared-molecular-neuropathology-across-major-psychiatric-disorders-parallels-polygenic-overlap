#5c_CompileRNAseqReplication

rm(list=ls())
library(reshape); library(ggplot2); library(gridExtra)

plot_pdf=T

pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  rho = cor(ds1,ds2,method="spearman")
  return(list(slope,intercept, rho))
}

setwd("~/Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap//")
array.asd = read.csv("./results/tables/Microarray_ASD_metaanalysis_092017.csv")
array.scz = read.csv("./results/tables/Microarray_SCZ_metaanalysis_092017.csv")
array.bd = read.csv("./results/tables/Microarray_BD_metaanalysis_092017.csv")

rnaseq.asd = read.csv("./results/tables/RNAseq_ASD_4region_sumstats.csv")
gvex = read.csv("./results/tables/RNAseq_SCZ_BD_GVEX.csv")
cmc = read.csv(w)

all_genes.array = intersect(array.asd$X, intersect(array.scz$X, array.bd$X))
all_genes.rnaseq = intersect(rnaseq.asd$X,intersect(gvex$X,cmc$X))
all_genes = intersect(all_genes.array, all_genes.rnaseq)

array.asd=array.asd[match(all_genes,array.asd$X),]
array.scz=array.scz[match(all_genes,array.scz$X),]
array.bd=array.bd[match(all_genes,array.bd$X),]
rnaseq.asd=rnaseq.asd[match(all_genes,rnaseq.asd$X),]
gvex=gvex[match(all_genes,gvex$X),]
cmc = cmc[match(all_genes, cmc$X),]

dat= data.frame(SCZ=gvex$SCZ.logFC, ASD=rnaseq.asd$All.logFC, BD=gvex$BD.logFC, dataset="BrainGVEX")
dat= rbind(dat, data.frame(SCZ=cmc$SCZ.logFC, ASD=rnaseq.asd$All.logFC, BD=cmc$BD.logFC, dataset="CommonMind"))
dat2= melt(dat,id=c(1,4))
dat2$value = as.numeric(dat2$value)
colnames(dat2)[3] = "comparison"
fit = pcreg(gvex$SCZ.logFC, rnaseq.asd$All.logFC); linreg = data.frame(dataset="BrainGVEX", comparison="ASD", slope=fit[[1]], intercept=fit[[2]], color="#F8766D", rho=fit[[3]])
fit = pcreg(gvex$SCZ.logFC, gvex$BD.logFC); linreg = rbind(linreg, data.frame(dataset="BrainGVEX", comparison="BD", slope=fit[[1]], intercept=fit[[2]], color="#00BA38", rho=fit[[3]]))
fit = pcreg(cmc$SCZ.logFC, rnaseq.asd$All.logFC); linreg = rbind(linreg,data.frame(dataset="CommonMind", comparison="ASD", slope=fit[[1]], intercept=fit[[2]], color="#F8766D",rho=fit[[3]]))
fit = pcreg(cmc$SCZ.logFC, cmc$BD.logFC); linreg = rbind(linreg,data.frame(dataset="CommonMind", comparison="BD", slope=fit[[1]], intercept=fit[[2]], color="#00BA38",rho=fit[[3]]))

TxSlope.rnaseq=  ggplot(dat2,aes(x=SCZ,y=value,color=comparison)) + geom_point(alpha=.5) + 
  geom_abline(slope=1, lty=2) + xlim(-1,1) + ylim(-1,1) + 
  geom_abline(data=linreg, aes(color=comparison, slope=slope, intercept=intercept)) +  facet_grid(.~dataset) +
  scale_colour_manual(values=c("#F8766D","#00BA38")) +
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) + ggtitle("RNAseq Transcriptome Severity") + geom_text(data=linreg, size=3, aes(x=1,y=-1, hjust=1, label=paste("slope=", signif(linreg$slope,2), "\nrho=", signif(linreg$rho,2))),position = position_fill(reverse = T))



#----Repeat with qSVA normalized data
rnaseq.asd.qsva = read.csv("./results/tables/RNAseq_ASD_4region_sumstats_qSVA.csv")
gvex.qsva = read.csv("./results/tables/RNAseq_SCZ_BD_GVEX.qsva.csv")
cmc.qsva = read.csv("./results/tables/RNAseq_SCZ_BD_CMC.qsva.csv")

rnaseq.asd.qsva=rnaseq.asd.qsva[match(all_genes,rnaseq.asd.qsva$X),]
gvex.qsva=gvex.qsva[match(all_genes,gvex.qsva$X),]
cmc.qsva = cmc.qsva[match(all_genes, cmc.qsva$X),]

datqsva= data.frame(SCZ=gvex.qsva$SCZ.logFC, ASD=rnaseq.asd.qsva$logFC, BD=gvex.qsva$BD.logFC, dataset="BrainGVEX")
datqsva= rbind(datqsva, data.frame(SCZ=cmc.qsva$SCZ.logFC, ASD=rnaseq.asd.qsva$logFC, BD=cmc.qsva$BD.logFC, dataset="CommonMind"))
datqsva2= melt(datqsva,id=c(1,4))
datqsva2$value = as.numeric(datqsva2$value)
colnames(datqsva2)[3] = "comparison"

fitqsva = pcreg(gvex.qsva$SCZ.logFC, rnaseq.asd.qsva$logFC); linreg.qsva = data.frame(dataset="BrainGVEX", comparison="ASD", slope=fitqsva[[1]], intercept=fitqsva[[2]], color="#F8766D", rho=fitqsva[[3]])
fitqsva = pcreg(gvex.qsva$SCZ.logFC, gvex.qsva$BD.logFC); linreg.qsva = rbind(linreg.qsva, data.frame(dataset="BrainGVEX", comparison="BD", slope=fitqsva[[1]], intercept=fitqsva[[2]], color="#00BA38", rho=fitqsva[[3]]))
fitqsva = pcreg(cmc.qsva$SCZ.logFC, rnaseq.asd.qsva$logFC); linreg.qsva = rbind(linreg.qsva,data.frame(dataset="CommonMind", comparison="ASD", slope=fitqsva[[1]], intercept=fitqsva[[2]], color="#F8766D",rho=fitqsva[[3]]))
fitqsva = pcreg(cmc.qsva$SCZ.logFC, cmc.qsva$BD.logFC); linreg.qsva = rbind(linreg.qsva,data.frame(dataset="CommonMind", comparison="BD", slope=fitqsva[[1]], intercept=fitqsva[[2]], color="#00BA38",rho=fitqsva[[3]]))

TxSlope.rnaseq.qsva=  ggplot(datqsva2,aes(x=SCZ,y=value,color=comparison)) + geom_point(alpha=.5) + 
  geom_abline(slope=1, lty=2) + xlim(-1,1) + ylim(-1,1) + 
  geom_abline(data=linreg.qsva, aes(color=comparison, slope=slope, intercept=intercept)) +  facet_grid(.~dataset) +
  scale_colour_manual(values=c("#F8766D","#00BA38")) +
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) + ggtitle("RNAseq Transcriptome Severity (qSVA)") + 
  geom_text(data=linreg.qsva, size=3, aes(x=1,y=-1, hjust=1, label=paste("slope=", signif(linreg.qsva$slope,2), "\nrho=", signif(linreg.qsva$rho,2))),position = position_fill(reverse = T))






# Individual Datasets
sig.asd = array.asd$fdr<.05
sig.scz = array.scz$fdr<.05
sig.bd = array.bd$fdr<.05

to_plot = rbind(data.frame(Microarray=array.scz$beta[sig.scz], RNAseq=gvex$SCZ.logFC[sig.scz], Group="SCZ-GVEX", Dx="SCZ", dataset="BrainGVEX"),
                data.frame(Microarray=array.asd$beta[sig.asd], RNAseq=rnaseq.asd$All.logFC[sig.asd], Group="ASD", Dx="ASD", dataset='ASD-pancortical'),
                data.frame(Microarray=array.bd$beta[sig.bd], RNAseq=gvex$BD.logFC[sig.bd], Group="BD-GVEX", Dx="BD", dataset="BrainGVEX"),
                data.frame(Microarray=array.bd$beta[sig.bd], RNAseq=cmc$BD.logFC[sig.bd], Group="BD-CMC", Dx="BD", dataset="CommonMind"),
                data.frame(Microarray=array.scz$beta[sig.scz], RNAseq=cmc$SCZ.logFC[sig.scz], Group="SCZ-CMC", Dx="SCZ", dataset="CommonMind"))
#to_plot$Group = factor(to_plot$Group, levels=c("ASD","SCZ-GVEX", "BD-GVEX", "SCZ-CMC","BD-CMC"))
to_plot$Dx = factor(to_plot$Dx, levels=c("ASD", 'SCZ', 'BD'))

c = cor.test(array.asd$beta[sig.asd], rnaseq.asd$All.logFC[sig.asd],method="spearman")
datLabel = data.frame(Dx="ASD", dataset='ASD-pancortical', rho=c$estimate, p=c$p.value)
c = cor.test(array.scz$beta[sig.scz], gvex$SCZ.logFC[sig.scz],method="spearman")  #0.85
datLabel = rbind(datLabel,data.frame(Dx="SCZ", dataset='BrainGVEX', rho=c$estimate, p=c$p.value))
c=cor.test(array.bd$beta[sig.bd], gvex$BD.logFC[sig.bd],method="spearman")   #0.78
datLabel = rbind(datLabel,data.frame(Dx="BD", dataset='BrainGVEX', rho=c$estimate, p=c$p.value))
c= cor.test(array.scz$beta[sig.scz], cmc$SCZ.logFC[sig.scz],method="spearman")   #0.65
datLabel = rbind(datLabel,data.frame(Dx="SCZ", dataset='CommonMind', rho=c$estimate, p=c$p.value))
c= cor.test(array.bd$beta[sig.bd], cmc$BD.logFC[sig.bd],method="spearman")   #0.56
datLabel = rbind(datLabel,data.frame(Dx="BD", dataset='CommonMind', rho=c$estimate, p=c$p.value))
datLabel$Dx = factor(datLabel$Dx, levels=c("ASD", "SCZ", "BD"))
datLabel$dataset = factor(datLabel$dataset)
  
RNAseq.rep=  ggplot(to_plot,aes(x=Microarray,y=RNAseq,color=dataset)) + geom_point(alpha=.5) + 
  geom_smooth(method="lm",fullrange=T) + geom_abline(slope=1, intercept = 0, lty=2) +
  xlab("Microarray log2FC") + ylab("RNAseq log2FC") + coord_fixed(ratio=1) +
  ggtitle(" RNAseq Replication")+ facet_grid(.~Dx) +  
  geom_text(data=datLabel, size=3, aes(x=1,y=-1,label=paste0("rho=", signif(rho,2))),position = position_fill(reverse = T)) 
            
if(plot_pdf) pdf(file="./results/figures/Manuscript//FigS7-RNAseq Replication.pdf", width=10,height=10)
grid.arrange(grobs=list(TxSlope.rnaseq, TxSlope.rnaseq.qsva, RNAseq.rep),ncol=1)
if(plot_pdf) dev.off()



#Inidividual Table
source("./code/00_scripts/fisher_overlap.R")
replicationTable=data.frame()
fisher = data.frame(t(ORA(all_genes[sign(array.asd$beta)==sign(rnaseq.asd$All.logFC) & rnaseq.asd$All.P.Value<0.05], all_genes[sig.asd], all_genes,all_genes)))

replicationTable = rbind(replicationTable, data.frame(Group="ASD", RNAseq="ASD-pancortical", Array.numDGE=sum(na.omit(sig.asd)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

fisher = data.frame(t(ORA(all_genes[sign(array.scz$beta)==sign(gvex$SCZ.logFC) & gvex$SCZ.P.Value<0.05], all_genes[sig.scz], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="SCZ", RNAseq="BrainGVEX", Array.numDGE=sum(na.omit(sig.scz)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

fisher = data.frame(t(ORA(all_genes[sign(array.bd$beta)==sign(gvex$BD.logFC) & gvex$BD.P.Value<0.05], all_genes[sig.bd], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="BD", RNAseq="BrainGVEX", Array.numDGE=sum(na.omit(sig.bd)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))


fisher = data.frame(t(ORA(all_genes[sign(array.scz$beta)==sign(cmc$SCZ.logFC) & cmc$SCZ.P.Value<0.05], all_genes[sig.scz], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="SCZ", RNAseq="CMC", Array.numDGE=sum(na.omit(sig.scz)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

fisher = data.frame(t(ORA(all_genes[sign(array.bd$beta)==sign(cmc$BD.logFC) & cmc$BD.P.Value<0.05], all_genes[sig.bd], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="BD", RNAseq="CMC", Array.numDGE=sum(na.omit(sig.bd)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

replicationTable

