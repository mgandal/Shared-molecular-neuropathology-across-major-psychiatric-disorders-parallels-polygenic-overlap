#3e_network_moduleTrait.R

rm(list=ls()); options(stringsAsFactors = F)
#install.packages("pSI"); library(pSI); library(pSI.data)
library(WGCNA);library(ggplot2); 
library(reshape); library(nlme)

PLOT_PDF = T

rootdir = "~//Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap//"; setwd(rootdir)
load("./working_data//NetworkAnalysis//finalizedNetwork_092017.RData")

all_colors = unique(colors)
all_colors = all_colors[!grepl("grey",all_colors)]
all_genes = colors
names(all_genes) = rownames(datExpr)

#Step 1 - Calculate module-trait P values, beta, and SEM
moduleTraitP = matrix(NA,nrow=length(all_colors),ncol=12)
colnames(moduleTraitP) = c("ANOVA", "ASD", "BD", "ETOH", "MDD", "SCZ", "Sex", "Age", "PMI", "pH", "RIN", "RNAdeg")
rownames(moduleTraitP) = all_colors
moduleTraitB = moduleTraitSE = moduleTraitP

for (m in all_colors) {
  me_name = paste("ME", m, sep="")
  me = MEs$eigengenes[[me_name]]
  i = which(m == rownames(moduleTraitP))
  
  mixedmodel = lme(me ~ Group + Age + Sex, data = datMeta, random = ~1|Subject)
  moduleTraitP[i,"ANOVA"]=anova(mixedmodel)["Group","p-value"]
  mixedmodel = summary(mixedmodel)$tTable
  for(var in c("ASD", "ETOH", "BD", "SCZ", "MDD", "Age", "Sex")) {
    moduleTraitP[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 5]
    moduleTraitB[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 1]
    moduleTraitSE[i,var] = mixedmodel[grep(var, rownames(mixedmodel)), 2]
  }
  
  for(cov in c("Sex", "Age", "PMI", "pH", "RIN", "RNAdeg")) {
    mod = summary(lm(me ~ datMeta[,cov]))
    moduleTraitP[i,cov] = mod$coefficients[2, 4]
    moduleTraitB[i,cov] = mod$coefficients[2, 1]
    moduleTraitSE[i,cov] = mod$coefficients[2, 2]
  }
}

moduleTraitP.fdr = p.adjust(moduleTraitP, "fdr")
dim(moduleTraitP.fdr) = dim(moduleTraitP);dimnames(moduleTraitP.fdr) = dimnames(moduleTraitP);

out = cbind(moduleTraitB, moduleTraitP, moduleTraitP.fdr)
colnames(out)[1:12] = paste("Beta.", colnames(out)[1:12], sep="")
colnames(out)[13:24] = paste("P.", colnames(out)[13:24], sep="")
colnames(out)[25:36] = paste("FDR.", colnames(out)[25:36], sep="")
write.csv(file= "./results/tables/Manuscript/TableS3-modTraitB.csv",out)


row_idx = which(rownames(moduleTraitB) %in% c("green", "turquoise", "salmon", "purple","greenyellow", "tan","yellow", "blue"))
col_idx = c(2:6)


bpdata = melt(moduleTraitB[row_idx,col_idx])
semdata = melt(moduleTraitSE[row_idx,col_idx])
pdata = melt(moduleTraitP.fdr[row_idx,col_idx])
bpdata$sem = semdata$value
bpdata$p = pdata$value
bpdata$p.symbol = ""
bpdata$p.symbol[bpdata$p<0.1] = "#"
bpdata$p.symbol[bpdata$p<0.05] = "*"
bpdata$p.symbol[bpdata$p<0.01] = "**"
bpdata$p.symbol[bpdata$p<0.001] = "***"

bpdata$X1 = factor(bpdata$X1, levels=c("turquoise", "green", "purple", "salmon","yellow", "greenyellow","tan", "blue"))

Fig3C.byModule=ggplot(bpdata, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol))+ facet_wrap(~X1,ncol=1, scales="free") + 
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(bpdata$X1)) + 
  labs(y="beta", x="") + 
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(angle=45, hjust=1))

FigS21.SEMbyModule=ggplot(bpdata, aes(x=X2, y=as.numeric(sem),fill=X1,group=X1, label=p.symbol))+ facet_wrap(~X1,ncol=4, scales="free_x") + 
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  theme_minimal() + scale_fill_manual(name="Group",values=levels(bpdata$X1)) + 
  labs(y="Std Error of Beta", x="") + 
  #geom_text(color="red",size=4,aes(y=sem*1.1), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(angle=45, hjust=1))
FigS21.SEMbyModule

Fig3C.byDisease= ggplot(bpdata, aes(x=X1, y=value,fill=X1,group=X1, label=p.symbol))+ facet_wrap(~X2,ncol=1) + 
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(bpdata$X1)) + 
  labs(y="beta", x="") + 
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(angle=45, hjust=1))

if(PLOT_PDF) ggsave(Fig3C.byDisease, filename="./results/figures/Manuscript/Fig3C.pdf",width = 2.5,height=10)


dat= bpdata[bpdata$X1 %in% c("green","purple", "salmon", "turquoise"),]
neuron=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Neuron Modules") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(neuron, filename="./figures/WGCNA/Modules.Neuron.pdf", width=6,height=3)


dat= bpdata[bpdata$X1 %in% c("yellow"),]
astro=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Astrocyte Module") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(astro, filename="./figures/WGCNA/Modules.Astro.pdf", width=4,height=3)

dat= bpdata[bpdata$X1 %in% c("greenyellow"),]
microglia=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Microglial Module") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(microglia, filename="./figures/WGCNA/Modules.Microglia.pdf", width=4,height=3)

dat= bpdata[bpdata$X1 %in% c("tan"),]
endo=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Endothelial Module") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(endo, filename="./figures/WGCNA/Modules.Endo.pdf", width=4,height=3)

dat= bpdata[bpdata$X1 %in% c("blue"),]
bluemod=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Blue Module") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(bluemod, filename="./figures/WGCNA/Modules.Blue.pdf", width=4,height=3)

