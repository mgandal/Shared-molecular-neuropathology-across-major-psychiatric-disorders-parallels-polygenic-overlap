#4b_rareVariantEnrich.R

rm(list=ls()); options(stringsAsFactors = F)
setwd("~/Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap//")
source("./code/00_scripts/fisher_overlap.R")

library(ggplot2); library(reshape)
library(biomaRt)

#Load finalized network modules 
load("./working_data//NetworkAnalysis/finalizedNetwork_092017.RData")

#Load genes compiled from supplement of WES studies in SCZ, ASD, ID, and control/siblings
#Compiled from:
#"De Rubieus 2014" "deLigt"          "Fromer 2014"     "Girard 2011"     "Gulsuner 2013"   "Iossifov 2012"  
#"Iossifov 2014"   "Neale 2012"      "O'Roak 2012"     "Rauch 2012"      "Sanders 2012"    "Xu 2011"    
compilation = read.csv("./raw_data/RareVariants//WES - deNovo Compilation 120915.csv")
table(tolower(compilation$Type))

mut.nonsyn = (compilation$Type=="LOF") | grepl("missense", tolower(compilation$Type)) | grepl("shift", tolower(compilation$Type)) | grepl("nonsense", tolower(compilation$Type)) & !grepl("no-frame-shift", tolower(compilation$Type))
mut.silent = grepl("synonymous", tolower(compilation$Type)) | grepl("silent", tolower(compilation$Type))

gene.lists = vector(mode="list")
asd.idx = compilation$Population=="ASD"
scz.idx = (toupper(compilation$Population)=="SCHIZOPHRENIA")
control.idx = grepl("sibling",tolower(compilation$Population))

gene.lists[[1]] = unique(compilation$Gene[asd.idx & mut.nonsyn])
gene.lists = c(gene.lists, list(unique(compilation$Gene[scz.idx & mut.nonsyn])))
gene.lists = c(gene.lists, list(unique(compilation$Gene[control.idx & mut.nonsyn])))
gene.lists = c(gene.lists, list(unique(compilation$Gene[asd.idx & mut.silent])))
gene.lists = c(gene.lists, list(unique(compilation$Gene[scz.idx & mut.silent])))
gene.lists = c(gene.lists, list(unique(compilation$Gene[control.idx & mut.silent])))
gene.lists = c(gene.lists, list(read.delim("./raw_data/RareVariants/ASD.cnv.txt",head=F)$V1))
gene.lists = c(gene.lists, list(read.delim("./raw_data/RareVariants/SCZ.cnv.txt",head=F)$V1))


names(gene.lists)= c("ASD.nonsynonymous", "SCZ.nonsynonymous", "CTL.nonsynonymous", "ASD.silent", "SCZ.silent", "CTL.silent",  "ASD.cnv", "SCZ.cnv")
for(i in 1:length(gene.lists)) {
  l = length(gene.lists[[i]])
  names(gene.lists)[i] = paste(names(gene.lists)[i], "\n(n=",l,")",sep="")
}



table.p = matrix(NA,nrow=14,ncol=length(gene.lists))
rownames(table.p) = unique(colors);
colnames(table.p) = names(gene.lists) 
table.or = logit.p = logit.or = table.p

for(i in 1:ncol(table.or)) {
  for(j in 1:dim(table.or)[1]) {
    col = rownames(table.or)[j]
    module.genes = datProbes$external_gene_id[colors==col]
    
    binaryMat = as.data.frame(cbind(as.numeric(colors==col), datProbes$size_kb, datProbes$cds_length, 0))
    colnames(binaryMat) = c("Module", "GeneLen", "CDSlen", "GeneSet")
    idx = match(gene.lists[[i]], datProbes$external_gene_id)
    binaryMat$GeneSet[idx]=1
    
    if(grepl("cnv",tolower(colnames(table.or)[i]))) {
      glm.out <- glm(binaryMat$GeneSet~binaryMat$Module+ binaryMat$GeneLen, family=binomial())
    } else {
      glm.out <- glm(binaryMat$GeneSet~binaryMat$Module+ binaryMat$CDSlen, family=binomial())
    }
    summary(glm.out)
    
    logit.or[j,i] = exp(coefficients(glm.out)[2])  #Calculate odds ratio from logistic regression
    logit.p[j,i] = summary(glm.out)$coefficients[2,4]
  }
}

table.p.fdr = p.adjust(logit.p,method="fdr")
dim(table.p.fdr) = dim(logit.p)
dimnames(table.p.fdr) = dimnames(logit.p)

to_plot = as.data.frame(table.p.fdr)
to_plot[logit.or<1]=1
to_plot = -log10(to_plot)
to_plot$module = rownames(to_plot)
to_keep = rownames(to_plot) %in% c("tan","blue","yellow","purple","turquoise","green","greenyellow","salmon")

to_plot = melt(to_plot[to_keep,])

p=ggplot(to_plot,aes(x=variable,y=value,fill=module)) + #ggtitle("Compilation\nLogit: Module ~ GeneSet + gene_len") +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  scale_fill_manual(values=sort(unique(to_plot$module)))+
  geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +theme_classic()+
  theme(axis.text.x=element_text(angle=50, size=10, hjust=1))
p

ggsave("./results/figures/Manuscript/Fig4B.pdf",p, width = 6.43, height=4.25)

