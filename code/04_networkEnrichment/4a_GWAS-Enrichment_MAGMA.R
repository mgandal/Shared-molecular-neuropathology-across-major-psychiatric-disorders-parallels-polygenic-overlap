#4A_GWAS-Enrichment.R
## This script uses MAGMA to calculate gene-level significance values from 
## GWAS summary statistics and then assesses module enrichment

rm(list=ls()); options(stringsAsFactors = F)
library(WGCNA); library(biomaRt); library(reshape); library(ggplot2)
setwd("~/Github/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/")


load("working_data/NetworkAnalysis/finalizedNetwork_092017.RData")
MEs = moduleEigengenes(t(datExpr), colors,softPower = wgcna_parameters$powers)
modulekME = signedKME(t(datExpr),MEs$eigengenes,corFnc = "bicor")


##Run MAGMA
if(FALSE) {
  #Download gencode v19 annotation: gencode.v19.annotation.gtf
  #ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
  gencode = read.delim("./gencode.v19.annotation.gtf", skip=5, head=F,stringsAsFactors = F)
  gencode = gencode[gencode$V3=="gene",]
  gencode$ensg= substr(gencode$V9, 9,23)
  gencode$chr = gsub("chr", "", gencode$V1)
  out = data.frame(ENSG=gencode$ensg, CHR=gencode$chr, START=gencode$V4, STOP=gencode$V5)
  write.table(out, file="gencode.v19.genes.out", quote=F, row.names = F, col.names = F)
  
  # Run MAGMA 
  # #!/bin/bash
  # GENES=/hp_shares/mgandal/datasets/refGenome/GRCh37/gencode.v19.genes.out
  # SNPS=/hp_shares/mgandal/datasets/MAGMA/Ref/g1000.bim
  # MAGMA=~/bin/magma
  # 
  # $MAGMA --annotate --gene-loc $GENES --snp-loc $SNPS --out /hp_shares/mgandal/datasets/MAGMA/Ref/hg19_gencodeV19_magma
  # 
  # ANNOTATION=/hp_shares/mgandal/datasets/MAGMA/Ref/hg19_gencodeV19_magma.genes.annot
  # BFILE=/hp_shares/mgandal/datasets/MAGMA/Ref/g1000
  # CDMODS=/hp_shares/mgandal/datasets/MAGMA/CrossDisorder/CDmods_ensg_ForMagma.txt
  # OUTDIR=/hp_shares/mgandal/datasets/MAGMA/CrossDisorder
  # 
  # $MAGMA --bfile $BFILE --gene-annot $ANNOTATION --set-annot $CDMODS col=1,2\
  #   --pval /hp_shares/mgandal/datasets/GWAS/ASD.PGC.2015/ASD.PGC.2015_euro.all.25Mar2015_hg19.txt N=10610 \
  #   --out $OUTDIR/ASD.PGC.2015
  # 
  # $MAGMA --bfile $BFILE --gene-annot $ANNOTATION --set-annot $CDMODS col=1,2 \
  # --pval /hp_shares/mgandal/datasets/GWAS/SCZ.PGC.2014/scz2.snp.results.txt N=82315 \
  # --out $OUTDIR/SCZ.PGC.2014
  # 
  # $MAGMA --bfile $BFILE --gene-annot $ANNOTATION --set-annot $CDMODS col=1,2 \
  # --pval /hp_shares/mgandal/datasets/GWAS/BP.PGC.2012/BP.PGC.2012_pgc.bip.full.2012-04_hg18.txt N=63766 \
  # --out $OUTDIR/BD.PGC.2012
  # 
  # $MAGMA --bfile $BFILE --gene-annot $ANNOTATION --set-annot $CDMODS col=1,2 \
  # --pval /hp_shares/mgandal/datasets/GWAS/MDD.PGC.2012/MDD.PGC.2012_full.2012-04_hg18.txt N=18759 \
  # --out $OUTDIR/MDD.PGC.2012
  # 
  # 
  # $MAGMA --bfile $BFILE --gene-annot $ANNOTATION --set-annot $CDMODS col=1,2 \
  # --pval /hp_shares/mgandal/datasets/GWAS/ETOH.AlcGen.2011/ETOH.gwas.txt N=23347 \
  # --out $OUTDIR/AAD.AlcGen.2011
  # 
  # $MAGMA --bfile $BFILE --gene-annot $ANNOTATION --set-annot $CDMODS col=1,2 \
  # --pval /hp_shares/mgandal/datasets/GWAS/IBD.2015/EUR.IBD.gwas_info03_filtered.assoc N=23347 \
  # --out $OUTDIR/IBD.IDDRC.2015
}


d = dir(path="./results/tables/MAGMA/", pattern="[.]genes.out")
gwas = vector(mode="list")
for(i in 1:length(d)) {
  gwas[[i]] = read.table(paste("./results/tables/MAGMA/", d[[i]],sep=""),head=T)
}
names(gwas) = gsub(".genes.out", "", d)

#Match and replace with ENSG
genes = union(gwas[[1]]$GENE, gwas[[2]]$GENE)
for(i in 3:length(gwas)) genes = union(genes, gwas[[i]]$GENE)

for(i in 1:length(gwas)){
  idx = match(gwas[[i]]$GENE, datProbes$ensembl_gene_id)
  gwas[[i]]$gene_name = datProbes$external_gene_id[idx]
}


#Calculate spearman's correlation between gene module membership and GWAS gene significance
table.kme.cor.p =table.kme.cor.r<- matrix(NA,nrow=length(unique(colors)),ncol=length(gwas))
rownames(table.kme.cor.r) = rownames(table.kme.cor.p) = unique(colors)
colnames(table.kme.cor.r) = colnames(table.kme.cor.p) = names(gwas) 

for(m in unique(colors)) {
  for(i in 1:length(gwas)) {
    col = paste("kME", m, sep="")
    genes = intersect(gwas[[i]]$GENE, rownames(datExpr))
    x  = -log10(gwas[[i]]$P[match(genes, gwas[[i]]$GENE)])
    y = modulekME[match(genes,rownames(modulekME)), col]
    cor = cor.test(x,y,method="spearman")
    table.kme.cor.r[m,i] = cor$estimate
    table.kme.cor.p[m,i] = cor$p.value
  }
}

table.kme.cor.p.fdr = p.adjust(table.kme.cor.p, method="fdr")
dim(table.kme.cor.p.fdr) = dim(table.kme.cor.p);  dimnames(table.kme.cor.p.fdr) = dimnames(table.kme.cor.p)

d = -log10(table.kme.cor.p.fdr) * sign(table.kme.cor.r) 
labeledHeatmap(d,textMatrix = signif(table.kme.cor.r,1), xLabels = colnames(d), yLabels = rownames(d),invertColors = T, colors = blueWhiteRed(1000), main="GWAS - kME correlation", cex.text = 0.6)

dat = as.data.frame(table.kme.cor.p.fdr*sign(table.kme.cor.r))[c("tan","blue","yellow","purple","turquoise","green","greenyellow","salmon"),]
dat[dat<0]=1 #Only look for positive enrichment
dat = -log10(dat)
dat$module = rownames(dat)
dat2 = melt(dat)
dat2$variable=as.character(dat2$variable)

p=ggplot(melt(dat),aes(x=variable,y=value,fill=module)) + 
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  scale_fill_manual(values=sort(unique(dat2$module)))+ theme_classic() +
  geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
  theme(axis.text.x=element_text(angle=50, size=10, hjust=1))
p

ggsave(p, filename = "./results/figures/Manuscript/Fig4A.pdf",width=6,height=4)
 