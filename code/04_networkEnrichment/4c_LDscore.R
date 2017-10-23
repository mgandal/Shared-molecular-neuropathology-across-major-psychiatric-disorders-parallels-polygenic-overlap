#4c_LDscore.R

rm(list=ls()); options(stringsAsFactors = F)
setwd("~/Github//Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/")
library(ggplot2)

if(FALSE) {
  #RUN LDSC from command line:
  #--------------------------
  # ldsc=~/bin/ldsc/ldsc.py
  # 
  # ANNOT1=/hp_shares/mgandal/datasets/LDscore/1000G_Phase1/baseline/baseline.
  # ANNOT2=/hp_shares/mgandal/datasets/LDscore/CrossDisorderLDSC/CDmod/cd.mod.
  # OUTDIR=/hp_shares/mgandal/datasets/LDscore/CrossDisorderLDSC/
  # 
  # $ldsc --ref-ld-chr $ANNOT2,$ANNOT1 --w-ld-chr /hp_shares/mgandal/datasets/LDscore/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr /hp_shares/mgandal/datasets/LDscore/1000G_Phase1/1000G_frq/1000G.mac5eur. --print-coefficients \
  # --h2 /hp_shares/mgandal/datasets/GWAS/SCZ.PGC.2014/SCZ.PGC.2014.sumstats.gz \
  # --pop-prev 0.01 --samp-prev 0.43 \
  # --out /hp_shares/mgandal/datasets/LDscore/CrossDisorderLDSC/SCZ.PGC.2014
#   
#   $ldsc --ref-ld-chr $ANNOT2,$ANNOT1 --w-ld-chr /hp_shares/mgandal/datasets/LDscore/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr /hp_shares/mgandal/datasets/LDscore/1000G_Phase1/1000G_frq/1000G.mac5eur. --print-coefficients \
#   --h2 /hp_shares/mgandal/datasets/GWAS/BP.PGC.2012/BD.PGC.2012.sumstats.gz \
#   --pop-prev 0.01 --samp-prev 0.45 \
#   --out $OUTDIR/BD.PGC.2012
#   
#   $ldsc --ref-ld-chr $ANNOT2,$ANNOT1 --w-ld-chr /hp_shares/mgandal/datasets/LDscore/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr /hp_shares/mgandal/datasets/LDscore/1000G_Phase1/1000G_frq/1000G.mac5eur. --print-coefficients \
#   --h2 /hp_shares/mgandal/datasets/GWAS/ASD.PGC.2015/*.sumstats.gz \
#   --pop-prev 0.01 --samp-prev 0.5 \
#   --out $OUTDIR/ASD.PGC.2015
#   
#   $ldsc --ref-ld-chr $ANNOT2,$ANNOT1 --w-ld-chr /hp_shares/mgandal/datasets/LDscore/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr /hp_shares/mgandal/datasets/LDscore/1000G_Phase1/1000G_frq/1000G.mac5eur. --print-coefficients \
#   --h2 /hp_shares/mgandal/datasets/GWAS/MDD.PGC.2012/*.sumstats.gz \
#   --pop-prev 0.13 --samp-prev 0.5 \
#   --out $OUTDIR/MDD.PGC.2012
# #   
#   $ldsc --ref-ld-chr $ANNOT2,$ANNOT1 --w-ld-chr /hp_shares/mgandal/datasets/LDscore/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr /hp_shares/mgandal/datasets/LDscore/1000G_Phase1/1000G_frq/1000G.mac5eur. --print-coefficients \
#   --h2 /hp_shares/mgandal/datasets/GWAS/ETOH.AlcGen.2011/*.sumstats.gz \
#   --out $OUTDIR/AAD.AlcGen.2011
#   
#   $ldsc --ref-ld-chr $ANNOT2,$ANNOT1 --w-ld-chr /hp_shares/mgandal/datasets/LDscore/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr /hp_shares/mgandal/datasets/LDscore/1000G_Phase1/1000G_frq/1000G.mac5eur. --print-coefficients \
#   --h2 /hp_shares/mgandal/datasets/GWAS/IBD.2015/*.sumstats.gz \
#   --pop-prev 0.004 --samp-prev 0.35 \
#   --out $OUTDIR/IBD.2015
#   
}

d = dir(path="./results/tables/LDSC/", pattern=".results")

ldsc = data.frame()

for(i in 1:length(d)) {
  dat = read.table(paste("./results/tables//LDSC/", d[[i]], sep=""),head=T)
  idx = grep("cd_", dat$Category)
  dat = dat[idx,]
  dat$GWAS = gsub(".results", "", d[[i]],)
  log = unlist(strsplit(gsub("Total Liability scale h2: ", "", system(paste("grep Liability ./results/tables/LDSC/", dat$GWAS[[1]], ".log", sep=""),intern=T))," "))
                
  
  dat$TotalH2 = as.numeric(log[1])
  dat$TotalH2_sem = as.numeric(gsub("[)]","", gsub("[(]","", log[2])))
  ldsc = rbind(ldsc, dat)
  
}


ldsc$Module = gsub("cd_", "", gsub("L2_0", "", ldsc$Category))
ldsc$FDR = p.adjust(ldsc$Enrichment_p,method="fdr")
ldsc$p.symbol=""
ldsc$p.symbol[ldsc$FDR<0.05]= "*"
ldsc$p.symbol[ldsc$FDR<0.01]= "**"
ldsc$p.symbol[ldsc$FDR<0.001]= "***"
ldsc$p.symbol[ldsc$Enrichment<1] = ""

SNPh2 = ldsc[ldsc$Module=="black",]
Fig4C = ggplot(SNPh2, aes(x=GWAS, y=TotalH2)) + geom_bar(stat="identity") + 
  geom_errorbar(aes(ymin=(TotalH2 - TotalH2_sem), ymax=(TotalH2 + TotalH2_sem)), position=position_dodge(.9), size=0.3,width=.3, color="black") +
  labs(x="",y="SNP Heritability (h2)") + ggtitle("Total SNP-Based\nHeritability") + theme_classic() +
  theme(axis.text.x=element_text(angle=50, size=10, hjust=1),legend.position = "none")
ggsave(filename = "./results/figures/Manuscript/Fig4C.pdf",plot = Fig4C,width = 3,height=4)


ldsc = ldsc[ldsc$Module %in% c("green", "purple", "turquoise", "salmon", "yellow", "tan", "greenyellow", "blue"),]

write.csv(file="./results/tables/Manuscript/DataTableS5-LDSC.csv", ldsc)

Fig4D = ggplot(ldsc,aes(x=GWAS,y=Prop._h2,fill=Module,label=p.symbol)) +
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  geom_errorbar(aes(ymin=(Prop._h2 - Prop._h2_std_error), ymax=(Prop._h2 + Prop._h2_std_error)), position=position_dodge(.9), size=0.3,width=.3) +
  scale_fill_manual(values=sort(unique(ldsc$Module)))+ ggtitle("Partitioned Heritability") + 
  labs(x="",y="Proportion of Heritability (h2)") +theme_classic()+
  geom_text(color="red",size=4,aes(y=Prop._h2+ sign(Prop._h2)*Prop._h2_std_error + sign(Prop._h2)*.005), position=position_dodge(.9))  + 
  theme(axis.text.x=element_text(angle=50, size=10, hjust=1),legend.position = "none")
ggsave(filename = "./results/figures/Manuscript/Fig4D.pdf",plot = Fig4D,width = 8,height=4)
