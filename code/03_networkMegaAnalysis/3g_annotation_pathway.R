#3g_Annotation_Pathway.R

rm(list=ls()); options(stringsAsFactors = F)
library(pSI); library(gProfileR); library(gplots); library(biomaRt); library(WGCNA)

load("./working_data/NetworkAnalysis/finalizedNetwork_092017.RData")

eigmat = MEs$eigengenes; colnames(eigmat) = gsub("ME","",colnames(eigmat))
kME = signedKME(t(datExpr), MEs$eigengenes); colnames(kME) = gsub("kME", "", colnames(kME))


#Write module genes for GO elite
for(c in unique(colors)) {
  mod = datProbes$ensembl_gene_id[colors==c]
  write.table(file = paste("./results/tables/GeneSets/modules/input/", c, ".txt", sep=""), data.frame(Gene=mod, Id="En"), row.names = F, quote=F)
}
write.table(file="./results/tables/GeneSets/modules/background/background.txt", data.frame(Gene=datProbes$ensembl_gene_id, Id="En"), row.names=F, quote=F)



## Calculate GO enrichment of Disease-associated modules
resultsGO = data.frame()
for (m in unique(colors)) {
  me_name = paste("ME", m, sep="")
  i = which(m == unique(colors))  
  
  ## GO Enrichment
  query = rownames(datExpr)[colors==m]
  go = gprofiler(query, organism="hsapiens", ordered_query = F, significant = T, exclude_iea = F, region_query = F,max_p_value = 1,correction_method = "fdr", custom_bg = rownames(datExpr),
                 max_set_size = 2000,hier_filtering = "moderate", domain_size = "annotated", 
                 numeric_ns = "", include_graph = F,src_filter = c("GO", "KEGG"))
  go = go[order(go$p.value),]
  
  num_to_plot=10
  par(oma=c(0,10,0,0))
  bp = barplot(-log10(as.numeric(go$p.value[num_to_plot:1])), main=paste("\n\n\n\n\n",m, "- GO, KEGG"), horiz=T, yaxt='n', col=m, xlab='-log10(FDR)\n',cex.main=0.7, cex.axis = .7)
  axis(2,at=bp,labels=go$term.name[num_to_plot:1],tick=FALSE,las=2,cex.axis=.8);
  abline(v=-log10(0.01), col="red", lwd=2,lty=2)
  
  resultsGO = rbind(resultsGO, cbind(m, go))
}
write.csv(file="./results/tables/GeneSets/gProfiler_resultsGO.092017.csv", resultsGO)


resultsGOsub2 = data.frame()
for(c in sort(c("yellow", "green", "tan", "turquoise", "salmon", "purple", "greenyellow", "blue"))) {
  idx = which(resultsGO$m==c) 
  resultsGOsub2 = rbind(resultsGO[idx[2:1],],resultsGOsub2)  
}


pdf("./results/figures/Manuscript/Fig3E-Module-GO.pdf",width = 4,height=3)
par(oma=c(0,6,0,0), mar=c(5,4,2,2))
bp = barplot(-log10(resultsGOsub2$p.value), main="Pathway Enrichment", horiz=T, yaxt='n', col=resultsGOsub2$m, xlab='-log10(p)',cex.main=0.7, cex.axis = .7, border=NA)
axis(2,at=bp,labels=resultsGOsub2$term.name,tick=FALSE,las=2,cex.axis=0.6);
abline(v=-log10(0.01), col="black", lwd=2,lty=2)
dev.off()

pdf("./results/figures/WGCNA/Modules.Neuron_GO.pdf",width = 4,height=2.5)
resultsGOsub3 = resultsGOsub2[resultsGOsub2$m %in% c("green", "turquoise", "salmon", "purple"),]
par(oma=c(0,6,0,0), mar=c(5,4,2,2))
bp = barplot(-log10(resultsGOsub3$p.value),  main="", horiz=T, yaxt='n', col=resultsGOsub3$m, xlab='-log10(FDR)',cex.main=0.5, cex.axis = .7, cex.names = .2, border="black")
axis(2,at=bp,labels=resultsGOsub3$term.name,tick=FALSE,las=2,cex.axis=0.6);
abline(v=-log10(0.01), col="black", lwd=2,lty=2)
dev.off()

pdf("./results/figures/WGCNA/Modules.GO.pdf",width = 4,height=2)
for(m in c("blue", "yellow", "greenyellow", "tan")) {
  resultsGOsub3 = resultsGOsub2[resultsGOsub2$m %in% m,]
  par(oma=c(0,6,0,0), mar=c(5,4,2,2))
  bp = barplot(-log10(resultsGOsub3$p.value),  main="", horiz=T, yaxt='n', col=resultsGOsub3$m, xlab='-log10(FDR)',cex.main=0.5, cex.axis = .7, cex.names = .2, border="black")
  axis(2,at=bp,labels=resultsGOsub3$term.name,tick=FALSE,las=2,cex.axis=0.6);
  abline(v=-log10(0.01), col="black", lwd=2,lty=2)
}
dev.off()

##Download human transcriptome data from:
##Zhang, Y. et al. Purification and Characterization of Progenitor and Mature Human Astrocytes Reveals Transcriptional and Functional Differences with Mouse. Neuron 89, 37–53 (2016).
##Supplemental Table 3, "Human data only" sheet
##---
if(FALSE) {
  zhang.datExpr = read.csv("./raw_data/Annotation/zhang_GSE73721_mmc3-human.csv",skip=3,nrow=23223,head=F) 
  zhang.datMeta = data.frame(row.names=1:41,t(read.csv("./raw_data/Annotation/zhang_GSE73721_mmc3-human.csv",nrow=3,head=F)[,-1]))
  zhang.datMeta$CellType = NA
  zhang.datProbes = data.frame(symbol=zhang.datExpr$V1)
  zhang.datExpr = zhang.datExpr[,-1]
  colnames(zhang.datMeta) = c("X1", "Age", "Gender","CellType")
  zhang.datMeta$CellType[15:26] = "Astrocyte"
  zhang.datMeta$CellType[27] = "Neuron"
  zhang.datMeta$CellType[28:32] = "Oligo"
  zhang.datMeta$CellType[33:35] = "Microglia"
  zhang.datMeta$CellType[36:37] = "Endothelial"
  zhang.datMeta$CellType[38:41] = "WholeCortex"
  
  zhang.datExpr2 = data.frame(matrix(NA, nrow=nrow(zhang.datExpr), ncol=5)); colnames(zhang.datExpr2)=  c("Neuron", "Astrocyte", "Oligo", "Microglia","Endothelial")
  zhang.datExpr2$Neuron = zhang.datExpr[,which(zhang.datMeta$CellType=="Neuron")]
  for(cell in colnames(zhang.datExpr2)[2:5]) {
    zhang.datExpr2[,cell] = apply(zhang.datExpr[,which(zhang.datMeta$CellType==cell)],1,mean)  
  }
  
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  bm = getBM(attributes = c("ensembl_gene_id", "external_gene_id"), filters = "hgnc_symbol", values = zhang.datProbes$symbol, mart=ensembl)
  
  zhang.datProbes = data.frame(symbol=zhang.datProbes$symbol, ensg=bm$ensembl_gene_id[match(zhang.datProbes$symbol, bm$external_gene_id)])
  cr =collapseRows(zhang.datExpr2, rowGroup = zhang.datProbes$ensg, rowID=1:nrow(zhang.datExpr2))
  zhang.datExpr = cr$datETcollapsed
  rm(zhang.datExpr2, zhang.datProbes)
  save(file=".")
}


# Calculate Cell-Type specificity of modules
#Zhang using pSI
zhang.datExpr = read.csv("./raw_data/Annotations//datExpr.zhangHuman.avgForPSI.csv",row.names=1)[,-1]
set.seed(100)
pSI.output = specificity.index(pSI.in=zhang.datExpr,bts=100,p_max=.1, e_min=0.3); 
pSI.count(pSI.output)

cell.p.zhang = matrix(NA, 14,5);  rownames(cell.p.zhang) = unique(colors)
colnames(cell.p.zhang) = colnames(pSI.output)


for(mod in rownames(cell.p.zhang)) {
    f = fisher.iteration(pSI.output, rownames(datExpr)[colors==mod],p.adjust = F)
    cell.p.zhang[mod,] = f$`0.05 - nominal`
}

cell.p.zhang.fdr = p.adjust(cell.p.zhang,"fdr")
dim(cell.p.zhang.fdr) = dim(cell.p.zhang); dimnames(cell.p.zhang.fdr) = dimnames(cell.p.zhang);
to_plot = cell.p.zhang.fdr[c("yellow", "salmon", "green", "turquoise", "purple", "blue", "greenyellow", "tan"),]


dendro.col = as.dendrogram(hclust(as.dist(1-bicor(zhang.datExpr)), method="average"))
denro.row= as.dendrogram(hclust(as.dist(1-bicor(eigmat[,c("yellow", "salmon", "green", "turquoise", "purple", "blue", "greenyellow", "tan")])),method="average"))

pdf("./results/figures/Manuscript/Fig3F-CellType.pdf",width=6,height=5)
heatmap.2(-log10(to_plot),col=blueWhiteRed(1000,1)[500:1000],
          scale="none",trace="none",cexRow = 0.8,cexCol = .8, density.info = "none",
          colsep=0:7,rowsep=0:8,sepcolor="grey",sepwidth=c(0.02,0.02),
          srtCol=45,offsetRow=0,offsetCol=-0.5,
          Rowv=denro.row, Colv=dendro.col,
          key=T,key.xlab="-log10(P)", cellnote=signif(to_plot,1), notecex=.8, notecol="black",main="Enrichment")
dev.off()


out = cell.p.zhang.fdr
colnames(out) = paste("Enrichment.", colnames(out), ".FDR",sep="")
write.csv(out,file="./results/tables/Manuscript/TableS2-CellType.csv")



## Calculate TFBS enrichment of Disease-associated modules
resultsTFBS = data.frame()
for (m in unique(colors)[!grepl("grey",unique(colors))]) {
  me_name = paste("ME", m, sep="")
  i = which(m == unique(colors))  
  
  ## GO Enrichment
  query = rownames(datExpr)[colors==m]
  go = gprofiler(query, organism="hsapiens", ordered_query = F, significant = T, exclude_iea = F, region_query = F,max_p_value = 1,correction_method = "fdr", custom_bg = rownames(datExpr),
                 hier_filtering = "strong", domain_size = "annotated", numeric_ns = "", include_graph = F,src_filter = c("TF"))
  go = go[order(go$p.value),]
  
  resultsTFBS = rbind(resultsTFBS, cbind(m, go[1:min(nrow(go),20),]))
}
write.csv(resultsTFBS,file="./results/tables/Manuscript/TableS2-TFBS.csv")



## Identify hub genes transcription factors
hubGenes = data.frame(Module=NA, Gene=NA, Symbol=NA,Rank=NA)
for (m in unique(colors)) {
  mod = rownames(datExpr)[colors==m]
  mod = mod[order(kME[mod,m],decreasing = T)[1:20]]
  sym = datProbes[mod,"external_gene_id"]
  hubGenes=rbind(hubGenes, data.frame(Module=m, Gene=mod, Symbol = sym, Rank=1:20))
}
hubGenes = hubGenes[-1,]


ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
a=listAttributes(ensembl);f=listFilters(ensembl)
bm = getBM(attributes = c("ensembl_gene_id", "go_id", "go_linkage_type","goslim_goa_accession","goslim_goa_description"), filters = "ensembl_gene_id", values = hubGenes$Gene, mart=ensembl)
bm = bm[grep("transcription factor",bm$goslim_goa_description),]
hubGenes$TF = bm$goslim_goa_description[match(hubGenes$Gene, bm$ensembl_gene_id)]

write.csv(hubGenes,file="./results/tables/TableS2-HubGeneTFs.csv", row.names = F)


## eRNA overlap with co-expression modules
## Download eRNA modules from http://www.nature.com/neuro/journal/v18/n8/extref/nn.4063-S12.xlsx
eRNAnetwork=  read.csv("raw_data/Annotations//nn.4063-S12.csv")
source("./code/00_scripts/fisher_overlap.R")
table.p = matrix(NA, nrow=length(unique(colors)), ncol=19)
rownames(table.p) = unique(colors); 
colnames(table.p) = paste("M", 1:19,sep="")
table.or = table.p
hgnc = datProbes$external_gene_id
for (m in unique(colors)) {
  for(e in colnames(table.p)) {
    f = ORA(hgnc[colors==m],eRNAnetwork$Symbol[eRNAnetwork$Module.Label==e], hgnc, eRNAnetwork$Symbol)
    table.or[m,e] = as.numeric(f[[1]])
    table.p[m,e] = as.numeric(f[[2]])
  }
}

table.p.fdr = p.adjust(table.p,"fdr")
dim(table.p.fdr) = dim(table.p); dimnames(table.p.fdr) = dimnames(table.p)
table.p.fdr[table.or<1] = 1

to_plot = table.p.fdr[-which(rownames(table.p.fdr)=="grey"),]

pdf(file="./results/figures/Manuscript/FigS11.pdf",width=8,height=6)
heatmap.2(-log10(to_plot),col=blueWhiteRed(1000,1)[500:1000],
          scale="none",trace="none",cexRow = 0.8,cexCol = .8, density.info = "none",
          colsep=0:19,rowsep=0:13,sepcolor="grey",sepwidth=c(0.02,0.02),
          srtCol=45,offsetRow=0,offsetCol=-0.5, Colv=F, keysize=1.1, xlab="eRNA-gene coexpression module",
          key=T,key.xlab="-log10(P)", cellnote=signif(to_plot,1), notecex=.6, notecol="black",main="Brain-Specific Enhancer\nNetwork Enrichment")
dev.off()




#Calculate overlap with Winden modules
#Download winden modules from http://msb.embopress.org/content/5/1/291.long#sec-23
winden= read.csv("raw_data/Annotations//winden.csv")
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="feb2014.archive.ensembl.org")
a=listAttributes(ensembl);f=listFilters(ensembl)
bm1= getBM(attributes = c("affy_moe430a","ensembl_gene_id"),
           filters = "affy_moe430a", values = winden$Affymetrix.ID, mart=ensembl)
bm2 = getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), mart=ensembl)
bm = merge(bm1,bm2)
idx = match(winden$Affymetrix.ID, bm$affy_moe430a)
winden$ensg= bm$hsapiens_homolog_ensembl_gene[idx]

table.p = matrix(NA, nrow=length(unique(colors)), ncol=length(table(winden$Module.Assigned)))
rownames(table.p) = unique(colors); 
colnames(table.p) = na.omit(unique(winden$Module.Assigned))
table.or = table.p

for (m in unique(colors)) {
  for(e in colnames(table.p)) {
    f = ORA(datProbes$ensembl_gene_id[colors==m],winden$ensg[winden$Module.Assigned==e], datProbes$ensembl_gene_id, winden$ensg)
    table.or[m,e] = as.numeric(f[[1]])
    table.p[m,e] = as.numeric(f[[2]])
  }
}
table.p.fdr = p.adjust(table.p,"fdr")
dim(table.p.fdr) = dim(table.p); dimnames(table.p.fdr) = dimnames(table.p)
table.p.fdr[table.or<1] = 1

to_plot=data.frame(Module=rownames(table.p.fdr), Mitochondria="Synaptic", Enrichment=(table.or[,"turquoise"]), FDR = table.p.fdr[,"turquoise"])
to_plot=rbind(to_plot, data.frame(Module=rownames(table.p.fdr), Mitochondria="Non-Synaptic", Enrichment=(table.or[,"blue"]), FDR = table.p.fdr[,"blue"]))
to_plot = to_plot[to_plot$Module %in% c("yellow", "green", "tan", "turquoise", "salmon", "purple", "greenyellow", "blue"),]
to_plot$Module = factor(to_plot$Module)
to_plot$Mitochondria = factor(to_plot$Mitochondria)

g.mito=ggplot(to_plot, aes(x=Module, y=Enrichment,fill=Mitochondria, alpha=FDR<0.05)) + coord_flip() +
  geom_bar(stat="identity",position=position_dodge(),color="grey60") + 
  geom_hline(yintercept = 1,lty=2) + ylab("Enrichment (OR)") + xlab("Cross Disorder Module") +
  ggtitle("Enrichment for Synaptic vs Non-Synaptic Mitochondria", subtitle = "(Winden et al., Mol Syst Biol, 2009)") + 
  theme(plot.subtitle = element_text(hjust=0.5),plot.title = element_text(hjust=0.5))
 ggsave(g.mito, file="results/figures/WGCNA/Modules.mitochondria.pdf") 
 write.csv(to_plot,file="./results/tables/Manuscript/TableS2-Winden.csv",row.names = F)


