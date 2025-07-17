library(readr)
library("DESeq2")

#Load count matrix
ipsc_counts<-read.csv("ipsc_counts.csv",row.names = 1)

### Differential Expression analysis ###

#Prepare the data for DESEQ2
ipsc_counts<-ipsc_counts[2:22]
double_lines=c("Bima1", "Bima1", "Kolf2", "Kolf2", "Kolf3", "Kolf3", "Kucg2", "Kucg2", "Letw5", "Letw5", "Podx1", "Podx1", "Qolg1", "Qolg1", "Sojd3", "Sojd3", "Wibj2", "Wibj2", "Yoch6", "Yoch6", "Yoch6")
vars=c("Bima1r1", "Bima1r2", "Kolf2r1", "Kolf2r2", "Kolf3r1", "Kolf3r2", "Kucg2r1", "Kucg2r2", "Letw5r1", "Letw5r2", "Podx1r1", "Podx1r2", "Qolg1r1", "Qolg1r2", "Sojd3r1", "Sojd3r2", "Wibj2r1", "Wibj2r2", "Yoch6r1", "Yoch6r2", "Yoch6r3")
coldf=data.frame(vars,double_lines)

#DESeq2 call
dds <- DESeqDataSetFromMatrix(countData = ipsc_counts,
                               colData=coldf,
                               design=~double_lines)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
ddsr <- DESeq(dds)

#Output 1:Background genes (later used for GO analysis)
background_genes<-rownames(ipsc_counts[keep,])
write.csv(background_genes, "genes_expressed_ipsc.csv")

#Pvalues. 
#Retrieve the pairwise p-values. Adjust them for multiple comparison. 
adjpval_ipsc<-as.data.frame(annot_auto)

for (i in 1:9){
    for (j in (i+1):10){
        adjp<-as.data.frame(results(ddsr, contrast=c("double_lines",lines[i],lines[j])))[6]
        colnames(adjp)<-c(paste(lines[i],lines[j], sep= ""))
        adjp$geneid<-rownames(adjp)
        adjpval_ipsc<-merge(adjpval_ipsc,adjp,by="geneid")}}

adjpval_ipsc[is.na(adjpval_ipsc)]<-1
rownames(adjpval_ipsc)<-adjpval_ipsc[,1]
adjpval_ipsc<-adjpval_ipsc[,2:46]

a2pvals_ipsc<-t(apply(adjpval_ipsc,1,p.adjust,method="BH"))

#Log2 Fold Changes (pairwise).
l2fc_ipsc<-as.data.frame(annot_auto)
for (i in 1:9){
for (j in (i+1):10){
lfc<-as.data.frame(results(ddsr, contrast=c("double_lines",lines[i],lines[j])))[2]
colnames(lfc)<-c(paste(lines[i],lines[j], sep= ""))
lfc$geneid<-rownames(lfc)
l2fc_ipsc<-merge(l2fc_ipsc,lfc,by="geneid")}}

rownames(l2fc_ipsc)<-l2fc_ipsc$geneid
l2fc_ipsc<-l2fc_ipsc[2:46]
l2fc_ipsc<-abs(l2fc_ipsc)

#Filters to find the DE genes for the set thresholds (in the same pairwise comparison).
l2fc_ipsc_bin<-l2fc_ipsc>3
a2pvals_ipsc_bin<-a2pvals_ipsc<0.01
l2fc_ipsc_bin<-as.data.frame(l2fc_ipsc_bin)
a2pvals_ipsc_bin<-as.data.frame(a2pvals_ipsc_bin)

filter<-a2pvals_ipsc_bin*l2fc_ipsc_bin
filter_ns<-as.data.frame(apply(filter,1,sum))
colnames(filter_ns)<-c("sum_filters")

#Final result
filter_ns$geneid<-rownames(filter_ns)
ffilter_ns<-filter_ns[filter_ns$sum_filters>0,]

#Output 2: Set of Differentially Expressed genes
write.csv(ffilter_ns , "DE_hiPSC.csv")

### Gene Ontology analysis ###

#Preparation
background_genes<-as.data.frame(background_genes)

library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)

#GO enrichment analysis for Molecular Function
ego <- enrichGO(gene          = ffilter_ns$geneid,
                 universe      = background_genes$background_genes,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 keyType="ENSEMBL",
		readable      = TRUE,)

#Output 3: Molecular Function GO of DE genes (as a png image)
png(filename="MF_GO_DE.png",  width = 7, height = 4.5, units = 'in', res=300)
dotplot(ego, showCategory=5, font.size = 15)+ ggplot2::theme(legend.text=ggplot2::element_text(size=15),text = ggplot2::element_text(size=15))
dev.off()

#GO enrichment analysis for Biological Process
ego <- enrichGO(gene          = ffilter_ns$geneid,
                 universe      = background_genes$background_genes,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 keyType="ENSEMBL",
		readable      = TRUE,)

dotplot(ego, showCategory=5, font.size = 15)+ ggplot2::theme(legend.text=ggplot2::element_text(size=15),text = ggplot2::element_text(size=15))

#Output 3: Biological Process GO of DE genes (as a png image)
png(filename="BP_GO_DE.png",  width = 7, height = 4.5, units = 'in', res=300)
dotplot(ego, showCategory=5, font.size = 15)+ ggplot2::theme(legend.text=ggplot2::element_text(size=15),text = ggplot2::element_text(size=15))
dev.off()




