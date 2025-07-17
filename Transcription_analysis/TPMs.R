#Obtain exon sizes per gene and add them up
library(GenomicFeatures)

#Prior to run the script, download file from Ensembl
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.111.gff3",format="gff3")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
exons_per_gene=data.frame(exonic.gene.sizes)

#Output1: Exonic gene size file.
write.csv(exons_per_gene, "exonic_gene_sizes.csv")

#Load files. Cell lines and replicates. We've got a third Yoch6 replicate.
lines=c("Bima1", "Kolf2", "Kolf3", "Kucg2", "Letw5", "Podx1", "Qolg1", "Sojd3", "Wibj2", "Yoch6")

files<-{}
for (i in 1:10)
 {filename=paste(lines[i],"_rep1_merged_rmd.tsv", sep= "")
 files=c(files,filename)
 filename=paste(lines[i],"_rep2_merged_rmd.tsv", sep= "")
 files=c(files,filename)}

files=c(files,"Yoch6_rep3_merged_rmd.tsv")

#Read the files. First load the annotation of protein coding genes.

annot<-read.delim("HsGRCh38_Ens_gene_protein.gff3", header=F, comment.char="#")
annot_auto<-annot[!is.na(as.numeric(annot[,1])),]
annot_auto<-as.data.frame(annot_auto$V9)
colnames(annot_auto)<-c("geneid")
rownames(annot_auto)<-annot_auto$geneid

#Now the actual loading of the data.
ipsc_counts<-as.data.frame(annot_auto)
for (i in 1:10){
    df<-readr::read_tsv(files[2*i-1], col_names = FALSE, show_col_types = FALSE)
    colnames(df)<-c("geneid",paste(lines[i],"r1",sep= ""))
    ipsc_counts<-merge(ipsc_counts,df,by="geneid")
    df<-readr::read_tsv(files[2*i], col_names = FALSE, show_col_types = FALSE)
    colnames(df)<-c("geneid",paste(lines[i],"r2",sep= ""))
    ipsc_counts<-merge(ipsc_counts,df,by="geneid")}
                                                                                
df<-readr::read_tsv(files[21], col_names = FALSE, show_col_types = FALSE)
colnames(df)<-c("geneid",paste(lines[10],"r3",sep= ""))                       
ipsc_counts<-merge(ipsc_counts,df,by="geneid")
rownames(ipsc_counts)<-ipsc_counts[,1]

#Output2: Count matrix for each replicate
write.csv(ipsc_counts , "ipsc_counts.csv")

#Recover the exonic lengths as estimate of transcript length
exons_gene<-read.csv("exonic_gene_sizes.csv")
colnames(exons_gene)<-c("geneid","size")

ipsc_counts<-merge(ipsc_counts,exons_gene,by="geneid")

#The actual TPM computation.

counts_gene_PKB<-ipsc_counts[,2:22]*1e3/ipsc_counts[,23]
TPM_gene_ipsc<-t(t(counts_gene_PKB)*1e6/colSums(counts_gene_PKB))
TPM_gene_ipsc<-as.data.frame(TPM_gene_ipsc)
row.names(TPM_gene_ipsc)<-ipsc_counts[,1]
TPM_gene_ipsc$geneid<-rownames(TPM_gene_ipsc)

#Output2: Matrix of TPMs per line and replicate
write.csv(TPM_gene_ipsc , "TPMS_hiPSC.csv")
