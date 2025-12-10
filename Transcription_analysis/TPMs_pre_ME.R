exons_gene<-read.csv("exonic_gene_sizes.csv")
lines=c("Bima1", "Kolf2", "Kolf3", "Kucg2", "Letw5", "Podx1", "Qolg1", "Sojd3", "Wibj2", "Yoch6")
files<-{}
for (i in 1:10)
{filename=paste(lines[i],"_meso_rep1_smerged_rmd_stranded_reverse.tsv", sep= "")
files=c(files,filename)
filename=paste(lines[i],"_meso_rep2_smerged_rmd_stranded_reverse.tsv", sep= "")
files=c(files,filename)}
annot<-read.delim("HsGRCh38_Ens_gene_protein.gff3", header=F, comment.char="#")
annot_auto<-annot[!is.na(as.numeric(annot[,1])),]
annot_auto<-as.data.frame(annot_auto$V9)
colnames(annot_auto)<-c("geneid")
rownames(annot_auto)<-annot_auto$geneid
ipsc_counts<-as.data.frame(annot_auto)

for (i in 1:10){
df<-readr::read_tsv(files[2*i-1], col_names = FALSE, show_col_types = FALSE)
colnames(df)<-c("geneid",paste(lines[i],"r1",sep= ""))
ipsc_counts<-merge(ipsc_counts,df,by="geneid")
df<-readr::read_tsv(files[2*i], col_names = FALSE, show_col_types = FALSE)
colnames(df)<-c("geneid",paste(lines[i],"r2",sep= ""))
ipsc_counts<-merge(ipsc_counts,df,by="geneid")}

rownames(ipsc_counts)<-ipsc_counts[,1]
colnames(exons_gene)<-c("geneid","size")
ipsc_counts<-merge(ipsc_counts,exons_gene,by="geneid")
counts_gene_PKB<-ipsc_counts[,2:21]*1e3/ipsc_counts[,22]
TPM_gene_ipsc<-t(t(counts_gene_PKB)*1e6/colSums(counts_gene_PKB))
TPM_gene_ipsc<-as.data.frame(TPM_gene_ipsc)
row.names(TPM_gene_ipsc)<-ipsc_counts[,1]
TPM_gene_ipsc$geneid<-rownames(TPM_gene_ipsc)
write.csv(TPM_gene_ipsc , "TPMS_pre-meso.csv")
