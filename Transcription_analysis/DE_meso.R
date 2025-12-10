double_lines=c("Bima1", "Bima1", "Kolf2", "Kolf2", "Kolf3", "Kolf3", "Kucg2", "Kucg2", "Letw5", "Letw5", "Podx1", "Podx1", "Qolg1", "Qolg1", "Sojd3", "Sojd3", "Wibj2", "Wibj2", "Yoch6", "Yoch6")
lines=c("Bima1", "Kolf2", "Kolf3", "Kucg2", "Letw5", "Podx1", "Qolg1", "Sojd3", "Wibj2", "Yoch6")

#Load files into single matrix
library(readr)
library("DESeq2")

annot<-read.delim("HsGRCh38_Ens_gene_protein.gff3", header=F, comment.char="#")
annot_auto<-annot[!is.na(as.numeric(annot[,1])),]
annot_auto<-as.data.frame(annot_auto$V9)
colnames(annot_auto)<-c("geneid")
rownames(annot_auto)<-annot_auto$geneid

files<-{}
for (i in 1:10)
{filename=paste(lines[i],"_meso_rep1_smerged_rmd_stranded_reverse.tsv", sep= "")
files=c(files,filename)
filename=paste(lines[i],"_meso_rep2_smerged_rmd_stranded_reverse.tsv", sep= "")
files=c(files,filename)}

meso_counts<-as.data.frame(annot_auto)
for (i in 1:10){
df<-readr::read_tsv(files[2*i-1], col_names = FALSE, show_col_types = FALSE)
colnames(df)<-c("geneid",paste(lines[i],"r1",sep= ""))
meso_counts<-merge(meso_counts,df,by="geneid")
df<-readr::read_tsv(files[2*i], col_names = FALSE, show_col_types = FALSE)
colnames(df)<-c("geneid",paste(lines[i],"r2",sep= ""))
meso_counts<-merge(meso_counts,df,by="geneid")}

rownames(meso_counts)<-meso_counts[,1]
meso_counts<-meso_counts[,2:21]


vars=c("Bima1r1", "Bima1r2", "Kolf2r1", "Kolf2r2", "Kolf3r1", "Kolf3r2", "Kucg2r1", "Kucg2r2", "Letw5r1", "Letw5r2", "Podx1r1", "Podx1r2", "Qolg1r1", "Qolg1r2", "Sojd3r1", "Sojd3r2", "Wibj2r1", "Wibj2r2", "Yoch6r1", "Yoch6r2")

coldf=data.frame(vars,double_lines)

dds <- DESeqDataSetFromMatrix(countData = meso_counts,
colData=coldf,
design=~double_lines)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
ddsr <- DESeq(dds)

adjpval_meso<-as.data.frame(annot_auto)

for (i in 1:9){
    for (j in (i+1):10){
        adjp<-as.data.frame(results(ddsr, contrast=c("double_lines",lines[i],lines[j])))[6]
        colnames(adjp)<-c(paste(lines[i],lines[j], sep= ""))
        adjp$geneid<-rownames(adjp)
        adjpval_meso<-merge(adjpval_meso,adjp,by="geneid")}}

adjpval_meso[is.na(adjpval_meso)]<-1

rownames(adjpval_meso)<-adjpval_meso[,1]
adjpval_meso<-adjpval_meso[,2:46]

a2pvals_meso<-t(apply(adjpval_meso,1,p.adjust,method="BH"))

l2fc_meso<-as.data.frame(annot_auto)
for (i in 1:9){
for (j in (i+1):10){
lfc<-as.data.frame(results(ddsr, contrast=c("double_lines",lines[i],lines[j])))[2]
colnames(lfc)<-c(paste(lines[i],lines[j], sep= ""))
lfc$geneid<-rownames(lfc)
l2fc_meso<-merge(l2fc_meso,lfc,by="geneid")}}

rownames(l2fc_meso)<-l2fc_meso$geneid
l2fc_meso<-l2fc_meso[2:46]
l2fc_meso<-abs(l2fc_meso)

#Filters
l2fc_meso_bin<-l2fc_meso>3
a2pvals_meso_bin<-a2pvals_meso<0.01
l2fc_meso_bin<-as.data.frame(l2fc_meso_bin)
a2pvals_meso_bin<-as.data.frame(a2pvals_meso_bin)

filter<-a2pvals_meso_bin*l2fc_meso_bin
filter_ns<-as.data.frame(apply(filter,1,sum))
colnames(filter_ns)<-c("sum_filters")

#Final result
filter_ns$geneid<-rownames(filter_ns)
ffilter_ns<-filter_ns[filter_ns$sum_filters>0,]

write.csv(ffilter_ns , "DE_meso.csv")
