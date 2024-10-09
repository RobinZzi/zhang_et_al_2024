
fpkm <- fread("TCGA-LIHC.htseq_fpkm.tsv.gz")
count <- fread("TCGA-LIHC.htseq_counts.tsv.gz")

count$Ensembl_ID <- str_sub(count$Ensembl_ID,1,15)
count <- merge(biomart_table[c(2,1,4)],by = "Ensembl_ID",count)


count<- data.frame(avereps(count,ID=count$Symbol))

count_filt <- count[!is.na(count$Symbol),]
rownames(count_filt) <- count_filt$Symbol
colnames(count_filt) <- colnames(count)