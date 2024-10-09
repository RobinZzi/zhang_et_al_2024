rm(list=ls())
setwd("~/bismark/samtools_depth_agg/sum")
save.image("CNV_norm.Rdata")


nt <- as.data.frame(fread("nt.txt"))
pt1 <- as.data.frame(fread("pt1.txt"))
pt2 <- as.data.frame(fread("pt2.txt"))
pt3 <- as.data.frame(fread("pt3.txt"))


common_chr <- Reduce(intersect,list(nt$V1,pt1$V1,pt2$V1,pt3$V1))

nt_sub <- subset(nt,subset=V1 %in% common_chr)
nt_sub$V1 <- NULL
nt_sub <- as.data.frame(lapply(nt_sub,as.numeric))
row.names(pt1_sub) <- common_chr

pt1_sub <- subset(pt1,subset=V1 %in% common_chr)
pt1_sub$V1 <- NULL
pt1_sub <- as.data.frame(lapply(pt1_sub,as.numeric))
row.names(pt1_sub) <- common_chr

pt2_sub <- subset(pt2,subset=V1 %in% common_chr)
pt2_sub$V1 <- NULL
pt2_sub <- as.data.frame(lapply(pt2_sub,as.numeric))
row.names(pt2_sub) <- common_chr

pt3_sub <- subset(pt3,subset=V1 %in% common_chr)
pt3_sub$V1 <- NULL
pt3_sub <- as.data.frame(lapply(pt3_sub,as.numeric))
row.names(pt3_sub) <- common_chr



big_mtx <- as.data.frame(cbind(nt_sub,pt1_sub,pt2_sub,pt3_sub))

sumdepth <- colSums(big_mtx)
sumco <- mean(sumdepth)/sumdepth
big_mtx_norm <- sweep(big_mtx,2,sumco, `*`)  
sumdepth_norm <- colSums(big_mtx_norm)

nt_norm <- select(big_mtx_norm,1:9)

base_line <-as.data.frame(cbind(common_chr,rowMeans(nt_norm)))
colnames(base_line) <- c("bin_name","level")
base_line$level <- as.numeric(base_line$level)

pt_norm <- big_mtx_norm[,10:97]

pt_norm <- sweep(pt_norm,1,base_line$level, `/`)  


pt_norm_t <- as.data.frame(t(pt_norm))

col_anno <- as.data.frame(str_split_fixed(common_chr,"_",2))
col_anno[,3:4] <- str_split_fixed(col_anno$V1,"r",2)
col_anno$V2 <- as.numeric(col_anno$V2)
col_anno$V4 <- as.numeric(col_anno$V4)
col_anno <- mutate(col_anno,id = 100*V4+V2)
row.names(col_anno) <- common_chr
col_anno <- col_anno[order(col_anno$id),]
row_names_factor <- factor(rownames(col_anno)) 
col_anno$bin_name <- rownames(col_anno)
levels(row_names_factor) <- col_anno$bin_name


col_anno_2 <- as.data.frame(dplyr::select(col_anno,1))
colnames(col_anno_2) <- c("chr")

pt_norm_t_sort <- dplyr::select(pt_norm_t,levels(row_names_factor))

col_anno_2$chr <- factor(col_anno_2$chr,levels = c("chr1","chr2","chr3","chr4","chr5",
                                                   "chr6","chr7","chr8","chr9","chr10",
                                                   "chr11","chr12","chr13","chr14","chr15",
                                                   "chr16","chr17","chr18","chr19","chr20",
                                                   "chr21","chr22","chrX","chrY"))

row_anno <- as.data.frame(str_split_fixed(row.names(pt_norm_t_sort),"_",2))
row.names(row_anno) <- row.names(pt_norm_t_sort)
row_anno <- as.data.frame(dplyr::select(row_anno,1))
colnames(row_anno) <- c("sample")


rowname_sort <- c(colnames(pt3_sub),colnames(pt2_sub),colnames(pt1_sub))

ann_colors =list(
  chr=c('chr1'='#E5D2DD','chr2'='#53A85F','chr3'= '#F1BB72','chr4'='#F3B1A0','chr5'='#D6E7A3',
        'chr6'='#57C3F3','chr7'='#476D87','chr8'='#E95C59','chr9'= '#E59CC4','chr10'= '#AB3282', 
        'chr11'='#23452F', 'chr12'='#BD956A', 'chr13'='#8C549C', 'chr14'='#585658', 'chr15'='#9FA3A8',
        'chr16'='#E0D4CA', 'chr17'='#5F3D69','chr18'= '#C5DEBA','chr19'= '#58A4C3','chr20'='#E4C755', 
        'chr21'='#F7F398', 'chr22'='#AA9A59','chrX'= '#E63863', 'chrY'='#E39A35'),
  sample=c('pt1'='#8C549C','pt2'='#D6E7A3','pt3'='#E59CC4')
)



pheatmap(pt_norm_t_sort[rowname_sort,],
         scale = 'row',cluster_rows = F,cluster_cols = F,angle_col = 315,
         breaks = seq(-2, 2, length.out = 100),
         show_colnames = F,show_rownames = F,
         annotation_col = col_anno_2,annotation_row = row_anno,annotation_colors = ann_colors,
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
)


pt1_norm_t <- pt_norm_t_sort[colnames(pt1_sub),]
pt1_norm_avg <- colMeans(pt1_norm_t)
pt2_norm_t <- pt_norm_t_sort[colnames(pt2_sub),]
pt2_norm_avg <- colMeans(pt2_norm_t)
pt3_norm_t <- pt_norm_t_sort[colnames(pt3_sub),]
pt3_norm_avg <- colMeans(pt3_norm_t)


pt_norm_t_avg <- as.data.frame(rbind(pt3_norm_avg,pt2_norm_avg,pt1_norm_avg))
row.names(pt_norm_t_avg) <- c("pt3","pt2","pt1")
pheatmap(pt_norm_t_avg,
         scale = 'row',cluster_rows = T,cluster_cols = F,angle_col = 315,
         breaks = seq(-2, 2, length.out = 100),
         show_colnames = F,show_rownames = T,
         annotation_col = col_anno_2,annotation_colors = ann_colors,
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
)



hcc2_dist = dist(pt_norm_t_avg, method = "euclidean")
hclust_hcc2 = hclust(hcc2_dist, method = "average")
plot(hclust_hcc2)
