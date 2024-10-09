meth <- fread("TCGA-LIHC.methylation450.tsv.gz",header = T)
anno <- fread("illuminaMethyl450_hg38_GDC")



meth_filt <- na.omit(meth)



bigmeth <- cbind(anno,meth_filt)
bigmeth[,5] <- round(bigmeth[,3]/100000)
bigmeth <- tidyr::unite(bigmeth,"id",chrom,'Composite Element REF',
                        sep="_",remove=TRUE)
bigmeth <- aggregate(bigmeth[,3:432],by=list(id=bigmeth$id),FUN=mean)



bigmeth  <- lapply(bigmeth, function(z) {
  mtx <- cbind(c(NA, head(z, -1)), z, c(tail(z, -1), NA))
  mtx[is.na(mtx[,2]) & rowSums(is.na(mtx)) > 1,] <- NA
  out <- ifelse(is.na(mtx[,2]), rowMeans(mtx, na.rm = TRUE), mtx[,2])
  out[is.nan(out)] <- NA
  out
})
bigmeth <- na.omit(bigmeth)

TCGA_pheatmap <- pheatmap(bigmeth,
                          show_rownames = F, show_colnames = F,
                          cluster_rows = T, cluster_cols = T,
                          clustering_method = "average",
                          color = colorRampPalette(c("#3f72af", "#fcefee", "#d72323"))(20),
                          treeheight_row = 100,
                          treeheight_col = 20,
                          fontsize          = 12)


TCGA_order_row = TCGA_pheatmap$tree_row$order
TCGA_order_col = TCGA_pheatmap$tree_col$order
sort_TAGA_data = data.frame(bigmeth[TCGA_order_row,TCGA_order_col])
sort_TAGA_data_mean <- as.data.frame(colMeans(sort_TAGA_data))
colnames(sort_TAGA_data_mean) <- 'mean'



hypogroup_sm <- subset(sort_TAGA_data2_mean,subset = mean <0.4)
hypergroup_sm <- subset(sort_TAGA_data2_mean,subset = mean >0.55)