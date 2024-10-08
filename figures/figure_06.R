#Figure_6A
cnvkit.py heatmap *_sorted.cns -d

#Figure_6B
pheatmap(hcc2_sum_subset,
         show_rownames = F,show_colnames = T,cluster_cols = T,cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,border_color = NA,
         treeheight_col = 10)


hcc2_sum_mtx <- t(hcc2_sum_subset)
hcc2_sum_dist_snv = dist(hcc2_sum_mtx, method = "euclidean")
hclust_hcc2_sum_snv = hclust(hcc2_sum_dist_snv, method = "complete")
plot(hclust_hcc2_sum_snv)



#Figure_6C
infercnv::plot_cnv(
  infercnv_obj,
  out_dir = paste0("result/infercnv/",condition,"/",sample),
  custom_color_pal = color.palette(c("#4a74a4", "#f5f6f7", "#b11a2b"),c(2, 2)),
  plot_chr_scale = T,
  output_filename = "infercnv_replot",
  output_format = "png"
)

#Figure_6D
pheatmap(hcc2_norm_t_avg,
         scale = 'row',cluster_rows = T,cluster_cols = F,angle_col = 315,
         breaks = seq(-2, 2, length.out = 100),
         show_colnames = F,show_rownames = T,
         annotation_col = col_anno_2,annotation_colors = ann_colors,
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
)



hcc2_dist_cnv = dist(pt_norm_t_avg, method = "euclidean")
hclust_hcc2_cnv = hclust(hcc2_dist_cnv, method = "average")
plot(hclust_hcc2_cnv)


#Figure_6E
pheatmap(hcc2_psu_removena,
         annotation_col = big_colanno_new_psu[,2:3],annotation_row = pmd_anno_psu,
         annotation_colors = big_ann_colors_new,
         show_rownames = F,show_colnames = T,
         clustering_method = "mcquitty",
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 10,
         angle_col = 45,
         fontsize_col= 15,
         cellwidth=80)
hcc2_psu_removena_t <- as.data.frame(t(hcc2_psu_removena))
hcc2_psu_dist = dist(hcc2_psu_removena_t, method = "euclidean")
hclust_hcc2_psu = hclust(hcc2_psu_dist, method = "average")
plot(hclust_hcc2_psu)

#Figure_6F
DimPlot(hcc2_hpc_sub, label = F,group.by = "orig.ident",cols=c('#8C549C','#D6E7A3','#E59CC4'))
hcc2_hpc_PT1 <- subset(hcc2_hpc , subset = orig.ident == "PT1")
hcc2_PT1_distance <- as.data.frame(Embeddings(object = hcc2_hpc_PT1[["umap"]]))
hcc2_hpc_PT2 <- subset(hcc2_hpc , subset = orig.ident == "PT2")
hcc2_PT2_distance <- as.data.frame(Embeddings(object = hcc2_hpc_PT2[["umap"]]))
hcc2_hpc_PT3 <- subset(hcc2_hpc , subset = orig.ident == "PT3")
hcc2_PT3_distance <- as.data.frame(Embeddings(object = hcc2_hpc_PT3[["umap"]]))


hcc2_distance_umap <- as.data.frame(rbind(colMeans(hcc2_PT1_distance),colMeans(hcc2_PT2_distance),colMeans(hcc2_PT3_distance)))
row.names(hcc2_distance_umap) <- c("PT1","PT2","PT3")


dist_hcc2_umap = dist(hcc2_distance, method = "euclidean")
hclust_dist_hcc2_umap = hclust(dist_hcc2_umap, method = "complete")
plot(hclust_dist_hcc2_umap)