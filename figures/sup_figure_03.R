#Sup_Fig_3A
pheatmap(big_psu_meth,
         annotation_col = big_colanno_new_psu[,2:3],annotation_row = pmd_anno_psu,
         annotation_colors = big_ann_colors_new,
         show_rownames = F,show_colnames = T,
         clustering_method = "mcquitty",
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 0,
         angle_col = 45,
         fontsize_col= 8)


#Sup_Fig_3B
DimPlot(bigmeth_seurat, label = F,group.by = "origin",
        cols = c('HCC2'='#E5D2DD','HCC3'='#BD956A','HCC6'='#F1BB72','HCC7'='#8C549C','HCC8'='#D6E7A3','HCC9'='#E59CC4'))
DimPlot(bigmeth_seurat, label = F,group.by = "tissue",
        cols=c("Normal Tissue"="#3f72af","Tumor Tissue"="#d72323"))


#Sup_Fig_3C
pheatmap("sample"_base_removena,
         annotation_col = "sample"_colanno,annotation_row = "sample"_anno,
         annotation_colors = "sample"_ann_colors,
         show_rownames = F,show_colnames = F,
         clustering_method = "complete",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15)
