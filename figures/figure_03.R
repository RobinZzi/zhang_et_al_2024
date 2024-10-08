#Fig_3_A
pheatmap(big_meth,
         annotation_col = big_colanno,annotation_row = pmd_anno,
         annotation_colors = big_ann_colors_new,
         show_rownames = F,show_colnames = F,
         clustering_method = "mcquitty",
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize_col= 15)

ggplot(sorted_bigmeth_colmeans)+
  geom_point(aes(x=sample,y=hmd_level),color='#f9ed69',size=1)+
  geom_line(aes(x=factor(sample),y=hmd_level,group=1),size = 1,color='#f9ed69')+
  geom_line(aes(x=factor(sample),y=pmd_level,group=1),size = 1,color='#3f72af')+
  geom_line(aes(x=factor(sample),y=mean_level,group=1),size = 1,color='#2ca260')+
  geom_point(aes(x=sample,y=pmd_level),color='#3f72af',size=1)+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.title.x = element_text(size = 0),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 0,color="black",angle = 45),axis.text.y = element_text(size = 10,color="black"))



#Fig_3_B
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

ggplot(bigmeth_colmeans,aes(x=patient, y=mean_meth,fill=origin))+
  geom_flat_violin(scale = "width",trim = F)+coord_flip()+geom_jitter(width = 0.1,size=0.5)


#Fig_3_C

pheatmap(sort_tcga_data,
         show_rownames = F,show_colnames = F,
         annotation_col = colanno,annotation_row = pmd_anno_psu,cluster_cols = T,cluster_rows = T,
         clustering_method = "mcquitty",annotation_colors = big_ann_colors_new,
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 1,
         cellheight = 0.03)

pheatmap(sort_trioseq_data,
         show_rownames = F,show_colnames = F,
         annotation_col = colanno,annotation_row = pmd_anno_psu,cluster_cols = F,cluster_rows =T,
         clustering_method = "mcquitty",annotation_colors = big_ann_colors_new,
         clustering_distance_rows = "euclidean",
         color = colorRampPalette(c("#4a74a4", "#f5f6f7", "#b11a2b"))(100),
         treeheight_row = 0,
         treeheight_col = 1,
         fontsize_col= 8,
         angle_col = 45,
         cellwidth = 10,
         cellheight = 0.03)

ggplot(meanlevel_sum[28:457,])+
  geom_point(aes(x=sample,y=hmd_level),color='#f9ed69',size=2)+
  geom_point(aes(x=sample,y=pmd_level),color='#3f72af',size=2)+
  geom_point(aes(x=sample,y=mean_level),color='#2ca260',size=2)+
  geom_line(aes(x=factor(sample),y=hmd_level,group=1),size = 1.5,color='#f9ed69')+
  geom_line(aes(x=factor(sample),y=pmd_level,group=1),size = 1.5,color='#3f72af')+
  geom_line(aes(x=factor(sample),y=mean_level,group=1),size = 1.5,color='#2ca260')+
  theme_bw()+ylim(0.2,0.8)+
  theme(panel.grid = element_blank(),axis.title.x = element_text(size = 0),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 0,color="black",angle = 45),axis.text.y = element_text(size = 10,color="black"))

ggplot(meanlevel_sum[1:27,])+
  geom_point(aes(x=sample,y=hmd_level),color='#f9ed69',size=2)+
  geom_point(aes(x=sample,y=pmd_level),color='#3f72af',size=2)+
  geom_point(aes(x=sample,y=mean_level),color='#2ca260',size=2)+
  geom_line(aes(x=factor(sample),y=hmd_level,group=1),size = 1.5,color='#f9ed69')+
  geom_line(aes(x=factor(sample),y=pmd_level,group=1),size = 1.5,color='#3f72af')+
  geom_line(aes(x=factor(sample),y=mean_level,group=1),size = 1.5,color='#2ca260')+
  theme_bw()+ylim(0.2,0.8)+
  theme(panel.grid = element_blank(),axis.title.x = element_text(size = 0),axis.title.y = element_text(size = 14),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 0,color="black",angle = 45),axis.text.y = element_text(size = 10,color="black"))


#Fig_3_D
ggplot(he_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=12,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))+
  stat_compare_means(label = "p.signif")
ggplot(eu_sum,aes(x=hm_type,y=length,fill=region))+
  geom_boxplot(outlier.colour = NA)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=12,face = "bold"),
        axis.title.x=element_text(size=0),
        axis.title.y=element_text(size=13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12))+
  stat_compare_means(label = "p.signif")


