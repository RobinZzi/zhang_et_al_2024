library(ggpubr)

hcc1_av <- AverageExpression(HCC1_HPC,group.by = "orig.ident",assays = 'RNA')
hcc1_av <- hcc1_av[[1]]
hcc1_cg=names(tail(sort(apply(hcc1_av, 1, sd)),1000))
pheatmap::pheatmap(cor(hcc1_av[hcc1_cg,],method = 'spearman'))




hcc2_av <- AverageExpression(HCC2_HPC,group.by = "orig.ident",assays = 'RNA')
hcc2_av <- hcc2_av[[1]]
hcc2_cg=names(tail(sort(apply(hcc2_av, 1, sd)),1000))
pheatmap::pheatmap(cor(hcc2_av[hcc2_cg,],method = 'spearman'))




hcc3_av <- AverageExpression(HCC3_HPC,group.by = "orig.ident",assays = 'RNA')
hcc3_av <- hcc3_av[[1]]
hcc3_cg=names(tail(sort(apply(hcc3_av, 1, sd)),1000))
pheatmap::pheatmap(cor(hcc3_av[hcc3_cg,],method = 'spearman'))





hcc1_all <- as.data.frame(GetAssayData(object = HCC1_HPC, slot = "data"))
hcc1_top <- hcc1_all[top3000_gene,]
hcc1_top <- na.omit(hcc1_top)
hcc1_cor <- as.data.frame(cor(hcc1_top,method = 'pearson'))

hcc1_anno <- as.data.frame(HCC1_HPC$orig.ident)
hcc1_anno2 <- as.data.frame(HCC1_HPC$satellite_region)
colnames(hcc1_anno) <- "origin"
colnames(hcc1_anno2) <- "origin"
pheatmap(hcc1_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc1_anno,annotation_col = hcc1_anno,
         cluster_rows = F,cluster_cols = F)
pheatmap(hcc1_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc1_anno2,annotation_col = hcc1_anno2,
         cluster_rows = F,cluster_cols = F)
pheatmap(hcc1_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc1_anno2,annotation_col = hcc1_anno2,
         cluster_rows = T,cluster_cols = T)

table(HCC1_HPC$orig.ident)

hcc1_pt2_pt4_mean <- mean(unlist(hcc1_cor[149:1393,1394:2180]))
hcc1_pt2_pt5_mean <- mean(unlist(hcc1_cor[149:1393,2181:2522]))
hcc1_pt4_pt5_mean <- mean(unlist(hcc1_cor[1394:2180,2181:2522]))

hcc1_pt2_pt2_mean <- mean(unlist(hcc1_cor[149:1393,149:1393]))
hcc1_pt4_pt4_mean <- mean(unlist(hcc1_cor[1394:2180,1394:2180]))
hcc1_pt5_pt5_mean <- mean(unlist(hcc1_cor[2181:2522,2181:2522]))

hcc1_pt2_pt4 <- as.data.frame(unlist(hcc1_cor[149:1393,1394:2180]))
colnames(hcc1_pt2_pt4) <- "cor"
hcc1_pt2_pt4$pair <- "hcc1_pt2_pt4"
hcc1_pt2_pt4$type <- "cross"
hcc1_pt2_pt4$patient <- "hcc1"
hcc1_pt2_pt5 <- as.data.frame(unlist(hcc1_cor[149:1393,2181:2522]))
colnames(hcc1_pt2_pt5) <- "cor"
hcc1_pt2_pt5$pair <- "hcc1_pt2_pt5"
hcc1_pt2_pt5$type <- "cross"
hcc1_pt2_pt5$patient <- "hcc1"
hcc1_pt4_pt5 <- as.data.frame(unlist(hcc1_cor[1394:2180,2181:2522]))
colnames(hcc1_pt4_pt5) <- "cor"
hcc1_pt4_pt5$pair <- "hcc1_pt4_pt5"
hcc1_pt4_pt5$type <- "cross"
hcc1_pt4_pt5$patient <- "hcc1"

hcc1_pt2_pt2 <- as.data.frame(unlist(hcc1_cor[149:1393,149:1393]))
colnames(hcc1_pt2_pt2) <- "cor"
hcc1_pt2_pt2$pair <- "hcc1_pt2_pt2"
hcc1_pt2_pt2$type <- "self"
hcc1_pt2_pt2$patient <- "hcc1"
hcc1_pt4_pt4 <- as.data.frame(unlist(hcc1_cor[1394:2180,1394:2180]))
colnames(hcc1_pt4_pt4) <- "cor"
hcc1_pt4_pt4$pair <- "hcc1_pt4_pt4"
hcc1_pt4_pt4$type <- "self"
hcc1_pt4_pt4$patient <- "hcc1"
hcc1_pt5_pt5 <- as.data.frame(unlist(hcc1_cor[2181:2522,2181:2522]))
colnames(hcc1_pt5_pt5) <- "cor"
hcc1_pt5_pt5$pair <- "hcc1_pt5_pt5"
hcc1_pt5_pt5$type <- "self"
hcc1_pt5_pt5$patient <- "hcc1"

hcc1_cor_all <- as.data.frame(rbind(hcc1_pt2_pt4,hcc1_pt2_pt5,hcc1_pt4_pt5,
                                    hcc1_pt2_pt2,hcc1_pt4_pt4,hcc1_pt5_pt5))

ggplot(hcc1_cor_all,aes(type,cor,fill=type))+
  geom_boxplot(width=0.5)+
  scale_fill_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())




hcc2_all <- as.data.frame(GetAssayData(object = HCC2_HPC, slot = "data"))
hcc2_top <- hcc2_all[_gene,]
hcc2_top <- na.omit(hcc2_top)
hcc2_cor <- as.data.frame(cor(hcc2_top,method = 'pearson'))
hcc2_anno <- as.data.frame(HCC2_HPC$orig.ident)
colnames(hcc2_anno) <- "origin"
pheatmap(hcc2_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc2_anno,annotation_col = hcc2_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC2_HPC$orig.ident)

hcc2_pt1_pt2_mean <- mean(unlist(hcc2_cor[74:1723,1624:2871]))
hcc2_pt1_pt3_mean <- mean(unlist(hcc2_cor[74:1623,2872:3169]))
hcc2_pt1_pt4_mean <- mean(unlist(hcc2_cor[74:1623,3170:3255]))
hcc2_pt2_pt3_mean <- mean(unlist(hcc2_cor[1624:2871,2872:3169]))
hcc2_pt2_pt4_mean <- mean(unlist(hcc2_cor[1624:2871,3170:3255]))
hcc2_pt3_pt4_mean <- mean(unlist(hcc2_cor[2872:3169,3170:3255]))


hcc2_pt1_pt1_mean <- mean(unlist(hcc2_cor[74:1623,74:1623]))
hcc2_pt2_pt2_mean <- mean(unlist(hcc2_cor[1624:2871,1624:2871]))
hcc2_pt3_pt3_mean <- mean(unlist(hcc2_cor[2872:3169,2872:3169]))
hcc2_pt4_pt4_mean <- mean(unlist(hcc2_cor[3170:3255,3170:3255]))





hcc2_pt1_pt2 <- as.data.frame(unlist(hcc2_cor[74:1723,1624:2871]))
colnames(hcc2_pt1_pt2) <- "cor"
hcc2_pt1_pt2$pair <- "hcc2_pt1_pt2"
hcc2_pt1_pt2$type <- "cross"
hcc2_pt1_pt2$patient <- "hcc2"

hcc2_pt1_pt3 <- as.data.frame(unlist(hcc2_cor[74:1623,2872:3169]))
colnames(hcc2_pt1_pt3) <- "cor"
hcc2_pt1_pt3$pair <- "hcc2_pt1_pt3"
hcc2_pt1_pt3$type <- "cross"
hcc2_pt1_pt3$patient <- "hcc2"

hcc2_pt1_pt4 <- as.data.frame(unlist(hcc2_cor[74:1623,3170:3255]))
colnames(hcc2_pt1_pt4) <- "cor"
hcc2_pt1_pt4$pair <- "hcc2_pt1_pt4"
hcc2_pt1_pt4$type <- "cross"
hcc2_pt1_pt4$patient <- "hcc2"

hcc2_pt2_pt3 <- as.data.frame(unlist(hcc2_cor[1624:2871,2872:3169]))
colnames(hcc2_pt2_pt3) <- "cor"
hcc2_pt2_pt3$pair <- "hcc2_pt2_pt3"
hcc2_pt2_pt3$type <- "cross"
hcc2_pt2_pt3$patient <- "hcc2"

hcc2_pt2_pt4 <- as.data.frame(unlist(hcc2_cor[1624:2871,3170:3255]))
colnames(hcc2_pt2_pt4) <- "cor"
hcc2_pt2_pt4$pair <- "hcc2_pt2_pt4"
hcc2_pt2_pt4$type <- "cross"
hcc2_pt2_pt4$patient <- "hcc2"

hcc2_pt3_pt4 <- as.data.frame(unlist(hcc2_cor[2872:3169,3170:3255]))
colnames(hcc2_pt3_pt4) <- "cor"
hcc2_pt3_pt4$pair <- "hcc2_pt3_pt4"
hcc2_pt3_pt4$type <- "cross"
hcc2_pt3_pt4$patient <- "hcc2"


hcc2_pt1_pt1 <- as.data.frame(unlist(hcc2_cor[74:1623,74:1623]))
colnames(hcc2_pt1_pt1) <- "cor"
hcc2_pt1_pt1$pair <- "hcc2_pt1_pt1"
hcc2_pt1_pt1$type <- "self"
hcc2_pt1_pt1$patient <- "hcc2"

hcc2_pt2_pt2 <- as.data.frame(unlist(hcc2_cor[1624:2871,1624:2871]))
colnames(hcc2_pt2_pt2) <- "cor"
hcc2_pt2_pt2$pair <- "hcc2_pt2_pt2"
hcc2_pt2_pt2$type <- "self"
hcc2_pt2_pt2$patient <- "hcc2"

hcc2_pt3_pt3 <- as.data.frame(unlist(hcc2_cor[2872:3169,2872:3169]))
colnames(hcc2_pt3_pt3) <- "cor"
hcc2_pt3_pt3$pair <- "hcc2_pt3_pt3"
hcc2_pt3_pt3$type <- "self"
hcc2_pt3_pt3$patient <- "hcc2"

hcc2_pt4_pt4 <- as.data.frame(unlist(hcc2_cor[3170:3255,3170:3255]))
colnames(hcc2_pt4_pt4) <- "cor"
hcc2_pt4_pt4$pair <- "hcc2_pt4_pt4"
hcc2_pt4_pt4$type <- "self"
hcc2_pt4_pt4$patient <- "hcc2"

hcc2_cor_all <- as.data.frame(rbind(hcc2_pt1_pt2,hcc2_pt1_pt3,hcc2_pt1_pt4,
                                    hcc2_pt2_pt3,hcc2_pt2_pt4,hcc2_pt3_pt4,
                                    hcc2_pt1_pt1,hcc2_pt2_pt2,hcc2_pt3_pt3,
                                    hcc2_pt4_pt4))

ggplot(hcc2_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

hcc3_all <- as.data.frame(GetAssayData(object = HCC3_HPC, slot = "data"))
hcc3_top <- hcc3_all[_gene,]
hcc3_top <- na.omit(hcc3_top)
hcc3_cor <- as.data.frame(cor(hcc3_top,method = 'pearson'))
hcc3_anno <- as.data.frame(HCC3_HPC$orig.ident)
colnames(hcc3_anno) <- "origin"
pheatmap(hcc3_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc3_anno,annotation_col = hcc3_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC3_HPC$orig.ident)

hcc3_pt1_pt2_mean <- mean(unlist(hcc3_cor[13:156,157:245]))
hcc3_pt1_pt3_mean <- mean(unlist(hcc3_cor[13:156,245:248]))
hcc3_pt2_pt3_mean <- mean(unlist(hcc3_cor[157:245,245:248]))

hcc3_pt1_pt1_mean <- mean(unlist(hcc3_cor[13:156,13:156]))
hcc3_pt2_pt2_mean <- mean(unlist(hcc3_cor[157:245,157:245]))
hcc3_pt3_pt3_mean <- mean(unlist(hcc3_cor[245:248,245:248]))

hcc3_pt1_pt2 <- as.data.frame(unlist(hcc3_cor[13:156,157:245]))
colnames(hcc3_pt1_pt2) <- "cor"
hcc3_pt1_pt2$pair <- "hcc3_pt1_pt2"
hcc3_pt1_pt2$type <- "cross"
hcc3_pt1_pt2$patient <- "hcc3"
hcc3_pt1_pt3 <- as.data.frame(unlist(hcc3_cor[13:156,245:248]))
colnames(hcc3_pt1_pt3) <- "cor"
hcc3_pt1_pt3$pair <- "hcc3_pt1_pt3"
hcc3_pt1_pt3$type <- "cross"
hcc3_pt1_pt3$patient <- "hcc3"
hcc3_pt2_pt3 <- as.data.frame(unlist(hcc3_cor[157:245,245:248]))
colnames(hcc3_pt2_pt3) <- "cor"
hcc3_pt2_pt3$pair <- "hcc3_pt2_pt3"
hcc3_pt2_pt3$type <- "cross"
hcc3_pt2_pt3$patient <- "hcc3"

hcc3_pt1_pt1 <- as.data.frame(unlist(hcc3_cor[13:156,13:156]))
colnames(hcc3_pt1_pt1) <- "cor"
hcc3_pt1_pt1$pair <- "hcc3_pt1_pt1"
hcc3_pt1_pt1$type <- "self"
hcc3_pt1_pt1$patient <- "hcc3"
hcc3_pt2_pt2 <- as.data.frame(unlist(hcc3_cor[157:245,157:245]))
colnames(hcc3_pt2_pt2) <- "cor"
hcc3_pt2_pt2$pair <- "hcc3_pt2_pt2"
hcc3_pt2_pt2$type <- "self"
hcc3_pt2_pt2$patient <- "hcc3"
hcc3_pt3_pt3 <- as.data.frame(unlist(hcc3_cor[245:248,245:248]))
colnames(hcc3_pt3_pt3) <- "cor"
hcc3_pt3_pt3$pair <- "hcc3_pt3_pt3"
hcc3_pt3_pt3$type <- "self"
hcc3_pt3_pt3$patient <- "hcc3"

hcc3_cor_all <- as.data.frame(rbind(hcc3_pt1_pt2,hcc3_pt1_pt3,hcc3_pt2_pt3,
                                    hcc3_pt1_pt1,hcc3_pt2_pt2,hcc3_pt3_pt3))


ggplot(hcc3_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(hcc3_cor_all,aes(type,cor,color=type))+
  geom_violin(width=0.5)+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


hcc4_all <- as.data.frame(GetAssayData(object = HCC4_HPC, slot = "data"))
hcc4_top <- hcc4_all[_gene,]
hcc4_top <- na.omit(hcc4_top)
hcc4_cor <- as.data.frame(cor(hcc4_top,method = 'pearson'))
hcc4_anno <- as.data.frame(HCC4_HPC$orig.ident)
colnames(hcc4_anno) <- "origin"
pheatmap(hcc4_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc4_anno,annotation_col = hcc4_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC4_HPC$orig.ident)

hcc4_pt1_pt2_mean <- mean(unlist(hcc4_cor[3:16,17:61]))
hcc4_pt1_pt3_mean <- mean(unlist(hcc4_cor[3:16,62:69]))
hcc4_pt2_pt3_mean <- mean(unlist(hcc4_cor[17:61,62:69]))


hcc4_pt1_pt1_mean <- mean(unlist(hcc4_cor[3:16,3:16]))
hcc4_pt2_pt2_mean <- mean(unlist(hcc4_cor[17:61,17:61]))
hcc4_pt3_pt3_mean <- mean(unlist(hcc4_cor[62:69,62:69]))


hcc4_pt1_pt2 <- as.data.frame(unlist(hcc4_cor[3:16,17:61]))
colnames(hcc4_pt1_pt2) <- "cor"
hcc4_pt1_pt2$pair <- "hcc4_pt1_pt2"
hcc4_pt1_pt2$type <- "cross"
hcc4_pt1_pt2$patient <- "hcc4"
hcc4_pt1_pt3 <- as.data.frame(unlist(hcc4_cor[3:16,62:69]))
colnames(hcc4_pt1_pt3) <- "cor"
hcc4_pt1_pt3$pair <- "hcc4_pt1_pt3"
hcc4_pt1_pt3$type <- "cross"
hcc4_pt1_pt3$patient <- "hcc4"
hcc4_pt2_pt3 <- as.data.frame(unlist(hcc4_cor[17:61,62:69]))
colnames(hcc4_pt2_pt3) <- "cor"
hcc4_pt2_pt3$pair <- "hcc4_pt2_pt3"
hcc4_pt2_pt3$type <- "cross"
hcc4_pt2_pt3$patient <- "hcc4"

hcc4_pt1_pt1 <- as.data.frame(unlist(hcc4_cor[3:16,3:16]))
colnames(hcc4_pt1_pt1) <- "cor"
hcc4_pt1_pt1$pair <- "hcc4_pt1_pt1"
hcc4_pt1_pt1$type <- "self"
hcc4_pt1_pt1$patient <- "hcc4"
hcc4_pt2_pt2 <- as.data.frame(unlist(hcc4_cor[17:61,17:61]))
colnames(hcc4_pt2_pt2) <- "cor"
hcc4_pt2_pt2$pair <- "hcc4_pt2_pt2"
hcc4_pt2_pt2$type <- "self"
hcc4_pt2_pt2$patient <- "hcc4"
hcc4_pt3_pt3 <- as.data.frame(unlist(hcc4_cor[62:69,62:69]))
colnames(hcc4_pt3_pt3) <- "cor"
hcc4_pt3_pt3$pair <- "hcc4_pt3_pt3"
hcc4_pt3_pt3$type <- "self"
hcc4_pt3_pt3$patient <- "hcc4"


hcc4_cor_all <- as.data.frame(rbind(hcc4_pt1_pt2,hcc4_pt1_pt3,hcc4_pt2_pt3,
                                    hcc4_pt1_pt1,hcc4_pt2_pt2,hcc4_pt3_pt3))


ggplot(hcc4_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

hcc5_all <- as.data.frame(GetAssayData(object = HCC5_HPC, slot = "data"))
hcc5_top <- hcc5_all[_gene,]
hcc5_top <- na.omit(hcc5_top)
hcc5_cor <- as.data.frame(cor(hcc5_top,method = 'pearson'))
hcc5_anno <- as.data.frame(HCC5_HPC$orig.ident)
colnames(hcc5_anno) <- "origin"
pheatmap(hcc5_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc5_anno,annotation_col = hcc5_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC5_HPC$orig.ident)

hcc5_pt1_pt2_mean <- mean(unlist(hcc5_cor[78:84,85:103]))
hcc5_pt1_pt3_mean <- mean(unlist(hcc5_cor[78:84,104:144]))
hcc5_pt1_pt4_mean <- mean(unlist(hcc5_cor[78:84,145:169]))
hcc5_pt1_pt5_mean <- mean(unlist(hcc5_cor[78:84,170:181]))
hcc5_pt2_pt3_mean <- mean(unlist(hcc5_cor[85:103,104:144]))
hcc5_pt2_pt4_mean <- mean(unlist(hcc5_cor[85:103,145:169]))
hcc5_pt2_pt5_mean <- mean(unlist(hcc5_cor[85:103,170:181]))
hcc5_pt3_pt4_mean <- mean(unlist(hcc5_cor[104:144,145:169]))
hcc5_pt3_pt5_mean <- mean(unlist(hcc5_cor[104:144,170:181]))
hcc5_pt4_pt5_mean <- mean(unlist(hcc5_cor[145:169,170:181]))


hcc5_pt1_pt1_mean <- mean(unlist(hcc5_cor[78:84,78:84]))
hcc5_pt2_pt2_mean <- mean(unlist(hcc5_cor[85:103,85:103]))
hcc5_pt3_pt3_mean <- mean(unlist(hcc5_cor[104:144,104:144]))
hcc5_pt4_pt4_mean <- mean(unlist(hcc5_cor[145:169,145:169]))
hcc5_pt5_pt5_mean <- mean(unlist(hcc5_cor[170:181,170:181]))



hcc5_pt1_pt2 <- as.data.frame(unlist(hcc5_cor[78:84,85:103]))
colnames(hcc5_pt1_pt2) <- "cor"
hcc5_pt1_pt2$pair <- "hcc5_pt1_pt2"
hcc5_pt1_pt2$type <- "cross"
hcc5_pt1_pt2$patient <- "hcc4"
hcc5_pt1_pt3 <- as.data.frame(unlist(hcc5_cor[78:84,104:144]))
colnames(hcc5_pt1_pt3) <- "cor"
hcc5_pt1_pt3$pair <- "hcc5_pt1_pt3"
hcc5_pt1_pt3$type <- "cross"
hcc5_pt1_pt3$patient <- "hcc5"
hcc5_pt1_pt4 <- as.data.frame(unlist(hcc5_cor[78:84,145:169]))
colnames(hcc5_pt1_pt4) <- "cor"
hcc5_pt1_pt4$pair <- "hcc5_pt1_pt4"
hcc5_pt1_pt4$type <- "cross"
hcc5_pt1_pt4$patient <- "hcc5"
hcc5_pt1_pt5 <- as.data.frame(unlist(hcc5_cor[78:84,170:181]))
colnames(hcc5_pt1_pt5) <- "cor"
hcc5_pt1_pt5$pair <- "hcc4_pt1_pt2"
hcc5_pt1_pt5$type <- "cross"
hcc5_pt1_pt5$patient <- "hcc5"
hcc5_pt2_pt3 <- as.data.frame(unlist(hcc5_cor[85:103,104:144]))
colnames(hcc5_pt2_pt3) <- "cor"
hcc5_pt2_pt3$pair <- "hcc5_pt2_pt3"
hcc5_pt2_pt3$type <- "cross"
hcc5_pt2_pt3$patient <- "hcc5"
hcc5_pt2_pt4 <- as.data.frame(unlist(hcc5_cor[85:103,145:169]))
colnames(hcc5_pt2_pt4) <- "cor"
hcc5_pt2_pt4$pair <- "hcc5_pt2_pt4"
hcc5_pt2_pt4$type <- "cross"
hcc5_pt2_pt4$patient <- "hcc5"
hcc5_pt2_pt5 <- as.data.frame(unlist(hcc5_cor[85:103,170:181]))
colnames(hcc5_pt2_pt5) <- "cor"
hcc5_pt2_pt5$pair <- "hcc5_pt2_pt5"
hcc5_pt2_pt5$type <- "cross"
hcc5_pt2_pt5$patient <- "hcc5"
hcc5_pt3_pt4 <- as.data.frame(unlist(hcc5_cor[104:144,145:169]))
colnames(hcc5_pt3_pt4) <- "cor"
hcc5_pt3_pt4$pair <- "hcc5_pt3_pt4" 
hcc5_pt3_pt4$type <- "cross"
hcc5_pt3_pt4$patient <- "hcc5"
hcc5_pt3_pt5 <- as.data.frame(unlist(hcc5_cor[104:144,170:181]))
colnames(hcc5_pt3_pt5) <- "cor"
hcc5_pt3_pt5$pair <- "hcc5_pt3_pt5"
hcc5_pt3_pt5$type <- "cross"
hcc5_pt3_pt5$patient <- "hcc5"
hcc5_pt4_pt5 <- as.data.frame(unlist(hcc5_cor[145:169,170:181]))
colnames(hcc5_pt4_pt5) <- "cor"
hcc5_pt4_pt5$pair <- "hcc5_pt4_pt5"
hcc5_pt4_pt5$type <- "cross"
hcc5_pt4_pt5$patient <- "hcc5"


hcc5_pt1_pt1 <- as.data.frame(unlist(hcc5_cor[78:84,78:84]))
colnames(hcc5_pt1_pt1) <- "cor"
hcc5_pt1_pt1$pair <- "hcc5_pt1_pt1"
hcc5_pt1_pt1$type <- "self"
hcc5_pt1_pt1$patient <- "hcc5"
hcc5_pt2_pt2 <- as.data.frame(unlist(hcc5_cor[85:103,85:103]))
colnames(hcc5_pt2_pt2) <- "cor"
hcc5_pt2_pt2$pair <- "hcc5_pt2_pt2"
hcc5_pt2_pt2$type <- "self"
hcc5_pt2_pt2$patient <- "hcc5"
hcc5_pt3_pt3 <- as.data.frame(unlist(hcc5_cor[104:144,104:144]))
colnames(hcc5_pt3_pt3) <- "cor"
hcc5_pt3_pt3$pair <- "hcc5_pt3_pt3"
hcc5_pt3_pt3$type <- "self"
hcc5_pt3_pt3$patient <- "hcc5"
hcc5_pt4_pt4 <- as.data.frame(unlist(hcc5_cor[145:169,145:169]))
colnames(hcc5_pt4_pt4) <- "cor"
hcc5_pt4_pt4$pair <- "hcc5_pt4_pt4"
hcc5_pt4_pt4$type <- "self"
hcc5_pt4_pt4$patient <- "hcc5"
hcc5_pt5_pt5 <- as.data.frame(unlist(hcc5_cor[170:181,170:181]))
colnames(hcc5_pt5_pt5) <- "cor"
hcc5_pt5_pt5$pair <- "hcc5_pt5_pt5"
hcc5_pt5_pt5$type <- "self"
hcc5_pt5_pt5$patient <- "hcc5"


hcc5_cor_all <- as.data.frame(rbind(hcc5_pt1_pt2, hcc5_pt1_pt3, hcc5_pt1_pt4, hcc5_pt1_pt5, hcc5_pt2_pt3, hcc5_pt2_pt4,
                                    hcc5_pt2_pt5, hcc5_pt3_pt4, hcc5_pt3_pt5, hcc5_pt4_pt5, hcc5_pt1_pt1, hcc5_pt2_pt2,
                                    hcc5_pt3_pt3, hcc5_pt4_pt4, hcc5_pt5_pt5))


ggplot(hcc5_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())



hcc6_all <- as.data.frame(GetAssayData(object = HCC6_HPC, slot = "data"))
hcc6_top <- hcc6_all[_gene,]
hcc6_top <- na.omit(hcc6_top)
hcc6_cor <- as.data.frame(cor(hcc6_top,method = 'pearson'))
hcc6_anno <- as.data.frame(HCC6_HPC$orig.ident)
colnames(hcc6_anno) <- "origin"
hcc6_anno2 <- as.data.frame(HCC6_HPC$satellite_region)
colnames(hcc6_anno2) <- "origin"
pheatmap(hcc6_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc6_anno,annotation_col = hcc6_anno,
         cluster_rows = F,cluster_cols = F)
table(HCC6_HPC$orig.ident)

hcc6_pt1_pt2_mean <- mean(unlist(hcc6_cor[19:115,116:169]))
hcc6_pt1_pt4_mean <- mean(unlist(hcc6_cor[19:115,170:298]))
hcc6_pt1_pt5_mean <- mean(unlist(hcc6_cor[19:115,299:429]))
hcc6_pt1_pt6_mean <- mean(unlist(hcc6_cor[19:115,430:434]))
hcc6_pt2_pt4_mean <- mean(unlist(hcc6_cor[116:169,170:298]))
hcc6_pt2_pt5_mean <- mean(unlist(hcc6_cor[116:169,299:429]))
hcc6_pt2_pt6_mean <- mean(unlist(hcc6_cor[116:169,104:144]))
hcc6_pt4_pt5_mean <- mean(unlist(hcc6_cor[104:144,170:298]))
hcc6_pt4_pt6_mean <- mean(unlist(hcc6_cor[104:144,430:434]))
hcc6_pt5_pt6_mean <- mean(unlist(hcc6_cor[170:298,430:434]))



hcc6_pt1_pt1_mean <- mean(unlist(hcc6_cor[19:115,19:115]))
hcc6_pt2_pt2_mean <- mean(unlist(hcc6_cor[116:169,116:169]))
hcc6_pt4_pt4_mean <- mean(unlist(hcc6_cor[104:144,104:144]))
hcc6_pt5_pt5_mean <- mean(unlist(hcc6_cor[170:298,170:298]))
hcc6_pt6_pt6_mean <- mean(unlist(hcc6_cor[430:434,430:434]))




hcc6_pt1_pt2 <- as.data.frame(unlist(hcc6_cor[19:115,116:169]))
colnames(hcc6_pt1_pt2) <- "cor"
hcc6_pt1_pt2$pair <- "hcc6_pt1_pt2"
hcc6_pt1_pt2$type <- "cross"
hcc6_pt1_pt2$patient <- "hcc6"
hcc6_pt1_pt4 <- as.data.frame(unlist(hcc6_cor[19:115,170:298]))
colnames(hcc6_pt1_pt4) <- "cor"
hcc6_pt1_pt4$pair <- "hcc6_pt1_pt4"
hcc6_pt1_pt4$type <- "cross"
hcc6_pt1_pt4$patient <- "hcc6"
hcc6_pt1_pt5 <- as.data.frame(unlist(hcc6_cor[19:115,299:429]))
colnames(hcc6_pt1_pt5) <- "cor"
hcc6_pt1_pt5$pair <- "hcc6_pt1_pt5"
hcc6_pt1_pt5$type <- "cross"
hcc6_pt1_pt5$patient <- "hcc6"
hcc6_pt1_pt6 <- as.data.frame(unlist(hcc6_cor[19:115,430:434]))
colnames(hcc6_pt1_pt6) <- "cor"
hcc6_pt1_pt6$pair <- "hcc6_pt1_pt6"
hcc6_pt1_pt6$type <- "cross"
hcc6_pt1_pt6$patient <- "hcc6"
hcc6_pt2_pt4 <- as.data.frame(unlist(hcc6_cor[116:169,170:298]))
colnames(hcc6_pt2_pt4) <- "cor"
hcc6_pt2_pt4$pair <- "hcc6_pt2_pt4"
hcc6_pt2_pt4$type <- "cross"
hcc6_pt2_pt4$patient <- "hcc6"
hcc6_pt2_pt5 <- as.data.frame(unlist(hcc6_cor[116:169,299:429]))
colnames(hcc6_pt2_pt5) <- "cor"
hcc6_pt2_pt5$pair <- "hcc6_pt2_pt5"
hcc6_pt2_pt5$type <- "cross"
hcc6_pt2_pt5$patient <- "hcc6"
hcc6_pt2_pt6 <- as.data.frame(unlist(hcc6_cor[116:169,104:144]))
colnames(hcc6_pt2_pt6) <- "cor"
hcc6_pt2_pt6$pair <- "hcc6_pt2_pt6"
hcc6_pt2_pt6$type <- "cross"
hcc6_pt2_pt6$patient <- "hcc6"
hcc6_pt4_pt5 <- as.data.frame(unlist(hcc6_cor[104:144,170:298]))
colnames(hcc6_pt4_pt5) <- "cor"
hcc6_pt4_pt5$pair <- "hcc6_pt4_pt5"
hcc6_pt4_pt5$type <- "cross"
hcc6_pt4_pt5$patient <- "hcc6"
hcc6_pt4_pt6 <- as.data.frame(unlist(hcc6_cor[104:144,430:434]))
colnames(hcc6_pt4_pt6) <- "cor"
hcc6_pt4_pt6$pair <- "hcc6_pt4_pt6"
hcc6_pt4_pt6$type <- "cross"
hcc6_pt4_pt6$patient <- "hcc6"
hcc6_pt5_pt6 <- as.data.frame(unlist(hcc6_cor[170:298,430:434]))
colnames(hcc6_pt5_pt6) <- "cor"
hcc6_pt5_pt6$pair <- "hcc6_pt5_pt6"
hcc6_pt5_pt6$type <- "cross"
hcc6_pt5_pt6$patient <- "hcc6"



hcc6_pt1_pt1 <- as.data.frame(unlist(hcc6_cor[19:115,19:115]))
colnames(hcc6_pt1_pt1) <- "cor"
hcc6_pt1_pt1$pair <- "hcc6_pt1_pt1"
hcc6_pt1_pt1$type <- "self"
hcc6_pt1_pt1$patient <- "hcc6"
hcc6_pt2_pt2 <- as.data.frame(unlist(hcc6_cor[116:169,116:169]))
colnames(hcc6_pt2_pt2) <- "cor"
hcc6_pt2_pt2$pair <- "hcc6_pt2_pt2"
hcc6_pt2_pt2$type <- "self"
hcc6_pt2_pt2$patient <- "hcc6"
hcc6_pt4_pt4 <- as.data.frame(unlist(hcc6_cor[104:144,104:144]))
colnames(hcc6_pt4_pt4) <- "cor"
hcc6_pt4_pt4$pair <- "hcc6_pt4_pt4"
hcc6_pt4_pt4$type <- "self"
hcc6_pt4_pt4$patient <- "hcc6"
hcc6_pt5_pt5 <- as.data.frame(unlist(hcc6_cor[170:298,170:298]))
colnames(hcc6_pt5_pt5) <- "cor"
hcc6_pt5_pt5$pair <- "hcc6_pt5_pt5"
hcc6_pt5_pt5$type <- "cross"
hcc6_pt5_pt5$patient <- "hcc6"
hcc6_pt6_pt6 <- as.data.frame(unlist(hcc6_cor[430:434,430:434]))
colnames(hcc6_pt6_pt6) <- "cor"
hcc6_pt6_pt6$pair <- "hcc6_pt6_pt6"
hcc6_pt6_pt6$type <- "self"
hcc6_pt6_pt6$patient <- "hcc6"

hcc6_cor_all <- as.data.frame(rbind(hcc6_pt1_pt2, hcc6_pt1_pt4, hcc6_pt1_pt5, hcc6_pt1_pt6, hcc6_pt2_pt4, hcc6_pt2_pt5,
                                    hcc6_pt2_pt6, hcc6_pt4_pt5, hcc6_pt4_pt6, hcc6_pt5_pt6, hcc6_pt1_pt1, hcc6_pt2_pt2,
                                    hcc6_pt4_pt4, hcc6_pt5_pt5, hcc6_pt6_pt6))





ggplot(hcc6_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


hcc7_all <- as.data.frame(GetAssayData(object = HCC7_HPC, slot = "data"))
hcc7_top <- hcc7_all[_gene,]
hcc7_top <- na.omit(hcc7_top)
hcc7_cor <- as.data.frame(cor(hcc7_top,method = 'pearson'))
hcc7_anno <- as.data.frame(HCC7_HPC$orig.ident)
colnames(hcc7_anno) <- "origin"
pheatmap(hcc7_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc7_anno,annotation_col = hcc7_anno,
         cluster_rows = F,cluster_cols = F)
table(HCC7_HPC$orig.ident)

hcc7_pt1_pt2_mean <- mean(unlist(hcc7_cor[11:42,43:58]))
hcc7_pt1_pt3_mean <- mean(unlist(hcc7_cor[11:42,59:116]))
hcc7_pt1_pt4_mean <- mean(unlist(hcc7_cor[11:42,117:159]))
hcc7_pt1_pt5_mean <- mean(unlist(hcc7_cor[11:42,160:174]))
hcc7_pt1_pt6_mean <- mean(unlist(hcc7_cor[11:42,175:187]))
hcc7_pt2_pt3_mean <- mean(unlist(hcc7_cor[43:58,59:116]))
hcc7_pt2_pt4_mean <- mean(unlist(hcc7_cor[43:58,117:159]))
hcc7_pt2_pt5_mean <- mean(unlist(hcc7_cor[43:58,160:174]))
hcc7_pt2_pt6_mean <- mean(unlist(hcc7_cor[43:58,175:187]))
hcc7_pt3_pt4_mean <- mean(unlist(hcc7_cor[59:116,117:159]))
hcc7_pt3_pt5_mean <- mean(unlist(hcc7_cor[59:116,160:174]))
hcc7_pt3_pt6_mean <- mean(unlist(hcc7_cor[59:116,175:187]))
hcc7_pt4_pt5_mean <- mean(unlist(hcc7_cor[117:159,160:174]))
hcc7_pt4_pt6_mean <- mean(unlist(hcc7_cor[117:159,175:187]))
hcc7_pt5_pt6_mean <- mean(unlist(hcc7_cor[160:174,175:187]))


hcc7_pt1_pt1_mean <- mean(unlist(hcc7_cor[11:42,11:42]))
hcc7_pt2_pt2_mean <- mean(unlist(hcc7_cor[43:58,43:58]))
hcc7_pt3_pt3_mean <- mean(unlist(hcc7_cor[59:116,59:116]))
hcc7_pt4_pt4_mean <- mean(unlist(hcc7_cor[117:159,117:159]))
hcc7_pt5_pt5_mean <- mean(unlist(hcc7_cor[160:174,160:174]))
hcc7_pt6_pt6_mean <- mean(unlist(hcc7_cor[175:187,175:187]))



hcc7_pt1_pt2 <- as.data.frame(unlist(hcc7_cor[11:42,43:58]))
colnames(hcc7_pt1_pt2) <- "cor"
hcc7_pt1_pt2$pair <- "hcc7_pt1_pt2"
hcc7_pt1_pt2$type <- "cross"
hcc7_pt1_pt2$patient <- "hcc7"
hcc7_pt1_pt3 <- as.data.frame(unlist(hcc7_cor[11:42,59:116]))
colnames(hcc7_pt1_pt3) <- "cor"
hcc7_pt1_pt3$pair <- "hcc7_pt1_pt3"
hcc7_pt1_pt3$type <- "cross"
hcc7_pt1_pt3$patient <- "hcc7"
hcc7_pt1_pt4 <- as.data.frame(unlist(hcc7_cor[11:42,117:159]))
colnames(hcc7_pt1_pt4) <- "cor"
hcc7_pt1_pt4$pair <- "hcc7_pt1_pt4"
hcc7_pt1_pt4$type <- "cross"
hcc7_pt1_pt4$patient <- "hcc7"
hcc7_pt1_pt5 <- as.data.frame(unlist(hcc7_cor[11:42,160:174]))
colnames(hcc7_pt1_pt5) <- "cor"
hcc7_pt1_pt5$pair <- "hcc7_pt1_pt5"
hcc7_pt1_pt5$type <- "cross"
hcc7_pt1_pt5$patient <- "hcc7"
hcc7_pt1_pt6 <- as.data.frame(unlist(hcc7_cor[11:42,175:187]))
colnames(hcc7_pt1_pt6) <- "cor"
hcc7_pt1_pt6$pair <- "hcc7_pt1_pt6"
hcc7_pt1_pt6$type <- "cross"
hcc7_pt1_pt6$patient <- "hcc7"
hcc7_pt2_pt3 <- as.data.frame(unlist(hcc7_cor[43:58,59:116]))
colnames(hcc7_pt2_pt3) <- "cor"
hcc7_pt2_pt3$pair <- "hcc7_pt2_pt3"
hcc7_pt2_pt3$type <- "cross"
hcc7_pt2_pt3$patient <- "hcc7"
hcc7_pt2_pt4 <- as.data.frame(unlist(hcc7_cor[43:58,117:159]))
colnames(hcc7_pt2_pt4) <- "cor"
hcc7_pt2_pt4$pair <- "hcc7_pt2_pt4"
hcc7_pt2_pt4$type <- "cross"
hcc7_pt2_pt4$patient <- "hcc7"
hcc7_pt2_pt5 <- as.data.frame(unlist(hcc7_cor[43:58,160:174]))
colnames(hcc7_pt2_pt5) <- "cor"
hcc7_pt2_pt5$pair <- "hcc7_pt2_pt5"
hcc7_pt2_pt5$type <- "cross"
hcc7_pt2_pt5$patient <- "hcc7"
hcc7_pt2_pt6 <- as.data.frame(unlist(hcc7_cor[43:58,175:187]))
colnames(hcc7_pt2_pt6) <- "cor"
hcc7_pt2_pt6$pair <- "hcc7_pt2_pt6"
hcc7_pt2_pt6$type <- "cross"
hcc7_pt2_pt6$patient <- "hcc7"
hcc7_pt3_pt4 <- as.data.frame(unlist(hcc7_cor[59:116,117:159]))
colnames(hcc7_pt3_pt4) <- "cor"
hcc7_pt3_pt4$pair <- "hcc7_pt3_pt4"
hcc7_pt3_pt4$type <- "cross"
hcc7_pt3_pt4$patient <- "hcc7"
hcc7_pt3_pt5 <- as.data.frame(unlist(hcc7_cor[59:116,160:174]))
colnames(hcc7_pt3_pt5) <- "cor"
hcc7_pt3_pt5$pair <- "hcc7_pt3_pt5"
hcc7_pt3_pt5$type <- "cross"
hcc7_pt3_pt5$patient <- "hcc7"
hcc7_pt3_pt6 <- as.data.frame(unlist(hcc7_cor[59:116,175:187]))
colnames(hcc7_pt3_pt6) <- "cor"
hcc7_pt3_pt6$pair <- "hcc7_pt3_pt6"
hcc7_pt3_pt6$type <- "cross"
hcc7_pt3_pt6$patient <- "hcc7"
hcc7_pt4_pt5 <- as.data.frame(unlist(hcc7_cor[117:159,160:174]))
colnames(hcc7_pt4_pt5) <- "cor"
hcc7_pt4_pt5$pair <- "hcc7_pt4_pt5"
hcc7_pt4_pt5$type <- "cross"
hcc7_pt4_pt5$patient <- "hcc7"
hcc7_pt4_pt6 <- as.data.frame(unlist(hcc7_cor[117:159,175:187]))
colnames(hcc7_pt4_pt6) <- "cor"
hcc7_pt4_pt6$pair <- "hcc7_pt4_pt6"
hcc7_pt4_pt6$type <- "cross"
hcc7_pt4_pt6$patient <- "hcc7"
hcc7_pt5_pt6 <- as.data.frame(unlist(hcc7_cor[160:174,175:187]))
colnames(hcc7_pt5_pt6) <- "cor"
hcc7_pt5_pt6$pair <- "hcc7_pt5_pt6"
hcc7_pt5_pt6$type <- "cross"
hcc7_pt5_pt6$patient <- "hcc7"



hcc7_pt1_pt1 <- as.data.frame(unlist(hcc7_cor[11:42,11:42]))
colnames(hcc7_pt1_pt1) <- "cor"
hcc7_pt1_pt1$pair <- "hcc7_pt1_pt1"
hcc7_pt1_pt1$type <- "self"
hcc7_pt1_pt1$patient <- "hcc7"
hcc7_pt2_pt2 <- as.data.frame(unlist(hcc7_cor[43:58,43:58]))
colnames(hcc7_pt2_pt2) <- "cor"
hcc7_pt2_pt2$pair <- "hcc7_pt2_pt2"
hcc7_pt2_pt2$type <- "self"
hcc7_pt2_pt2$patient <- "hcc7"
hcc7_pt3_pt3 <- as.data.frame(unlist(hcc7_cor[59:116,59:116]))
colnames(hcc7_pt3_pt3) <- "cor"
hcc7_pt3_pt3$pair <- "hcc7_pt3_pt3"
hcc7_pt3_pt3$type <- "self"
hcc7_pt3_pt3$patient <- "hcc7"
hcc7_pt4_pt4 <- as.data.frame(unlist(hcc7_cor[117:159,117:159]))
colnames(hcc7_pt4_pt4) <- "cor"
hcc7_pt4_pt4$pair <- "hcc7_pt4_pt4"
hcc7_pt4_pt4$type <- "self"
hcc7_pt4_pt4$patient <- "hcc7"
hcc7_pt5_pt5 <- as.data.frame(unlist(hcc7_cor[160:174,160:174]))
colnames(hcc7_pt5_pt5) <- "cor"
hcc7_pt5_pt5$pair <- "hcc7_pt5_pt5"
hcc7_pt5_pt5$type <- "self"
hcc7_pt5_pt5$patient <- "hcc7"
hcc7_pt6_pt6 <- as.data.frame(unlist(hcc7_cor[175:187,175:187]))
colnames(hcc7_pt6_pt6) <- "cor"
hcc7_pt6_pt6$pair <- "hcc7_pt6_pt6"
hcc7_pt6_pt6$type <- "self"
hcc7_pt6_pt6$patient <- "hcc7"


hcc7_cor_all <- as.data.frame(rbind(hcc7_pt1_pt2, hcc7_pt1_pt3, hcc7_pt1_pt4, hcc7_pt1_pt5, hcc7_pt1_pt6, hcc7_pt2_pt3,
                                    hcc7_pt2_pt4, hcc7_pt2_pt5, hcc7_pt2_pt6, hcc7_pt3_pt4, hcc7_pt3_pt5, hcc7_pt3_pt6,
                                    hcc7_pt4_pt5, hcc7_pt4_pt6, hcc7_pt5_pt6,
                                    hcc7_pt1_pt1, hcc7_pt2_pt2, hcc7_pt3_pt3, hcc7_pt4_pt4, hcc7_pt5_pt5, hcc7_pt6_pt6))


ggplot(hcc7_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

hcc8_all <- as.data.frame(GetAssayData(object = HCC8_HPC, slot = "data"))
hcc8_top <- hcc8_all[_gene,]
hcc8_top <- na.omit(hcc8_top)
hcc8_cor <- as.data.frame(cor(hcc8_top,method = 'pearson'))
hcc8_anno <- as.data.frame(HCC8_HPC$orig.ident)
colnames(hcc8_anno) <- "origin"
pheatmap(hcc8_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc8_anno,annotation_col = hcc8_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC8_HPC$orig.ident)

hcc8_pt1_pt2_mean <- mean(unlist(hcc8_cor[78:1145,1146:2083]))
hcc8_pt1_pt4_mean <- mean(unlist(hcc8_cor[78:1145,2084:2114]))
hcc8_pt2_pt4_mean <- mean(unlist(hcc8_cor[1146:2083,2084:2114]))



hcc8_pt1_pt1_mean <- mean(unlist(hcc8_cor[78:1145,78:1145]))
hcc8_pt2_pt2_mean <- mean(unlist(hcc8_cor[1146:2083,1146:2083]))
hcc8_pt4_pt4_mean <- mean(unlist(hcc8_cor[2084:2114,2084:2114]))




hcc8_pt1_pt2 <- as.data.frame(unlist(hcc8_cor[78:1145,1146:2083]))
colnames(hcc8_pt1_pt2) <- "cor"
hcc8_pt1_pt2$pair <- "hcc8_pt1_pt2"
hcc8_pt1_pt2$type <- "cross"
hcc8_pt1_pt2$patient <- "hcc8"
hcc8_pt1_pt4 <- as.data.frame(unlist(hcc8_cor[78:1145,2084:2114]))
colnames(hcc8_pt1_pt4) <- "cor"
hcc8_pt1_pt4$pair <- "hcc8_pt1_pt4"
hcc8_pt1_pt4$type <- "cross"
hcc8_pt1_pt4$patient <- "hcc8"
hcc8_pt2_pt4 <- as.data.frame(unlist(hcc8_cor[1146:2083,2084:2114]))
colnames(hcc8_pt2_pt4) <- "cor"
hcc8_pt2_pt4$pair <- "hcc8_pt2_pt4"
hcc8_pt2_pt4$type <- "cross"
hcc8_pt2_pt4$patient <- "hcc8"


hcc8_pt1_pt1 <- as.data.frame(unlist(hcc8_cor[78:1145,78:1145]))
colnames(hcc8_pt1_pt1) <- "cor"
hcc8_pt1_pt1$pair <- "hcc8_pt1_pt1"
hcc8_pt1_pt1$type <- "self"
hcc8_pt1_pt1$patient <- "hcc8"
hcc8_pt2_pt2 <- as.data.frame(unlist(hcc8_cor[1146:2083,1146:2083]))
colnames(hcc8_pt2_pt2) <- "cor"
hcc8_pt2_pt2$pair <- "hcc8_pt2_pt2"
hcc8_pt2_pt2$type <- "self"
hcc8_pt2_pt2$patient <- "hcc8"
hcc8_pt4_pt4 <- as.data.frame(unlist(hcc8_cor[2084:2114,2084:2114]))
colnames(hcc8_pt4_pt4) <- "cor"
hcc8_pt4_pt4$pair <- "hcc8_pt4_pt4"
hcc8_pt4_pt4$type <- "self"
hcc8_pt4_pt4$patient <- "hcc8"




hcc8_cor_all <- as.data.frame(rbind(hcc8_pt1_pt2,hcc8_pt1_pt4,hcc8_pt2_pt4,
                                    hcc8_pt1_pt1,hcc8_pt2_pt2,hcc8_pt4_pt4))


ggplot(hcc8_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


hcc9_all <- as.data.frame(GetAssayData(object = HCC9_HPC, slot = "data"))
hcc9_top <- hcc9_all[_gene,]
hcc9_top <- na.omit(hcc9_top)
hcc9_cor <- as.data.frame(cor(hcc9_top,method = 'pearson'))
hcc9_anno <- as.data.frame(HCC9_HPC$orig.ident)
colnames(hcc9_anno) <- "origin"
pheatmap(hcc9_cor,show_rownames = F,show_colnames = F,
         annotation_row = hcc9_anno,annotation_col = hcc9_anno,
         cluster_rows = F,cluster_cols = F)

table(HCC9_HPC$orig.ident)

hcc9_pt1_pt3_mean <- mean(unlist(hcc9_cor[164:494,495:733]))
hcc9_pt1_pt4_mean <- mean(unlist(hcc9_cor[164:494,734:1148]))
hcc9_pt3_pt4_mean <- mean(unlist(hcc9_cor[495:733,734:1148]))


hcc9_pt1_pt1_mean <- mean(unlist(hcc9_cor[164:494,164:494]))
hcc9_pt3_pt3_mean <- mean(unlist(hcc9_cor[495:733,495:733]))
hcc9_pt4_pt4_mean <- mean(unlist(hcc9_cor[734:1148,734:1148]))


hcc9_pt1_pt3 <- as.data.frame(unlist(hcc9_cor[164:494,495:733]))
colnames(hcc9_pt1_pt3) <- "cor"
hcc9_pt1_pt3$pair <- "hcc9_pt1_pt3"
hcc9_pt1_pt3$type <- "cross"
hcc9_pt1_pt3$patient <- "hcc9"
hcc9_pt1_pt4 <- as.data.frame(unlist(hcc9_cor[164:494,734:1148]))
colnames(hcc9_pt1_pt4) <- "cor"
hcc9_pt1_pt4$pair <- "hcc9_pt1_pt4"
hcc9_pt1_pt4$type <- "cross"
hcc9_pt1_pt4$patient <- "hcc9"
hcc9_pt3_pt4 <- as.data.frame(unlist(hcc9_cor[495:733,734:1148]))
colnames(hcc9_pt3_pt4) <- "cor"
hcc9_pt3_pt4$pair <- "hcc9_pt3_pt4"
hcc9_pt3_pt4$type <- "cross"
hcc9_pt3_pt4$patient <- "hcc9"


hcc9_pt1_pt1 <- as.data.frame(unlist(hcc9_cor[164:494,164:494]))
colnames(hcc9_pt1_pt1) <- "cor"
hcc9_pt1_pt1$pair <- "hcc9_pt1_pt1"
hcc9_pt1_pt1$type <- "self"
hcc9_pt1_pt1$patient <- "hcc9"
hcc9_pt3_pt3 <- as.data.frame(unlist(hcc9_cor[495:733,495:733]))
colnames(hcc9_pt3_pt3) <- "cor"
hcc9_pt3_pt3$pair <- "hcc9_pt3_pt3"
hcc9_pt3_pt3$type <- "self"
hcc9_pt3_pt3$patient <- "hcc9"
hcc9_pt4_pt4 <- as.data.frame(unlist(hcc9_cor[734:1148,734:1148]))
colnames(hcc9_pt4_pt4) <- "cor"
hcc9_pt4_pt4$pair <- "hcc9_pt4_pt4"
hcc9_pt4_pt4$type <- "self"
hcc9_pt4_pt4$patient <- "hcc9"




hcc9_cor_all <- as.data.frame(rbind(hcc9_pt1_pt3,hcc9_pt1_pt4,hcc9_pt3_pt4,
                                    hcc9_pt1_pt1,hcc9_pt3_pt3,hcc9_pt4_pt4))


ggplot(hcc9_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())



hcc1_cor_mean <- as.data.frame(rbind(hcc1_pt2_pt4_mean,hcc1_pt2_pt5_mean ,hcc1_pt4_pt5_mean,
                                     hcc1_pt2_pt2_mean,hcc1_pt4_pt4_mean,hcc1_pt5_pt5_mean))
colnames(hcc1_cor_mean) <- "cor"
hcc1_cor_mean$comb <- row.names(hcc1_cor_mean)
hcc1_cor_mean$sample <- "hcc1"
hcc1_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))





hcc2_cor_mean <- as.data.frame(rbind(hcc2_pt1_pt2_mean,hcc2_pt1_pt3_mean,hcc2_pt1_pt4_mean,
                                     hcc2_pt2_pt3_mean,hcc2_pt2_pt4_mean,hcc2_pt3_pt4_mean,
                                     hcc2_pt1_pt1_mean,hcc2_pt2_pt2_mean,hcc2_pt3_pt3_mean,
                                     hcc2_pt4_pt4_mean))
colnames(hcc2_cor_mean) <- "cor"
hcc2_cor_mean$comb <- row.names(hcc2_cor_mean)
hcc2_cor_mean$sample <- "hcc2"
hcc2_cor_mean$comb_type <- c(rep("cross",6),rep("self",4))





hcc3_cor_mean <- as.data.frame(rbind(hcc3_pt1_pt2_mean,hcc3_pt1_pt3_mean,hcc3_pt2_pt3_mean,
                                     hcc3_pt1_pt1_mean,hcc3_pt2_pt2_mean,hcc3_pt3_pt3_mean))
colnames(hcc3_cor_mean) <- "cor"
hcc3_cor_mean$comb <- row.names(hcc3_cor_mean)
hcc3_cor_mean$sample <- "hcc3"
hcc3_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))




hcc4_cor_mean <- as.data.frame(rbind(hcc4_pt1_pt2_mean,hcc4_pt1_pt3_mean,hcc4_pt2_pt3_mean,
                                     hcc4_pt1_pt1_mean,hcc4_pt2_pt2_mean,hcc4_pt3_pt3_mean))
colnames(hcc4_cor_mean) <- "cor"
hcc4_cor_mean$comb <- row.names(hcc4_cor_mean)
hcc4_cor_mean$sample <- "hcc4"
hcc4_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))






hcc5_cor_mean <- as.data.frame(rbind(hcc5_pt1_pt2_mean,hcc5_pt1_pt3_mean,hcc5_pt1_pt4_mean,hcc5_pt1_pt5_mean,
                                     hcc5_pt2_pt3_mean,hcc5_pt2_pt4_mean,hcc5_pt2_pt5_mean,hcc5_pt3_pt4_mean,
                                     hcc5_pt3_pt5_mean,hcc5_pt4_pt5_mean,
                                     hcc5_pt1_pt1_mean,hcc5_pt2_pt2_mean,hcc5_pt3_pt3_mean,hcc5_pt4_pt4_mean,
                                     hcc5_pt5_pt5_mean))
colnames(hcc5_cor_mean) <- "cor"
hcc5_cor_mean$comb <- row.names(hcc5_cor_mean)
hcc5_cor_mean$sample <- "hcc5"
hcc5_cor_mean$comb_type <- c(rep("cross",10),rep("self",5))





hcc6_cor_mean <- as.data.frame(rbind(hcc6_pt1_pt2_mean,hcc6_pt1_pt6_mean,hcc6_pt1_pt4_mean,hcc6_pt1_pt5_mean,
                                     hcc6_pt2_pt6_mean,hcc6_pt2_pt4_mean,hcc6_pt2_pt5_mean,hcc6_pt4_pt6_mean,
                                     hcc6_pt5_pt6_mean,hcc6_pt4_pt5_mean,
                                     hcc6_pt1_pt1_mean,hcc6_pt2_pt2_mean,hcc6_pt6_pt6_mean,hcc6_pt4_pt4_mean,
                                     hcc6_pt5_pt5_mean))
colnames(hcc6_cor_mean) <- "cor"
hcc6_cor_mean$comb <- row.names(hcc6_cor_mean)
hcc6_cor_mean$sample <- "hcc6"
hcc6_cor_mean$comb_type <- c(rep("cross",10),rep("self",5))



hcc7_cor_mean <- as.data.frame(rbind(hcc7_pt1_pt2_mean,hcc7_pt1_pt3_mean,hcc7_pt1_pt4_mean,hcc7_pt1_pt5_mean,
                                     hcc7_pt1_pt6_mean,hcc7_pt2_pt3_mean,hcc7_pt2_pt4_mean,hcc7_pt2_pt5_mean,
                                     hcc7_pt2_pt6_mean,hcc7_pt3_pt4_mean,hcc7_pt3_pt5_mean,hcc7_pt3_pt6_mean,
                                     hcc7_pt4_pt5_mean,hcc7_pt4_pt6_mean,hcc7_pt5_pt6_mean,
                                     hcc7_pt1_pt1_mean,hcc7_pt2_pt2_mean,hcc7_pt3_pt3_mean,hcc7_pt4_pt4_mean,
                                     hcc7_pt5_pt5_mean,hcc7_pt6_pt6_mean))
colnames(hcc7_cor_mean) <- "cor"
hcc7_cor_mean$comb <- row.names(hcc7_cor_mean)
hcc7_cor_mean$sample <- "hcc7"
hcc7_cor_mean$comb_type <- c(rep("cross",15),rep("self",6))


hcc8_cor_mean <- as.data.frame(rbind(hcc8_pt1_pt2_mean,hcc8_pt1_pt4_mean,hcc8_pt2_pt4_mean,
                                     hcc8_pt1_pt1_mean,hcc8_pt2_pt2_mean,hcc8_pt4_pt4_mean))
colnames(hcc8_cor_mean) <- "cor"
hcc8_cor_mean$comb <- row.names(hcc8_cor_mean)
hcc8_cor_mean$sample <- "hcc8"
hcc8_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))


hcc9_cor_mean <- as.data.frame(rbind(hcc9_pt1_pt4_mean,hcc9_pt1_pt3_mean,hcc9_pt3_pt4_mean,
                                     hcc9_pt1_pt1_mean,hcc9_pt4_pt4_mean,hcc9_pt3_pt3_mean))
colnames(hcc9_cor_mean) <- "cor"
hcc9_cor_mean$comb <- row.names(hcc9_cor_mean)
hcc9_cor_mean$sample <- "hcc9"
hcc9_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))




sc_cor_sum <- rbind(hcc1_cor_mean,hcc2_cor_mean,hcc3_cor_mean,
                    hcc4_cor_mean,hcc5_cor_mean,hcc6_cor_mean,
                    hcc7_cor_mean,hcc8_cor_mean,hcc9_cor_mean)


ggplot(sc_cor_sum,aes(x=sample,y=cor,fill=comb_type,color=comb_type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")


VlnPlot(HPC,features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,split.by = "lib.method",group.by = "patient")






ggplot(sc_cor_sum,aes(comb_type,cor,color=comb_type))+
  geom_boxplot(width=0.5)+
  geom_point()+
  geom_jitter(width = 0.25)+
  facet_grid(~sample)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


sc_cor_all <- rbind(hcc1_cor_all,hcc2_cor_all,hcc3_cor_all,
                    hcc4_cor_all,hcc5_cor_all,hcc6_cor_all,
                    hcc7_cor_all,hcc8_cor_all,hcc9_cor_all)

ggplot(sc_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~patient)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(sc_cor_all,aes(type,cor,color=type))+
  geom_violin(width=0.5,outlier.size=0)+
  facet_grid(~patient)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

hcc1_cross_mean <- mean(subset(hcc1_cor_all,subset=type=="cross")$cor)
hcc1_self_mean <- mean(subset(hcc1_cor_all,subset=type=="self")$cor)
hcc2_cross_mean <- mean(subset(hcc2_cor_all,subset=type=="cross")$cor)
hcc2_self_mean <- mean(subset(hcc2_cor_all,subset=type=="self")$cor)
hcc3_cross_mean <- mean(subset(hcc3_cor_all,subset=type=="cross")$cor)
hcc3_self_mean <- mean(subset(hcc3_cor_all,subset=type=="self")$cor)
hcc4_cross_mean <- mean(subset(hcc4_cor_all,subset=type=="cross")$cor)
hcc4_self_mean <- mean(subset(hcc4_cor_all,subset=type=="self")$cor)
hcc5_cross_mean <- mean(subset(hcc5_cor_all,subset=type=="cross")$cor)
hcc5_self_mean <- mean(subset(hcc5_cor_all,subset=type=="self")$cor)
hcc6_cross_mean <- mean(subset(hcc6_cor_all,subset=type=="cross")$cor)
hcc6_self_mean <- mean(subset(hcc6_cor_all,subset=type=="self")$cor)
hcc7_cross_mean <- mean(subset(hcc7_cor_all,subset=type=="cross")$cor)
hcc7_self_mean <- mean(subset(hcc7_cor_all,subset=type=="self")$cor)
hcc8_cross_mean <- mean(subset(hcc8_cor_all,subset=type=="cross")$cor)
hcc8_self_mean <- mean(subset(hcc8_cor_all,subset=type=="self")$cor)
hcc9_cross_mean <- mean(subset(hcc9_cor_all,subset=type=="cross")$cor)
hcc9_self_mean <- mean(subset(hcc9_cor_all,subset=type=="self")$cor)

hcc1_mean <- as.data.frame(abs(hcc1_self_mean-hcc1_cross_mean))
colnames(hcc1_mean) <- "cor_delta"
hcc1_mean$patient <- "hcc1"
hcc1_mean$type <- "CMN"


hcc2_mean <- as.data.frame(abs(hcc2_self_mean-hcc2_cross_mean))
colnames(hcc2_mean) <- "cor_delta"
hcc2_mean$patient <- "hcc2"
hcc2_mean$type <- "CMN"

hcc3_mean <- as.data.frame(abs(hcc3_self_mean-hcc3_cross_mean))
colnames(hcc3_mean) <- "cor_delta"
hcc3_mean$patient <- "hcc3"
hcc3_mean$type <- "SN"


hcc4_mean <- as.data.frame(abs(hcc4_self_mean-hcc4_cross_mean))
colnames(hcc4_mean) <- "cor_delta"
hcc4_mean$patient <- "hcc4"
hcc4_mean$type <- "CMN"


hcc5_mean <- as.data.frame(abs(hcc5_self_mean-hcc5_cross_mean))
colnames(hcc5_mean) <- "cor_delta"
hcc5_mean$patient <- "hcc5"
hcc5_mean$type <- "SN"

hcc6_mean <- as.data.frame(abs(hcc6_self_mean-hcc6_cross_mean))
colnames(hcc6_mean) <- "cor_delta"
hcc6_mean$patient <- "hcc6"
hcc6_mean$type <- "CMN"


hcc7_mean <- as.data.frame(abs(hcc7_self_mean-hcc7_cross_mean))
colnames(hcc7_mean) <- "cor_delta"
hcc7_mean$patient <- "hcc7"
hcc7_mean$type <- "SN"


hcc8_mean <- as.data.frame(abs(hcc8_self_mean-hcc8_cross_mean))
colnames(hcc8_mean) <- "cor_delta"
hcc8_mean$patient <- "hcc8"
hcc8_mean$type <- "CMN"



hcc9_mean <- as.data.frame(abs(hcc9_self_mean-hcc9_cross_mean))
colnames(hcc9_mean) <- "cor_delta"
hcc9_mean$patient <- "hcc9"
hcc9_mean$type <- "CMN"

sc_mean_sum <- rbind(hcc1_mean,hcc2_mean,hcc3_mean,
                     hcc4_mean,hcc5_mean,hcc6_mean,
                     hcc7_mean,hcc8_mean,hcc9_mean)


ggplot(sc_mean_sum,aes(type,cor_delta,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  geom_jitter(width=0.25)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("CMN","SN")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =0.2)+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())




ggplot(sc_mean_sum,aes(type,cor_delta,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  geom_jitter(width=0.25)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("CMN","SN")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =0.2)+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())



ggplot(sc_cor_all,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~patient)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

sc_cor_cmn <-  subset(sc_cor_all,subset=patient %in% c("hcc1","hcc2","hcc4","hcc6","hcc8","hcc9"))
sc_cor_sn <-  subset(sc_cor_all,subset=patient %in% c("hcc3","hcc5","hcc7"))

ggplot(sc_cor_cmn,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~patient)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(sc_cor_sn,aes(type,cor,fill=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~patient)+
  scale_fill_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

HPC <- FindVariableFeatures(HPC, selection.method = "vst", nfeatures = 3000)
top3000_gene <- head(VariableFeatures(HPC), 3000)





hcc1_pt2 <- subset(HCC1_HPC, subset = orig.ident=="PT2")
hcc1_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc1_pt2, slot = "data"))
hcc1_pt2_top_mtx <- hcc1_pt2_mtx[_gene,]
hcc1_pt2_top_mtx <- na.omit(hcc1_pt2_top_mtx)
hcc1_pt2_cor <- as.data.frame(cor(hcc1_pt2_top_mtx,method = 'pearson'))


for(i in 1:nrow(hcc1_pt2_cor)){
  vector <- as.vector(hcc1_pt2_cor[i,1:(nrow(hcc1_pt2_cor)-i+1)])
  if(i == 1){
    hcc1_pt2_cor_nu <- vector
  }
  else{
    hcc1_pt2_cor_nu <- c(hcc1_pt2_cor_nu,vector)
  }
}



hcc1_pt4 <- subset(HCC1_HPC, subset = orig.ident=="PT4")
hcc1_pt4_mtx <-  as.data.frame(GetAssayData(object = hcc1_pt4, slot = "data"))
hcc1_pt4_top_mtx <- hcc1_pt4_mtx[top3000_gene,]
hcc1_pt4_top_mtx <- na.omit(hcc1_pt4_top_mtx)
hcc1_pt4_cor <- as.data.frame(cor(hcc1_pt4_top_mtx,method = 'pearson'))
for(i in 1:nrow(hcc1_pt4_cor)){
  vector <- as.vector(hcc1_pt4_cor[i,1:(nrow(hcc1_pt4_cor)-i+1)])
  if(i == 1){
    hcc1_pt4_cor_nu <- vector
  }
  else{
    hcc1_pt4_cor_nu <- c(hcc1_pt4_cor_nu,vector)
  }
}

hcc1_pt5 <- subset(HCC1_HPC, subset = orig.ident=="PT5")
hcc1_pt5_mtx <-  as.data.frame(GetAssayData(object = hcc1_pt5, slot = "data"))
hcc1_pt5_top_mtx <- hcc1_pt5_mtx[top3000_gene,]
hcc1_pt5_top_mtx <- na.omit(hcc1_pt5_top_mtx)
hcc1_pt5_cor <- as.data.frame(cor(hcc1_pt5_top_mtx,method = 'pearson'))
for(i in 1:nrow(hcc1_pt5_cor)){
  vector <- as.vector(hcc1_pt5_cor[i,1:(nrow(hcc1_pt5_cor)-i+1)])
  if(i == 1){
    hcc1_pt5_cor_nu <- vector
  }
  else{
    hcc1_pt5_cor_nu <- c(hcc1_pt5_cor_nu,vector)
  }
}


hcc1_primary_cor3 <- as.data.frame(c(hcc1_pt4_cor_nu,hcc1_pt2_cor_nu))
hcc1_primary_cor3 <-as.data.frame(t(hcc1_primary_cor3))
colnames(hcc1_primary_cor3) <- "cor"
hcc1_primary_cor3$type <- "primary"
hcc1_satellite_cor3 <- as.data.frame(unlist(unique(hcc1_pt5_cor_nu)))
colnames(hcc1_satellite_cor3) <- "cor"
hcc1_satellite_cor3$type <- "satellite"
hcc1_ps_info2 <- rbind(hcc1_primary_cor3,hcc1_satellite_cor3)



ggplot(data=hcc1_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC1")







hcc3_pt1 <- subset(HCC3_HPC, subset = orig.ident=="PT1")
hcc3_pt1_mtx <-  as.data.frame(GetAssayData(object = hcc3_pt1, slot = "data"))
hcc3_pt1_top_mtx <- hcc3_pt1_mtx[top3000_gene,]
hcc3_pt1_top_mtx <- na.omit(hcc3_pt1_top_mtx)
hcc3_pt1_cor <- as.data.frame(cor(hcc3_pt1_top_mtx,method = 'pearson'))
hcc3_pt1_cor_nu <- unique(as.vector(as.matrix(hcc3_pt1_cor)))
for(i in 1:nrow(hcc3_pt1_cor)){
  vector <- as.vector(hcc3_pt1_cor[i,1:(nrow(hcc3_pt1_cor)-i+1)])
  if(i == 1){
    hcc3_pt1_cor_nu <- vector
  }
  else{
    hcc3_pt1_cor_nu <- c(hcc3_pt1_cor_nu,vector)
  }
}

hcc3_pt2 <- subset(HCC3_HPC, subset = orig.ident=="PT2")
hcc3_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc3_pt2, slot = "data"))
hcc3_pt2_top_mtx <- hcc3_pt2_mtx[top3000_gene,]
hcc3_pt2_top_mtx <- na.omit(hcc3_pt2_top_mtx)
hcc3_pt2_cor <- as.data.frame(cor(hcc3_pt2_top_mtx,method = 'pearson'))
hcc3_pt2_cor_nu <- unique(as.vector(as.matrix(hcc3_pt2_cor)))
for(i in 1:nrow(hcc3_pt2_cor)){
  vector <- as.vector(hcc3_pt2_cor[i,1:(nrow(hcc3_pt2_cor)-i+1)])
  if(i == 1){
    hcc3_pt2_cor_nu <- vector
  }
  else{
    hcc3_pt2_cor_nu <- c(hcc3_pt2_cor_nu,vector)
  }
}

hcc3_pt3 <- subset(HCC3_HPC, subset = orig.ident=="PT3")
hcc3_pt3_mtx <-  as.data.frame(GetAssayData(object = hcc3_pt3, slot = "data"))
hcc3_pt3_top_mtx <- hcc3_pt3_mtx[top3000_gene,]
hcc3_pt3_top_mtx <- na.omit(hcc3_pt3_top_mtx)
hcc3_pt3_cor <- as.data.frame(cor(hcc3_pt3_top_mtx,method = 'pearson'))
hcc3_pt3_cor_nu <- unique(as.vector(as.matrix(hcc3_pt3_cor)))
for(i in 1:nrow(hcc3_pt3_cor)){
  vector <- as.vector(hcc3_pt3_cor[i,1:(nrow(hcc3_pt3_cor)-i+1)])
  if(i == 1){
    hcc3_pt3_cor_nu <- vector
  }
  else{
    hcc3_pt3_cor_nu <- c(hcc3_pt3_cor_nu,vector)
  }
}

hcc3_primary_cor3 <- as.data.frame(c(hcc3_pt1_cor_nu,hcc3_pt2_cor_nu))
hcc3_primary_cor3 <-as.data.frame(t(hcc3_primary_cor3))
colnames(hcc3_primary_cor3) <- "cor"
hcc3_primary_cor3$type <- "primary"
hcc3_satellite_cor3 <- as.data.frame(unlist(unique(hcc3_pt3_cor_nu)))
colnames(hcc3_satellite_cor3) <- "cor"
hcc3_satellite_cor3$type <- "satellite"
hcc3_ps_info2 <- rbind(hcc3_primary_cor3,hcc3_satellite_cor3)
ggplot(data=hcc3_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC3")


######hcc5-pt1/pt2/pt3/pt4####pt5#####
hcc5_pt1 <- subset(HCC5_HPC, subset = orig.ident=="PT1")
hcc5_pt1_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt1, slot = "data"))
hcc5_pt1_top_mtx <- hcc5_pt1_mtx[top3000_gene,]
hcc5_pt1_top_mtx <- na.omit(hcc5_pt1_top_mtx)
hcc5_pt1_cor <- as.data.frame(cor(hcc5_pt1_top_mtx,method = 'pearson'))
hcc5_pt1_cor_nu <- unique(as.vector(as.matrix(hcc5_pt1_cor)))
for(i in 1:nrow(hcc5_pt1_cor)){
  vector <- as.vector(hcc5_pt1_cor[i,1:(nrow(hcc5_pt1_cor)-i+1)])
  if(i == 1){
    hcc5_pt1_cor_nu <- vector
  }
  else{
    hcc5_pt1_cor_nu <- c(hcc5_pt1_cor_nu,vector)
  }
}




hcc5_pt2 <- subset(HCC5_HPC, subset = orig.ident=="PT2")
hcc5_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt2, slot = "data"))
hcc5_pt2_top_mtx <- hcc5_pt2_mtx[top3000_gene,]
hcc5_pt2_top_mtx <- na.omit(hcc5_pt2_top_mtx)
hcc5_pt2_cor <- as.data.frame(cor(hcc5_pt2_top_mtx,method = 'pearson'))
hcc5_pt2_cor_nu <- unique(as.vector(as.matrix(hcc5_pt2_cor)))
for(i in 1:nrow(hcc5_pt2_cor)){
  vector <- as.vector(hcc5_pt2_cor[i,1:(nrow(hcc5_pt2_cor)-i+1)])
  if(i == 1){
    hcc5_pt2_cor_nu <- vector
  }
  else{
    hcc5_pt2_cor_nu <- c(hcc5_pt2_cor_nu,vector)
  }
}






hcc5_pt3 <- subset(HCC5_HPC, subset = orig.ident=="PT3")
hcc5_pt3_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt3, slot = "data"))
hcc5_pt3_top_mtx <- hcc5_pt3_mtx[top3000_gene,]
hcc5_pt3_top_mtx <- na.omit(hcc5_pt3_top_mtx)
hcc5_pt3_cor <- as.data.frame(cor(hcc5_pt3_top_mtx,method = 'pearson'))
hcc5_pt3_cor_nu <- unique(as.vector(as.matrix(hcc5_pt3_cor)))
for(i in 1:nrow(hcc5_pt3_cor)){
  vector <- as.vector(hcc5_pt3_cor[i,1:(nrow(hcc5_pt3_cor)-i+1)])
  if(i == 1){
    hcc5_pt3_cor_nu <- vector
  }
  else{
    hcc5_pt3_cor_nu <- c(hcc5_pt3_cor_nu,vector)
  }
}




hcc5_pt4 <- subset(HCC5_HPC, subset = orig.ident=="PT4")
hcc5_pt4_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt4, slot = "data"))
hcc5_pt4_top_mtx <- hcc5_pt4_mtx[top3000_gene,]
hcc5_pt4_top_mtx <- na.omit(hcc5_pt4_top_mtx)
hcc5_pt4_cor <- as.data.frame(cor(hcc5_pt4_top_mtx,method = 'pearson'))
hcc5_pt4_cor_nu <- unique(as.vector(as.matrix(hcc5_pt4_cor)))
for(i in 1:nrow(hcc5_pt4_cor)){
  vector <- as.vector(hcc5_pt4_cor[i,1:(nrow(hcc5_pt4_cor)-i+1)])
  if(i == 1){
    hcc5_pt4_cor_nu <- vector
  }
  else{
    hcc5_pt4_cor_nu <- c(hcc5_pt4_cor_nu,vector)
  }
}



hcc5_pt5 <- subset(HCC5_HPC, subset = orig.ident=="PT5")
hcc5_pt5_mtx <-  as.data.frame(GetAssayData(object = hcc5_pt5, slot = "data"))
hcc5_pt5_top_mtx <- hcc5_pt5_mtx[top3000_gene,]
hcc5_pt5_top_mtx <- na.omit(hcc5_pt5_top_mtx)
hcc5_pt5_cor <- as.data.frame(cor(hcc5_pt5_top_mtx,method = 'pearson'))
hcc5_pt5_cor_nu <- unique(as.vector(as.matrix(hcc5_pt5_cor)))
for(i in 1:nrow(hcc5_pt5_cor)){
  vector <- as.vector(hcc5_pt5_cor[i,1:(nrow(hcc5_pt5_cor)-i+1)])
  if(i == 1){
    hcc5_pt5_cor_nu <- vector
  }
  else{
    hcc5_pt5_cor_nu <- c(hcc5_pt5_cor_nu,vector)
  }
}



hcc5_primary_cor3 <- as.data.frame(c(hcc5_pt1_cor_nu,hcc5_pt2_cor_nu,hcc5_pt3_cor_nu,hcc5_pt4_cor_nu))
hcc5_primary_cor3 <-as.data.frame(t(hcc5_primary_cor3))
colnames(hcc5_primary_cor3) <- "cor"
hcc5_primary_cor3$type <- "primary"
hcc5_satellite_cor3 <- as.data.frame(unlist(unique(hcc5_pt5_cor_nu)))
colnames(hcc5_satellite_cor3) <- "cor"
hcc5_satellite_cor3$type <- "satellite"
hcc5_ps_info2 <- rbind(hcc5_primary_cor3,hcc5_satellite_cor3)
ggplot(data=hcc5_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC5")

######hcc6-pt1/pt2####pt4/pt5/pt6#####
hcc6_pt1 <- subset(HCC6_HPC, subset = orig.ident=="PT1")
hcc6_pt1_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt1, slot = "data"))
hcc6_pt1_top_mtx <- hcc6_pt1_mtx[top3000_gene,]
hcc6_pt1_top_mtx <- na.omit(hcc6_pt1_top_mtx)
hcc6_pt1_cor <- as.data.frame(cor(hcc6_pt1_top_mtx,method = 'pearson'))
hcc6_pt1_cor_nu <- unique(as.vector(as.matrix(hcc6_pt1_cor)))
for(i in 1:nrow(hcc6_pt1_cor)){
  vector <- as.vector(hcc6_pt1_cor[i,1:(nrow(hcc6_pt1_cor)-i+1)])
  if(i == 1){
    hcc6_pt1_cor_nu <- vector
  }
  else{
    hcc6_pt1_cor_nu <- c(hcc6_pt1_cor_nu,vector)
  }
}




hcc6_pt2 <- subset(HCC6_HPC, subset = orig.ident=="PT2")
hcc6_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt2, slot = "data"))
hcc6_pt2_top_mtx <- hcc6_pt2_mtx[top3000_gene,]
hcc6_pt2_top_mtx <- na.omit(hcc6_pt2_top_mtx)
hcc6_pt2_cor <- as.data.frame(cor(hcc6_pt2_top_mtx,method = 'pearson'))
hcc6_pt2_cor_nu <- unique(as.vector(as.matrix(hcc6_pt2_cor)))
for(i in 1:nrow(hcc6_pt2_cor)){
  vector <- as.vector(hcc6_pt2_cor[i,1:(nrow(hcc6_pt2_cor)-i+1)])
  if(i == 1){
    hcc6_pt2_cor_nu <- vector
  }
  else{
    hcc6_pt2_cor_nu <- c(hcc6_pt2_cor_nu,vector)
  }
}





hcc6_pt4 <- subset(HCC6_HPC, subset = orig.ident=="PT4")
hcc6_pt4_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt4, slot = "data"))
hcc6_pt4_top_mtx <- hcc6_pt4_mtx[top3000_gene,]
hcc6_pt4_top_mtx <- na.omit(hcc6_pt4_top_mtx)
hcc6_pt4_cor <- as.data.frame(cor(hcc6_pt4_top_mtx,method = 'pearson'))
hcc6_pt4_cor_nu <- unique(as.vector(as.matrix(hcc6_pt4_cor)))
for(i in 1:nrow(hcc6_pt4_cor)){
  vector <- as.vector(hcc6_pt4_cor[i,1:(nrow(hcc6_pt4_cor)-i+1)])
  if(i == 1){
    hcc6_pt4_cor_nu <- vector
  }
  else{
    hcc6_pt4_cor_nu <- c(hcc6_pt4_cor_nu,vector)
  }
}




hcc6_pt5 <- subset(HCC6_HPC, subset = orig.ident=="PT5")
hcc6_pt5_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt5, slot = "data"))
hcc6_pt5_top_mtx <- hcc6_pt5_mtx[top3000_gene,]
hcc6_pt5_top_mtx <- na.omit(hcc6_pt5_top_mtx)
hcc6_pt5_cor <- as.data.frame(cor(hcc6_pt5_top_mtx,method = 'pearson'))
hcc6_pt5_cor_nu <- unique(as.vector(as.matrix(hcc6_pt5_cor)))
for(i in 1:nrow(hcc6_pt5_cor)){
  vector <- as.vector(hcc6_pt5_cor[i,1:(nrow(hcc6_pt5_cor)-i+1)])
  if(i == 1){
    hcc6_pt5_cor_nu <- vector
  }
  else{
    hcc6_pt5_cor_nu <- c(hcc6_pt5_cor_nu,vector)
  }
}



hcc6_pt6 <- subset(HCC6_HPC, subset = orig.ident=="PT6")
hcc6_pt6_mtx <-  as.data.frame(GetAssayData(object = hcc6_pt6, slot = "data"))
hcc6_pt6_top_mtx <- hcc6_pt6_mtx[top3000_gene,]
hcc6_pt6_top_mtx <- na.omit(hcc6_pt6_top_mtx)
hcc6_pt6_cor <- as.data.frame(cor(hcc6_pt6_top_mtx,method = 'pearson'))
hcc6_pt6_cor_nu <- unique(as.vector(as.matrix(hcc6_pt6_cor)))
for(i in 1:nrow(hcc6_pt6_cor)){
  vector <- as.vector(hcc6_pt6_cor[i,1:(nrow(hcc6_pt6_cor)-i+1)])
  if(i == 1){
    hcc6_pt6_cor_nu <- vector
  }
  else{
    hcc6_pt6_cor_nu <- c(hcc6_pt6_cor_nu,vector)
  }
}






hcc6_primary_cor3 <- as.data.frame(c(hcc6_pt1_cor_nu,hcc6_pt2_cor_nu))
hcc6_primary_cor3 <-as.data.frame(t(hcc6_primary_cor3))
colnames(hcc6_primary_cor3) <- "cor"
hcc6_primary_cor3$type <- "primary"
hcc6_satellite_cor3 <- as.data.frame(c(hcc6_pt4_cor_nu,hcc6_pt5_cor_nu,hcc6_pt6_cor_nu))
hcc6_satellite_cor3 <-as.data.frame(t(hcc6_satellite_cor3))
colnames(hcc6_satellite_cor3) <- "cor"
hcc6_satellite_cor3$type <- "satellite"
hcc6_ps_info2 <- rbind(hcc6_primary_cor3,hcc6_satellite_cor3)
ggplot(data=hcc6_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC6")





######hcc7-pt1/pt2/pt3/pt4/pt5###pt6##

hcc7_pt1 <- subset(HCC7_HPC, subset = orig.ident=="PT1")
hcc7_pt1_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt1, slot = "data"))
hcc7_pt1_top_mtx <- hcc7_pt1_mtx[top3000_gene,]
hcc7_pt1_top_mtx <- na.omit(hcc7_pt1_top_mtx)
hcc7_pt1_cor <- as.data.frame(cor(hcc7_pt1_top_mtx,method = 'pearson'))
hcc7_pt1_cor_nu <- unique(as.vector(as.matrix(hcc7_pt1_cor)))
for(i in 1:nrow(hcc7_pt1_cor)){
  vector <- as.vector(hcc7_pt1_cor[i,1:(nrow(hcc7_pt1_cor)-i+1)])
  if(i == 1){
    hcc7_pt1_cor_nu <- vector
  }
  else{
    hcc7_pt1_cor_nu <- c(hcc7_pt1_cor_nu,vector)
  }
}



hcc7_pt2 <- subset(HCC7_HPC, subset = orig.ident=="PT2")
hcc7_pt2_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt2, slot = "data"))
hcc7_pt2_top_mtx <- hcc7_pt2_mtx[top3000_gene,]
hcc7_pt2_top_mtx <- na.omit(hcc7_pt2_top_mtx)
hcc7_pt2_cor <- as.data.frame(cor(hcc7_pt2_top_mtx,method = 'pearson'))
hcc7_pt2_cor_nu <- unique(as.vector(as.matrix(hcc7_pt2_cor)))
for(i in 1:nrow(hcc7_pt2_cor)){
  vector <- as.vector(hcc7_pt2_cor[i,1:(nrow(hcc7_pt2_cor)-i+1)])
  if(i == 1){
    hcc7_pt2_cor_nu <- vector
  }
  else{
    hcc7_pt2_cor_nu <- c(hcc7_pt2_cor_nu,vector)
  }
}



hcc7_pt3 <- subset(HCC7_HPC, subset = orig.ident=="PT3")
hcc7_pt3_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt3, slot = "data"))
hcc7_pt3_top_mtx <- hcc7_pt3_mtx[top3000_gene,]
hcc7_pt3_top_mtx <- na.omit(hcc7_pt3_top_mtx)
hcc7_pt3_cor <- as.data.frame(cor(hcc7_pt3_top_mtx,method = 'pearson'))
hcc7_pt3_cor_nu <- unique(as.vector(as.matrix(hcc7_pt3_cor)))
for(i in 1:nrow(hcc7_pt3_cor)){
  vector <- as.vector(hcc7_pt3_cor[i,1:(nrow(hcc7_pt3_cor)-i+1)])
  if(i == 1){
    hcc7_pt3_cor_nu <- vector
  }
  else{
    hcc7_pt3_cor_nu <- c(hcc7_pt3_cor_nu,vector)
  }
}




hcc7_pt4 <- subset(HCC7_HPC, subset = orig.ident=="PT4")
hcc7_pt4_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt4, slot = "data"))
hcc7_pt4_top_mtx <- hcc7_pt4_mtx[top3000_gene,]
hcc7_pt4_top_mtx <- na.omit(hcc7_pt4_top_mtx)
hcc7_pt4_cor <- as.data.frame(cor(hcc7_pt4_top_mtx,method = 'pearson'))
hcc7_pt4_cor_nu <- unique(as.vector(as.matrix(hcc7_pt4_cor)))
for(i in 1:nrow(hcc7_pt4_cor)){
  vector <- as.vector(hcc7_pt4_cor[i,1:(nrow(hcc7_pt4_cor)-i+1)])
  if(i == 1){
    hcc7_pt4_cor_nu <- vector
  }
  else{
    hcc7_pt4_cor_nu <- c(hcc7_pt4_cor_nu,vector)
  }
}




hcc7_pt5 <- subset(HCC7_HPC, subset = orig.ident=="PT5")
hcc7_pt5_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt5, slot = "data"))
hcc7_pt5_top_mtx <- hcc7_pt5_mtx[top3000_gene,]
hcc7_pt5_top_mtx <- na.omit(hcc7_pt5_top_mtx)
hcc7_pt5_cor <- as.data.frame(cor(hcc7_pt5_top_mtx,method = 'pearson'))
hcc7_pt5_cor_nu <- unique(as.vector(as.matrix(hcc7_pt5_cor)))
for(i in 1:nrow(hcc7_pt5_cor)){
  vector <- as.vector(hcc7_pt5_cor[i,1:(nrow(hcc7_pt5_cor)-i+1)])
  if(i == 1){
    hcc7_pt5_cor_nu <- vector
  }
  else{
    hcc7_pt5_cor_nu <- c(hcc7_pt5_cor_nu,vector)
  }
}






hcc7_pt6 <- subset(HCC7_HPC, subset = orig.ident=="PT6")
hcc7_pt6_mtx <-  as.data.frame(GetAssayData(object = hcc7_pt6, slot = "data"))
hcc7_pt6_top_mtx <- hcc7_pt6_mtx[top3000_gene,]
hcc7_pt6_top_mtx <- na.omit(hcc7_pt6_top_mtx)
hcc7_pt6_cor <- as.data.frame(cor(hcc7_pt6_top_mtx,method = 'pearson'))
hcc7_pt6_cor_nu <- unique(as.vector(as.matrix(hcc7_pt6_cor)))
for(i in 1:nrow(hcc7_pt6_cor)){
  vector <- as.vector(hcc7_pt6_cor[i,1:(nrow(hcc7_pt6_cor)-i+1)])
  if(i == 1){
    hcc7_pt6_cor_nu <- vector
  }
  else{
    hcc7_pt6_cor_nu <- c(hcc7_pt6_cor_nu,vector)
  }
}






hcc7_primary_cor3 <- as.data.frame(c(hcc7_pt1_cor_nu,hcc7_pt2_cor_nu,hcc7_pt3_cor_nu,hcc7_pt4_cor_nu,hcc7_pt5_cor_nu))
hcc7_primary_cor3 <-as.data.frame(t(hcc7_primary_cor3))
colnames(hcc7_primary_cor3) <- "cor"
hcc7_primary_cor3$type <- "primary"
hcc7_satellite_cor3 <- as.data.frame(unlist(unique(hcc7_pt6_cor_nu)))
colnames(hcc7_satellite_cor3) <- "cor"
hcc7_satellite_cor3$type <- "satellite"
hcc7_ps_info2 <- rbind(hcc7_primary_cor3,hcc7_satellite_cor3)
ggplot(data=hcc7_ps_info2,aes(x=cor,color=type))+
  geom_density(aes(fill = type), alpha=0.4)+
  theme_bw()+
  labs(y="Frequency",x="Cor",title = "HCC7")

hcc1_ps_info2$sample <- "hcc1"
hcc3_ps_info2$sample <- "hcc3"
hcc5_ps_info2$sample <- "hcc5"
hcc6_ps_info2$sample <- "hcc6"
hcc7_ps_info2$sample <- "hcc7"

sta_prim_sum <- rbind(hcc1_ps_info2,hcc3_ps_info2,hcc5_ps_info2,
                      hcc6_ps_info2,hcc7_ps_info2)
ggplot(sta_prim_sum,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("primary","satellite")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1)+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())
ggplot(sta_prim_sum,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~sample)+
  scale_color_manual(values =c('#8BABD3','#D7B0B0'))+
  stat_compare_means(comparisons = list(c("primary","satellite")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())




sta_prim_sum_cmn <- subset(sta_prim_sum,sample %in% c("hcc3","hcc5","hcc7"))
sta_prim_sum_sn <- subset(sta_prim_sum,sample %in% c("hcc1","hcc6"))

ggplot(sta_prim_sum_cmn,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~sample)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("primary","satellite")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


ggplot(sta_prim_sum_sn,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~sample)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("primary","satellite")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())



sta_prim_sum_cmn$nd_type <- 'CMN'
sta_prim_sum_sn$nd_type <- 'SN'

sta_prim_sum_snd <- rbind(sta_prim_sum_cmn,sta_prim_sum_sn)

sta_prim_sum_snd_pr <- subset(sta_prim_sum_snd,subset= type=='primary')
sta_prim_sum_snd_st <- subset(sta_prim_sum_snd,subset= type=='satellite')

ggplot(sta_prim_sum_snd_pr,aes(nd_type,cor,color=nd_type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("SN","CMN")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(sta_prim_sum_snd_st,aes(nd_type,cor,color=nd_type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("SN","CMN")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())




ggplot(sta_prim_sum_snd,aes(nd_type,cor,color=nd_type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("SN","CMN")),
                     method = "wilcox.test",label = "p.signif",
                     label.y =1 )+theme_bw() +
  facet_grid(~type)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())
