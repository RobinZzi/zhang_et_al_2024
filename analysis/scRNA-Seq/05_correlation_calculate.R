library(ggpubr)

hcc1_av <- AverageExpression(HCC1_HPC,group.by = "orig.ident",assays = 'RNA')
hcc1_av <- hcc1_av[[1]]
hcc1_cg=names(tail(sort(apply(hcc1_av, 1, sd)),1000))
pheatmap::pheatmap(cor(hcc1_av[hcc1_cg,],method = 'spearman'))



hcc1_all <- as.data.frame(GetAssayData(object = HCC1_HPC, slot = "data"))
hcc1_top <- hcc1_all[top3000_gene,]
hcc1_top <- na.omit(hcc1_top)
hcc1_cor <- as.data.frame(cor(hcc1_top,method = 'pearson'))

hcc1_anno <- as.data.frame(HCC1_HPC$orig.ident)
hcc1_anno2 <- as.data.frame(HCC1_HPC$satellite_region)
colnames(hcc1_anno) <- "origin"
colnames(hcc1_anno2) <- "origin"

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




hcc1_cor_mean <- as.data.frame(rbind(hcc1_pt2_pt4_mean,hcc1_pt2_pt5_mean ,hcc1_pt4_pt5_mean,
                                     hcc1_pt2_pt2_mean,hcc1_pt4_pt4_mean,hcc1_pt5_pt5_mean))
colnames(hcc1_cor_mean) <- "cor"
hcc1_cor_mean$comb <- row.names(hcc1_cor_mean)
hcc1_cor_mean$sample <- "hcc1"
hcc1_cor_mean$comb_type <- c(rep("cross",3),rep("self",3))



sc_cor_sum <- rbind(hcc1_cor_mean,hcc2_cor_mean,hcc3_cor_mean,
                    hcc4_cor_mean,hcc5_cor_mean,hcc6_cor_mean,
                    hcc7_cor_mean,hcc8_cor_mean,hcc9_cor_mean)


ggplot(sc_cor_sum,aes(x=sample,y=cor,fill=comb_type,color=comb_type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")





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

hcc3_cross_mean <- mean(subset(hcc3_cor_all,subset=type=="cross")$cor)
hcc3_self_mean <- mean(subset(hcc3_cor_all,subset=type=="self")$cor)


hcc1_mean <- as.data.frame(abs(hcc1_self_mean-hcc1_cross_mean))
colnames(hcc1_mean) <- "cor_delta"
hcc1_mean$patient <- "hcc1"
hcc1_mean$type <- "CMN"


hcc3_mean <- as.data.frame(abs(hcc3_self_mean-hcc3_cross_mean))
colnames(hcc3_mean) <- "cor_delta"
hcc3_mean$patient <- "hcc3"
hcc3_mean$type <- "SN"


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
hcc1_pt2_top_mtx <- hcc1_pt2_mtx[top3000_gene,]
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
