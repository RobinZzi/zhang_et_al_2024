#Sup_Fig_2_A

cellchat_color_cmn <- c("#E5D2DD","#E95C59","#57C3F3","#585658","#E59CC4","#53A85F", "#D6E7A3","#F1BB72","#E0D4CA",
                        "#23452F","#476D87","#AB3282","#BD956A","#F3B1A0")
netVisual_circle(cellchat_hcc6_pt2@net$count,
                 color.use = cellchat_color_cmn,
                 vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

cellchat_color_sn <- c("#E95C59","#57C3F3","#E59CC4","#53A85F", "#D6E7A3","#F1BB72","#E0D4CA",
                       "#23452F","#476D87","#AB3282","#BD956A")
netVisual_circle(cellchat_hcc3_pt2@net$count,
                 color.use = cellchat_color_sn,
                 vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")


#Sup_Fig_2_B
ggplot(sc_cor_cmn,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  facet_grid(~patient)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test", 
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(sc_cor_sn,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  facet_grid(~patient)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("cross","self")),
                     method = "wilcox.test",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())



#Sup_Fig_2_C
ggplot(sta_prim_sum_snd,aes(nd_type,cor,color=nd_type))+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("SN","CMN")),
                     method = "wilcox.test", 
                     label.y =1 )+theme_bw() +
  facet_grid(~type)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


#Sup_Fig_2_D
ggplot(diff_gene_count,aes(Patient,Num,fill=Type))+
  geom_bar(width=0.5,stat = 'identity',position = position_dodge(0.5))+
  scale_fill_manual(values =c('#b11a2b','#4a74a4'))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",size=12,angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank())


#Sup_Fig_2_E
sta_down_Venn <- list(hcc1 = row.names(hcc1_sta_markers_sig_down), hcc3 = row.names(hcc3_sta_markers_sig_down), hcc5 = row.names(hcc5_sta_markers_sig_down),
                      hcc6 = row.names(hcc6_sta_markers_sig_down),hcc7 = row.names(hcc7_sta_markers_sig_down))

sta_up_Venn <- list(hcc1 = row.names(hcc1_sta_markers_sig_up), hcc3 = row.names(hcc3_sta_markers_sig_up), hcc5 = row.names(hcc5_sta_markers_sig_up),
                    hcc6 = row.names(hcc6_sta_markers_sig_up),hcc7 = row.names(hcc7_sta_markers_sig_up))

venn.diagram(sta_down_Venn, filename = 'sta_down.png', imagetype = 'png',fontfamily = 'serif',height = 5000, 
             width = 7000, 
             fill = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'), cat.cex = 0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'),cex = 2.5 )

venn.diagram(sta_up_Venn, filename = 'sta_up.png', imagetype = 'png',fontfamily = 'serif',height = 5000, 
             width = 7000, 
             fill = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'), cat.cex =0, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948','#ff7f57','#FFD4E7'),cex = 2.5 )