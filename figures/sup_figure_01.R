my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3',
               '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A',
               '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
               '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35',
               '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD',
               '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')


#Sup_Fig_1_A
DimPlot(bigseu,group.by = "seurat_clusters",label = T,cols=my37colors)+NoLegend()


#Sup_Fig_1_B
DimPlot(bigseu,group.by = "big_cell_cluster2_info",label=F,cols=my37colors[3:6])
FeaturePlot(bigseu,features = c("PTPRC","AMBP","PECAM1","COL1A2"))


#Sup_Fig_1_C
DimPlot(hcc2, label = F,group.by = "lib.method",cols= c("#e9e4d4","#57C3F3","#E95C59"),shuffle = T)
DimPlot(hcc2, label = T,group.by = "new.ident")+ NoLegend()

#Sup_Fig_1_D
VlnPlot(bigseu, features = c("nFeature_RNA", "nCount_RNA"), group.by = "lib.method",cols = my37colors,pt.size = 0)

#Sup_Fig_1_E
ggplot(prop_result_ptnt_merge_rm89,aes(x=type,y=proportion,color=type))+
  geom_boxplot(outlier.size = 0.1,outlier.alpha = 1,outlier.colour = "white")+
  geom_jitter(width = 0.1,shape = 20,size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,0.2)+
  scale_color_manual(values=c("pt"="#d72323","nt"="#3f72af"))+
  stat_compare_means(comparisons = list(c("pt","nt")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 0.15 )+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())


#Sup_Fig_1_F
ggplot(cell.prop_immune_merge,aes(sample,proportion,fill=origin))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  scale_fill_manual(values=my36colors)+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 0),legend.text=element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10,color="black",angle = 45,hjust = 1),axis.text.y = element_text(size = 10,color="black"))


#Sup_Fig_1_G
immuneprop_plot <- function(data) {
  ggplot(data,aes(sample,proportion,fill=origin))+
    geom_bar(stat="identity",position="fill")+
    ggtitle("")+
    scale_fill_manual(values=my36colors)+
    theme_bw()+
    theme(axis.ticks.length=unit(0.5,'cm'))+
    guides(fill=guide_legend(title=NULL))+
    labs(x=paste(data[1,4]))+
    theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 0),legend.text=element_text(size = 12))+
    theme(axis.text.x = element_text(size = 10,color="black",angle = 45,hjust = 1),axis.text.y = element_text(size = 10,color="black"))+
    NoLegend()
}

immuneprop_plot_list <- lapply(immune_merge_data_list, immuneprop_plot)

do.call(grid.arrange, c(immuneprop_plot_list, ncol = 5))