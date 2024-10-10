immune_merge <- subset(bigseu, subset = new.ident %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                         "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                         "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                         "PlasmaB cell","Neutrophil"))
nt_samples <- unique(immune_merge_nt@meta.data$patient_pt)

immune_merge_nt_result <- data.frame()
immune_merge_nt <- subset(immune_merge, subset = group %in% "NT")
immune_merge_pt <- subset(immune_merge, subset = group %in% "P")
immune_merge_st <- subset(immune_merge, subset = group %in% "S")

for(s in nt_samples){
  cells_of_sample <- subset(immune_merge_nt@meta.data, patient_pt == s)
  prop <- prop.table(table(cells_of_sample$new.ident))
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'nt',
    'patient' = unique(cells_of_sample$patient)
  )
  immune_merge_nt_result <- rbind(immune_merge_nt_result, prop_df)
}


immune_merge_pt_result <- data.frame()
pt_samples <- unique(immune_merge_pt@meta.data$patient_pt)

for(s in pt_samples){
  cells_of_sample <- subset(immune_merge_pt@meta.data, patient_pt == s)
  prop <- prop.table(table(cells_of_sample$new.ident))
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'pt',
    'patient' = unique(cells_of_sample$patient)
  )
  immune_merge_pt_result <- rbind(immune_merge_pt_result, prop_df)
}



immune_merge_st_result <- data.frame()
st_samples <- unique(immune_merge_st@meta.data$patient_pt)

for(s in st_samples){
  cells_of_sample <- subset(immune_merge_st@meta.data, patient_pt == s)
  prop <- prop.table(table(cells_of_sample$new.ident))
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'st',
    'patient' = unique(cells_of_sample$patient)
  )
  immune_merge_st_result <- rbind(immune_merge_st_result, prop_df)
}


prop_result_ptst_merge <- rbind(immune_merge_pt_result,immune_merge_st_result)
prop_result_ptnt_merge <- rbind(immune_merge_pt_result,immune_merge_nt_result)

prop_result_ptst_merge <- subset(prop_result_ptst_merge, subset = celltype %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                                                  "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                                                  "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                                                  "PlasmaB cell","Neutrophil"))
prop_result_ptnt_merge <- subset(prop_result_ptnt_merge, subset = celltype %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                                                  "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                                                  "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                                                  "PlasmaB cell","Neutrophil"))
prop_result_ptst_merge$celltype <-factor(prop_result_ptst_merge$celltype, levels = c( "CD8+ exhausted", "CD4+ memory", "CD4+ Treg","CD8+ cytotoxic",
                                                                                      "B cell", "CD8+ memory", "Dendritic cell","NK",
                                                                                      "Neutrophil","Macrophage", "PlasmaB cell","Proliferative T"))

prop_result_ptnt_merge_rm89  <- subset(prop_result_ptnt_merge,subset= patient %in% c('HCC1','HCC2','HCC3','HCC4','HCC5','HCC6','HCC7'))
ggplot(prop_result_ptnt_merge,aes(x=type,y=proportion,fill=type))+
  geom_boxplot(outlier.size = 0.1)+
  geom_jitter(width = 0.1,shape = 21, colour = "black",size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,0.4)+
  stat_compare_means(comparisons = list(c("pt","nt")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 0.32)


prop_result_ptnt_merge_rm89$celltype <-factor(prop_result_ptnt_merge_rm89$celltype, levels = c( "CD8+ exhausted", "CD4+ memory", "CD4+ Treg","CD8+ cytotoxic",
                                                                                                "B cell", "CD8+ memory", "Dendritic cell","NK",
                                                                                                "Neutrophil","Macrophage", "PlasmaB cell","Proliferative T"))


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





immune_merge_cmn <- subset(immune_merge, subset = cntype %in% "CMN")
immune_merge_sn <- subset(immune_merge, subset = cntype %in% "SN")

immune_merge_cmn_result <- data.frame()
cmn_samples <- unique(immune_merge_cmn@meta.data$patient_pt)

for(s in cmn_samples){
  cells_of_sample <- subset(immune_merge_cmn@meta.data, patient_pt == s)
  prop <- prop.table(table(cells_of_sample$new.ident))
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'cmn'
  )
  
  immune_merge_cmn_result <- rbind(immune_merge_cmn_result, prop_df)
}


immune_merge_sn_result <- data.frame()
sn_samples <- unique(immune_merge_sn@meta.data$patient_pt)

for(s in sn_samples){
  cells_of_sample <- subset(immune_merge_sn@meta.data, patient_pt == s)
  prop <- prop.table(table(cells_of_sample$new.ident))
  prop_df <- data.frame(
    'sample' = s, 
    'celltype' = names(prop),
    'proportion' = as.numeric(prop),
    'type' = 'sn'
  )
  
  immune_merge_sn_result <- rbind(immune_merge_sn_result, prop_df)
}


immune_prop_cmnsn_result_merge <- rbind(immune_merge_sn_result,immune_merge_cmn_result)


immune_prop_cmnsn_result_merge_rm89 <- subset(immune_prop_cmnsn_result_merge,subset = sample %in% prop_result_ptnt_merge_rm89$sample  )
immune_prop_cmnsn_result_merge_rm89 <- subset(immune_prop_cmnsn_result_merge_rm89,subset = type =='pt')





ggplot(immune_prop_cmnsn_result_merge_rm89,aes(x=type,y=proportion,fill=type))+
  geom_boxplot(outlier.size = 0)+
  geom_jitter(width = 0.1,shape = 21, colour = "black",size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,0.2)+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "wilcox.test",label = "p.signif",
                     label.y = 0.13)
















cell.prop_immune_merge <- as.data.frame(prop.table(table(immune_merge$patient_pt,immune_merge$new.ident)))

colnames(cell.prop_immune_merge) <- c("HCC","origin","proportion")

cell.prop_immune_merge <- subset(cell.prop_immune_merge, subset = origin %in% c("Macrophage","CD8+ exhausted","CD4+ memory",
                                                                                "CD8+ memory","Dendritic cell","CD4+ Treg",
                                                                                "NK","B cell","CD8+ cytotoxic","Proliferative T",
                                                                                "PlasmaB cell","Neutrophil"))



cell.prop_immune_merge[,'patient'] <- str_split_fixed(cell.prop_immune_merge$HCC,"_",2)[,1]
cell.prop_immune_merge[,'sample'] <- str_split_fixed(cell.prop_immune_merge$HCC,"_",2)[,2]

hcc1_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC1')
hcc2_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC2')
hcc3_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC3')
hcc4_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC4')
hcc5_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC5')
hcc6_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC6')
hcc7_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC7')
hcc8_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC8')
hcc9_cell.prop_immune_merge <- subset(cell.prop_immune_merge,subset = patient =='HCC9')

immune_merge_data_list <- list(hcc1_cell.prop_immune_merge, hcc2_cell.prop_immune_merge, hcc3_cell.prop_immune_merge,
                               hcc4_cell.prop_immune_merge, hcc5_cell.prop_immune_merge, hcc6_cell.prop_immune_merge,
                               hcc7_cell.prop_immune_merge, hcc8_cell.prop_immune_merge, hcc9_cell.prop_immune_merge)



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


