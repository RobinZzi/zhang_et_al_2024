
immune_merge_cmn_remove_mp <- subset(immune_merge_cmn, subset = group != "NT")
immune_merge_cmn_remove_mp <-  subset(immune_merge_cmn_remove_mp, subset = new.ident != "Macrophage")
immune_merge_cmn_remove_mp_result2 <- data.frame()
cmn_samples_remove_mp <- unique(immune_merge_cmn_remove_mp@meta.data$patient_pt)

cmn_samples_iflist <- list()
for(s in cmn_samples_remove_mp){
  cells_of_sample <- subset(immune_merge_cmn_remove_mp,subset=patient_pt==s)
  counts_of_sample <- GetAssayData(cells_of_sample, assay = "RNA", slot = "data")
  identity <- subset(cells_of_sample@meta.data, select = "new.ident")
  identity = droplevels(identity, exclude = setdiff(levels(identity),unique(identity)))
  cellchat <- createCellChat(object = counts_of_sample, meta = identity,  group.by = "new.ident")
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  cellchat <- cellchat %>% 
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    projectData(PPI.human) %>%
    computeCommunProb(raw.use = TRUE) %>%
    filterCommunication(min.cells = 3) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
  prop_df <- data.frame(
    'sample' = s, 
    'interactions' = sum(cellchat@net[["count"]]),
    'weight' = sum(cellchat@net[["weight"]]),
    'type' = 'cmn'
  )
  p <- rankNet(cellchat, mode = "single", stacked = F,do.stat = FALSE)
  cmn_samples_iflist[[s]] <- p$data
  immune_merge_cmn_remove_mp_result2 <- rbind(immune_merge_cmn_remove_mp_result2, prop_df)
}







immune_merge_sn_remove_mp <- subset(immune_merge_sn, subset = group != "NT")
immune_merge_sn_remove_mp <-  subset(immune_merge_sn_remove_mp, subset = new.ident != "Macrophage")
immune_merge_sn_remove_mp_result2 <- data.frame()
sn_samples_remove_mp <- unique(immune_merge_sn_remove_mp@meta.data$patient_pt)

sn_samples_iflist <- list()
for(s in sn_samples_remove_mp){
  cells_of_sample <- subset(immune_merge_sn_remove_mp,subset=patient_pt==s)
  counts_of_sample <- GetAssayData(cells_of_sample, assay = "RNA", slot = "data")
  identity <- subset(cells_of_sample@meta.data, select = "new.ident")
  identity = droplevels(identity, exclude = setdiff(levels(identity),unique(identity)))
  cellchat <- createCellChat(object = counts_of_sample, meta = identity,  group.by = "new.ident")
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  cellchat <- cellchat %>% 
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    projectData(PPI.human) %>%
    computeCommunProb(raw.use = TRUE) %>%
    filterCommunication(min.cells = 3) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
  prop_df <- data.frame(
    'sample' = s, 
    'interactions' = sum(cellchat@net[["count"]]),
    'weight' = sum(cellchat@net[["weight"]]),
    'type' = 'sn'
  )
  p <- rankNet(cellchat, mode = "single", stacked = F,do.stat = FALSE)
  sn_samples_iflist[[s]] <- p$data
  immune_merge_sn_remove_mp_result2 <- rbind(immune_merge_sn_remove_mp_result2, prop_df)
}





cmnsn_cellchat_sum_remove_mp <- rbind(immune_merge_cmn_remove_mp_result2,immune_merge_sn_remove_mp_result2)
ggplot(cmnsn_cellchat_sum_remove_mp,aes(x=type,y=interactions,fill=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  labs(x = "type", 
       y= "interaction",
       title="10x included")+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y =500 )+theme_bw()



cmnsn_cellchat_sum_remove_mp_drop <- subset(cmnsn_cellchat_sum_remove_mp,subset=sample%in%dropseq_list)

ggplot(cmnsn_cellchat_sum_remove_mp_drop,aes(x=type,y=interactions,fill=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  labs(x = "type", 
       y= "interaction",
       title="10x excluded")+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y =200 )+theme_bw()


p <- rankNet(cmn_cellchat, mode = "single", stacked = F,do.stat = FALSE)



cmn_pathway <- data.frame()
length(cmn_samples_iflist)
for (i in 1:length(cmn_samples_iflist)) {
  data <- cmn_samples_iflist[[i]]
  data[,'pt'] <- names(cmn_samples_iflist[i])
  data[,'type'] <- 'cmn'
  cmn_pathway <- rbind(cmn_pathway,data)
}
sn_pathway <- data.frame()
length(sn_samples_iflist)
for (i in 1:length(sn_samples_iflist)) {
  data <- cmn_samples_iflist[[i]]
  data[,'pt'] <- names(sn_samples_iflist[i])
  data[,'type'] <- 'sn'
  sn_pathway <- rbind(sn_pathway,data)
}

pathway_sum_cmn_sn <- rbind(sn_pathway,cmn_pathway)



