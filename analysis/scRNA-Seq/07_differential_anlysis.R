markers <- FindMarkers(seurat,ident.1 = "group1",ident.2 = "group2",group.by = "group_type" )
markers_sig <- subset(markers,subset = p_val_adj < 0.05)
markers_sig_up <- subset(markers_sig,subset = avg_log2FC > 0)
markers_sig_down <- subset(markers_sig,subset = avg_log2FC < 0)

up_go <- enrichGO(gene  = row.names(markers_sig_up),
                           OrgDb      = org.Hs.eg.db,
                           keyType    = 'SYMBOL',
                           ont        = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
up_go <- as.data.frame(up_go@result)
up_go [,"logp"] <- -log10(up_go$pvalue)
up_enrichment_fold <- apply(up_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
up_go$EF <- up_enrichment_fold

#Draw Bubble Plot
ggplot(data = up_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ")

down_go <- enrichGO(gene  = row.names(markers_sig_down),
                             OrgDb      = org.Hs.eg.db,
                             keyType    = 'SYMBOL',
                             ont        = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)
down_go <- as.data.frame(down_go@result)
down_go [,"logp"] <- -log10(down_go$pvalue)
down_enrichment_fold <- apply(down_go,1,function(x){ 
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"])) 
  enrichment_fold=round(GeneRatio/BgRatio,2) 
  enrichment_fold 
})
down_go$EF <- down_enrichment_fold

#Draw Bubble Plot
ggplot(data = down_go[1:20,],aes(y=reorder(Description,LogFDR),x=LogFDR))+
  geom_point(aes(color=LogFDR, size=EF))+
  scale_fill_gradient(expression(EF),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-LogFDR", 
       y= " ")