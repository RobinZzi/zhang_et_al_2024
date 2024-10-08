#Fig_4A,B,C
DimPlot(tx_HPC,group.by = "patient",cols = c('HCC2'='#fb9795','HCC8'='#00ADC4','HCC9'='#4C1A72'))
DimPlot(tx_HPC,group.by = "sample")
DimPlot(tx_HPC,group.by = "group",cols = c('P'='#E95C59','NT'='#57C3F3'))


#Fig_4D-E
hcc_up_venn <- list(hcc2_uplist$symbol,hcc8_uplist$x,hcc9_uplist$gene)
venn.diagram(hcc_up_venn, filename = 'up.png', imagetype = 'png', ,category.names = c("hcc2" , "hcc8" , "hcc9"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')

tcga_up_venn <- list(TCGA_neg_re$symbol,hcc2_uplist$x,TCGA_diff_hypo_up$gene)


venn.diagram(tcga_up_venn, filename = 'up.png', imagetype = 'png', ,category.names = c("tcga_neg" , "hypo_up" , "tcga_up"),
             fill = c('#4D157D', '#84C7DB', '#C8D948'), alpha = 0.50, 
             cat.col = c('#4D157D', '#84C7DB', '#C8D948'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('#4D157D', '#84C7DB', '#C8D948'), cex = 1.5, fontfamily = 'serif')

#Fig_4G
ggscatterstats(data = meth_fpkm_comb, 
               y =GADD45A,
               x = mean_meth,
               centrality.para = "mean"
)