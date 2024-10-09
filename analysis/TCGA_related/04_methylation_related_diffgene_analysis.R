hypergroup_sm_tm <- subset(hypergroup_sm,subset=row.names(hypergroup_sm) %in% tumorid_p)
count_order_hypo1 <- dplyr::select(count,row.names(hypogroup_sm))
count_order_hyper <- dplyr::select(count,row.names(hypergroup_sm_tm))

count_deg_data <- data.frame(cbind(count_order_hypo1,count_order_hyper))
count_deg_data <- as.data.frame(lapply(count_deg_data,as.numeric))
count_deg_Group <- factor(c(rep("hypo",61),rep("hyper",86)),levels=c('hypo','hyper'))

#edgeR 
dge <- DGEList(counts=count_deg_data,group=count_deg_Group)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge) 
design <- model.matrix(~0+count_deg_Group)
rownames(design)<-colnames(dge)
colnames(design)<-levels(count_deg_Group)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast=c(1,-1)) 
DEG=topTags(fit2, n=nrow(exp))
DEG=as.data.frame(DEG)
logFC_t=0
P.Value_t = 0.05
k1 = (DEG$FDR < P.Value_t)&(DEG$logFC < -logFC_t)
k2 = (DEG$FDR < P.Value_t)&(DEG$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
DEG$change <- change
DEG_edgeR <- DEG
DEG_edgeR[,'gene'] <- row.names(DEG_edgeR)

# down stable     up 
# 3080 34655   2366 


#wilcox-test
cpm_order_data <- edgeR::cpm(count_deg_data)
cpm_order_tr <- as.data.frame(t(cpm_order_data))
colnames(cpm_order_tr) <- count_deg_genes$`rownames(count_foridtransfer)`
cpm_order_tr[,"sample_type"] <- as.factor(c(rep("hypo",107),rep("hyper",150)))
cpm_order_tr <- dplyr::select(cpm_order_tr,"sample_type",everything())
lihc_meth  <- numeric(6)
lihc_methdiff_base <- as.data.frame(lihc_meth)
for(i in 2: length(cpm_order_tr)){
  test_diff <- wilcox.test(cpm_order_tr[,i] ~ cpm_order_tr$sample_type)
  test_diff_less <- wilcox.test(cpm_order_tr[,i] ~ cpm_order_tr$sample_type,alternative = "less")
  test_diff_more <- wilcox.test(cpm_order_tr[,i] ~ cpm_order_tr$sample_type,alternative = "greater")
  hypomean <- mean(log2(cpm_order_tr_hypo[,names(cpm_order_tr)[i]]+1))
  hypermean <- mean(log2(cpm_order_tr_hyper[,names(cpm_order_tr)[i]]+1))
  logFC <- hypomean-hypermean
  lihc_meth[1] <- test_diff$p.value
  lihc_meth[2] <- test_diff_less$p.value
  lihc_meth[3] <- test_diff_more$p.value
  lihc_meth[4] <- hypomean
  lihc_meth[5] <- hypermean
  lihc_meth[6] <- logFC
  lihc_methdiff_base <- cbind(lihc_methdiff_base,lihc_meth)
  names(lihc_methdiff_base)[i] <- names(cpm_order_tr)[i]
}

lihc_methdiff_result <- as.data.frame(t(lihc_methdiff_base))
colnames(lihc_methdiff_result) <- c("diff","hypo_more","hypo_less","hypomean","hypermean","log2FC")
lihc_methdiff_result <- lihc_methdiff_result[-1,]
lihc_methdiff_result$gene <- row.names(lihc_methdiff_result)
lihc_methdiff_result$diff_padj <- p.adjust(as.vector(lihc_methdiff_result$diff), "fdr", n = length(lihc_methdiff_result$diff))
lihc_methdiff_result$hypo_more_padj <- p.adjust(as.vector(lihc_methdiff_result$hypo_more), "fdr", n = length(lihc_methdiff_result$hypo_more))
lihc_methdiff_result$hypo_less_padj <- p.adjust(as.vector(lihc_methdiff_result$hypo_less), "fdr", n = length(lihc_methdiff_result$hypo_less))

lihc_methdiff_result$change <- ifelse(lihc_methdiff_result$hypo_more_padj < 0.05,"up",ifelse(lihc_methdiff_result$hypo_less_padj < 0.05,"down","stable"))

# down stable     up 
#4898  29198   6005  