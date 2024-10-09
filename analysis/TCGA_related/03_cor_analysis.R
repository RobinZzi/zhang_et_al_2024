#calculate mean methylation level and combine the result with fpkm mtrix
meth_mean <- as.data.frame(as.numeric(meth))
meth_id <-colnames(meth_mean)
exp_id <- colnames(fpkm)
sampleid <- intersect(exp_id,meth_id)
tumorid <- grep(pattern = '-01',x=sampleid,value = T)
normalid <- grep(pattern = '-11',x=sampleid,value = T)
comb <- cbind(meth_mean,fpkm)
tumorcomb <- comb[tumorid,]

for (i in 1:length(colnames(tumorcomb))){
  
  test <- cor.test(as.numeric(tumorcomb[,i]),y,method="spearman")
  
  cor_data_df[i,2] <- test$estimate
  
  cor_data_df[i,3] <- test$p.value
  
}

#take genes whose fdr below 0.05 when calculating correlation
cor_data_sig <- cor_data_df %>% 
  filter(fdr < 0.05)

#order the table by correlation
cor_data_sig <- cor_data_sig[order(cor_data_sig$correlation),]

#take top 1000 correlated-genes in neg/pos 
cor_data_df_neg_top <- cor_data_sig[1:1000,]
cor_data_df_pos_top <- cor_data_sig[9187:10186,] 

cor_data_df_neg_top_go <- as.data.frame(cor_data_df_neg_top_go@result)
cor_data_df_neg_top_go [,"LogFDR"] <- -log10(cor_data_df_neg_top_go$p.adjust)       

cor_data_df <- mutate(cor_data_df,state = case_when(fdr < 0.05 & correlation > 0 ~ 'pos',
                                                    fdr < 0.05 & correlation < 0 ~ 'neg'))
#neg  pos 
#5790 5111        