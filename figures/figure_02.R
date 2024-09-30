
#Sup_Fig_2_A
ggplot(immune_prop_cmnsn_result_merge_rm89,aes(x=type,y=proportion,color=type))+
  geom_boxplot(outlier.size = 0,outlier.alpha = 1,outlier.color = "white")+
  geom_jitter(width = 0.1,shape = 20,size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,1)+
  scale_color_manual(values=c("cmn"="#d72323","sn"="#3f72af"))+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",label = "p.signif",
                     label.y = 0.8 )+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())


#Sup_Fig_2_B
