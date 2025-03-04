#Fig_2_A
ggplot(immune_prop_cmnsn_result_merge_rm89,aes(x=type,y=proportion,color=type))+
  geom_boxplot(outlier.size = 0,outlier.alpha = 1,outlier.color = "white")+
  geom_jitter(width = 0.1,shape = 20,size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,1)+
  scale_color_manual(values=c("cmn"="#d72323","sn"="#3f72af"))+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "wilcox.test",
                     label.y = 0.9 )+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())


#Fig_2_B
ggplot(cmnsn_cellchat_sum_remove_mp_drop,aes(x=type,y=interactions,fill=type))+
  geom_boxplot()+geom_jitter(width = 0.1,shape = 21, colour = "black")+
  labs(x = "type", 
       y= "interaction",
       title="10x excluded")+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "wilcox.test",
                     label.y =200 )+theme_bw()


#Fig_2_C
ggplot(pathway_sum_cmn_sn_sig,aes(x=type,y=contribution.scaled,color=type))+
  geom_boxplot(outlier.size = 0,outlier.colour = "white")+geom_jitter(width = 0.1,size=0.5,shape = 20)+
  labs(x = "type", 
       y= "interaction",
       title="10x excluded")+
  scale_color_manual(values=c("sn"="#3f72af","cmn"="#d72323"))+
  stat_compare_means(comparisons = list(c("cmn","sn")),
                     method = "t.test",paired = F,
                     label.y =0.7)+theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=15,face="bold",angle = 45,hjust = 1,color = 'black'))+
  facet_grid(~name)+
  ylim(0,0.8)


#Fig_2_D
ggplot(prop_result_ptst_merge,aes(x=type,y=proportion,color=type))+
  geom_boxplot(outlier.size = 0,outlier.alpha = 1,outlier.color = "white")+
  geom_jitter(width = 0.1,shape = 20,size=0.5)+
  facet_wrap(~ celltype)+
  ylim(0,0.2)+
  scale_color_manual(values=c("pt"="#d72323","st"="#3f72af"))+
  stat_compare_means(comparisons = list(c("pt","st")),
                    method = "wilcox.test",
                   label.y = 0.5 )+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())


#Fig_2_E
ggplot(sta_prim_sum_cmn,aes(type,cor,color=type))+
  geom_boxplot(width=0.5,outlier.size=0)+
  facet_grid(~sample)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("primary","satellite")),
                     method = "wilcox.test",
                     label.y =1 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  
        axis.text.x = element_text(face="bold",angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())


#Fig_2_F
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


#Fig_2_G
ggplot(data = pt_up_sum_go[1:20,],aes(y=reorder(Description,logp),x=logp))+
  geom_point(aes(color=logp, size=GeneRatio))+
  scale_fill_gradient(expression(GeneRatio),low="blue",high="red")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 14, face = "bold"),
        axis.text=element_text(size=14,face = "bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))+
  labs(x = "-Logp", 
       y= " ",
       title="pt_up")