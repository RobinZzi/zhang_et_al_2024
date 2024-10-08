#Sup_Fig_%A
ggplot(G6_OE_meansum,aes(Group,HMD_Methy,color=Group))+
  geom_boxplot(width=0.5)+
  geom_jitter(width = 0.1,shape = 20,size=2)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("V","OE")),
                     method = "t.test",label = "p.signif",
                     label.y =0.62 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",size=12,angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(GA_OE_meansum,aes(Group,HMD_Methy,color=Group))+
  geom_boxplot(width=0.5)+
  geom_jitter(width = 0.1,shape = 20,size=2)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("V","OE")),
                     method = "t.test",label = "p.signif",
                     label.y =0.62 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",size=12,angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(G6_SH_meansum,aes(Group,HMD_Methy,color=Group))+
  geom_boxplot(width=0.5)+
  geom_jitter(width = 0.1,shape = 20,size=2)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("V","SH")),
                     method = "t.test",label = "p.signif",
                     label.y =0.62 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",size=12,angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())

ggplot(GA_CR_meansum,aes(Group,HMD_Methy,color=Group))+
  geom_boxplot(width=0.5)+
  geom_jitter(width = 0.1,shape = 20,size=2)+
  scale_color_manual(values =c('#b11a2b','#4a74a4'))+
  stat_compare_means(comparisons = list(c("V","CR")),
                     method = "t.test",label = "p.signif",
                     label.y =0.62 )+theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),  ##去掉背景网格
        axis.text.x = element_text(face="bold",size=12,angle = 45,hjust = 1,color = 'black'),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.direction = "vertical",
        legend.title =element_blank())