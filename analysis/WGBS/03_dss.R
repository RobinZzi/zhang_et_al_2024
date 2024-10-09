library(DSS)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggplot2)
setwd("~/dss")
rm(list = ls())



setwd("~/bismark_cov")
c1 <- fread("C1.bismark.cov.gz")
c2 <- fread("C2.bismark.cov.gz")
c3 <- fread("C3.bismark.cov.gz")
e1 <- fread("OE1.bismark.cov.gz")
e2 <- fread("OE2.bismark.cov.gz")
e3 <- fread("OE3.bismark.cov.gz")


c1$total <- c1$V5+c1$V6
c1_con <- dplyr::select(c1,1,2,7,5)
colnames(c1_con) <- c("chr","pos","N","X")

c2$total <- c2$V5+c2$V6
c2_con <- dplyr::select(c2,1,2,7,5)
colnames(c2_con) <- c("chr","pos","N","X")

c3$total <- c3$V5+c3$V6
c3_con <- dplyr::select(c3,1,2,7,5)
colnames(c3_con) <- c("chr","pos","N","X")

e1$total <- e1$V5+e1$V6
e1_con <- dplyr::select(e1,1,2,7,5)
colnames(e1_con) <- c("chr","pos","N","X")

e2$total <- e2$V5+e2$V6
e2_con <- dplyr::select(e2,1,2,7,5)
colnames(e2_con) <- c("chr","pos","N","X")

e3$total <- e3$V5+e3$V6
e3_con <- dplyr::select(e3,1,2,7,5)
colnames(e3_con) <- c("chr","pos","N","X")

BSobj <- makeBSseqData( list(c1_con, c2_con,c3_con, e1_con, e2_con,e3_con),
                             c("C1","C2","C3", "E1", "E2","E3") )

dmlTest.sm <- DMLtest(BSobj, group1=c("C1","C2","C3"), group2=c("E1","E2","E3"), smoothing=TRUE)




dmrs <- callDMR(dmlTest.sm, p.threshold=0.01,delt=0.2,minlen = 50, minCG = 3, dis.merge = 100)

dmrs <- mutate(dmrs,state = case_when(diff.Methy<0 ~ "Up-regulated", # 上调
                                                diff.Methy>0 ~ "Down-regulated", # 下调
                                                TRUE ~ "Unchanged"))


dmrs$region <- paste(dmrs$chr,dmrs$start,sep = "_")

dmls <- callDML(dmlTest.sm, p.threshold=0.01)
dmls <- mutate(dmls,state = case_when(diff<0 ~ "Up-regulated", # 上调
                                                diff>0 ~ "Down-regulated", # 下调
                                                TRUE ~ "Unchanged"))


table(dmrs$state)
showOneDMR(dmrs['1477',], BSobj)
df_color <- c("#1f76b6", "#ff7d0e")
dmrs_prop <- as.data.frame(t(table(dmrs$state)))
ggplot(dmrs_prop, aes(x = 1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity", position = position_fill(reverse = T), width = 1, show.legend = F)+
  coord_polar(theta = "y", direction = -1, start = pi*1.5)+  #1为原bar从下至上顺时针
  labs(x = NULL, y = NULL, title = "Pie Chart of Vehicle Class - Bad")+
  theme_bw()+
  theme(aspect.ratio = 1/1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))

table(dmls$state)
dmls_prop <- as.data.frame(t(table(dmls$state)))
ggplot(dmls_prop, aes(x = 1, y = Freq, fill = Var2))+
  geom_bar(stat = "identity", position = position_fill(reverse = T), width = 1, show.legend = F)+
  coord_polar(theta = "y", direction = -1, start = pi*1.5)+  #1为原bar从下至上顺时针
  labs(x = NULL, y = NULL, title = "Pie Chart of Vehicle Class - Bad")+
  theme_bw()+
  theme(aspect.ratio = 1/1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))



write.table(dmrs, "dmrs.bedGraph",sep = "\t",quote=F,row.names = F,col.names = F)
