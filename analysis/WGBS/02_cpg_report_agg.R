library(data.table)
library(stringr)
rm(list=ls())

chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
          "chr16","chr17","chr18","chr19","chr20","chr21","chr22")

setwd("/PATH/cpg_report")

fs <- list.files()

if(T){
  for (i in 1:length(fs)) {
    setwd("/PATH/cpg_report")
    data <- fread(paste(fs[i]))
    data <- subset(data,data$V1 %in% chrs)
    colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
    data[,3] <- round(data[,2]/10000)
    data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
    data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
    sample_id <-str_split_fixed(fs[i],".CpG_report.txt",2)[1]
    data[,"count"] <- data[,2]+data[,3]
    data[,paste(sample_id)] <- data[,2]/(data[,2]+data[,3])
    setwd("/PATH/methy_mtx/10kbin")
    write.table(data,paste(sample_id,"_10kbin.txt",sep = ""))
  }
  
  
  
  for (i in 1:length(fs)) {
    setwd("/PATH/cpg_report")
    data <- fread(paste(fs[i]))
    data <- subset(data,data$V1 %in% chrs)
    colnames(data) <- c("chr","site","sum","meth","dmeth","t1","t2")
    data[,3] <- round(data[,2]/100000)
    data <- tidyr::unite(data,"id",chr,sum,sep="_",remove=TRUE)
    data <- aggregate(data[,3:4],by=list(id=data$id),FUN=sum)
    sample_id <-str_split_fixed(fs[i],".CpG_report.txt",2)[1]
    data[,"count"] <- data[,2]+data[,3]
    data[,paste(sample_id)] <- data[,2]/(data[,2]+data[,3])
    setwd("/PATH/methy_mtx/100kbin")
    write.table(data,paste(sample_id,"_100kbin.txt",sep = ""))
  }
}



