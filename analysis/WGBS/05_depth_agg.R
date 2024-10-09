#10-Mb genomic windows
rm(list=ls())
data_path <- "~/bismark/samtools_depth_result"
setwd(data_path)





folders <- list.files()

chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
          "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
for (i in 1:length(folders)) {
  setwd(paste(data_path,folders[i],sep = "/"))
  mtxs <- list.files()
  for(j in 1:length(mtxs)){
    sample_name <- str_split(mtxs[j],"_sorted.depth")[[1]][1]
    print(paste(sample_name,"start",sep = ":"))
    data <- fread(mtxs[j])
    data <- subset(data,subset=V1 %in% chrs)
    data <- mutate(data,bin_id=floor(V2/10000000+1))
    data <- mutate(data,bin_id=paste(V1,bin_id,sep = "_"))
    data <- aggregate(data$V3,by=list(data$bin_id),FUN=sum)
    print(paste(sample_name,"end",sep = ":"))
    write.table(data,paste(sample_name,"_10M",".txt",sep = ""))
  }
  
}

data_path <- "~/bismark/samtools_depth_agg"

setwd(data_path)
folders <- list.files()
for(i in 1:length(folders)){
  setwd(paste(data_path,folders[i],sep = "/"))
  combined_df <- data.frame()
  file_list <- list.files()
  for (j in 1:length(file_list)) {  
    sample_name <- tools::file_path_sans_ext(basename(file_list[j])) 
    sample_name <- strsplit(sample_name,split = "_R1")[[1]][1]
    df <- as.data.frame(fread(file_list[j]))
    row.names(df) <- df$Group.1
    df <- dplyr::select(df,3)
    colnames(df) <-sample_name
    df$rownames <- rownames(df) 
    if(j == 1){
      combined_df <- df
    }
    else{
      combined_df <- merge(combined_df,df,by="rownames",all=FALSE)
    }
  } 
  rownames(combined_df) <- combined_df$rownames  
  combined_df$rownames <- NULL  
  write.table(combined_df,paste(folders[i],".txt",sep = ""))
}





