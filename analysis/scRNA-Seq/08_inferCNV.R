set.seed(1)
library("infercnv")
library("Seurat")
library(AnnoProbe)
options(scipen = 100)
infercnv_run <- function(sample, condition){
  se <- readRDS(paste0("data/",condition,"/",sample,".RDS"))
  if(sample=="hcc8_strt_hpc" & condition == "strt"){
     se <- subset(se, cells=rownames(se@meta.data)[which(se$sample_pt!="HCC8_pt4")])   
   }
  dat <- GetAssayData(se,assay = "RNA",slot = "counts")
  dat <- as.data.frame(dat)
  geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
  geneInfor=geneInfor[c(1,4:6)] 
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  dat=dat[match(geneInfor[,1], rownames(dat)),] 
  rownames(geneInfor) <- geneInfor$SYMBOL
  geneInfor <- geneInfor[,-1]
  cluster <- se$sample_pt
  cluster_factor <- unique(cluster)
  normal_items <- cluster_factor[grepl("_NT|_nt", cluster_factor)]  
  meta <- subset(se@meta.data, select = c("sample_pt"))
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                      annotations_file=meta,
                                      delim="\t",
                                      gene_order_file=geneInfor,
                                      ref_group_name=normal_items)
  dir.create(paste0("result/infercnv/",condition,"/"))
  dir.create(paste0("result/infercnv/",condition,"/",sample))
  if(condition=="10x"){
    cutoff=0.1
  }else{
    cutoff=1
  }
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=cutoff, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0("result/infercnv/",condition,"/",sample), 
                               cluster_by_groups=T,  
                               denoise=T,     
                               num_threads=16,
                               plot_chr_scale = T,
                               k_nn=30,  # 1.9.1 default param
                               leiden_resolution = 1,  # 1.9.1 default param
                               leiden_method = "simple",  # 1.9.1 default param
                               leiden_function = "modularity",
                               # analysis_mode = "subclusters",
                               HMM=T) 
  infercnv_final_obj <- readRDS(paste0("result/infercnv/",condition,"/",sample,"/run.final.infercnv_obj"))
  cnv_matrix<-infercnv_final_obj@expr.data
  write.csv(cnv_matrix,paste0("result/infercnv/",condition,"/",sample,"/final.infercnv_matrix.csv"))
}

condition <- "10x"
samples <- c("hcc2_hpc","hcc8_hpc")
for(sample in samples){
  infercnv_run(sample,condition)
}

condition <- "strt"
samples <- c("hcc2_strt_hpc","hcc3_strt_hpc","hcc4_strt_hpc","hcc5_strt_hpc","hcc6_strt_hpc","hcc7_strt_hpc","hcc8_strt_hpc")
for(sample in samples){
  infercnv_run(sample,condition)
}