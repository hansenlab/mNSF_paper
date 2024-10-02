set=3

# ml conda_R/4.4
# R

# conda create -n  spatialPCA -- also works
# conda activate spatialPCA
# conda install conda-forge::r-devtools 

library(SpatialPCA)
# spatialPCA package pacakge: https://github.com/shangll123/SpatialPCA
# spatialPCA analysis code:https://github.com/shangll123/SpatialPCA_analysis_codes
# spatialPCA tutorial: https://lulushang.org/SpatialPCA_Tutorial/
############################## 
dir_out  = "/dcs04/hansen/data/ywang/revision/domain/"
setwd(dir_out)
source("/users/ywang/Hansen_projects/mNSF/revision/spatialPCA/cere/function.R")

dir_Data = "/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/"
dir_dlpfc_data = paste0(dir_Data,"DLPFC/")
############################## 
for(kdonor in 1:3){
  LIBD =readRDS(paste0(dir_dlpfc_data, "LIBD_sample", (kdonor*4 - 4 +set), "_spatialPcaObj.rds"))
  
  if(kdonor==1){
    genes_svg = rownames(LIBD@normalized_expr)
  }else{
    genes_svg = intersect(genes_svg, rownames(LIBD@normalized_expr))
  }
  
}
# 
# length(genes_svg)
# # 905
saveRDS(genes_svg, file = paste0("genes_svg_dlpfc_set",set,".rds"))
# genes_svg = readRDS("genes_svg_dlpfc.rds")
############################## use dlpfc data with the same genes selected in mnsf paper

list_count = list()
list_location  = list()

for(kdonor in 1:3){
  ksample = (kdonor-1)*4+ set
  
  
  ############################## 
  load(paste0(dir_dlpfc_data, "LIBD_sample",ksample,".RData"))
  
  colnames(count_sub) = rownames(xy_coords) = paste0("Sample",kdonor,"_",colnames(count_sub))
  
  
  list_count[[kdonor]] = count_sub[genes_svg, ]
  list_location[[kdonor]] = data.matrix(xy_coords[ ,])
  
  
}

# head(out_SpatialPCA$SpatialPC_list$Sample1)
out_SpatialPCA = SpatialPCA_Multiple_Sample(count_list=list_count, location_list=list_location, 
                                            gene.type="hvg",sparkversion="spark",
                                            numCores_spark=1,gene.number=length(genes_svg), 
                                            customGenelist=NULL, min.loctions = 0, 
                                            min.features=0,bandwidth_common=0.1,SpatialPCnum_=10)

saveRDS(out_SpatialPCA, file = paste0("out_SpatialPCA_3samples_svgs_Set",set,"_L10.rds"))

