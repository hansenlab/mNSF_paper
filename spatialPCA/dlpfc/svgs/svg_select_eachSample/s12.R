kdonor=1

############################## 
# ml conda_R/4.4
# R

# conda create -n  spatialPCA -- also works
# conda activate spatialPCA
# conda install conda-forge::r-devtools 

############################## 
library(SpatialPCA)
# spatialPCA package pacakge: https://github.com/shangll123/SpatialPCA
# spatialPCA analysis code:https://github.com/shangll123/SpatialPCA_analysis_codes
# spatialPCA tutorial: https://lulushang.org/SpatialPCA_Tutorial/

############################## 
dir_out  = "/dcs04/hansen/data/ywang/revision/domain/"
setwd(dir_out)
library(matrixStats)

dir_Data = "/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/"
dir_dlpfc_data = paste0(dir_Data,"DLPFC/")
############################## use dlpfc data with the same genes selected in mnsf paper

for(ksample in 12){
  ############################## 
  load(paste0(dir_dlpfc_data, "LIBD_sample",ksample,".RData"))

  
  colnames(count_sub) = rownames(xy_coords) = paste0("Sample",ksample,"_",colnames(count_sub))
  

  ############################## 
  # location matrix: n x 2, count matrix: g x n.
  
  LIBD = CreateSpatialPCAObject(counts=count_sub, location=data.matrix(xy_coords), 
                                project = "SpatialPCA",
                                gene.type="spatial",sparkversion="spark",numCores_spark=5,
                                gene.number=3000, customGenelist=NULL, min.loctions = 20,
                                min.features=20)
  
  
  saveRDS(LIBD, file = paste0("LIBD_sample", ksample, "_spatialPcaObj.rds"))
  


  
  
}



