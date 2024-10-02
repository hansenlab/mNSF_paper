

############################## 
library(SpatialPCA)

############################## 
dir_out  = "/dcs04/hansen/data/ywang/revision/domain/"
setwd(dir_out)

dir_Data = "/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/"
dir_dlpfc_data = paste0(dir_Data,"DLPFC/")
############################## use dlpfc data with the same genes selected in mnsf paper
genes_svg = readRDS("genes_svg_dlpfc.rds")
############################## 
for(set in 2:4){
  genes_svg= readRDS(paste0("genes_svg_dlpfc_set",set,".rds"))

  for(kdonor in 1:3){
    ksample = (kdonor-1)*4+ set
    
    print(ksample)
    load(paste0(dir_dlpfc_data, "LIBD_sample",ksample,".RData"))
    
    Y = t(data.matrix(count_sub[genes_svg,]))
    write.csv(Y,row.names = F,#col.names = F,
              file=paste0("Y_ksample",ksample,"_svg.csv"))
    write.csv(data.matrix(xy_coords),row.names = F,#col.names = F,
              file=paste0("X_ksample",ksample,"_svg.csv"))
    
  }
  
}


############################## 


