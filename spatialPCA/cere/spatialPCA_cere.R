# ml conda_R/4.4
# R

# conda create -n  spatialPCA -- also works
# conda activate spatialPCA
# conda install conda-forge::r-devtools 

library(SpatialPCA)
# spatialPCA package pacakge: https://github.com/shangll123/SpatialPCA
# spatialPCA analysis code:https://github.com/shangll123/SpatialPCA_analysis_codes
# spatialPCA tutorial: https://lulushang.org/SpatialPCA_Tutorial/

source("/users/ywang/Hansen_projects/mNSF/revision/spatialPCA/cere/function.R")
############################## 
dir_out  = "/dcs04/hansen/data/ywang/revision/domain/"
setwd(dir_out)


############################## use dlpfc data with the same genes selected in mnsf paper
dir_Data = "/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/"
# dir_dlpfc_data = paste0(dir_Data,"DLPFC/")
dir_SlideseqCerebellum_data = paste0(dir_Data,"SlideseqCerebellum/")
list_count = list()
list_location  = list()

# load cere slide-seq data
kdonor = 1 
X=read.csv(paste0(dir_SlideseqCerebellum_data, 'X.csv'))
Y=read.csv(paste0(dir_SlideseqCerebellum_data, 'Y.csv'))
rownames(Y) = rownames(X) = paste0("Sample",kdonor,"_",rownames(Y))
list_count[[kdonor]] = t(Y[,])
list_location[[kdonor]] = data.matrix(X)[,]

kdonor = 2
Y=read.csv(('/dcs04/hansen/data/ywang/ST/data_10X_ST//mouse_Sagittal/put/Y_features_sele_sample3_v2_460genes.csv')) #_500genes
X=read.csv(('/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/X_sample3.csv'))
rownames(Y) = rownames(X) = paste0("Sample",kdonor,"_",rownames(Y))

list_count[[kdonor]] = t(Y[,])
list_location[[kdonor]] = data.matrix(X)[,]
# (colVars(data.matrix(X)))
# x       y 
# 3855219 5895616 
############################## 
# list_count[[1]]=list_count[[2]]
# list_location[[1]] = list_location[[2]]

# count=1
# (out_SpatialPCA$SpatialPC_list)
# :vacant
# 
out_SpatialPCA = SpatialPCA_Multiple_Sample(count_list=list_count, location_list=list_location, 
                                            gene.type="hvg",sparkversion="spark",
                           numCores_spark=5, gene.number=460, customGenelist=NULL, 
                           min.loctions = 0, 
                           min.features=0, bandwidth_common=0.1,SpatialPCnum_=10)

saveRDS(out_SpatialPCA, file = "out_SpatialPCA_2samples_cere_L10.rds")
  


