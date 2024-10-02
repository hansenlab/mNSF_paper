# ml conda_R/4.4
# R

# conda create -n  spatialPCA -- also works
# conda activate spatialPCA
# conda install conda-forge::r-devtools 

library(SpatialPCA)
# spatialPCA package pacakge: https://github.com/shangll123/SpatialPCA
# spatialPCA analysis code:https://github.com/shangll123/SpatialPCA_analysis_codes
# spatialPCA tutorial: https://lulushang.org/SpatialPCA_Tutorial/

dir_out  = "/dcs04/hansen/data/ywang/revision/domain/"
setwd(dir_out)

# ######################### 
# ######################### load spatialPCA results and data
# ######################### 
# ########### specify spatialPCA data locations
dir_Data = "/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/"
# dir_dlpfc_data = paste0(dir_Data,"DLPFC/")
dir_SlideseqCerebellum_data = paste0(dir_Data,"SlideseqCerebellum/")
# dir_SlideseqV2Hippocampus_data = paste0(dir_Data,"SlideseqV2Hippocampus/")
# # dir_BreastTumor_data = paste0(dir_Data,"BreastTumor/")
# # dir_Vizgen_data = paste0(dir_Data,"Vizgen/")
# 

############ load spatialPCA  data for dir_SlideseqCerebellum_data
# # load dir_SlideseqCerebellum_data data
load(paste0(dir_SlideseqCerebellum_data,"slideseq.rds") )
print(dim(sp_count)) # The count matrix
# [1] 17729 25551

print(dim(location)) # The location matrixxc
# [1] 25551     2

########### use the genes selected by mouse saggital section data
# load data
counts_mouseSaggital=read.csv(('/dcs04/hansen/data/ywang/ST/data_10X_ST//mouse_Sagittal/put/Y_features_sele_sample3_v2_500genes.csv')) #_500genes
# list_subdir[[1]]="V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix"
# list_subdir[[2]]="V1_Mouse_Brain_Sagittal_Anterior_Section_2_filtered_feature_bc_matrix"
# list_subdir[[3]]="V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix"
# list_subdir[[4]]="V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix"
# X=pd.read_csv(path.join(dpth,'/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/X_sample'+ str(ksample+1) +'.csv'))
genes_sele = intersect(colnames(counts_mouseSaggital), rownames(sp_count))
length(genes_sele  )
# 460

counts_mouseSaggital_genesShared = counts_mouseSaggital[,genes_sele]

write.csv(counts_mouseSaggital_genesShared,row.names = F,#col.names = F,
          file=paste0('/dcs04/hansen/data/ywang/ST/data_10X_ST//mouse_Sagittal/put/Y_features_sele_sample3_v2_460genes.csv'))

sum(colSums(data.matrix(sp_count)[genes_sele,])==0)
# 136

selec_spots = which(colSums(data.matrix(sp_count)[genes_sele,])>0)
write.csv(t(data.matrix(sp_count)[genes_sele,selec_spots]),row.names = F,#col.names = F,
          file=paste0(dir_SlideseqCerebellum_data, 'Y.csv'))

write.csv(location[selec_spots,],row.names = F,#col.names = F,
          file=paste0(dir_SlideseqCerebellum_data, 'X.csv'))

########### 
# 
########### 
library(matrixStats)
summary(colVars(data.matrix(sp_count)[genes_sele,]))
summary(colVars(t(counts_mouseSaggital_genesShared)[,]))

########### 
locs_mouseSaggital=read.csv(('/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out//X_sample3.csv')) #_500genes

summary(colVars(data.matrix(location)[,]))
summary(colVars(data.matrix(locs_mouseSaggital)[,]))
summary(colVars(data.matrix(locs_mouseSaggital)[,]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1219365 1317481 1415596 1415596 1513712 1611827 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3855219 4365318 4875417 4875417 5385517 5895616 


summary(colSums(data.matrix(sp_count)[genes_sele,]))
summary(colSums(t(counts_mouseSaggital_genesShared)[,]))
# > summary(colSums(data.matrix(sp_count)[genes_sele,]))
# summary(colSums(t(counts_mouseSaggital_genesShared)[,]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    7.00   16.00   26.63   33.00  717.00 
# Warning message:
#   In asMethod(object) :
#   sparse->dense coercion: allocating vector of size 3.4 GiB
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4    2212    3310    3782    4849   28044 

########### 

