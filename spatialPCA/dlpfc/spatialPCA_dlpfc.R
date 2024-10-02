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


############################## use dlpfc data with the same genes selected in mnsf paper

list_count = list()
list_location  = list()

for(kdonor in 1:3){
  ksample = (kdonor-1)*4+ 1
  
  X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))
  Y=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///Y_features_sele_sample',ksample,'_500genes.csv'))
  
  rownames(Y) = rownames(X) = paste0("donor",kdonor,"_",rownames(Y))
  list_count[[kdonor]] = t(Y)
  list_location[[kdonor]] = data.matrix(X)
  
  
}
# location matrix: n x 2, count matrix: g x n.
library(matrixStats)
(colVars(data.matrix(X)))
# pos_col_sample_tmp_ pos_row_sample_tmp_ 
# 9053.606           11359.037 

# Y=pd.read_csv(path.join('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data//Y_features_sele_sample'+str(ksample*4  + 1)+'_500genes.csv'))
out_SpatialPCA = SpatialPCA_Multiple_Sample(count_list=list_count, location_list=list_location, 
                                            gene.type="spatial",sparkversion="spark",
                                            # gene.type="hvg",sparkversion="spark",
                                            
                           numCores_spark=5,gene.number=500, customGenelist=NULL, min.loctions = 20, min.features=20,bandwidth_common=0.1)

saveRDS(out_SpatialPCA, file = "out_SpatialPCA_3samples.rds")
  
# ######################### 
# ######################### load spatialPCA results and data
# ######################### 
# ########### specify spatialPCA data locations
# dir_Data = "/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/"
# dir_dlpfc_data = paste0(dir_Data,"DLPFC/")
# dir_SlideseqCerebellum_data = paste0(dir_Data,"SlideseqCerebellum/")
# dir_SlideseqV2Hippocampus_data = paste0(dir_Data,"SlideseqV2Hippocampus/")
# # dir_BreastTumor_data = paste0(dir_Data,"BreastTumor/")
# # dir_Vizgen_data = paste0(dir_Data,"Vizgen/")
# 

# ########### load spatialPCA  data for dlpfc
# # load DLPFX data
# load(paste0(dir_dlpfc_data,"slideseq.rds") )
# # print(dim(sp_count)) # The count matrix
# # print(dim(location)) # The location matrix
# 
# ########### load spatialPCA result for dlpfc analysis
# load(paste0(dir_Data,"/Source_data_for_analysis/sourcedata/LIBD/SpatialPCA_LIBD_sample_2.RData") )
# ls(SpatialPCA_result)
# # SpatialPCs
# 
# 
# 
