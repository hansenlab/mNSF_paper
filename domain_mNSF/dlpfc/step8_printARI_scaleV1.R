
### load pacakges 
library(tidyverse)
library(ggplot2)
library(Matrix)
library(ggforce)#
library(cowplot)
library(RColorBrewer)
library(grid)
require(nnet)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

dir_processedData = "/dcs04/hansen/data/ywang/ST/DLPFC/processed_Data_dist005/"

################# 
dir_output="/dcs04/hansen/data/ywang/ST/DLPFC/threeSamples//"
# dir.create(dir_output)
setwd(dir_output)

####### load ARI of mNSF and spatialPCA on dlpfc 3-sample data
# process_data <- function()
# set1
ARI_mNSF = readRDS("ARI_vec_sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_mnsfrefined_svgs.rds")
ARI_spatialPCA = readRDS("ARI_vec_sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_spatialPCA_refined_svgs.rds")

df_ARI_mNSF = data.frame(ARI = ARI_mNSF, 
                         donor = (paste0("donor ", 1:3)), 
                         method = "mNSF",
                         set = "1")
df_ARI_spatialPCA = data.frame(ARI = ARI_spatialPCA, 
                         donor = (paste0("donor ", 1:3)), 
                         method = "spatialPCA",
                         set= "1")
df_ = rbind(df_ARI_mNSF, df_ARI_spatialPCA)

# set2-4
for( set in 2:4){
  ARI_mNSF = readRDS( paste0("ARI_vec_sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_mnsfrefined_scaleV1_3samples_034induced_set",set,".rds"))
  ARI_spatialPCA = readRDS( paste0("ARI_vec_sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_spatialPCA_refined_svgs_set",set,".rds"))
  # saveRDS(ARI_vec, file = paste0("ARI_vec_sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_mnsfrefined_scaleV2_3samples_034induced_set",set,".rds"))

  df_ARI_mNSF = data.frame(ARI = ARI_mNSF, 
                           donor = (paste0("donor ", 1:3)), 
                           method = "mNSF",
                           set= as.character(set))
  df_ARI_spatialPCA = data.frame(ARI = ARI_spatialPCA, 
                                 donor = (paste0("donor ", 1:3)), 
                                 method = "spatialPCA",
                                 set= as.character(set))
  df_ = rbind(df_, df_ARI_mNSF, df_ARI_spatialPCA)
  
}

df_$set = paste0("set ",df_$set)

#######  make plot
pdf("ARI_domainPred_dlpfc_svgs_3samples_3sets_ini_scaleV1.pdf", height = 2, width = 10)
ggplot(df_, aes(y=ARI, x= donor,
                fill=method,
                col=method))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_boxplot()+
  geom_point(size=3)+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~set,ncol=4)
dev.off()

