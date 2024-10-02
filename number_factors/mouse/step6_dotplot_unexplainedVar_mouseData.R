
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

# dir_processedData = "/dcs04/hansen/data/ywang/ST/DLPFC/processed_Data_dist005/"

dir_out = "/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/"
################# 
dir_output="/dcs04/hansen/data/ywang/ST/nfactor/"
# dir.create(dir_output)
setwd(dir_output)

####### load ARI of mNSF and spatialPCA on dlpfc 3-sample data
list_nfactor = c(8,16,24,32)
for(nfactor in list_nfactor){
  dev_unexplained = read.csv(paste0(dir_out,"vec_dev_500selectedFeatures_dev_interpolated_10percent_szMean_L",nfactor,".csv"))
  
  df_tmp = data.frame(dev_unexplained = dev_unexplained[,2], 
                     sample = (paste0("sample ", 1:4)), 
                     nfactor = nfactor)
  
  if(nfactor ==  list_nfactor[1]){
    df_dev_unexplained = df_tmp
  }else{
    df_dev_unexplained = rbind(df_dev_unexplained, df_tmp)
  }
}
# list_dev_c = readRDS("list_dev_c.rds")

# df_dev_unexplained = 
df_dev_unexplained$sample = factor (df_dev_unexplained$sample, levels = paste0("sample ",1:4))
# df_plot_dev = data.frame(dev = unlist(list_dev_c), 
#                          sample = rep(paste0("sample ", 1:12), 5), 
#                          nfactor = rep(c(2,4,6,8,10), each = 12))
# df_plot_dev$sample = factor (df_plot_dev$sample, levels = paste0("sample ",1:12))

#######  make plot
pdf("var_unexplained_nfactor_mouse_010induced.pdf", height = 5, width = 10)
ggplot(df_dev_unexplained, aes(y=dev_unexplained, x= nfactor,
                        # fill=sample,
                        col=sample))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
  geom_jitter(size=2,width= 0.25)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          text=element_text(size=21),
                                          
                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
