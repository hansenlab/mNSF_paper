# ml conda_R/4.3
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

dir_out = "/dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_keepAll_scTransform//"
################# 
dir_output="/dcs04/hansen/data/ywang/ST/disp/"
dir.create(dir_output)
setwd(dir_output)
# for disp in ['0001', '001', '01', 1, 10, 100]:# still need to run： 001, 01
####### load ARI of mNSF and spatialPCA on dlpfc 3-sample data
list_filenames = c('0001', '001', '01', 1, 10, 100)
list_disp = c(0.001, 0.01, 0.1, 1, 10, 100)
for(k in c(1:2,4:6)){
  disp= list_disp[k]
  disp_name = list_filenames[k]
  
  dev_unexplained = read.csv(paste0(dir_out,"vec_dev_500selectedFeatures_dev_interpolated_35percent_szMean_disp",disp_name,".csv"))
  
  df_tmp = data.frame(dev_unexplained = dev_unexplained[,2], 
                     sample = (paste0("sample ", 1:12)), 
                     disp = disp)
  
  if(k ==1 ){
    df_dev_unexplained = df_tmp
  }else{
    df_dev_unexplained = rbind(df_dev_unexplained, df_tmp)
  }
}

df_dev_unexplained$sample = factor (df_dev_unexplained$sample, levels = paste0("sample ",1:12))

#######  make plot
pdf("var_unexplained_disp_dlpfc.pdf", height = 5, width = 10)
ggplot(df_dev_unexplained, aes(y=dev_unexplained, x= log10(disp),
                        # fill=sample,
                        col=sample))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
  geom_jitter(size=2,width= 0.25)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
