
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

dir_dev = "/dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_keepAll_scTransform//"

################# 
dir_output="/dcs04/hansen/data/ywang/ST/nfactor/"
# dir.create(dir_output)
setwd(dir_output)
list_factor = c(4,8,12,16,20)

####### load ARI of mNSF and spatialPCA on dlpfc 3-sample data
for(nfactor in list_factor){
  dev_unexplained = read.csv(paste0(dir_dev,"vec_dev_500selectedFeatures_dev_interpolated_35percent_szMean_L",nfactor,"_010Induced.csv"))
  
  df_tmp = data.frame(dev_unexplained = dev_unexplained[,2], 
                      sample = (paste0("sample ", 1:12)), 
                      nfactor = nfactor)
  
  if(nfactor == list_factor[1]){
    df_dev_unexplained = df_tmp
  }else{
    df_dev_unexplained = rbind(df_dev_unexplained, df_tmp)
  }
}
# list_dev_c = readRDS("list_dev_c.rds")

# df_dev_unexplained = 
df_dev_unexplained$sample = factor (df_dev_unexplained$sample, levels = paste0("sample ",1:12))


#######  make plot
pdf("var_unexplained_nfactor_dlpfc_010Induced.pdf", height = 5, width = 10)
ggplot(df_dev_unexplained, aes(y=dev_unexplained, x= nfactor,
                               # fill=sample,
                               col=sample))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
  geom_jitter(size=2,width= 0.25)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          text=element_text(size=21),
                                          
                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
