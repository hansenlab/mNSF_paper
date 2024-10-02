
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
dir_output="/dcs04/hansen/data/ywang/ST/nfactor/"
# dir.create(dir_output)
setwd(dir_output)

####### load ARI of mNSF and spatialPCA on dlpfc 3-sample data
list_factor = c(4,8,12,16,20)

list_acc_c = readRDS("list_acc_c_010Induced.rds")
list_dev_c = readRDS("list_dev_c_010Induced.rds")

df_plot_acc = data.frame(acc = unlist(list_acc_c), 
                         sample = rep(paste0("sample ", 1:12), length(list_factor)), 
                         nfactor = rep(list_factor, each = 12))
df_plot_acc$sample = factor (df_plot_acc$sample, levels = paste0("sample ",1:12))
df_plot_dev = data.frame(dev = unlist(list_dev_c), 
                         sample = rep(paste0("sample ", 1:12), length(list_factor)), 
                         nfactor = rep(list_factor, each = 12))
df_plot_dev$sample = factor (df_plot_dev$sample, levels = paste0("sample ",1:12))

#######  make plot
pdf("acc_nfactor_010Induced.pdf", height = 5, width = 10)
ggplot(df_plot_acc, aes(y=acc, x= nfactor,
                # fill=sample,
                col=sample))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_boxplot()+
  geom_jitter(size=2,width= 0.25)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          text=element_text(size=21),
                                          
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()



pdf("dev_nfactor_010Induced.pdf", height = 5, width = 10)
ggplot(df_plot_dev, aes(y=dev, x= nfactor,
                        # fill=sample,
                        col=sample))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
  geom_jitter(size=2,width= 0.25)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          text=element_text(size=21),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
