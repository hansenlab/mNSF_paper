


### load packages 
library(tidyverse)
library(ggplot2)
library(Matrix)
# library(Rmisc)#
library(ggforce)#
# library(rjson)#
library(cowplot)
library(RColorBrewer)
library(grid)
# library(readbitmap)#
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


sample_=1
a_s1=read.csv(paste0("/dcs04/legacy-dcs01-hansen/hansen_lab1/ywang/ST/April8_2022_NSF/nsf-paper-main_9_ori/simulations/ggblocks_lr/results/sim_SS_factors_sample",sample_,".csv"))
a_s1=a_s1[,-1]
sample_=2
a_s2=read.csv(paste0("/dcs04/legacy-dcs01-hansen/hansen_lab1/ywang/ST/April8_2022_NSF/nsf-paper-main_9_ori/simulations/ggblocks_lr/results/sim_SS_factors_sample",sample_,".csv"))
a_s2=a_s2[,-1]
valueRange_ = range(c(range(a_s1),range(a_s2)))
range(a_s1[,2])
valueRange_
sum(a_s2)

for (sample_ in 1:2){
  
  # dim(a)
  
  a=read.csv(paste0("/dcs04/legacy-dcs01-hansen/hansen_lab1/ywang/ST/April8_2022_NSF/nsf-paper-main_9_ori/simulations/ggblocks_lr/results/sim_SS_factors_sample",sample_,".csv"))
  a=a[,-1]

  X=read.csv(paste0("/dcs04/legacy-dcs01-hansen/hansen_lab1/ywang/ST/May6_2022_sNMF_2samples/nsf-paper-main_2samples_LukasData_sample2_3_qsub_fasterVersion_sub500_v2_topleft_500genes_layer1to4//X",sample_,"_.csv"))
  X=X[,-1]
  
  plots_l=list()
  
  
  for (i in 1:4) {
    dim(a)
    df_tmp=data.frame(imagecol=X[,1],
                      imagerow=X[,2],
                      fill_tmp=a[,i])
    
    colors <- c('navyblue', 'skyblue1')
    plot_tmp =  ggplot(df_tmp,aes(x=imagecol,y=imagerow,fill=fill_tmp)) +
      # geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape = 22,  size = 3, stroke = NA)+
      coord_cartesian(expand=FALSE)+
      scale_fill_gradientn(colors = colors, limits = valueRange_)+
      xlab("") +
      ylab("") +
      ggtitle(paste0("M",i))+
      labs(fill = paste0(" "))+
      theme_set(theme_bw(base_size = 10))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none"
      )+ coord_equal()
    
    plots_l[[i]]=plot_tmp
    
    
  }
  
  
  
  # pdf(paste0("/dcs04/legacy-dcs01-hansen/hansen_lab1/ywang/ST/May6_2022_sNMF_2samples/nsf-paper-main_2samples_rotate_regularSizeData_SSlayer_modZ_nSamples_memorySaving_byBatches/plot__fullData_50iterations/*LdaFalse_L3.pdf /home/yi/Downloads//sample__",sample_,"_LdaFalse_L3.pdf"),height=3,width=18/6*4)
  pdf(paste0("/dcs04/hansen/data/ywang/ST/sim/mNSF_sim_SS_sample",sample_,"_final.pdf"),height=2.6,width=18/6*4.5)
  print(plot_grid(plotlist = plots_l,nrow=1))
  dev.off()
}

