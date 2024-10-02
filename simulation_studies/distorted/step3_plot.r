


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

setwd("/users/ywang/Hansen_projects/mNSF/revision/simulation/")

for(type in c("distorted","sizeDif", "noiseDif"))
{
  sample_=1
  a_s1=read.csv(paste0("sim_rotated_factors_sample1_",type,".csv"))
  a_s1=a_s1[,-1]
  X_s1 = read.csv(paste0("X_sample",1,"_",type,".csv"))
  X_s1 = X_s1[,-1]
  
  sample_=2
  a_s2=read.csv(paste0("sim_rotated_factors_sample2_",type,".csv"))
  a_s2=a_s2[,-1]
  X_s2 = read.csv(paste0("X_sample",2,"_",type,".csv"))
  X_s2 = X_s2[,-1]
  
  valueRange_ = range(c(as.numeric(data.matrix(a_s1)),as.numeric(data.matrix(a_s2))))
  coordX_Range_ = range(c(as.numeric(data.matrix(X_s1[,1])),
                          as.numeric(data.matrix(X_s2[,1]))))
  coordY_Range_ = range(c(as.numeric(data.matrix(X_s1[,2])),
                          as.numeric(data.matrix(X_s2[,2]))))
  
  for (sample_ in 1:2){
    a_true=read.csv(paste0("ggblocks_lr_1_", type, "_sample",sample_,".csv"))
    a_true=a_true[,-1]
    
    a=read.csv(paste0("sim_rotated_factors_sample",sample_,"_",type,".csv"))
    a=a[,-1]
    
    X = read.csv(paste0("X_sample",sample_,"_",type,".csv"))
    # X=read.csv(paste0("/dcs04/legacy-dcs01-hansen/hansen_lab1/ywang/ST/May6_2022_sNMF_2samples/nsf-paper-main_2samples_LukasData_sample2_3_qsub_fasterVersion_sub500_v2_topleft_500genes_layer1to4//X",sample_,"_.csv"))
    X=X[,-1]
    
    plots_l=list()
    
    
    #### plot mNSF factors
    for (i in 1:4) {
      dim(a)
      df_tmp=data.frame(imagecol=X[,1],
                        imagerow=X[,2],
                        fill_tmp=a[,i])
      
      colors <- c('navyblue', 'skyblue1')
      # size_ = ifelse(type=="noiseDif",3,6)
      size_=6
      plot_tmp =  ggplot(df_tmp,aes(x=imagecol,y=imagerow,fill=fill_tmp)) +
        # geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
        geom_point(shape = 22,  size = size_, stroke = NA)+
        coord_cartesian(expand=FALSE)+
        scale_fill_gradientn(colors = colors, limits = valueRange_)+
        xlab("") +
        ylab("") +
        xlim(coordX_Range_)+
        ylim(coordY_Range_)+
        ggtitle(paste0("M",i))+
        labs(fill = paste0(" "))+
        # theme_set(theme_bw(base_size = 10))+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              # axis.line = element_line(colour = "black"),
              axis.line =  element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = "none")+ coord_equal()
      
      plots_l[[i]]=plot_tmp
      
      
    }
    
    pdf(paste0("mNSF_sim_",type,"_sample",sample_,".pdf"),height=2.6,width=18/6*4.5)
    print(plot_grid(plotlist = plots_l,nrow=1))
    dev.off()
  }
  
  
  
}
