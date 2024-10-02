
### load pacakges 
library(tidyverse)
library(ggplot2)
library(Matrix)
library(ggforce)#
library(cowplot)
library(RColorBrewer)
library(grid)
require(nnet)

#######
myPalette_ <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# https://ggplot2-book.org/scale-colour.html
myPalette = scale_fill_brewer(palette = "Set2")

#######
setwd('/dcs04/hansen/data/ywang/revision/SlideseqCerebellum/')
dir_Data = "/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/"
dir_SlideseqCerebellum_data = paste0(dir_Data,"SlideseqCerebellum/")
#######


X_visium=read.csv(('/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/X_sample3.csv'))
X_visium[,1] = -X_visium[,1]
X_xenium=read.csv(paste0(dir_SlideseqCerebellum_data, 'X.csv'))
X_xenium = X_xenium[,c(2,1)]

list_X = list()
list_X[[1]] = X_visium
list_X[[2]] = X_xenium



#######
make_plot<-function(X,factor_mat,range_perFactor_,
                    samplename_){
  plots_l = list()
  for (i in 1:ncol(factor_mat)) {
    # dim(a)
    df_tmp=data.frame( imagerow=as.numeric(X[,1]),
                       imagecol=as.numeric(X[,2]),
                       fill_tmp=factor_mat[,i])

    plot_tmp =  ggplot(df_tmp,aes(x=imagecol,y=imagerow,fill=fill_tmp)) +
      # geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape = 21, colour = "black", size = 1, stroke = NA)+
      coord_cartesian(expand=FALSE)+
      scale_fill_gradientn(colours = myPalette_(100), limits=range_perFactor_[[i]])+
      
      xlab("") +
      ylab("") +
      ggtitle(paste0(samplename_,", mNSF factor ",i))+
      labs(fill = paste0(" "))+
      theme_set(theme_bw(base_size = 10))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    plots_l[[i]]=plot_tmp
    
  }

  plots_l
}

# Factors_df.to_csv(path.join("factors_sample"+str(sample1)+"_"+str(sample1+1)+"_s"+str(k+1)+"_L3.csv"))
range_perFactor = list()
for(ksample in 1:11){
  factor_mat_tmp=read.csv(paste0("factors_nb_szMean_sample_s",ksample,"_L10_saveMem.csv"),header=T)
  
  # factor_mat_tmp=read.csv(paste0("factors_nb_szMean_sample_s",ksample,"_L10_2samples_0.10InducedSlideSeq.csv"),header=T)
  factor_mat_tmp=factor_mat_tmp[,-1]
  for(l in 1:ncol(factor_mat_tmp)){
    if(ksample==1){ range_perFactor[[l]]=range(factor_mat_tmp[,l])}else{
      range_perFactor[[l]]=range(c(range_perFactor[[l]], factor_mat_tmp[,l]))
    }
    
  }
}


pdf(paste0("LdaFalse_L10_SlideSeq_memSave_plotAllOnePlot.pdf"),height=3,width=18/6*30/2)
# for(kpair in 1:6){
for(k in 1:2){
    ksample =  k 
      if(k==1){
        a_=read.csv(paste0("factors_nb_szMean_sample_s",ksample,"_L10_saveMem.csv"),header=T)
      }

    if(k==2){
      for(ksplit in 1:10){
        a_split=read.csv(paste0("factors_nb_szMean_sample_s",(ksample-1)+ksplit,"_L10_saveMem.csv"),header=T)
        if(ksplit ==1){a_ = a_split}else{
          a_ = rbind(a_, a_split)
        }
      }
    }
    
    a_ = a_[,-1]
    
    p1= make_plot(list_X[[k]], a_, range_perFactor, paste0("sample ",ksample))
    print(plot_grid(plotlist = p1,nrow=1))
  
}
dev.off()
