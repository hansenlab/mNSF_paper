
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

###############################

dir_out  = "/dcs04/hansen/data/ywang/revision/domain/"
setwd(dir_out)
# saveRDS(out_SpatialPCA, file = "out_SpatialPCA_2samples_cere.rds")
out_SpatialPCA = readRDS("out_SpatialPCA_2samples_cere_L10.rds")
# LS()
###############################
# setwd('/dcs04/hansen/data/ywang/revision/SlideseqCerebellum/')
dir_Data = "/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/"
dir_SlideseqCerebellum_data = paste0(dir_Data,"SlideseqCerebellum/")
#######

X_visium=read.csv(('/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/X_sample3.csv'))
X_visium[,1] = -X_visium[,1]
X_xenium=read.csv(paste0(dir_SlideseqCerebellum_data, 'X.csv'))
X_xenium = X_xenium[,c(2,1)]

list_X = list( X_xenium,X_visium)

# out_SpatialPCA_2samples_cere.rds
####### domain
sample_vec = c()
for(k in 1:2){
  # ksample =  k 

  a_ = t(out_SpatialPCA$SpatialPC_list[[k]][,])
  
  # if(k==1){
  #   set.seed(111)
  #   rd_ = sample(1:nrow(a_), 5000)
  #   a_ = a_[rd_,]
  # }
  if(k==1){a_combined = a_}else{
    a_combined = rbind(a_combined, a_)
  }
  sample_vec = c(sample_vec,rep(k,nrow(a_)))
}

set.seed(111)
clusterlabel= kmeans((a_combined),5)$cluster

# clusterlabel= louvain_clustering(clusternum=8,latent_dat=t(a_combined),knearest=round(sqrt(dim(a_combined)[2])) )
saveRDS(clusterlabel,file = "clusterlabel_L10_2samples_5clu_spatialPCA_L10.rds")

clusterlabel = readRDS("clusterlabel_L10_2samples_5clu_spatialPCA_L10.rds")
make_plot<-function(pos,
                    layer_sample_tmp___,title_,group.colors, size_){
  # plots_l = list()
  nspot = nrow(pos)
  df_tmp=data.frame( imagerow=pos[,1],
                     imagecol= pos[,2],
                     # fill_tmp=factor_mat[,i],
                     layer=layer_sample_tmp___)
  
  p=ggplot(df_tmp,aes(x=imagecol,y=imagerow,fill=layer)) +
    geom_point(shape = 21, colour = "black", size = size_, stroke = NA)+
    coord_cartesian(expand=FALSE)+
    xlab("") +
    ylab("") +
    
    ggtitle(title_)+
    
    labs(fill = paste0(" "))+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    
    scale_fill_manual(values=group.colors) +
    coord_equal()
  p
}

#####
pdf(paste0("cere_domain_5clu_spatialPCA_L10.pdf"),height=5,width=5)
for(k in 1:2){
  

  # a_ = t(out_SpatialPCA$SpatialPC_list[[k]][1:10,])
  
  
  X = list_X[[k]]
  cbp = c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
                     "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
                     "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,
                     "#6B6ECF", "#7B4173" )
  
  layers_pred_sampleTmp = clusterlabel[sample_vec==k]
  p1= make_plot(X[,], as.character(layers_pred_sampleTmp), paste0("predicted layer by mNSF, sample ",k ),
                group.colors=cbp , size_ = ifelse(k==2,2,0.6)
                )
  print(p1 )

}
dev.off()




