
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
# for(k in 1:10){
#   st_ =    k*2500-2500
#   if(k!=10){
#     end_ =    k*2500
#     
#   }else{
#     end_ = ncol(X_xenium)
#   }
#   
#   # list_X = c(list_X, X_xenium[st_:end_,])
#   list_X[[k+1]] =  X_xenium[st_:end_,]
# }

# Factors_df.to_csv(path.join(dir_output,"factors_nb_szMean_sample_s"+str(k+1)+"_L10_saveMem.csv"))
####### domain
sample_vec = c()
for(k in 1:11){
  ksample =  k 
  a_=read.csv(paste0("factors_nb_szMean_sample_s",ksample,"_L10_saveMem.csv"),header=T)
  a_ = a_[,-1]
  # print(dim(a_))
  
  if(k==1){a_combined = a_}else{
    a_combined = rbind(a_combined, a_)
  }
  
  ksample = ifelse(k==1,1,2)
  sample_vec = c(sample_vec,rep(ksample,nrow(a_)))
}


# library(WeightedCluster)
# library(wskm)
library(FactoClass)

# FactoClass
set.seed(111)
# clusterlabel= kmeans(scale(a_combined),7)$cluster
vec_weights = array(dim=length(sample_vec))
vec_weights[sample_vec==1] = 1/sum(sample_vec==1)
vec_weights[sample_vec==2] = 1/sum(sample_vec==2)

clusterlabel = kmeansW(scale(a_combined), 7, weight = vec_weights)$cluster

saveRDS(clusterlabel,file = "clusterlabel_L10_2samples_saveMem_7clu.rds")



library(SpatialPCA)
for(k in 1:2){
  print(k)
  X = list_X[[k]][sample_vec]
  layers_pred_refined = refine_cluster_10x(as.character(clusterlabel[sample_vec==k]),X,shape="square")
  saveRDS(layers_pred_refined,file = "clusterlabel_L10_2samples_saveMem_7clu_refined.rds")

}

# clusterlabel= louvain_clustering(clusternum=8,latent_dat=t(a_combined),knearest=round(sqrt(dim(a_combined)[2])) )
# clusterlabel = readRDS("clusterlabel_L10_2samples_saveMem_7clu.rds")
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



# library(SpatialPCA)
pdf(paste0("cere_domain_7clu_saveMem.pdf"),height=5,width=5)
for(k in 1:2){
  
  X = list_X[[k]]
  
  
  cbp = c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
                     "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
                     "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,
                     "#6B6ECF", "#7B4173" )
  layers_pred_sampleTmp = clusterlabel[sample_vec==k]
  
  
  # layers_pred_sampleTmp=refine_cluster_10x(as.character(layers_pred_sampleTmp),X,shape="square")
  
  
  p1= make_plot(X[,], as.character(layers_pred_sampleTmp), paste0("predicted layer by mNSF, sample ",k ),
                group.colors=cbp  ,
                size_ = ifelse(k==1,2,0.6))
  
  
  print(p1 )

}
dev.off()




