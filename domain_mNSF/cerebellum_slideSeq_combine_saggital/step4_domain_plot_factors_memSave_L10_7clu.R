
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
# list_X[[2]] = X_xenium
for(k in 1:10){
  st_ =    k*2500-2500
  if(k!=10){
    end_ =    k*2500

  }else{
    end_ = ncol(X_xenium)
  }

  # list_X = c(list_X, X_xenium[st_:end_,])
  list_X[[k+1]] =  X_xenium[st_:end_,]
}

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

# set.seed(111)
# library(WeightedCluster)
# library(wskm)
# library(FactoClass)

# FactoClass
# clusterlabel= kmeans(scale(a_combined),7)$cluster
vec_weights = array(dim=length(sample_vec))
vec_weights[sample_vec==1] = 1/sum(sample_vec==1)
vec_weights[sample_vec==2] = 1/sum(sample_vec==2)

###
# library(igraph)
# # g <- make_graph("Zachary")
# create_graph_from_adj_matrix <- function(adj_matrix) {
#   g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag=F,weighted = TRUE)
#   return(g)
# }

# graph_=create_graph_from_adj_matrix(data.matrix(dist((scale(a_combined[,])))))
# graph_=simplify(graph_)
# louvain_clusters <- cluster_louvain(graph_,resolution=1.05)
# print(louvain_clusters)

# table(louvain_clusters$membership)
# Get the membership vector (which node belongs to which community)
# membership <- membership(louvain_clusters)
# print(membership)
# library(SpatialPCA)
# clusterlabel= louvain_clustering(clusternum=5,latent_dat=t(scale(a_combined)),
                                 # knearest=round(sqrt(dim(a_combined)[2])) )
# 
# # clusterlabel = kmeans(scale(a_combined), 5)$cluster
# ######
# library(Seurat)
# obj_ = CreateSeuratObject(t(100000*(a_combined)))
# obj_ = NormalizeData(obj_)
# obj_ = ScaleData(obj_)
# obj_ = FindVariableFeatures(obj_)
# obj_ = RunPCA(obj_)
# # obj_[['pca']][1,] = 
# # obj_[["RNA"]]$data = t(scale(a_combined))
# 
# mat <- (scale(a_combined))[,1:9]
# colnames(mat) = colnames(obj_[['pca']]@cell.embeddings)
# rownames(mat) = rownames(obj_[['pca']]@cell.embeddings)
# 
# obj_[['pca']]@cell.embeddings  = mat
# obj_ = FindNeighbors(obj_, k.param=50, dims=9)
# obj_ = FindClusters(obj_,res=.1)
# clusterlabel=(Idents(obj_))
# table(Idents(obj_))

######

library(SWKM)
set.seed(111)
clusterlabel = kmeans.weight(scale(a_combined), 5, weight =vec_weights)$cluster

# clusterlabel = kmeans.weight(scale(a_combined), 5, weight = rep(1,length(vec_weights)))$cluster
# clusterlabel = (clusterlabel$cluster)
clusterlabel = as.vector(clusterlabel)

saveRDS(clusterlabel,file = "clusterlabel_L10_2samples_saveMem_5clu.rds")



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
pdf(paste0("cere_domain_5clu_saveMem.pdf"),height=5,width=5)
for(k in 1:2){

  X = list_X[[k]]
  
  cbp = c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
                     "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
                     "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,
                     "#6B6ECF", "#7B4173" )
  layers_pred_sampleTmp = clusterlabel[sample_vec==k]
  p1= make_plot(X[,], as.character(layers_pred_sampleTmp), paste0("predicted layer by mNSF, sample ",k ),
                # group.colors=c("1" = "#FC8D62", "2" = "#FFD92F", "3" ="#A6D854", "4" = "#66C2A5", "5" = "#00A9FF",
                               # "6"="#8DA0CB", "7"="#E78AC3", "8" = "brown")   ,
                group.colors=cbp  ,
                size_ = ifelse(k==1,2,2))
  print(p1 )
  
}
dev.off()

#  


