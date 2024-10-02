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


################# 
dir_output="/dcs04/hansen/data/ywang/ST/data_10X_ST/mouse_Sagittal_spaceRanger1_1_0/out/"
setwd(dir_output)
####### load ARI of mNSF and spatialPCA on dlpfc 3-sample data
list_filenames = c('0001', '001', '01', 1, 10, 100)
list_disp = c(0.001, 0.01, 0.1, 1, 10, 100)

list_filenames = c('0001', '001', '01', 10, 100)
list_disp = c(0.001, 0.01, 0.1,  10, 100)


get_full_factor = function(kfactor_){
  for(disp_name in list_filenames){
    for(ksample in 1:4){
      
      if(disp_name=="001"){
        a_tmp=read.csv(paste0("factors_sample",
                              ksample,
                              "_500selectedFeatures_dev_interpolated_35percent_szMean.csv"),header=T)
        
      }else{
        a_tmp=read.csv(paste0("factors_sample",
                              ksample,
                              "_500selectedFeatures_dev_interpolated_35percent_szMean_disp",disp_name,".csv"),header=T)
        # _500selectedFeatures_dev_interpolated_35percent_szMean_disp1
      }
        a_tmp = a_tmp[,-1]
      
      if(ksample == 1){a_c = (a_tmp[,kfactor_])}else{
        a_c = c(a_c, (a_tmp[,kfactor_]))
      }
      
    }
    if(disp_name == list_filenames[[1]]){a_mat = data.matrix(a_c)}else{
      a_mat = cbind(a_mat, data.matrix(a_c))
    }
  }
  return(a_mat)
}


list_factor_mat = list()
for(kfactor in 1:20){
  print(kfactor)
  list_factor_mat[[kfactor]]= get_full_factor(kfactor)
}
for(kfactor in 1:20){
  print(kfactor)
  # list_factor_mat[[kfactor]]= get_full_factor(kfactor)
  colnames(list_factor_mat[[kfactor]]) = paste0("disp=",list_disp)
}


saveRDS(list_factor_mat, file = "list_factor_mat_mouse.rds")



#######  make plot
library(ggplot2)
library(GGally)
library(ggrastr)
my_custom_cor <- function(data, mapping, method = "pearson", ndp = 2, sz = 5, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  cor <- cor(x, y, method = method)
  cor_text <- sprintf("%.2f", cor)
  
  ggplot(data = data, mapping = mapping) +
    geom_point_rast(alpha = 0.3, raster.dpi = 300) +  # Rasterized points
    geom_text(x = mean(range(x)), y = mean(range(y)), 
              label = cor_text, size = sz, color = "red")
}
# Create the pair plot
dir.create("scatter_pairs_disp_mouse")

for(kfactor in 1:20){
  
  pdf(paste0("scatter_pairs_disp_mouse/scatter_pairs_disp_factor",kfactor,"_mouse.pdf"))
  print(kfactor)
  p <- ggpairs((data.frame(list_factor_mat[[kfactor]]))[,], 
               lower = list(continuous = my_custom_cor),
               diag = list(continuous = "densityDiag"),
               upper = list(continuous = "blank")) + ggtitle(paste0("factor ",kfactor))
  print(p)
  dev.off()
  
}







