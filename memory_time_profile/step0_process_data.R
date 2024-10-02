# dir_out="/users/ywang/Hansen_projects/scRNA/Oct5_2021_Lukas_data_more_Genes/out/"
# setwd(dir_out)

dir_processedData = "/dcs04/hansen/data/ywang/ST/DLPFC/processed_Data/"
setwd(dir_processedData)

# for( j in 1){
j=1
  print(j)
  
  counts_sample_tmp_ = read.csv(paste0(dir_processedData,"Y_alllGenes_2percentSparseFilter_sample",j,".csv"))

  print(dim(counts_sample_tmp_))
  
  # write.csv(X_allSpots,
            # file=paste0("X_allSpots_sample",j,".csv"),row.names = F)
  
  X_allSpots = read.csv(paste0(dir_processedData,"X_allSpots_sample",j,".csv"))
  print(dim(X_allSpots))
  # [1] 4226    2

# }
# >   print(dim(counts_sample_tmp_))
# [1]  4226 13122

write.csv(counts_sample_tmp_[1:500,1:500],row.names = F,#col.names = F,
            file=paste0("Y_small_forMemTimeProfile.csv"))
write.csv(X_allSpots[1:500,],row.names = F,#col.names = F,
          file=paste0("X_small_forMemTimeProfile.csv"))
######
write.csv(counts_sample_tmp_[],row.names = F,#col.names = F,
          file=paste0("Y_forMemTimeProfile.csv"))
write.csv(X_allSpots[,],row.names = F,#col.names = F,
          file=paste0("X_forMemTimeProfile.csv"))

###### 500 genes
write.csv(counts_sample_tmp_[,1:500],row.names = F,#col.names = F,
          file=paste0("Y_500genes_forMemTimeProfile.csv"))

###### 500 spots
write.csv(counts_sample_tmp_[1:500,],row.names = F,#col.names = F,
          file=paste0("Y_500spots_forMemTimeProfile.csv"))

