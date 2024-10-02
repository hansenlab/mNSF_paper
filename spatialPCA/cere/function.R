SpatialPCA_Multiple_Sample = function(count_list,location_list,gene.type="spatial",sparkversion="spark",numCores_spark=5,gene.number=3000, 
                                      customGenelist=NULL,min.loctions = 0, min.features=0,bandwidth_common=0.1,SpatialPCnum_=20){
  
  # create each SpatialPCA object, this step is essentially for spatial gene selection
  spatialpca_list = list()
  for(count in 1:length(count_list)){
    print(paste0("Creating SpatialPCA object for dataset ",count))
    spatialpca_list[[count]] <- CreateSpatialPCAObject(counts=count_list[[count]], location=location_list[[count]], project = "SpatialPCA",
                                                       gene.type=gene.type,sparkversion=sparkversion,numCores_spark=numCores_spark,gene.number=gene.number, 
                                                       customGenelist=customGenelist,min.loctions = min.loctions, min.features=min.features)
    print("6666")
    # print(spatialpca_list[[count]]@params$M )

  }
  
  
  print(dim(spatialpca_list[[1]]@normalized_expr) )
  print(dim(spatialpca_list[[2]]@normalized_expr) )
  # [1]  17 200
  # [1] 144 200
  
  require(Seurat)
  if(length(count_list) != length(location_list)) {
    stop("The number of count matrix should match with the number of location matrix. ")
  }# end fi
  
  # create seurat object for each dataset
  seurat_list = list()
  for(count in 1:length(count_list)){
    colnames(count_list[[count]]) = paste0("Sample",count,"_",colnames(count_list[[count]]))
    rownames(location_list[[count]]) = colnames(count_list[[count]])
    seurat_list[[count]] <- CreateSeuratObject(counts = count_list[[count]], project = "multiple")
  }
  
  # normalize and identify variable features for each dataset independently
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    print(paste0("Normalizing dataset"))
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list,nfeatures = 10000)
  # integration, find anchors
  MultipleSample.anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
  # this command creates an 'integrated' data assay
  MultipleSample.combined <- Seurat::IntegrateData(anchorset = MultipleSample.anchors)
  DefaultAssay(MultipleSample.combined) <- "integrated"
  # Run the standard workflow for visualization and clustering
  MultipleSample.combined <- Seurat::ScaleData(MultipleSample.combined, verbose = FALSE)
  # obtain integrated and normalized data
  integrated_data = MultipleSample.combined@assays$integrated@scale.data
  
  # print(spatialpca_list[[count]]@params$M )
  
  matched_data = list()
  matched_location = list()
  Kernal_mat = list()
  for(count in 1:length(count_list)){
    print(paste0("Creating SpatialPCA kernel for dataset ",count))
    match_gene = na.omit(match(rownames( spatialpca_list[[count]]@normalized_expr), 
                               rownames(integrated_data)))
    # match_spot = grep(paste0("1_",count),colnames(integrated_data))
    match_spot = grep(paste0("Sample",count),colnames(integrated_data)) #changed by YW
    
    matched_data[[count]] = integrated_data[match_gene,match_spot]
    matched_location[[count]] = spatialpca_list[[count]]@location
    spatialpca_list[[count]]@normalized_expr = t(scale(t(matched_data[[count]])))
    spatialpca_list[[count]]@location = scale(matched_location[[count]])
    spatialpca_list[[count]] = SpatialPCA_buildKernel(spatialpca_list[[count]], kerneltype="gaussian", bandwidth.set.by.user=bandwidth_common)
    Kernal_mat[[count]] = spatialpca_list[[count]]@kernelmat
  }
  

  dim(spatialpca_list[[1]]@params$expr)
  # [1]   7 200
  dim(spatialpca_list[[2]]@params$expr)
  # [1] 144 200
  common_genes = unique(unlist(lapply(matched_data, function(x){ row.names(x)})))
  MultipleSample_merge = spatialpca_list[[1]] # initialize using the first data, will replace with integrated data later
  MultipleSample_merge@normalized_expr = t(scale(t(integrated_data[match(common_genes, rownames(integrated_data)),])))
  MultipleSample_merge@kernelmat = Matrix::bdiag(Kernal_mat)
  MultipleSample_merge@kernelmat = as.matrix(MultipleSample_merge@kernelmat)
  # MultipleSample_merge@kernelmat = as(MultipleSample_merge@kernelmat, "sparseMatrix")
  # MultipleSample_merge@sparseKernel = TRUE
  MultipleSample_merge@params$expr = MultipleSample_merge@normalized_expr
  MultipleSample_merge@location = do.call("rbind",lapply(spatialpca_list, function(x){ x@location}))
   MultipleSample_merge = SpatialPCA_EstimateLoading(MultipleSample_merge,fast=TRUE,SpatialPCnum=SpatialPCnum_)

 
  MultipleSample_merge = SpatialPCA_SpatialPCs(MultipleSample_merge, fast=TRUE)
  colnames(MultipleSample_merge@SpatialPCs) = rownames(MultipleSample_merge@location)
  
  id_list = list()
  SpatialPC_list = list()
  Location_spatialpc_list = list()
  
  # print("Column names of MultipleSample_merge@SpatialPCs:")
  # print(colnames(MultipleSample_merge@SpatialPCs))
  print("Dimensions of MultipleSample_merge@SpatialPCs:")
  print(dim(MultipleSample_merge@SpatialPCs))
  
  for(sample in 1:length(count_list)){
    id_list[[sample]] = grep(paste0("Sample",sample), colnames(MultipleSample_merge@SpatialPCs))
    SpatialPC_list[[sample]] = MultipleSample_merge@SpatialPCs[,id_list[[sample]]]
    Location_spatialpc_list[[sample]] = spatialpca_list[[sample]]@location
    
    print(paste("Sample", sample))
    # print("Matching column indices:")
    # print(id_list[[sample]])
  }
  
  names(SpatialPC_list) = paste0("Sample",1:length(SpatialPC_list))
  names(Location_spatialpc_list) = paste0("Sample",1:length(Location_spatialpc_list))
  return(list("MultipleSample_SpatialPCA_object"=MultipleSample_merge, 
              "SpatialPC_list" = SpatialPC_list, 
              "Location_spatialpc_list"=Location_spatialpc_list))
}

# maxiter=300;initial_tau=1;fast=FALSE;eigenvecnum=NULL;SpatialPCnum=20
# object = MultipleSample_merge

# maxiter=3000 is changed by YW
SpatialPCA_EstimateLoading = function(object, maxiter=300,
                                      initial_tau=1,
                                      # initial_tau=2, #changed by YW
                                      fast=FALSE,eigenvecnum=NULL,SpatialPCnum=20){
  

  
  suppressMessages(require(RSpectra))
  set.seed(1234)
  param_ini=log(initial_tau)
  object@SpatialPCnum = SpatialPCnum
  object@fast = fast
  object@params$X = scale(object@location)
  object@params$n = dim(object@params$X)[1]
  object@params$p=dim(object@params$X)[2]
  
  if(is.null(object@covariate)){
    object@params$H = matrix(1, dim(object@params$X)[1],1)
    HH_inv=solve(t(object@params$H)%*%object@params$H,tol = 1e-40)
    HH = object@params$H%*%HH_inv%*%t(object@params$H)
    object@params$M=diag(object@params$n)-HH #####!!!!!!!!!!!!!!!!!!!!!
  

    ### inserted by YW for debugging
    print("Starting SpatialPCA_EstimateLoading")
    print(paste("n =", object@params$n))
    print(paste("p =", object@params$p))
    
    
    print("Dimensions of object@params$expr:")
    print(dim(object@params$expr))
    print("Dimensions of object@params$M:")
    print(dim(object@params$M))
    
    # object@params$tr_YMY=sum(diag(object@params$expr%*%object@params$M%*%t(object@params$expr)))
    #######
    ####
    object@params$tr_YMY=sum(diag(object@params$expr%*%object@params$M%*%t(object@params$expr)))
    object@params$YM = object@params$expr%*%object@params$M
    object@params$q=1
  }else{
    object@params$q = dim(object@covariate)[2]+1
    object@params$H = matrix(0, object@params$n,object@params$q)
    object@params$H[,1]=1
    object@params$H[,2:object@params$q] = object@covariate
    HH_inv=solve(t(object@params$H)%*%object@params$H,tol = 1e-40)
    HH=object@params$H%*%HH_inv%*%t(object@params$H)
    
    ### added by YW for debugging
    print("Creating H matrix")
    print(dim(object@params$H))
    
    print("Calculating HH")
    print(dim(HH))
    
    print(paste("n =", object@params$n))
    print(paste("p =", object@params$p))
    print(paste("Is n NULL?", is.null(object@params$n)))
    print(paste("Is p NULL?", is.null(object@params$p)))
    
    print("Checking H matrix:")
    print(paste("Dimensions of H:", paste(dim(object@params$H), collapse="x")))
    print(paste("Is H NULL?", is.null(object@params$H)))
    
    print("Creating M matrix:")
    if(is.null(object@params$H)) {
      print("H is NULL, using default")
      object@params$H = matrix(1, object@params$n, 1)
    }
    HH_inv = try(solve(t(object@params$H) %*% object@params$H, tol = 1e-40))
    print(paste("Is HH_inv an error?", inherits(HH_inv, "try-error")))
    if(!inherits(HH_inv, "try-error")) {
      HH = object@params$H %*% HH_inv %*% t(object@params$H)
      object@params$M = diag(object@params$n) - HH
      print(paste("Dimensions of M:", paste(dim(object@params$M), collapse="x")))
    } else {
      print("Error in creating HH_inv")
    }
    print("Calculating M")
    object@params$M=diag(object@params$n)-HH #####!!!!!!!!!!!!!!!!!!!!!
    print(dim(object@params$M))
    ### 
    object@params$tr_YMY=sum(diag(object@params$expr%*%object@params$M%*%t(object@params$expr)))
    object@params$YM = object@params$expr%*%object@params$M
  }
  
  
  if(fast==FALSE){
    object@fast=fast
    print("Eigen decomposition on kernel matrix!")
    eigen_res = eigen(object@kernelmat)
    object@params$delta = eigen_res$values
    object@params$U = eigen_res$vectors
    print("Using all eigenvectors and eigenvalues in the Kernel matrix!")
  }else{
    object@fast=fast
    if(!is.null(eigenvecnum)){
      print("Eigen decomposition on kernel matrix!")
      object@eigenvecnum=eigenvecnum
      if(object@sparseKernel==TRUE){
        eigen_res = eigs_sym(object@kernelmat, k=object@eigenvecnum)
        object@params$delta = eigen_res$values
        object@params$U = eigen_res$vectors
      }else{
        eigen_res = eigs_sym(object@kernelmat, k=object@eigenvecnum, which = "LM")
        object@params$delta = eigen_res$values
        object@params$U = eigen_res$vectors
      }
      
      print("Low rank approximation!")
      print(paste0("Using user selected top ",object@eigenvecnum," eigenvectors and eigenvalues in the Kernel matrix!"))
    }else if(object@params$n>5000){
      print("Eigen decomposition on kernel matrix!")
      if(object@sparseKernel==TRUE){
        eigen_res = eigs_sym(object@kernelmat, k=20)
        object@params$delta = eigen_res$values
        object@params$U = eigen_res$vectors
      }else{
        eigen_res = eigs_sym(object@kernelmat, k=20, which = "LM")
        object@params$delta = eigen_res$values
        object@params$U = eigen_res$vectors
      }
      print("Low rank approximation!")
      print("Large sample, using top 20 eigenvectors and eigenvalues in the Kernel matrix!")
    }else{
      eigen_res = eigen(object@kernelmat)
      delta_all = eigen_res$values
      U_all = eigen_res$vectors
      ind = which(cumsum(delta_all/length(delta_all))>0.9)[1]
      print("Low rank approximation!")
      print(paste0("Small sample, using top ",ind," eigenvectors and eigenvalues in the Kernel matrix!"))
      object@params$delta = delta_all[1:ind]
      object@params$U = U_all[,1:ind]
      rm(U_all)
    }
  }
  
  
  object@params$MYt = object@params$M %*% t(object@params$expr)
  object@params$YMMYt = object@params$YM %*% object@params$MYt
  object@params$YMU = object@params$YM %*% object@params$U
  object@params$Xt = t(object@params$H)
  object@params$XtU = object@params$Xt %*% object@params$U
  object@params$Ut = t(object@params$U)
  object@params$UtX = object@params$Ut %*% object@params$H
  object@params$YMX = object@params$YM %*% object@params$H
  object@params$UtU = object@params$Ut %*% object@params$U
  object@params$XtX = object@params$Xt %*% object@params$H
  object@params$SpatialPCnum = SpatialPCnum
  
  
  # optim_result =try(optim(param_ini, SpatialPCA_estimate_parameter,params=object@params,
  #                         control = list(maxit = maxiter), lower = -10, upper = 10,method="Brent"),silent=T) #chaged by YW
  optimize_with_fallback <- function(param_ini, objective_func, params, maxiter = 1000, lower = -10, upper = 10) {
    tryCatch({
      # Attempt optimization
      result <- optim(param_ini, objective_func, params = params,
                      control = list(maxit = maxiter), 
                      lower = lower, upper = upper, method = "Brent")
      
      # Check if optimization converged
      if (result$convergence == 0) {
        return(list(par = result$par, success = TRUE))
      } else {
        warning("Optimization did not converge. Returning initial parameters.")
        return(list(par = param_ini, success = FALSE))
      }
    }, error = function(e) {
      warning(paste("Optimization failed with error:", e$message, "\nReturning initial parameters."))
      return(list(par = param_ini, success = FALSE))
    })
  }
  
  # Usage:
  optim_result <- optimize_with_fallback(param_ini, SpatialPCA_estimate_parameter, params = object@params)
  
  # Check the result
  if (optim_result$success) {
    # cat("Optimization succeeded. Optimal parameters:", optim_result$par, "\n")
    cat("Optimization succeeded. Optimal parameters: \n")
    
  } else {
    cat("Optimization failed. Using initial parameters: \n")
  }
  
  ######################
  ######################
  
  # Use the parameters (either optimized or initial)
  final_params <- optim_result$par
  
  object@tau = exp(optim_result$par)
  k = dim(object@params$expr)[1]
  n = dim(object@params$expr)[2]
  q=object@params$q
  tauD_UtU_inv = solve(object@tau*diag(object@params$delta) + object@params$UtU, tol = 1e-40)
  YMU_tauD_UtU_inv_Ut = object@params$YMU %*% tauD_UtU_inv %*% object@params$Ut
  YMU_tauD_UtU_inv_UtX = YMU_tauD_UtU_inv_Ut %*% object@params$H
  XtU_inv_UtX = object@params$XtU %*% tauD_UtU_inv %*% object@params$UtX
  left = object@params$YMX - YMU_tauD_UtU_inv_UtX
  right = t(left)
  middle = solve(-XtU_inv_UtX, tol = 1e-40)
  G_each = object@params$YMMYt - YMU_tauD_UtU_inv_Ut %*% object@params$MYt - left %*% middle %*% right
  object@W = eigs_sym(G_each, k=SpatialPCnum, which = "LM")$vectors
  object@sigma2_0 = as.numeric((object@params$tr_YMY+F_funct_sameG(object@W,G_each))/(k*(n-q)))
  
  rm(eigen_res)
  rm(tauD_UtU_inv)
  rm(YMU_tauD_UtU_inv_Ut)
  rm(YMU_tauD_UtU_inv_UtX)
  rm(XtU_inv_UtX)
  rm(left)
  rm(right)
  rm(middle)
  rm(G_each)
  gc()
  
  return(object)
}

spatialPCA_estimate_parameter = function(param_ini, params){
  # suppressMessages(require(RSpectra))
  set.seed(1234)
  tau=exp(param_ini[1])
  k = dim(params$expr)[1]
  n = dim(params$expr)[2]
  q=params$q
  PCnum=params$SpatialPCnum
  tauD_UtU_inv = solve(tau*diag(params$delta) + params$UtU, tol = 1e-40)
  YMU_tauD_UtU_inv_Ut = params$YMU %*% tauD_UtU_inv %*% params$Ut
  YMU_tauD_UtU_inv_UtX = YMU_tauD_UtU_inv_Ut %*% params$H
  XtU_inv_UtX = params$XtU %*% tauD_UtU_inv %*% params$UtX
  left = params$YMX - YMU_tauD_UtU_inv_UtX
  right = t(left)
  middle = solve(-XtU_inv_UtX, tol = 1e-40)
  G_each = params$YMMYt - YMU_tauD_UtU_inv_Ut %*% params$MYt - left %*% middle %*% right
  log_det_tauK_I = determinant(1/tau*diag(1/params$delta)+ params$UtU, logarithm=TRUE)$modulus[1] + determinant(tau*diag(params$delta), logarithm=TRUE)$modulus[1]
  Xt_invmiddle_X = params$XtX - params$XtU %*% solve(params$UtU + 1/tau *diag( 1/params$delta) , tol = 1e-40) %*% params$UtX
  log_det_Xt_inv_X = determinant(Xt_invmiddle_X, logarithm=TRUE)$modulus[1]
  sum_det=0
  sum_det=sum_det+(0.5*log_det_tauK_I+0.5*log_det_Xt_inv_X  )*PCnum
  
  rm(tauD_UtU_inv)
  rm(YMU_tauD_UtU_inv_Ut)
  rm(YMU_tauD_UtU_inv_UtX)
  rm(XtU_inv_UtX)
  rm(left)
  rm(middle)
  rm(right)
  rm(Xt_invmiddle_X)
  gc()
  
  W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors
  -(-sum_det -(k*(n-q))/2*log(params$tr_YMY+F_funct_sameG(W_est_here,G_each)))
}

F_funct_sameG = function(X,G){ # G is a matrix
  return_val=0
  for(i in 1: dim(X)[2]){
    return_val=return_val+t(X[,i])%*%G%*%X[,i]
  }
  -return_val
}