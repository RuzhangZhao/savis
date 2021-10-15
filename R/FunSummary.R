#' savis
#'
#' savis: single-cell RNAseq adaptive visualiztaion
#'
#' @details This function argument to the function
#'
#' @param expr_matrix The expression matrix: gene(feature) as row; cell(sample) as column.
#' @param assay_for_var_features 
#' @param is_count_matrix Whether expr_matrix is count matrix or normalized version. If the expression matrix is count matrix, normalization will be performed. Default is TRUE.
#' @param npcs The number of principle components will be computed. Default is 20.
#' @param nfeatures The number of highly variable genes will be selected. Default is 2000.
#' @param distance_metric The default is "euclidean". 
#' @param cluster_method The default is "louvain". User can choose from c("louvain","spectral"). But "louvain" performs much better.
#' @param resolution The resolution for The default is 0.5.
#' @param resolution_sub The default is 0.
#' @param memory_save The default is FALSE. This function will take some storage to temporarily save the data from memory. Don't worry. SAVIS will soon delete it!!! Also, SAVIS uses unique name for the temporary data to keep everything safe!!!
#' @param adaptive The default is FALSE.
#' @param max_stratification The default is 3.
#' @param scale_factor_separation The default is 3.
#' @param process_min_size The default is NULL.
#' @param process_min_count The default is NULL.
#' @param run_adaUMAP The default is TRUE.
#' @param adjust_UMAP The default is TRUE.
#' @param adjust_method The default is "all". Select from c("umap","mds").
#' @param adjust_rotate The default is TRUE.
#' @param shrink_distance The default is TRUE.
#' @param check_differential The default is FALSE.
#' @param verbose The default is TRUE.
#' @param show_cluster The default is FALSE.
#' @param return_cluster The default is FALSE.
#' @param verbose_more The default is FALSE.
#' @param seed.use The default is 42L
#' 
#' 
#' 
#' 
#' 
#' 
#'
#' @return nothing useful
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA DefaultAssay
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom uwot umap
#' @importFrom MASS isoMDS
#' @importFrom cluster pam
#' @import ggplot2
#' @import RColorBrewer
#' @import dplyr
#' 
#' @export
#'
#' @examples
#' a<-1
#'
savis<-function(
  expr_matrix,
  is_count_matrix=TRUE,
  assay_for_var_features = "rawcount",
  npcs = 20,
  nfeatures = 2000,
  distance_metric = "euclidean",
  cluster_method = "louvain",
  resolution = 0.5,
  resolution_sub = 0,
  memory_save = FALSE,
  adaptive = TRUE,
  max_stratification = 3,
  scale_factor_separation =3,
  process_min_size = NULL,
  process_min_count = NULL,
  run_adaUMAP = TRUE,
  adjust_UMAP = TRUE,
  adjust_method = "all",
  adjust_rotate = TRUE,
  shrink_distance = TRUE,
  density_adjust = TRUE,
  density_adjust_via_global_umap = FALSE,
  global_umap_embedding = NULL,
  check_differential = FALSE,
  verbose = TRUE,
  show_cluster = FALSE,
  return_cluster = FALSE,
  return_combined_PC = FALSE,
  verbose_more = FALSE,
  compressed_storage = TRUE,
  seed.use = 42L
){
  if (!inherits(x = expr_matrix, 'Matrix')) {
    expr_matrix <- as(object = as.matrix(x = expr_matrix), Class = 'Matrix')
  }
  if (!inherits(x = expr_matrix, what = 'dgCMatrix')) {
    expr_matrix <- as(object = expr_matrix, Class = 'dgCMatrix')
  }
  # change the seed.use to be integer
  if(!is.integer(seed.use)){
    seed.use<-as.integer(seed.use)
  }
  if(max_stratification == 1){
    stop("Please directly use umap: savis 
      supports adaptive visualization for 
      max_stratification larger than 1.")
  }
  
  if (nrow(expr_matrix) < nfeatures){
    stop("nfeatures should be smaller than
      the number of features in expression
      matrix")
  }
  if (ncol(expr_matrix) < npcs){
    stop("npcs(number of PC) should be smaller
      than the number of samples in expression
      matrix")
  }
  if(is.null(rownames(expr_matrix)[1])){
    rownames(expr_matrix)<-c(1:nrow(expr_matrix))
  }
  if(is.null(colnames(expr_matrix)[1])){
    colnames(expr_matrix)<-c(1:ncol(expr_matrix))
  }else if(length(unique(colnames(expr_matrix)))<
      ncol(expr_matrix) ) {
    print("WARN: There are duplicated cell names! Make cell names unique by renaming!")
    colnames(expr_matrix)<-make.unique(colnames(expr_matrix))
  }else if(length(unique(rownames(expr_matrix)))<
      nrow(expr_matrix) ) {
    print("WARN: There are duplicated gene names! Make gene names unique by renaming!")
    rownames(expr_matrix)<-make.unique(rownames(expr_matrix))
  }
  if (!(assay_for_var_features %in% c("rawcount","normalizedcount"))){
    stop("Please select assay_for_var_features from c('rawcount','normalizedcount')")
  }
  if(!is_count_matrix & assay_for_var_features == "rawcount"){
    cat("assay_for_var_features is set to be 'rawcount'. \n Please use count matrix as input. \n Also, please set is_count_matrix to be TRUE. \n")
    stop()
  }
  if(is_count_matrix & assay_for_var_features == "normalizedcount"){
    if(verbose){
      cat('\n')
      print("Normalizing Expression Matrix...")
      pb <- txtProgressBar(min = 0, max = 20, style = 3, file = stderr())
    }
    
    expr_matrix<-NormalizeData(
      expr_matrix,
      verbose = verbose_more)
    if(verbose){
      cat('\n')
      print("Finding Variable Features...")
      setTxtProgressBar(pb = pb, value = 1)
    }
   
    expr_matrix_hvg <- ExpFindVariableFeatures(
      expr_matrix,
      verbose = verbose_more)
  }else if (is_count_matrix & assay_for_var_features == "rawcount"){
    if(verbose){
      cat('\n')
      print("Normalizing Expression Matrix...")
      pb <- txtProgressBar(min = 0, max = 20, style = 3, file = stderr())
    }
    if(verbose){
      cat('\n')
      print("Finding Variable Features...")
      setTxtProgressBar(pb = pb, value = 1)
    }
    expr_matrix_hvg <- FindVariableFeatures(
      expr_matrix,
      verbose = verbose_more)$vst.variance.standardized
    expr_matrix_process<-NormalizeData(
      expr_matrix,
      verbose = verbose_more)
  }else{
    if(verbose){
      pb <- txtProgressBar(min = 0, max = 20, style = 3, file = stderr())
    }
    if(verbose){
      cat('\n')
      print("Finding Variable Features...")
      setTxtProgressBar(pb = pb, value = 1)
    }
    expr_matrix_hvg <- ExpFindVariableFeatures(
      expr_matrix,
      verbose = verbose_more)
  }
  
  hvg<-savis_nth(x = expr_matrix_hvg,
    k = nfeatures)
  if(assay_for_var_features == "rawcount"){
    expr_matrix_process<-expr_matrix_process[hvg,]
  }else{
    expr_matrix_process<-expr_matrix[hvg,]
  }
 
  if(verbose){
    cat('\n')
    print("Scaling Expression Matrix...")
    setTxtProgressBar(pb = pb, value = 2)
  }
  
  expr_matrix_process <- ScaleData(
    expr_matrix_process,
    verbose = verbose_more)
  
  if(verbose){
    cat('\n')
    print("Calculating Global PCA...")
    setTxtProgressBar(pb = pb, value = 3)
  }
  suppressWarnings(expr_matrix_pca <- RunPCA(
    object = expr_matrix_process,
    features = rownames(expr_matrix_process),
    npcs = npcs,
    verbose = verbose_more)@cell.embeddings)
  rm(expr_matrix_process)
  expr_matrix_pca<-data.frame(expr_matrix_pca)
  expr_matrix_pca<-as.matrix(expr_matrix_pca)
  if(verbose){
    cat('\n')
    print("Doing Clustering...")
    setTxtProgressBar(pb = pb, value = 5)
  }
  cluster_label<-DoCluster(
    pc_embedding = expr_matrix_pca,
    method = cluster_method,
    resolution = resolution)$cluster
  global_cluster_label<-cluster_label
  # sorted unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  
  size_cluster<-c()
  for ( i in 1:N_label){
    size_cluster<-c(size_cluster,
      sum(cluster_label == label_index[i]))
  }
  if(verbose){
    cat('\n')
    print(paste0("Clustering Results:",
      length(unique(cluster_label)),
      " clusters."))
    setTxtProgressBar(pb = pb, value = 6)
  }
  if(verbose){
    if(show_cluster){
      cat('\n')
      print(paste0("Size of Cluster: ",size_cluster))
    }
  }
  
  if(verbose){
    cat('\n')
    print("Calculating Local PCA...")
    setTxtProgressBar(pb = pb, value = 8)
  }

  if(max_stratification == 2 | adaptive == FALSE){
    
    combined_embedding<-FormCombinedEmbedding(
      expr_matrix=expr_matrix,
      expr_matrix_pca=expr_matrix_pca,
      cluster_label=cluster_label,
      npcs=npcs,
      nfeatures =nfeatures,
      assay_for_var_features = assay_for_var_features,
      #center_method = center_method,
      scale_factor_separation=scale_factor_separation)
    rm(expr_matrix)
    adaptive<-FALSE
  }
  if(adaptive){
    if (!is.null(process_min_count)){
      process_min_size<-sort(size_cluster,decreasing = T)[process_min_count]
    }else{
      if (is.null(process_min_size)){
        process_min_size<-mean(size_cluster)
      }
    }
    if(verbose){
      cat('\n')
      print(paste0("Process_min_size: ",process_min_size))
      setTxtProgressBar(pb = pb, value = 9.5)
    }
    
    
    if(verbose){
      cat('\n')
      print("Exploring if clusters can be separated further...")
      setTxtProgressBar(pb = pb, value = 10)
    }
    
    umap_res<-FormAdaptiveCombineList(
      expr_matrix=expr_matrix,
      expr_matrix_pca=expr_matrix_pca,
      max_stratification=max_stratification,
      stratification_count=1,
      assay_for_var_features = assay_for_var_features,
      scale_factor_separation = scale_factor_separation,
      resolution=resolution_sub,
      cluster_method=cluster_method,
      npcs=npcs,
      nfeatures=nfeatures,
      process_min_size=process_min_size,
      do_cluster = FALSE,
      cluster_label = cluster_label,
      check_differential = check_differential,
      verbose = verbose_more)
    rm(expr_matrix)
    cluster_label<-umap_res$cluster_label
    if(is.null(dim(cluster_label)[1])){
      combined_embedding<-data.frame(
        "Layer2Cluster"=cluster_label,
        umap_res$combined_embedding)
      if(ncol(umap_res$combined_embedding)!=(2*npcs)){
        stop("label and combined embedding size do not match")
      }
      metric_count<-1
    }else if (ncol(cluster_label) == 2) {
      colnames(cluster_label)<-paste0("Layer",
        2:(ncol(cluster_label)+1),"Cluster")
      combined_embedding<-data.frame(
        cluster_label,
        umap_res$combined_embedding)
      if(ncol(umap_res$combined_embedding)!=(3*npcs)){
        stop("label and combined embedding size do not match")
      }
      metric_count<-2
    }else if (ncol(cluster_label) == 3){
      colnames(cluster_label)<-paste0("Layer",
        2:(ncol(cluster_label)+1),"Cluster")
      combined_embedding<-data.frame(
        cluster_label,
        umap_res$combined_embedding)
      if(ncol(umap_res$combined_embedding)!=(4*npcs)){
        stop("label and combined embedding size do not match")
      }
      metric_count<-3
    }else{
      colnames(cluster_label)<-paste0("Layer",
        2:(ncol(cluster_label)+1),"Cluster")
      combined_embedding<-data.frame("Num_Layer"=(ncol(cluster_label)+1),
        cluster_label,
        umap_res$combined_embedding)
      metric_count<-4
    }
    rm(umap_res)
  }else{
    combined_embedding<-data.frame(
      "cluster_label"=cluster_label,
      combined_embedding)
    metric_count <- 1
  }
  
  
  combined_embedding_list<-list()
  if(compressed_storage){
    if(metric_count > 1){
      #combined_embedding_list[[1]]<-cluster_label
      ncol_cluster<-ncol(cluster_label)
      start_col<-ncol_cluster+1
      for ( i in 2:ncol_cluster){
        if(sum(cluster_label[,i]== -1) > 0){
          start_col<-i
          break
        }
      }
      if(start_col <= ncol_cluster){
        combined_embedding_list[[1]]<-combined_embedding[,1:(ncol_cluster+start_col*npcs)]
        for ( i in start_col:ncol_cluster){
          index_i<-which(cluster_label[,i]!= -1)
          combined_embedding_list[[i]]<-combined_embedding[index_i,(ncol_cluster+i*npcs+1):(ncol_cluster+(i+1)*npcs)]
        } 
      }
    }
  }
  
  if(!run_adaUMAP){
    if(return_cluster){
      if(length(combined_embedding_list)>0){
        newList<-list("combined_embedding"=combined_embedding_list,
          "cluster_label"=cluster_label) 
      }else{
        newList<-list("combined_embedding"=combined_embedding,
          "cluster_label"=cluster_label)
      }
    }else{
      if(length(combined_embedding_list)>0){
        newList<-combined_embedding_list
      }else{
        newList<-combined_embedding
      }
    }
  }else{
    if(verbose){
      cat('\n')
      print("Running Adaptive UMAP...")
      setTxtProgressBar(pb = pb, value = 12)
    }
    #print(metric_count)
    umap_embedding<-RunAdaUMAP(
      X = combined_embedding,
      metric = distance_metric,
      metric_count = metric_count,
      seed.use = seed.use)
    if(adjust_UMAP){
      if(verbose){
        cat('\n')
        print("Adjusting UMAP...")
        setTxtProgressBar(pb = pb, value = 17)
      }
      expr_matrix_umap = NULL
      if(density_adjust_via_global_umap){
        if(is.null(global_umap_embedding)){
          expr_matrix_umap<-umap(
            X = expr_matrix_pca,
            a = 1.8956, 
            b = 0.8006, 
            metric = distance_metric
          )
        }else{
          expr_matrix_umap<-global_umap_embedding
        }
      }
      
      umap_embedding<<-umap_embedding
      expr_matrix_pca<<-expr_matrix_pca
      global_cluster_label<<-global_cluster_label
      
      umap_embedding<-adjustUMAP(
        pca_embedding = expr_matrix_pca,
        umap_embedding = umap_embedding,
        cluster_label = global_cluster_label,
        global_umap_embedding = expr_matrix_umap,
        adjust_method = adjust_method,
        density_adjust = density_adjust,
        shrink_distance = shrink_distance,
        rotate = adjust_rotate,
        seed.use = seed.use)
    }
    if(return_cluster){
      if(length(combined_embedding_list)>0){
        newList<-list("combined_embedding"=combined_embedding_list,
          "savis_embedding"=umap_embedding,
          "cluster_label"=cluster_label) 
      }else{
        newList<-list("combined_embedding"=combined_embedding,
          "savis_embedding"=umap_embedding,
          "cluster_label"=cluster_label)
      }
      #newList<-list("savis_embedding"=umap_embedding,
       # "cluster_label"=cluster_label)
    }else{
      newList<-umap_embedding
    }
  }
  
  
  
  if(verbose){
    cat('\n')
    print("Finished...")
    setTxtProgressBar(pb = pb, value = 20)
  }
  return(newList)
}


#' RunPreSAVIS
#'
#' savis: single-cell RNAseq adaptive visualiztaion
#'
#' @details This function argument to the function
#'
#' @param object sdsd
#' @param assay_for_var_features sds
#' @param distance_metric The default is "euclidean". 
#' @param cluster_method The default is "louvain". User can choose from c("louvain","spectral"). But "louvain" performs much better.
#' @param resolution The resolution for The default is 0.5.
#' @param resolution_sub The default is 0.
#' @param max_stratification The default is 3.
#' @param scale_factor_separation The default is 3.
#' @param process_min_size The default is NULL.
#' @param process_min_count The default is NULL.
#' @param check_differential The default is FALSE.
#' @param verbose The default is TRUE.
#' @param show_cluster The default is FALSE.
#' @param return_cluster The default is FALSE.
#' @param verbose_more The default is FALSE.
#' 
#' 
#' 
#' 
#' 
#' 
#'
#' @return nothing useful
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures VariableFeatures ScaleData RunPCA DefaultAssay
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom methods new representation setClass as
#' @export
#'
#' @examples
#' a<-1
#'
RunPreSAVIS<-function(
  object,
  assay_for_var_features = "rawcount",
  distance_metric = "euclidean",
  cluster_method = "louvain",
  resolution = 0.5,
  resolution_sub = 0,
  adaptive = TRUE,
  max_stratification = 3,
  scale_factor_separation =3,
  process_min_size = NULL,
  process_min_count = NULL,
  check_differential = FALSE,
  verbose = TRUE,
  show_cluster = FALSE,
  return_cluster = FALSE,
  verbose_more = FALSE
){
  default_assay<-DefaultAssay(object)
  print(paste0("PreSAVIS is based on the default assay: ",default_assay))
  if(verbose){
    pb <- txtProgressBar(min = 0, max = 10, style = 3, file = stderr())
  }
  nfeatures<-length(VariableFeatures(object))
  npcs<-ncol(object@reductions$pca@cell.embeddings)
  
  if(max_stratification == 1){
    stop("Please directly use umap: savis 
      supports adaptive visualization for 
      max_stratification larger than 1.")
  }
  
  if(verbose){
    cat('\n')
    print("Global PCs are captured from SeuratObject...")
    setTxtProgressBar(pb = pb, value = 1)
  }
  expr_matrix_pca <- object@reductions$pca@cell.embeddings
  expr_matrix_pca<-as.matrix(expr_matrix_pca)
  if(verbose){
    cat('\n')
    print("Doing Clustering...")
    setTxtProgressBar(pb = pb, value = 2)
  }
  cluster_label<-DoCluster(
    pc_embedding = expr_matrix_pca,
    method = cluster_method,
    resolution = resolution)$cluster
  global_cluster_label<-cluster_label
  # sorted unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  
  size_cluster<-c()
  for ( i in 1:N_label){
    size_cluster<-c(size_cluster,
      sum(cluster_label == label_index[i]))
  }
  if(verbose){
    cat('\n')
    print(paste0("Clustering Results:",
      length(unique(cluster_label)),
      " clusters."))
    setTxtProgressBar(pb = pb, value = 3)
  }
  if(verbose){
    if(show_cluster){
      cat('\n')
      print(paste0("Size of Cluster: ",size_cluster))
    }
  }
  
  if(verbose){
    cat('\n')
    print("Calculating Local PCA...")
    setTxtProgressBar(pb = pb, value = 4)
  }
  if(max_stratification == 2 | adaptive == FALSE){
    
    if(assay_for_var_features == "rawcount"){
      combined_embedding<-FormCombinedEmbedding(
        expr_matrix=object@assays[[default_assay]]@counts,
        expr_matrix_pca=expr_matrix_pca,
        cluster_label=cluster_label,
        npcs=npcs,
        nfeatures =nfeatures,
        assay_for_var_features = "rawcount",
        #center_method = center_method,
        scale_factor_separation=scale_factor_separation)
    }else if(assay_for_var_features == "normalizedcount"){
      combined_embedding<-FormCombinedEmbedding(
        expr_matrix=object@assays[[default_assay]]@data,
        expr_matrix_pca=expr_matrix_pca,
        cluster_label=cluster_label,
        npcs=npcs,
        nfeatures =nfeatures,
        assay_for_var_features = "normalizedcount",
        #center_method = center_method,
        scale_factor_separation=scale_factor_separation)
    }else{
      stop("Please select assay_for_var_features from c('rawcount','normalizedcount')")
    }
    
    adaptive<-FALSE
  }
  if(adaptive){
    if (!is.null(process_min_count)){
      process_min_size<-sort(size_cluster,decreasing = T)[process_min_count]
    }else{
      if (is.null(process_min_size)){
        process_min_size<-mean(size_cluster)
      }
    }
    if(verbose){
      cat('\n')
      print(paste0("Process_min_size: ",process_min_size))
      setTxtProgressBar(pb = pb, value = 5)
    }
    
    
    if(verbose){
      cat('\n')
      print("Exploring if clusters can be separated further...")
      setTxtProgressBar(pb = pb, value = 7)
    }
    
    if(assay_for_var_features == "rawcount"){
      umap_res<-FormAdaptiveCombineList(
        expr_matrix=object@assays[[default_assay]]@counts,
        expr_matrix_pca=expr_matrix_pca,
        max_stratification=max_stratification,
        stratification_count=1,
        scale_factor_separation = scale_factor_separation,
        resolution=resolution_sub,
        cluster_method=cluster_method,
        npcs=npcs,
        nfeatures=nfeatures,
        process_min_size=process_min_size,
        assay_for_var_features = "rawcount",
        do_cluster = FALSE,
        cluster_label = cluster_label,
        check_differential = check_differential,
        verbose = verbose_more)
    }else if(assay_for_var_features == "normalizedcount"){
      umap_res<-FormAdaptiveCombineList(
        expr_matrix=object@assays[[default_assay]]@data,
        expr_matrix_pca=expr_matrix_pca,
        max_stratification=max_stratification,
        stratification_count=1,
        scale_factor_separation = scale_factor_separation,
        resolution=resolution_sub,
        cluster_method=cluster_method,
        npcs=npcs,
        nfeatures=nfeatures,
        process_min_size=process_min_size,
        assay_for_var_features = "normalizedcount",
        do_cluster = FALSE,
        cluster_label = cluster_label,
        check_differential = check_differential,
        verbose = verbose_more)
    }else{
      stop("Please select assay_for_var_features from c('rawcount','normalizedcount')")
    }
    cluster_label<-umap_res$cluster_label
    if(is.null(dim(cluster_label)[1])){
      combined_embedding<-data.frame(
        "Layer2Cluster"=cluster_label,
        umap_res$combined_embedding)
      if(ncol(umap_res$combined_embedding)!=(2*npcs)){
        stop("label and combined embedding size do not match")
      }
      metric_count<-1
    }else if (ncol(cluster_label) == 2) {
      colnames(cluster_label)<-paste0("Layer",
        2:(ncol(cluster_label)+1),"Cluster")
      combined_embedding<-data.frame(
        cluster_label,
        umap_res$combined_embedding)
      if(ncol(umap_res$combined_embedding)!=(3*npcs)){
        stop("label and combined embedding size do not match")
      }
      metric_count<-2
    }else if (ncol(cluster_label) == 3){
      colnames(cluster_label)<-paste0("Layer",
        2:(ncol(cluster_label)+1),"Cluster")
      combined_embedding<-data.frame(
        cluster_label,
        umap_res$combined_embedding)
      if(ncol(umap_res$combined_embedding)!=(4*npcs)){
        stop("label and combined embedding size do not match")
      }
      metric_count<-3
    }else{
      colnames(cluster_label)<-paste0("Layer",
        2:(ncol(cluster_label)+1),"Cluster")
      combined_embedding<-data.frame("Num_Layer"=(ncol(cluster_label)+1),
        cluster_label,
        umap_res$combined_embedding)
      metric_count<-4
    }
    rm(umap_res)
  }else{
    combined_embedding<-data.frame(
      "cluster_label"=cluster_label,
      combined_embedding)
    metric_count <- 1
  }
  
  setClass("savis_class", representation(
    combined_embedding = "matrix", 
    cluster_label = "numeric",
    global_cluster_label = "numeric",
    savis_embedding = "matrix",
    distance_metric = "character",
    metric_count = "numeric"))
  savis_pre <- new("savis_class", 
    combined_embedding = as.matrix(combined_embedding), 
    cluster_label = cluster_label,
    global_cluster_label = global_cluster_label,
    distance_metric = distance_metric,
    metric_count = metric_count)
  
  object@reductions$savis<-savis_pre
  
  if(verbose){
    cat('\n')
    print("Finished...")
    setTxtProgressBar(pb = pb, value = 10)
  }
  return(object)
}






#' RunSAVIS
#'
#' RunSAVIS
#'
#' @details This function argument to the function
#'
#' @param object sds
#' @param adjust_UMAP = TRUE,
#' @param adjust_method = "all",
#' @param adjust_rotate = TRUE,
#' @param shrink_distance = TRUE,
#' @param density_adjust = TRUE,
#' @param verbose = TRUE,
#' @param seed.use = 42L
#' 
#' 
#' 
#' 
#' 
#'
#' @return nothing useful
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA DefaultAssay
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' 
#' @examples
#' a<-1
#' 
RunSAVIS<-function(
  object,
  adjust_UMAP = TRUE,
  adjust_method = "all",
  adjust_rotate = TRUE,
  shrink_distance = TRUE,
  density_adjust = TRUE,
  verbose = TRUE,
  seed.use = 42L
){
  # change the seed.use to be integer
  if(!is.integer(seed.use)){
    seed.use<-as.integer(seed.use)
  }
  default_assay<-DefaultAssay(object)
  if(is.null(object@reductions$savis)){
    stop("Please apply RunPreSAVIS before RunSAVIS")
  }
  if(verbose){
    pb <- txtProgressBar(min = 0, max = 10, style = 3, file = stderr())
  }
  if(!is.null(object@reductions$umap)){
    density_adjust_via_global_umap<-TRUE
  }else{
    density_adjust_via_global_umap<-FALSE
  }
  if(verbose){
    cat('\n')
    print("Running SAVIS with Adaptive Settings...")
    setTxtProgressBar(pb = pb, value = 1)
  }
  #print(metric_count)
  savis_embedding<-RunAdaUMAP(
    X = object@reductions$savis@combined_embedding,
    metric = object@reductions$savis@distance_metric,
    metric_count = object@reductions$savis@metric_count,
    seed.use = seed.use)
  if(adjust_UMAP){
    if(verbose){
      cat('\n')
      print("Adjusting SAVIS...")
      setTxtProgressBar(pb = pb, value = 8)
    }
    expr_matrix_umap = NULL
    if(density_adjust_via_global_umap){
      expr_matrix_umap<-object@reductions$umap@cell.embeddings
    }
    expr_matrix_pca<-object@reductions$pca@cell.embeddings
    savis_embedding<-adjustUMAP(
      pca_embedding = expr_matrix_pca,
      umap_embedding = savis_embedding,
      cluster_label = object@reductions$savis@global_cluster_label,
      global_umap_embedding = expr_matrix_umap,
      adjust_method = adjust_method,
      density_adjust = density_adjust,
      shrink_distance = shrink_distance,
      rotate = adjust_rotate,
      seed.use = seed.use)
  }
  object@reductions$savis@savis_embedding<-savis_embedding
  if(verbose){
    cat('\n')
    print("Finished...")
    setTxtProgressBar(pb = pb, value = 10)
  }
  return(object)
}


savis_nth<- function(x, k) {
  if(sum(is.na(x))>0){
    x[is.na(x)]<-min(x[!is.na(x)])-0.1
  }
  ## might have problem when k is too large for nan case
  p <- length(x) - k
  if(p < 0){
    stop("savis_nth: input k too larger") 
  }else if(p == 0){
    res<-1:length(x)
  }else{
    xp <- base::sort(x, partial=p)[p]
    res<-which(x > xp)
  }
  res
}

#' @importFrom Seurat FindVariableFeatures
#' @importFrom cli print.boxx
ExpFindVariableFeatures<-function(expr_matrix,verbose=F){
  if (!inherits(x = expr_matrix, 'Matrix')) {
    expr_matrix <- as(object = as.matrix(x = expr_matrix), Class = 'Matrix')
  }
  if (!inherits(x = expr_matrix, what = 'dgCMatrix')) {
    expr_matrix <- as(object = expr_matrix, Class = 'dgCMatrix')
  }
  expr_matrix@x<-exp(expr_matrix@x)-1
  FindVariableFeatures(
    expr_matrix,
    verbose = verbose)$vst.variance.standardized
}

#' RunAdaUMAP
#'
#' Use Adaptive Distance Metric to Run UMAP
#'
#' @details This function argument to the function
#'
#' @param X sds
#' @param metric = 'euclidean',
#' @param metric_count = 1,
#' @param py_envir = globalenv(),
#' @param n.neighbors = 30L,
#' @param n.components = 2L,
#' @param n.epochs = NULL,
#' @param learning.rate = 1.0,
#' @param min.dist = 0.3,
#' @param spread = 1.0,
#' @param set.op.mix.ratio = 1.0,
#' @param local.connectivity = 1L,
#' @param repulsion.strength = 1,
#' @param negative.sample.rate = 5,
#' @param a = 1.8956, 
#' @param b = 0.8006,
#' @param uwot.sgd = FALSE,
#' @param seed.use = 42L,
#' @param metric.kwds = NULL,
#' @param angular.rp.forest = FALSE,
#' @param verbose = FALSE
#'
#' @return nothing useful
#'
#' @importFrom reticulate py_run_string import import_main py_get_attr py_module_available py_set_seed
#' @importFrom glue glue
#'
#' @examples
#' a<-1
#'
RunAdaUMAP<-function(
  X,
  metric = 'euclidean',
  metric_count = 1,
  py_envir = globalenv(),
  n.neighbors = 30L,
  n.components = 2L,
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = 1.8956, 
  b = 0.8006,
  uwot.sgd = FALSE,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  verbose = FALSE){
  if(inherits(X, "list")){
    X<-get_matrix_from_list(X)
  }
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip in command line(e.g. pip install umap-learn) or through reticulate package in R (e.g. reticulate::py_install('umap') )")
  }
  py_func_names<-c("adaptive_dist_grad",
    "adaptive_dist2_grad",
    "adaptive_dist3_grad",
    "adaptive_dist_general_grad")
  
  # source the python script into the main python module
  py_run_string(glue(
    "
import numba 
import numpy as np
import warnings
numba.set_num_threads(16)
from umap import distances as dist
py_metric='{metric}' 
py_dist = dist.named_distances_with_gradients[py_metric]
warnings.filterwarnings('ignore')
@numba.njit(fastmath=True)
def adaptive_dist_grad(x, y):
    result = 0.0
    npcs = int((len(x)-1)/2)
    if x[0] != y[0]:
        d,grad = py_dist(x[1:(npcs+1)],y[1:(npcs+1)])
    else:
        d,grad = py_dist(x[(npcs+1):(2*npcs+1)],y[(npcs+1):(2*npcs+1)])
    return d, grad

@numba.njit(fastmath=True)
def adaptive_dist2_grad(x, y):
    result = 0.0
    npcs = int((len(x)-2)/3)
    if x[0] != y[0]:
        d,grad = py_dist(x[2:(npcs+2)],y[2:(npcs+2)])
    else:
        if x[1] != y[1] or x[1] == -1:
            d,grad = py_dist(x[(npcs+2):(2*npcs+2)],y[(npcs+2):(2*npcs+2)])
        else:
            d,grad = py_dist(x[(2*npcs+2):(3*npcs+2)],y[(2*npcs+2):(3*npcs+2)])
    return d, grad 

@numba.njit(fastmath=True)
def adaptive_dist3_grad(x, y):
    
    result = 0.0
    npcs = int((len(x)-3)/4)
    if x[0] != y[0]:
        d,grad = py_dist(x[3:(npcs+3)],y[3:(npcs+3)])
    else:
        if x[1] != y[1] or x[1] == -1:
            d,grad = py_dist(x[(npcs+3):(2*npcs+3)],y[(npcs+3):(2*npcs+3)])
        else:
            if x[2] != y[2] or x[2] == -1:
                d,grad = py_dist(x[(2*npcs+3):(3*npcs+3)],y[(2*npcs+3):(3*npcs+3)])
            else:
                d,grad = py_dist(x[(3*npcs+3):(4*npcs+3)],y[(3*npcs+3):(4*npcs+3)])
    return d, grad

@numba.njit(fastmath=True)
def adaptive_dist_general_grad(x, y):
    result = 0.0
    num_layer = x[0]
    npcs = int((len(x)-num_layer)/num_layer)
    processed = False
    for layer in range(1,num_layer):
        if x[layer] != y[layer]:
            print(layer)
            d,grad = py_dist(x[((layer-1)*npcs+num_layer):(layer*npcs+num_layer)],y[((layer-1)*npcs+num_layer):(layer*npcs+num_layer)])
            processed = True
            break
    if not processed:
        d,grad=py_dist(x[((num_layer-1)*npcs+num_layer):(num_layer*npcs+num_layer)],y[((num_layer-1)*npcs+num_layer):(num_layer*npcs+num_layer)])
         
    return d, grad
")  
    ,local = FALSE, convert = TRUE)
  
  
  # copy objects from the main python module into the specified R environment
  py_main <- import_main(convert = TRUE)
  py_main_dict <- py_get_attr(py_main, "__dict__")
  
  Encoding(py_func_names) <- "UTF-8"
  #for (py_name in py_func_names){
  #  py_value <- py_main_dict[[py_name]]
  #  assign(py_name, py_value, envir = py_envir) 
  #}
  py_name<- py_func_names[metric_count]
  py_value <- py_main_dict[[py_name]]
  assign(py_name, py_value, envir = py_envir) 
  if (!is.null(x = seed.use)) {
    py_set_seed(seed = seed.use)
  }
  if (typeof(x = n.epochs) == "double") {
    n.epochs <- as.integer(x = n.epochs)
  }
  
  if (metric_count == 1){
    adaptive_metric<-adaptive_dist_grad
  }
  if (metric_count == 2){
    adaptive_metric<-adaptive_dist2_grad
  }
  if (metric_count == 3){
    adaptive_metric<-adaptive_dist3_grad
  }
  if (metric_count == 4){
    adaptive_metric<-adaptive_dist_general_grad
  }
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(x = n.neighbors),
    n_components = as.integer(x = n.components),
    metric = adaptive_metric,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    metric_kwds = metric.kwds,
    angular_rp_forest = angular.rp.forest,
    random_state=seed.use,
    verbose = verbose
  )
  umap_embedding<-umap$fit_transform(as.matrix(x = X))
  umap_embedding<-data.frame(umap_embedding)
  colnames(umap_embedding)<-paste0("SAVIS_",1:n.components)
  umap_embedding
}

### If the storage is compressed, it should be recovered 
## to be a matrix 
get_matrix_from_list<-function(
  combined_embedding_list){
  list_len<-length(combined_embedding_list)
  npcs<-ncol(combined_embedding_list[[list_len]])
  cluster_index<-which(substr(colnames(combined_embedding_list[[1]]),nchar(colnames(combined_embedding_list[[1]]))-6,nchar(colnames(combined_embedding_list[[1]])))
    == "Cluster")
  cluster_label<-combined_embedding_list[[1]][,cluster_index]
  ncol_cluster<-length(cluster_index)
  start_col<-0
  for ( i in 2:ncol_cluster){
    if(sum(cluster_label[,i]== -1) > 0){
      start_col<-i
      break
    }
  }
  combined_embedding<-combined_embedding_list[[1]]
  for ( i in 2:ncol_cluster){
    sub_PC_supp<-matrix(0,nrow = nrow(combined_embedding_list[[1]]),ncol=npcs)
    rownames(sub_PC_supp)<-rownames(combined_embedding_list[[1]])
    colnames(sub_PC_supp)<-paste0("Layer",(i+1),"PC",1:npcs)
    index_i<-which(cluster_label[,i]!= -1)
    sub_PC_supp[index_i,]<-as.matrix(combined_embedding_list[[i]])
    combined_embedding<-cbind(combined_embedding,sub_PC_supp)
  }
  combined_embedding
}


#' tsMDS
#'
#' two step MDS
#' 
#' @importFrom stats cmdscale dist
#' @importFrom Rfast Dist
#' @importFrom mize mize
#' @importFrom glue glue
#' @export
#'
#' @examples
#' a<-1
#'
tsMDS<-function(
  dist_full,
  main_index,
  dist_main=NULL){
  N<-nrow(dist_full)
  if(is.null(dist_main)){
    dist_main<-dist_full[main_index,main_index] 
  }
  main_initial<-cmdscale(dist(dist_main),k=2)
  ### First Step of two-step MDS  
  cost_fun <- function(R, D) {
    diff2 <- (R - D) ^ 2
    sum(diff2) * 0.5
  }
  
  cost_grad <- function(R, D, y) {
    K <- (R - D) / (D + 1.e-10)
    
    G <- matrix(nrow = nrow(y), ncol = ncol(y))
    
    for (i in 1:nrow(y)) {
      dyij <- sweep(-y, 2, -y[i, ])
      G[i, ] <- apply(dyij * K[, i], 2, sum)
    }
    
    as.vector(t(G)) * -2
  }
  
  mmds_fn <- function(par) {
    R <- dist_main
    y <- matrix(par, ncol = 2, byrow = TRUE)
    D <- Dist(y)
    #D <- as.matrix(parDist(y))
    cost_fun(R, D)
  }
  
  mmds_gr <- function(par) {
    R <- dist_main
    y <- matrix(par, ncol = 2, byrow = TRUE)
    D <- Dist(y)
    #D <- as.matrix(parDist(y))
    
    cost_grad(R, D, y)
  }
  
  initial_val_main<-c(t(main_initial))
  res_main <- mize(initial_val_main, list(fn = mmds_fn, gr = mmds_gr), 
    method = "L-BFGS", verbose = FALSE, 
    grad_tol = 1e-5, check_conv_every = 10)
  
  main_mds<-matrix(res_main$par, ncol = 2, byrow = TRUE)
  
  ### Second Step of two-step MDS
  remain_index<-c(1:N)[which(!c(1:N)%in%main_index)]
  if(length(remain_index) == 0){
    tsMDS_res<-main_mds
    return(tsMDS_res)
  }else if(length(remain_index) == 1){
    remain_initial<-c(0,0)
  }else{
    dist_remain<-dist_full[remain_index,remain_index]
    remain_initial<-cmdscale(dist(dist_remain),k=2) 
  }
    cost_fun <- function(R, D) {
      diff2 <- (R - D) ^ 2
      sum(diff2) * 0.5
    }
    cost_grad1 <- function(R, D, y1, y2) {
      K <- (R - D) / (D + 1.e-10)
      y<-rbind(y1,y2)
      G <- matrix(nrow = nrow(y)-nrow(y1), ncol = ncol(y))
      
      for (i in 1:(nrow(y)-nrow(y1))) {
        i1<-nrow(y1)+i
        dyij <- sweep(-y, 2, -y[i1, ])
        G[i, ] <- apply(dyij * K[, i1], 2, sum)
      }
      
      as.vector(t(G)) * -2
    }
    
    mmds_fn <- function(par) {
      R <- dist_full
      y1<- main_mds
      y2 <- matrix(par, ncol = 2, byrow = TRUE)
      y<-rbind(y1,y2)
      D <- Dist(y)
      #D <- as.matrix(parDist(y))
      cost_fun(R, D)
    }
    
    mmds_gr <- function(par) {
      R <- dist_full
      y1<- main_mds
      y2 <- matrix(par, ncol = 2, byrow = TRUE)
      y<-rbind(y1,y2)
      D <- Dist(y)
      #D <- as.matrix(parDist(y))
      cost_grad1(R, D, y1, y2)
    }
    initial_val_remain<-c(t(remain_initial))
    
    res_remain <- mize(initial_val_remain, list(fn = mmds_fn, gr = mmds_gr), 
      method = "L-BFGS", verbose = FALSE, 
      grad_tol = 1e-5, check_conv_every = 10)
    remain_mds<-matrix(res_remain$par, ncol = 2, byrow = TRUE)
    
    tsMDS_res<-rbind(main_mds,remain_mds)
  

  tsMDS_res
}


########## newly added part for UMAP adjust

rotation = function(x,y){
  u=x/sqrt(sum(x^2))
  
  v=y-sum(u*y)*u
  v=v/sqrt(sum(v^2))
  
  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
  
  sint=sqrt(1-cost^2);
  
  diag(length(x)) - u %*% t(u) - v %*% t(v) + 
    cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
}


Detect_edge<-function(
  whole,
  whole_mean,
  edge){
  rotation = function(x,y){
    u=x/sqrt(sum(x^2))
    
    v=y-sum(u*y)*u
    v=v/sqrt(sum(v^2))
    
    cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
    
    sint=sqrt(1-cost^2);
    
    diag(length(x)) - u %*% t(u) - v %*% t(v) + 
      cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
  }
  
  vec<-edge-whole_mean
  Rvec2xaxis<-rotation(as.numeric(vec),c(1,0))
  whole_rotated<-t(t(whole)-as.numeric(whole_mean))%*%Rvec2xaxis
  index_cone<-which(abs(whole_rotated[,2]/whole_rotated[,1])<0.5)
  index_cone_dist<-sapply(c(index_cone), function(i){
    sum(whole_rotated[i,]^2)
  })
  new_edge_index<-index_cone[which.max(index_cone_dist)]
  new_edge_index
}

#' tsMDS
#'
#' two step MDS
#' 
#' @importFrom pdist pdist
#'
#' @examples
#' a<-1
#'
Detect_farthest<-function(
  whole,
  whole_mean){
  a<-pdist(whole_mean,whole)
  which.max(a@dist)
}

#' ScaleFactor
#'
#' Combined PC embedding with scale factor for subPC
#'
#' @importFrom pdist pdist
#' @importFrom Rfast Dist
#' @importFrom stats dist as.dist
#' @importFrom uwot umap
#'
#'
#' @examples
#' a<-1
#'
#'
get_umap_embedding_adjust_umap<-function(
  pca_embedding,
  pca_center,
  pca_anchor_index,
  pca_dist,
  umap_embedding,
  N_label_,
  cluster_,
  label_index_,
  main_index = NULL,
  pca_dist_main=NULL,
  distance_metric = "euclidean",
  scale_factor = 0.9,
  rotate = TRUE,
  seed.use = 42
){
  Rotation2to1<-function(umap_center1,umap_center2,pos1,pos2){
    N_label_<-nrow(umap_center1)
    umap_center1_tmp<-t(t(umap_center1[-pos1,])-as.numeric(umap_center1[pos1,]))
    weight_1<-1/pdist(umap_center1_tmp,c(0,0))@dist
    weight_1<-weight_1/sum(weight_1)
    umap_center2_tmp<-t(t(umap_center2[-pos2,])-as.numeric(umap_center2[pos2,]))
    #weight_2<-1/pdist(umap_center2_tmp,c(0,0))@dist
    #weight_2<-weight_2/sum(weight_2)
    angles<-sapply(1:(N_label_-1), function(i){
      
      umap1<-umap_center1_tmp[i,]
      umap2<-umap_center2_tmp[i,]
      umap1<-umap1/sqrt(sum(umap1^2))
      umap2<-umap2/sqrt(sum(umap2^2))
      
      Rumap2toumap1<-rotation(umap2,umap1)
      Rumap2toumap1 <- pmax(Rumap2toumap1,-1)
      Rumap2toumap1 <- pmin(Rumap2toumap1,1)
      angle<-acos(Rumap2toumap1[1,1])
      
      if(Rumap2toumap1[2,1]>=0){
        angle<-acos(Rumap2toumap1[1,1])
      }else{
        angle<- -acos(Rumap2toumap1[1,1])
      }
      angle
    })
    #angle2to1<-mean(angles)
    #angle2to1<-median(angles)
    angle2to1<-sum(angles*weight_1)
    R2to1<-diag(cos(angle2to1),2)
    R2to1[2,1]<-sin(angle2to1)
    R2to1[1,2]<- -R2to1[2,1]
    R2to1
  }
  
  if(is.null(pca_dist_main)){
    pca_dist_main<-pca_dist[main_index,main_index]
  }
  set.seed(seed.use)
  main_umap_center <-
    umap(
      X = dist(pca_dist_main),
      n_neighbors = as.integer(x = length(main_index)-1),
      n_components = as.integer(x =2L),
      metric = distance_metric,
      learning_rate = 1.0,
      min_dist = 0.3,
      spread =  1.0,
      set_op_mix_ratio =  1.0,
      local_connectivity =  1.0,
      repulsion_strength = 1,
      negative_sample_rate = 5,
      fast_sgd = FALSE
    )
  colnames(main_umap_center)<-c("UMAP_1","UMAP_2")
  
  main_umap_center<-data.frame(main_umap_center)
  sf1<-(max(umap_embedding[,1])-min(umap_embedding[,1]))/(max(main_umap_center[,1]) -min(main_umap_center[,1]))
  sf2<-(max(umap_embedding[,2])-min(umap_embedding[,2]))/(max(main_umap_center[,2]) -min(main_umap_center[,2]))
  main_umap_center[,1]<-main_umap_center[,1]*sf1*scale_factor
  main_umap_center[,2]<-main_umap_center[,2]*sf2*scale_factor
  umap_embedding_mean<-t(sapply(1:N_label_, function(i){
    index_i<-which(cluster_ == label_index_[i])
    colMeans(as.matrix(umap_embedding[index_i,]))
  }))
  umap_embedding_adjust<-umap_embedding
  
  if (rotate){
    
    main_umap_center_flip<-main_umap_center
    main_umap_center_flip[,1]<-main_umap_center_flip[,1]*(-1)
    angle_var<-c()
    weight_sample<-c()
    angle_no_flip<-list()
    angle_flip<-list()
    for (i in 1:length(main_index)){
      
      weight_sample<-c(weight_sample,sum(cluster_ == label_index_[main_index[i]]))
      R2to1<-Rotation2to1(
        umap_center1 = main_umap_center,
        umap_center2 = pca_center[,c(1,2)],
        pos1 = i,
        pos2 = main_index[i])  
      
      R2to1_flip<-Rotation2to1(
        umap_center1 = main_umap_center_flip,
        umap_center2 = pca_center[,c(1,2)],
        pos1 = i,
        pos2 = main_index[i])  
      
      angle_vec<-sapply(1:length(pca_anchor_index[[i]]), function(j){
        y<-as.numeric(pca_embedding[pca_anchor_index[[i]][j],c(1,2)]-pca_center[main_index[i],c(1,2)])%*%R2to1
        y<-as.numeric(y)
        y<-y/sqrt(sum(y^2))
        x<-umap_embedding[pca_anchor_index[[i]][j],]-umap_embedding_mean[main_index[i],]
        x<-as.numeric(x)
        y<-y/sqrt(sum(y^2))
        x<-x/sqrt(sum(x^2))
        
        Rx2y <- rotation(x,y)
        Rx2y <- pmax(Rx2y,-1)
        Rx2y <- pmin(Rx2y,1)
        
        if(Rx2y[2,1]>=0){
          angle<-acos(Rx2y[1,1])
        }else{
          angle<- - acos(Rx2y[1,1])
        }
        angle
      })
      angle_vec_flip<-sapply(1:length(pca_anchor_index[[i]]), function(j){
        y<-as.numeric(pca_embedding[pca_anchor_index[[i]][j],c(1,2)]-pca_center[main_index[i],c(1,2)])%*%R2to1_flip
        y<-as.numeric(y)
        y<-y/sqrt(sum(y^2))
        x<-umap_embedding[pca_anchor_index[[i]][j],]-umap_embedding_mean[main_index[i],]
        x<-as.numeric(x)
        y<-y/sqrt(sum(y^2))
        x<-x/sqrt(sum(x^2))
        
        Rx2y <- rotation(x,y)
        
        
        if(Rx2y[2,1]>=0){
          angle<- acos(Rx2y[1,1])
        }else{
          angle<- -1*acos(Rx2y[1,1])
        }
        angle
      })
      
      angle_vec<-angle_vec[angle_vec>min(angle_vec)]
      angle_vec<-angle_vec[angle_vec<max(angle_vec)]
      angle_vec_flip<-angle_vec_flip[angle_vec_flip>min(angle_vec_flip)]
      angle_vec_flip<-angle_vec_flip[angle_vec_flip<max(angle_vec_flip)]
      angle_no_flip[[i]]<- angle_vec
      angle_flip[[i]] <- angle_vec_flip
      angle_var<-rbind(angle_var,(c(var(angle_vec),var(angle_vec_flip))))
    }
    weight_sample<-weight_sample/sum(weight_sample)
    if (sum(angle_var[,1]*weight_sample)>=sum(angle_var[,2]*weight_sample)){
      # use flip
      angle_vec<-angle_flip
    }else{
      # use no flip
      angle_vec<-angle_no_flip
    }
    
    for ( i in 1:length(main_index)){
      anglex2y<-mean(angle_vec[[i]])
      #print(anglex2y)
      Rx2y<-diag(cos(anglex2y),2)
      Rx2y[2,1]<-sin(anglex2y)
      Rx2y[1,2]<- -Rx2y[2,1]
      index_i<-which(cluster_ == label_index_[main_index[i]])
      umap_embedding_adjust[index_i,]<-t(t(t(t(umap_embedding[index_i,])-as.numeric(umap_embedding_mean[main_index[i],]))%*%Rx2y)+as.numeric(main_umap_center[i,]))
    }
    
  }else{
    
    for(i in 1:length(main_index)){
      index_i<-which(cluster_ == label_index_[main_index[i]])
      umap_embedding_adjust[index_i,]<-t(t(umap_embedding[index_i,])-as.numeric(umap_embedding_mean[main_index[i],]-main_umap_center[i,]))
    }
  }
  newList<-list("main_umap_center"=main_umap_center,
    "umap_embedding_adjust"=umap_embedding_adjust)
  return(newList)
}




#' ScaleFactor
#'
#' Combined PC embedding with scale factor for subPC
#'
#' @importFrom pdist pdist
#' @importFrom stats dist as.dist
#' @importFrom uwot umap
#'
#'
#' @examples
#' a<-1
#'
#'
adjustUMAP_via_umap<-function(
  pca_embedding,
  umap_embedding,
  cluster_label,
  global_umap_embedding = NULL,
  distance_metric = "euclidean",
  scale_factor = 0.9,
  rotate = TRUE,
  density_adjust = TRUE,
  seed.use = 42,
  min_size = 100,
  maxit_push = NULL
){
 if(!is.matrix(pca_embedding)){
    pca_embedding<-as.matrix(pca_embedding)
 }
  if(!is.matrix(umap_embedding)){
    umap_embedding<-as.matrix(umap_embedding)
  }
  if(is.null(rownames(umap_embedding)[1])){
    rownames(umap_embedding)<-1:nrow(umap_embedding)
  }
  # This is for the clustering results
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  
  size_cluster<-c()
  for ( i in 1:N_label){
    size_cluster<-c(size_cluster,
      sum(cluster_label == label_index[i]))
  }
  snn_<- FindNeighbors(object = umap_embedding,
    nn.method = "rann",
    verbose = F)$snn
  cluster_ <- FindClusters(snn_,
    resolution = 0,
    verbose = F)[[1]]
  cluster_ <- as.numeric(as.character(cluster_))
  
  N_label_<-length(unique(cluster_))
  label_index_<-sort(unique(cluster_))
  if(N_label_ <= 2){
    return(umap_embedding)
  }
  cluster_size_<-sapply(1:N_label_, function(i){
    sum(cluster_==label_index_[i])
  })
  N_sample<-nrow(pca_embedding)
  cutoff_main_remain<-0.01*N_sample 
  main_index<-c(1:N_label_)[which(cluster_size_ > cutoff_main_remain)]
  remain_index<-c(1:N_label_)[which(!c(1:N_label_)%in%main_index)]  
  length_large_main_index<-length(main_index)
  cutoff_small_size_cluster<-mean(size_cluster)
  small_size_cluster_index<-which(size_cluster < cutoff_small_size_cluster)
  large_size_cluster_index<-which(size_cluster >= cutoff_small_size_cluster)
  
  cluster_remain_index_collection<-list()
  
  for(i in large_size_cluster_index){
    index_i<-which(cluster_label == label_index[i])
    cluster_i_collect<-which(label_index_%in%unique(cluster_[index_i]))
    cluster_i_collect<-intersect(cluster_i_collect,remain_index)
    cluster_remain_index_collection[[i]]<-cluster_i_collect
  }
  
  for(i in small_size_cluster_index){
    index_i<-which(cluster_label == label_index[i])
    cluster_i_collect<-which(label_index_%in%unique(cluster_[index_i]))
    cluster_i_sizes<-sapply(cluster_i_collect, function(j){
      sum(cluster_[index_i] == label_index_[j])
    }) 
    if(max(cluster_i_sizes)>cutoff_main_remain){
      cluster_i_collect<-intersect(cluster_i_collect,remain_index)
    }else{
      tmp<- cluster_i_collect[which.max(cluster_i_sizes)]
      main_index<-c(main_index,tmp)
      remain_index<-remain_index[remain_index!=tmp]
      cluster_i_collect<-intersect(cluster_i_collect,remain_index) 
    }
    cluster_remain_index_collection[[i]]<-cluster_i_collect
  }
  main_index<-main_index[order(main_index)]
  
  if (density_adjust & !is.null(global_umap_embedding)){
    prop_density<-sapply(1:length(main_index), function(j){
      i<-main_index[j]
      index_i<-which(cluster_ == label_index_[i])
      set.seed(seed.use)
      sample_index_i<-sample(index_i,min(min_size,length(index_i)) )
      sample_global_dist<-Dist(global_umap_embedding[sample_index_i,])
      #sample_global_dist<-as.matrix(parDist(global_umap_embedding[sample_index_i,]))
      sample_local_dist<-Dist(umap_embedding[sample_index_i,])
      #sample_local_dist<-as.matrix(parDist(umap_embedding[sample_index_i,]))
      mean(c(sample_global_dist))/mean(c(sample_local_dist))
    })
    for(j in 1:length(main_index)){
      i<-main_index[j]
      index_i<-which(cluster_ == label_index_[i])
      cur_umap<-umap_embedding[index_i,]
      umap_embedding[index_i,]<-t((t(cur_umap)-as.numeric(colMeans(cur_umap)))*min(3,prop_density[j])+as.numeric(colMeans(cur_umap)))
    }
    for(j in 1:length(remain_index)){
      i<-remain_index[j]
      index_i<-which(cluster_ == label_index_[i])
      cur_umap<-umap_embedding[index_i,]
      umap_embedding[index_i,]<-t((t(cur_umap)-as.numeric(colMeans(cur_umap)))*0.5+as.numeric(colMeans(cur_umap)))
    }
  }else if(density_adjust & is.null(global_umap_embedding)){
    for(j in 1:length(main_index)){
      i<-main_index[j]
      index_i<-which(cluster_ == label_index_[i])
      cur_umap<-umap_embedding[index_i,]
      if(j <= length_large_main_index){
        cur_sf_here<-1.5
      }else{
        cur_sf_here<-3
      }
      
      umap_embedding[index_i,]<-t((t(cur_umap)-
          as.numeric(colMeans(cur_umap)))*cur_sf_here+
          as.numeric(colMeans(cur_umap)))
    }
    for(j in 1:length(remain_index)){
      i<-remain_index[j]
      index_i<-which(cluster_ == label_index_[i])
      cur_umap<-umap_embedding[index_i,]
      umap_embedding[index_i,]<-t((t(cur_umap)-as.numeric(colMeans(cur_umap)))*0.5+as.numeric(colMeans(cur_umap)))
    }
  }
  
  pca_center<-t(sapply(1:N_label_, function(i){
    index_i<-which(cluster_ == label_index_[i])
    colMeans(pca_embedding[index_i,])
  }))
  
  pca_anchor_index<-lapply(main_index, function(i){
    index_i<-which(cluster_ == label_index_[i])
    set.seed(seed.use)
    sample_index_i<-sample(index_i,min(min_size,length(index_i)) )
    sample_index_dist<-pdist(pca_embedding[sample_index_i,c(1,2)],pca_center[i,c(1,2)])@dist
    savis_nth(x=sample_index_dist,
      k = max(1,min(ceiling(min_size/5),length(index_i))))
  })
  
  pca_dist1<-Dist(pca_center)
  #pca_dist1<-as.matrix(parDist(pca_center))
  
  step1_res<-get_umap_embedding_adjust_umap(
    pca_embedding=pca_embedding,
    pca_center=pca_center,
    pca_anchor_index=pca_anchor_index,
    pca_dist=pca_dist1,
    pca_dist_main = NULL,
    main_index = main_index,
    distance_metric=distance_metric,
    umap_embedding=umap_embedding,
    N_label_=N_label_,
    cluster_ = cluster_,
    label_index_ = label_index_,
    scale_factor =scale_factor,
    rotate = rotate,
    seed.use = seed.use)
  main_umap_center<-step1_res$main_umap_center
  umap_embedding_adjust<-step1_res$umap_embedding_adjust
  umap_embedding_adjust_main<-umap_embedding_adjust[which(cluster_%in%label_index_[main_index]),]
  cluster_main<-cluster_[which(cluster_%in%label_index_[main_index])]
  for(i in 1:length(main_index)){
    cluster_main[which(cluster_main == label_index_[main_index[i]])]<- i-1
  }
  label_index_main<-sort(unique(cluster_main))
  snn_1<- FindNeighbors(object = umap_embedding_adjust_main,
    nn.method = "rann",
    verbose = F)$snn
  cluster_1 <- FindClusters(snn_1,
    resolution = 0,
    verbose = F)[[1]]
  cluster_1 <- as.numeric(as.character(cluster_1))
  
  N_label_1<-length(unique(cluster_1))
  if (N_label_1 < length(main_index)){
    ## First Adjustment
    label_index_1<-sort(unique(cluster_1))
    
    bad_index<-list()
    list_index<-0
    for ( i in 1:N_label_1){
      index_i1<-which(cluster_1 == label_index_1[i])
      cur_index_len<-length(unique(cluster_main[index_i1]))
      if(cur_index_len > 1){
        list_index<-list_index+1
        bad_index[[list_index]]<-c(unique(cluster_main[index_i1]))
      }
    }
    
    cur_iter<-0
    if(is.null(maxit_push)){
      maxit_push<-min(N_label_/3,10)
    }
    while (N_label_1 < length(main_index) & cur_iter < maxit_push) {
      cur_iter<-cur_iter+1
      for (i in 1:length(bad_index)){
        pos<-min(bad_index[[i]])
        other_pos<-bad_index[[i]][bad_index[[i]]>pos]
        pos<-pos+1
        other_pos<-other_pos+1
        dist_mat<-pdist(main_umap_center[pos,],main_umap_center)@dist
        target_distance<-min(dist_mat[dist_mat>max(dist_mat[other_pos])])
        target_distance<-(target_distance+dist_mat[other_pos])/2
        for (k in 1:length(other_pos)){
          cur_pos<-other_pos[k]
          cur_center<-main_umap_center[pos,]
          cur_arrow<-main_umap_center[cur_pos,]
          prop<-target_distance[k]/dist_mat[cur_pos]
          cur_arrow<-(cur_arrow-cur_center)*prop+cur_center
          index_i<-which(cluster_main == label_index_main[cur_pos])
          umap_embedding_adjust_main[index_i,]<-t(t(umap_embedding_adjust_main[index_i,])-as.numeric(main_umap_center[cur_pos,])+as.numeric(cur_arrow))
        }
      }
      
      snn_1<- FindNeighbors(object = umap_embedding_adjust_main,
        nn.method = "rann",
        verbose = F)$snn
      cluster_1 <- FindClusters(snn_1,
        resolution = 0,
        verbose = F)[[1]]
      cluster_1 <- as.numeric(as.character(cluster_1))
      N_label_1<-length(unique(cluster_1))
      label_index_1<-sort(unique(cluster_1))
      
      bad_index<-list()
      list_index<-0
      for ( i in 1:N_label_1){
        index_i1<-which(cluster_1 == label_index_1[i])
        cur_index_len<-length(unique(cluster_main[index_i1]))
        if(cur_index_len > 1){
          list_index<-list_index+1
          bad_index[[list_index]]<-c(unique(cluster_main[index_i1]))
        }
      }
      
    }
    
    if(N_label_1 < length(main_index)){
      #if(verbose){
      #  cat('\n')
      #  print("Rescaling SAVIS...")
      #  setTxtProgressBar(pb = pb, value = 19.5)
      #}
      ## Third Adjustment Rescale
      for (i in 1:length(bad_index)){
        pos<-min(bad_index[[i]])
        other_pos<-bad_index[[i]][bad_index[[i]]>pos]
        pos<-pos+1
        other_pos<-other_pos+1
        dist_mat<-pdist(main_umap_center[pos,],main_umap_center)@dist
        target_distance<-min(dist_mat[dist_mat>=max(dist_mat[other_pos])])
        re_sf<-target_distance/max(dist_mat[other_pos])
        
        umap_embedding_adjust<-umap_embedding_adjust*re_sf
        umap_embedding_adjust_main<-umap_embedding_adjust_main*re_sf
      }
    }
  }
  
  
  # Return to the whole umap_embedding,
  umap_embedding_adjust[which(cluster_%in%label_index_[main_index]),]<-umap_embedding_adjust_main
  
  ###############---Begin to deal with small clusters---#########
  
  unique_cluster_remain_index<-sort(unique(unlist(cluster_remain_index_collection)))
  unlist_cluster_remain_index<-unlist(cluster_remain_index_collection)
  unique_cluster_remain_index_count<-sapply(unique_cluster_remain_index, function(i){
    sum(unlist_cluster_remain_index == i)
  })
  duplicated_remain_index<-unique_cluster_remain_index[which(unique_cluster_remain_index_count > 1)]
  
  for( dup_index in duplicated_remain_index){
    index_dup<-which(cluster_ == label_index_[dup_index])
    location_of_cluster<-unique(cluster_label[index_dup])
    size_location_of_cluster<-sapply(location_of_cluster, function(i){
      sum(cluster_label[index_dup] == i)
    })
    location_of_cluster[which.max(size_location_of_cluster)]
    rm_index<-which(label_index%in%location_of_cluster)[-which.max(size_location_of_cluster)]
    for(i in rm_index){
      tmp<-cluster_remain_index_collection[[i]]
      tmp<-tmp[which(tmp!=dup_index)]
      cluster_remain_index_collection[[i]]<-tmp
    }
  }
  
  for(i in c(large_size_cluster_index,small_size_cluster_index)){
    cluster_i_collect<-cluster_remain_index_collection[[i]]
    if(length(cluster_i_collect)>0){
      ## Do something when there is thing in list
      index_i<-which(cluster_label == label_index[i])
      index_remain<-which(cluster_ %in% label_index_[cluster_i_collect])
      index_i_noremain<-setdiff(index_i,index_remain)
      
      cur_remain_pc_mean<-t(sapply(cluster_i_collect, function(cur_remain_index){
        index_cur_remain<-which(cluster_ == label_index_[cur_remain_index])
        colMeans(pca_embedding[index_cur_remain,])
      }))
      
      res_knn<-FNN::knn(train = pca_embedding[index_i_noremain,],
        test = cur_remain_pc_mean,cl = cluster_[index_i_noremain],k = 1)
      close_cluster_label_<-res_knn[1:length(cluster_i_collect)]
      close_index<-index_i_noremain[attr(res_knn,"nn.index")[,1]]
      
      #unique_close_cluster_label_<-unique(close_cluster_label_)
      #close_umap_mean<-t(sapply(unique_close_cluster_label_, function(unique_i){
      #  index_unique_close_cluster_label_<-which(cluster_==unique_i)
      #  colMeans(umap_embedding_adjust[index_unique_close_cluster_label_,])
      #}))
      
      # go back to umap, check if these selected points are at edge
      for(cur_ in 1:length(close_cluster_label_)){
        
        index_cur<-which(cluster_==close_cluster_label_[cur_])
        index_cur_remain<-which(cluster_==label_index_[cluster_i_collect[cur_]])
        whole_remain<-umap_embedding_adjust[index_cur_remain,]
        whole_remain_mean<-colMeans(whole_remain)
        farthest_remain_index<-Detect_farthest(whole_remain,whole_remain_mean)
        remain_vec1<-whole_remain[farthest_remain_index,]
        #cur_close_umap_mean<-close_umap_mean[
        #  which(unique_close_cluster_label_ 
        #    == close_cluster_label_[cur_close_index]),]
        cur_edge_umap<-umap_embedding_adjust[close_index[cur_],]
        whole<-umap_embedding_adjust[index_cur,]
        whole_mean<-colMeans(whole)
        cur_close_index<-Detect_edge(whole,whole_mean,cur_edge_umap)
        cur_close_vec1<-whole[cur_close_index,]
        extra_vec<-cur_close_vec1-whole_mean
        cur_close_vec2<-cur_close_vec1+0.2*extra_vec
        Rremain2cur<-rotation(as.numeric(remain_vec1 - whole_remain_mean),
          as.numeric(whole_mean-cur_close_vec1))
        whole_remain_adjust<-t(t(whole_remain)-whole_remain_mean)%*%Rremain2cur
        whole_remain_adjust<-t(t(whole_remain_adjust)-as.numeric(whole_remain_adjust[farthest_remain_index,])+as.numeric(cur_close_vec2))
        umap_embedding_adjust[index_cur_remain,]<-whole_remain_adjust
      }
    }
  }
  
  
  return(umap_embedding_adjust)
}


get_umap_embedding_adjust_tsMDS<-function(
  pca_embedding,
  pca_center,
  pca_anchor_index,
  pca_dist,
  umap_embedding,
  N_label,
  cluster_,
  label_index,
  adjust_method = "tsMDS",
  main_index = NULL,
  pca_dist_main=NULL,
  distance_metric = "euclidean",
  scale_factor = 0.9,
  rotate = TRUE,
  seed.use = 42
){
  Rotation2to1<-function(umap_center1,umap_center2,pos){
    N_label<-nrow(umap_center1)
    umap_center1_tmp<-t(t(umap_center1[-pos,])-as.numeric(umap_center1[pos,]))
    weight_1<-1/pdist(umap_center1_tmp,c(0,0))@dist
    weight_1<-weight_1/sum(weight_1)
    umap_center2_tmp<-t(t(umap_center2[-pos,])-as.numeric(umap_center2[pos,]))
    #weight_2<-1/pdist(umap_center2_tmp,c(0,0))@dist
    #weight_2<-weight_2/sum(weight_2)
    angles<-sapply(1:(N_label-1), function(i){
      
      umap1<-umap_center1_tmp[i,]
      umap2<-umap_center2_tmp[i,]
      umap1<-umap1/sqrt(sum(umap1^2))
      umap2<-umap2/sqrt(sum(umap2^2))
      
      Rumap2toumap1<-rotation(umap2,umap1)
      Rumap2toumap1 <- pmax(Rumap2toumap1,-1)
      Rumap2toumap1 <- pmin(Rumap2toumap1,1)
      angle<-acos(Rumap2toumap1[1,1])
      
      if(Rumap2toumap1[2,1]>=0){
        angle<-acos(Rumap2toumap1[1,1])
      }else{
        angle<- -acos(Rumap2toumap1[1,1])
      }
      angle
    })
    #angle2to1<-mean(angles)
    #angle2to1<-median(angles)
    angle2to1<-sum(angles*weight_1)
    R2to1<-diag(cos(angle2to1),2)
    R2to1[2,1]<-sin(angle2to1)
    R2to1[1,2]<- -R2to1[2,1]
    R2to1
  }
  
  if (adjust_method == "tsMDS"){
    set.seed(seed.use)
    umap_center<-tsMDS(dist_full = pca_dist,
      main_index = main_index,
      dist_main = pca_dist_main)
    colnames(umap_center)<-c("tsMDS_1","tsMDS_2")
  }else if (adjust_method == "isoMDS"){
    umap_center<-isoMDS(dist(pca_dist))$points 
    colnames(umap_center)<-c("isoMDS_1","isoMDS_2")
  }else{
    stop("wrong adjust method")
  }
  
  umap_center<-data.frame(umap_center)
  sf1<-(max(umap_embedding[,1])-min(umap_embedding[,1]))/(max(umap_center[,1]) -min(umap_center[,1]))
  sf2<-(max(umap_embedding[,2])-min(umap_embedding[,2]))/(max(umap_center[,2]) -min(umap_center[,2]))
  umap_center[,1]<-umap_center[,1]*sf1*scale_factor
  umap_center[,2]<-umap_center[,2]*sf2*scale_factor
  umap_embedding_mean<-t(sapply(1:N_label, function(i){
    index_i<-which(cluster_ == label_index[i])
    colMeans(as.matrix(umap_embedding[index_i,]))
  }))
  umap_embedding_adjust<-umap_embedding
  
  if (rotate){
    
    umap_center_flip<-umap_center
    umap_center_flip[,1]<-umap_center_flip[,1]*(-1)
    angle_var<-c()
    weight_sample<-c()
    angle_no_flip<-list()
    angle_flip<-list()
    for (i in 1:N_label){
      weight_sample<-c(weight_sample,sum(cluster_ == label_index[i]))
      R2to1<-Rotation2to1(
        umap_center1 = umap_center,
        umap_center2 = pca_center[,c(1,2)],
        pos = i)  
      
      R2to1_flip<-Rotation2to1(
        umap_center1 = umap_center_flip,
        umap_center2 = pca_center[,c(1,2)],
        pos = i) 
      
      angle_vec<-sapply(1:length(pca_anchor_index[[i]]), function(j){
        y<-as.numeric(pca_embedding[pca_anchor_index[[i]][j],c(1,2)]-pca_center[i,c(1,2)])%*%R2to1
        y<-as.numeric(y)
        y<-y/sqrt(sum(y^2))
        x<-umap_embedding[pca_anchor_index[[i]][j],]-umap_embedding_mean[i,]
        x<-as.numeric(x)
        y<-y/sqrt(sum(y^2))
        x<-x/sqrt(sum(x^2))
        
        Rx2y <- rotation(x,y)
        Rx2y <- pmax(Rx2y,-1)
        Rx2y <- pmin(Rx2y,1)
        
        if(Rx2y[2,1]>=0){
          i
          angle<-acos(Rx2y[1,1])
        }else{
          angle<- - acos(Rx2y[1,1])
        }
        angle
      })
      angle_vec_flip<-sapply(1:length(pca_anchor_index[[i]]), function(j){
        y<-as.numeric(pca_embedding[pca_anchor_index[[i]][j],c(1,2)]-pca_center[i,c(1,2)])%*%R2to1_flip
        y<-as.numeric(y)
        y<-y/sqrt(sum(y^2))
        x<-umap_embedding[pca_anchor_index[[i]][j],]-umap_embedding_mean[i,]
        x<-as.numeric(x)
        y<-y/sqrt(sum(y^2))
        x<-x/sqrt(sum(x^2))
        
        Rx2y <- rotation(x,y)
        
        
        if(Rx2y[2,1]>=0){
          angle<- acos(Rx2y[1,1])
        }else{
          angle<- -1*acos(Rx2y[1,1])
        }
        angle
      })
      
      angle_vec<-angle_vec[angle_vec>min(angle_vec)]
      angle_vec<-angle_vec[angle_vec<max(angle_vec)]
      angle_vec_flip<-angle_vec_flip[angle_vec_flip>min(angle_vec_flip)]
      angle_vec_flip<-angle_vec_flip[angle_vec_flip<max(angle_vec_flip)]
      angle_no_flip[[i]]<- angle_vec
      angle_flip[[i]] <- angle_vec_flip
      angle_var<-rbind(angle_var,(c(var(angle_vec),var(angle_vec_flip))))
    }
    weight_sample<-weight_sample/sum(weight_sample)
    if (sum(angle_var[,1]*weight_sample)>=sum(angle_var[,2]*weight_sample)){
      # use flip
      angle_vec<-angle_flip
    }else{
      # use no flip
      angle_vec<-angle_no_flip
    }
    
    for ( i in 1:N_label){
      anglex2y<-mean(angle_vec[[i]])
      #print(anglex2y)
      Rx2y<-diag(cos(anglex2y),2)
      Rx2y[2,1]<-sin(anglex2y)
      Rx2y[1,2]<- -Rx2y[2,1]
      index_i<-which(cluster_ == label_index[i])
      umap_embedding_adjust[index_i,]<-t(t(t(t(umap_embedding[index_i,])-as.numeric(umap_embedding_mean[i,]))%*%Rx2y)+as.numeric(umap_center[i,]))
    }
    
  }else{
    
    for(i in 1:N_label){
      index_i<-which(cluster_ == label_index[i])
      umap_embedding_adjust[index_i,]<-t(t(umap_embedding[index_i,])-as.numeric(umap_embedding_mean[i,]-umap_center[i,]))
    }
  }
  newList<-list("umap_center"=umap_center,
    "umap_embedding_adjust"=umap_embedding_adjust)
  return(newList)
}


#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom pdist pdist
adjustUMAP_via_tsMDS<-function(
  pca_embedding,
  umap_embedding,
  global_umap_embedding = NULL,
  distance_metric = "euclidean",
  scale_factor = 0.9,
  rotate = TRUE,
  density_adjust = TRUE,
  shrink_distance = TRUE,
  adjust_method = "tsMDS",
  shrink_factor = 0.2,
  seed.use = 42,
  min_size = 100,
  maxit_push = NULL
){
  if(!is.matrix(pca_embedding)){
    pca_embedding<-as.matrix(pca_embedding)
  }
  if(!is.matrix(umap_embedding)){
    umap_embedding<-as.matrix(umap_embedding)
  }
  if(is.null(rownames(umap_embedding)[1])){
    rownames(umap_embedding)<-1:nrow(umap_embedding)
  }
  snn_<- FindNeighbors(object = umap_embedding,
    nn.method = "rann",
    verbose = F)$snn
  cluster_ <- FindClusters(snn_,
    resolution = 0,
    verbose = F)[[1]]
  cluster_ <- as.numeric(as.character(cluster_))
  
  N_label<-length(unique(cluster_))
  label_index<-sort(unique(cluster_))
  if(N_label <= 2){
    return(umap_embedding)
  }
  cluster_size<-sapply(1:N_label, function(i){
    sum(cluster_==label_index[i])
  })
  N_sample<-nrow(pca_embedding)
  
  main_index<-c(1:N_label)[which(cluster_size > 0.015*N_sample)]
  if (density_adjust & !is.null(global_umap_embedding)){
    prop_density<-sapply(1:N_label, function(i){
      index_i<-which(cluster_ == label_index[i])
      set.seed(seed.use)
      sample_index_i<-sample(index_i,min(min_size,length(index_i)) )
      sample_global_dist<-Dist(global_umap_embedding[sample_index_i,])
      #sample_global_dist<-as.matrix(parDist(global_umap_embedding[sample_index_i,]))
      sample_local_dist<-Dist(umap_embedding[sample_index_i,])
      #sample_local_dist<-as.matrix(parDist(umap_embedding[sample_index_i,]))
      mean(c(sample_global_dist))/mean(c(sample_local_dist))
    })
    for(i in 1:N_label){
      index_i<-which(cluster_ == label_index[i])
      cur_umap<-umap_embedding[index_i,]
      umap_embedding[index_i,]<-t((t(cur_umap)-as.numeric(colMeans(cur_umap)))*min(3,prop_density[i])+as.numeric(colMeans(cur_umap)))
    }
  }else if (density_adjust & is.null(global_umap_embedding)){
    for(j in 1:length(main_index)){
      i<-main_index[j]
      index_i<-which(cluster_ == label_index[i])
      cur_umap<-umap_embedding[index_i,]
      cur_sf_here<-1.5
      umap_embedding[index_i,]<-t((t(cur_umap)-as.numeric(colMeans(cur_umap)))*cur_sf_here+as.numeric(colMeans(cur_umap)))
    }
  }
  
  pca_center<-t(sapply(1:N_label, function(i){
    index_i<-which(cluster_ == label_index[i])
    colMeans(as.matrix(pca_embedding[index_i,]))
  }))
  
  
  pca_anchor_index<-lapply(1:N_label, function(i){
    index_i<-which(cluster_ == label_index[i])
    set.seed(seed.use)
    sample_index_i<-sample(index_i,min(min_size,length(index_i)) )
    sample_index_dist<-pdist(pca_embedding[sample_index_i,c(1,2)],pca_center[i,c(1,2)])@dist
    savis_nth(x=sample_index_dist,
      k = max(1,min(ceiling(min_size/5),length(index_i))))
  })
  pca_dist1<-Dist(pca_center)
  #pca_dist1<-as.matrix(parDist(pca_center))
  
  pca_dist2<-pca_dist1
  pca_dist_main<-pca_dist1[main_index,main_index]
  
  if(shrink_distance){
    remain_index<-c(1:N_label)[which(!c(1:N_label)%in%main_index)]
    #prop_<-sqrt(exp(cluster_size/max(cluster_size))/
    #    max(exp(cluster_size/max(cluster_size))))
    prop_<-sqrt(exp(cluster_size/max(cluster_size))/
        mean(exp(cluster_size/max(cluster_size))[main_index]))
    for(i in remain_index){
      pca_dist1[remain_index,i]<-pca_dist1[remain_index,i]*prop_[i]
      pca_dist1[i,remain_index]<-pca_dist1[i,remain_index]*prop_[i]
      #pca_dist1[main_index,i]<-pca_dist1[main_index,i]*prop_[i]
      #pca_dist1[i,main_index]<-pca_dist1[i,main_index]*prop_[i]
      #pca_dist1[,i]<-pca_dist1[,i]*prop_[i]
      #pca_dist1[i,]<-pca_dist1[i,]*prop_[i]
    } 
    for(i in remain_index){
      x<-main_index[which.min(pca_dist1[main_index,i])]
      min_x<-min(pca_dist1[,i])*0.5
      pca_dist1[x,i]<-min_x
      pca_dist1[i,x]<-min_x
    }
  }
  if(shrink_distance){
    pam_res<-pam(x = pca_dist1,k = 2)
    
    var_1<-sum(pca_center[which(pam_res$clustering==1),]^2)/sum(pam_res$clustering==1)
    var_2<-sum(pca_center[which(pam_res$clustering==2),]^2)/sum(pam_res$clustering==2)
    
    if(var_1 >= var_2){
      clu1<-1
      clu2<-2
    }else{
      clu1<-2
      clu2<-1
    }
    index_clu1<-which(pam_res$clustering==clu1)
    index_clu2<-which(pam_res$clustering==clu2)
    for( i in index_clu1){
      closed_dist<-min(pca_dist1[i,index_clu2])
      pca_dist1[i,-i]<-pca_dist1[i,-i] - closed_dist*shrink_factor
      if(sum(pca_dist1[i,-i]<0)>0){
        index_neg<-which(pca_dist1[i,-i]<0)
        pca_dist1[i,index_neg]<-min(pca_dist1[i,which(pca_dist1[i,]>0)])
      }
      pca_dist1[,i]<-pca_dist1[i,]
    }
  }
  
  step1_res<-get_umap_embedding_adjust_tsMDS(
    pca_embedding=pca_embedding,
    pca_center=pca_center,
    pca_anchor_index=pca_anchor_index,
    pca_dist=pca_dist1,
    pca_dist_main = pca_dist_main,
    adjust_method = adjust_method,
    main_index = main_index,
    distance_metric=distance_metric,
    umap_embedding=umap_embedding,
    N_label=N_label,
    cluster_ = cluster_,
    label_index = label_index,
    scale_factor =scale_factor,
    rotate = rotate,
    seed.use = seed.use)
  umap_center<-step1_res$umap_center
  umap_embedding_adjust<-step1_res$umap_embedding_adjust
  snn_1<- FindNeighbors(object = umap_embedding_adjust,
    nn.method = "rann",
    verbose = F)$snn
  cluster_1 <- FindClusters(snn_1,
    resolution = 0,
    verbose = F)[[1]]
  cluster_1 <- as.numeric(as.character(cluster_1))
  
  N_label1<-length(unique(cluster_1))
  if (N_label1 < N_label){
    ## First Adjustment
    label_index1<-sort(unique(cluster_1))
    
    bad_index<-list()
    list_index<-0
    for ( i in 1:N_label1){
      index_i1<-which(cluster_1 == label_index1[i])
      cur_index_len<-length(unique(cluster_[index_i1]))
      if(cur_index_len > 1){
        list_index<-list_index+1
        bad_index[[list_index]]<-c(unique(cluster_[index_i1]))
      }
    }
    for (i in 1:length(bad_index)){
      
      pos<-min(bad_index[[i]])
      other_pos<-bad_index[[i]][bad_index[[i]]>pos]
      pos<-pos+1
      other_pos<-other_pos+1
      dist_vec<-pca_dist1[,pos]
      all_pos<- dist_vec <= max(dist_vec[other_pos])
      pca_dist1[all_pos,pos]<-min(dist_vec[dist_vec>=max(dist_vec[other_pos])])
      pca_dist1[pos,all_pos]<-pca_dist1[all_pos,pos]
    }
    pca_dist<-as.dist(pca_dist1)
    pca_dist1<-as.matrix(pca_dist)
    step2_res<-get_umap_embedding_adjust_tsMDS(
      pca_embedding=pca_embedding,
      pca_center=pca_center,
      pca_anchor_index=pca_anchor_index,
      pca_dist=pca_dist1,
      adjust_method = adjust_method,
      main_index = main_index,
      distance_metric=distance_metric,
      umap_embedding=umap_embedding,
      N_label=N_label,
      cluster_ = cluster_,
      label_index = label_index,
      scale_factor =scale_factor,
      rotate = rotate,
      seed.use = seed.use)
    umap_center<-step2_res$umap_center
    umap_embedding_adjust<-step2_res$umap_embedding_adjust
    snn_2<- FindNeighbors(object = umap_embedding_adjust,
      nn.method = "rann",
      verbose = F)$snn
    cluster_2 <- FindClusters(snn_2,
      resolution = 0,
      verbose = F)[[1]]
    cluster_2 <- as.numeric(as.character(cluster_2))
    
    N_label2<-length(unique(cluster_2))
    if(N_label2 < N_label){
      ## Second Adjustment Push Away
      label_index2<-sort(unique(cluster_2))
      
      bad_index<-list()
      list_index<-0
      for ( i in 1:N_label2){
        index_i2<-which(cluster_2 == label_index2[i])
        cur_index_len<-length(unique(cluster_[index_i2]))
        if(cur_index_len > 1){
          list_index<-list_index+1
          bad_index[[list_index]]<-c(unique(cluster_[index_i2]))
        }
      }
      N_label3<-N_label - 1
      cur_iter<-0
      if(is.null(maxit_push)){
        maxit_push<-max(10,N_label/3)
      }
      while (N_label3 < N_label & cur_iter < maxit_push) {
        cur_iter<-cur_iter+1
        for (i in 1:length(bad_index)){
          pos<-min(bad_index[[i]])
          other_pos<-bad_index[[i]][bad_index[[i]]>pos]
          pos<-pos+1
          other_pos<-other_pos+1
          dist_mat<-pdist(umap_center[pos,],umap_center)@dist
          target_distance<-min(dist_mat[dist_mat>max(dist_mat[other_pos])])
          target_distance<-(target_distance+dist_mat[other_pos])/2
          for (k in 1:length(other_pos)){
            cur_pos<-other_pos[k]
            cur_center<-umap_center[pos,]
            cur_arrow<-umap_center[cur_pos,]
            prop<-target_distance[k]/dist_mat[cur_pos]
            cur_arrow<-(cur_arrow-cur_center)*prop+cur_center
            index_i<-which(cluster_ == label_index[cur_pos])
            umap_embedding_adjust[index_i,]<-t(t(umap_embedding_adjust[index_i,])-as.numeric(umap_center[cur_pos,])+as.numeric(cur_arrow))
          }
        }
        
        snn_3<- FindNeighbors(object = umap_embedding_adjust,
          nn.method = "rann",
          verbose = F)$snn
        cluster_3 <- FindClusters(snn_3,
          resolution = 0,
          verbose = F)[[1]]
        cluster_3 <- as.numeric(as.character(cluster_3))
        N_label3<-length(unique(cluster_3))
        label_index3<-sort(unique(cluster_3))
        
        bad_index<-list()
        list_index<-0
        for ( i in 1:N_label3){
          index_i3<-which(cluster_3 == label_index3[i])
          cur_index_len<-length(unique(cluster_[index_i3]))
          if(cur_index_len > 1){
            list_index<-list_index+1
            bad_index[[list_index]]<-c(unique(cluster_[index_i3]))
          }
        }
        
      }
      if(N_label3 < N_label){
        #if(verbose){
        #  cat('\n')
        #  print("Rescaling SAVIS...")
        #  setTxtProgressBar(pb = pb, value = 19.5)
        #}
        ## Third Adjustment Rescale
        for (i in 1:length(bad_index)){
          pos<-min(bad_index[[i]])
          other_pos<-bad_index[[i]][bad_index[[i]]>pos]
          pos<-pos+1
          other_pos<-other_pos+1
          dist_mat<-pdist(umap_center[pos,],umap_center)@dist
          target_distance<-min(dist_mat[dist_mat>=max(dist_mat[other_pos])])
          re_sf<-target_distance/max(dist_mat[other_pos])
          
          umap_embedding_adjust<-umap_embedding_adjust*re_sf
        }
      }
      
    }
  }
  return(umap_embedding_adjust)
}







#' adjustUMAP
#'
#' Adjust UMAP to deal with distortion 
#'
#' @details amazing umap
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom pdist pdist
#' @importFrom uwot umap
#' @importFrom MASS isoMDS
#' @importFrom stats cmdscale var as.dist dist
#'
#' @export
#' 
#' @examples
#' a<-1
#'
#'
adjustUMAP<-function(
  pca_embedding,
  umap_embedding,
  cluster_label = NULL,
  global_umap_embedding = NULL,
  adjust_method = "umap",
  distance_metric = "euclidean",
  scale_factor = 0.9,
  rotate = TRUE,
  density_adjust = TRUE,
  shrink_distance = TRUE,
  shrink_factor = 0.2,
  seed.use = 42,
  min_size = 100,
  maxit_push = NULL
){
  if(adjust_method == "all"){
    umap_adjust<-adjustUMAP(
      pca_embedding=pca_embedding,
      umap_embedding=umap_embedding,
      cluster_label = cluster_label,
      global_umap_embedding=global_umap_embedding,
      adjust_method = "umap",
      distance_metric =distance_metric,
      scale_factor = scale_factor,
      rotate = rotate,
      density_adjust = density_adjust,
      #shrink_all_distance = shrink_all_distance,
      shrink_distance = shrink_distance,
      seed.use = seed.use,
      min_size = min_size,
      maxit_push = maxit_push
    )
    tsMDS_adjust<-adjustUMAP(
      pca_embedding=pca_embedding,
      umap_embedding=umap_embedding,
      global_umap_embedding=global_umap_embedding,
      adjust_method = "tsMDS",
      distance_metric =distance_metric,
      scale_factor = scale_factor,
      rotate = rotate,
      density_adjust = density_adjust,
      #shrink_all_distance = shrink_all_distance,
      shrink_distance = shrink_distance,
      seed.use = seed.use,
      min_size = min_size,
      maxit_push = maxit_push
    )
    newList<-list("umap" = umap_adjust,
      "tsMDS"=tsMDS_adjust)
    return(newList)
  }else if (adjust_method == "umap"| 
      adjust_method == "Umap"|
      adjust_method == "UMAP"){
    umap_adjust<-adjustUMAP_via_umap(
      pca_embedding = pca_embedding,
      umap_embedding = umap_embedding,
      cluster_label = cluster_label,
      global_umap_embedding = global_umap_embedding,
      distance_metric = distance_metric,
      scale_factor = scale_factor,
      rotate = rotate,
      density_adjust = density_adjust,
      seed.use = seed.use,
      min_size = min_size,
      maxit_push = maxit_push)
    return(umap_adjust)
  }else if(adjust_method == "MDS"|
      adjust_method == "Mds"|
      adjust_method == "mds"|
      adjust_method == "tsMDS"|
      adjust_method == "tsmds"|
      adjust_method == "TSMDS"){
    tsMDS_adjust<-adjustUMAP_via_tsMDS(
      pca_embedding = pca_embedding,
      umap_embedding = umap_embedding,
      global_umap_embedding = global_umap_embedding,
      distance_metric = distance_metric,
      scale_factor = scale_factor,
      rotate = rotate,
      density_adjust = density_adjust,
      shrink_distance = shrink_distance,
      adjust_method = "tsMDS",
      shrink_factor = shrink_factor,
      seed.use = seed.use,
      min_size = min_size,
      maxit_push = maxit_push)
    return(tsMDS_adjust)
  }else if (adjust_method == "isoMDS"|
      adjust_method == "ISOMDS"|
      adjust_method == "isomds"){
    tsMDS_adjust<-adjustUMAP_via_tsMDS(
      pca_embedding = pca_embedding,
      umap_embedding = umap_embedding,
      global_umap_embedding = global_umap_embedding,
      distance_metric = distance_metric,
      scale_factor = scale_factor,
      rotate = rotate,
      density_adjust = density_adjust,
      shrink_distance = shrink_distance,
      adjust_method = "isoMDS",
      shrink_factor = shrink_factor,
      seed.use = seed.use,
      min_size = min_size,
      maxit_push = maxit_push)
    return(tsMDS_adjust)
  }else{
    print("Wrong adjust method!")
    return(0)
  }
}









#' ScaleFactor
#'
#' Combined PC embedding with scale factor for subPC
#'
#' @details Adapts source code from: https://github.com/rstudio/reticulate/blob/master/R/source.R
#' @details Now source_python_ supports loading from http
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom pdist pdist
#' @importFrom stats dist as.dist
#'
#'
#' @examples
#' a<-1
#'
#'
ScaleFactor<-function(
  combined_embedding,
  cluster_label,
  npcs=NULL,
  center_method = "mean",
  scale_factor_separation =3,
  return_factor=TRUE){
  
  # sorted unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  if (is.null(npcs)){
    npcs<-ncol(combined_embedding)/2
    if (npcs %% 1 != 0){
      stop("combined_embedding combines global PC
        and local PC, column number of combined
        embedding should be even.")
    }
  }else{
    if (2*npcs != ncol(combined_embedding)){
      stop("column number of combined embedding
        is not equal to 2*npcs")
    }
  }
  
  if (center_method == "geometry"){
    
    cluster_center<-t(sapply(1:N_label, function(i){
      index_i<-which(cluster_label == label_index[i])
      tmp<-sapply(1:(2*npcs), function(j){
        (max(combined_embedding[index_i,j]) +
            min(combined_embedding[index_i,j]))/2
      })
      names(tmp)<-c(paste0("PC_",1:npcs),
        paste0("subPC_",1:npcs))
      tmp
    }))
    
  }else if (center_method == "mean"){
    # cluster center by mean
    cluster_center<-t(sapply(1:N_label, function(i){
      index_i<-which(cluster_label == label_index[i])
      colMeans(combined_embedding[index_i,1:(2*npcs)])
    }))
    colnames(cluster_center)<-c(paste0("PC_",1:npcs),
      paste0("subPC_",1:npcs))
    
  }else{
    stop("center_method is chosen from geometry or mean")
  }
  
  # distance matrix for cluster center
  cluster_center_dist<-Dist(cluster_center[,1:npcs])
  #cluster_center_dist<-as.matrix(parDist(cluster_center[,1:npcs]))
  diag(cluster_center_dist)<-NA
  
  # pdist is used here, which is fast
  cluster_semi_d<-sapply(1:N_label, function(i){
    index_i<-which(cluster_label == label_index[i])
    max((pdist(cluster_center[i,(npcs+1):(2*npcs)],
      combined_embedding[index_i,(npcs+1):(2*npcs)]))@dist)
  })
  
  scale_factor<-sapply(1:N_label, function(i){
    min(cluster_center_dist[i,],na.rm = T)/scale_factor_separation/cluster_semi_d[i]
  })
  
  for ( i in 1:N_label){
    index_i<-which(cluster_label == label_index[i])
    combined_embedding[index_i,(1+npcs):(2*npcs)]<-
      combined_embedding[index_i,(1+npcs):(2*npcs)]*scale_factor[i]
  }
  if(return_factor){
    newList<-list("combined_embedding"=combined_embedding,
      "scale_factor"=scale_factor)
    return(newList)
  }
  combined_embedding
}

#' DoCluster
#'
#' Do cluster for PC embedding
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom Spectrum Spectrum
#'
#'
#' @examples
#' a<-1
#'
#'
DoCluster<-function(
  pc_embedding,
  method = "louvain",
  resolution = 0.5,
  subcluster = FALSE,
  fixk = NULL,
  seed.use = 42,
  verbose = FALSE,
  verbose_more = FALSE){
  if (verbose_more){
    verbose <- TRUE
  }
  set.seed(seed.use)
  if (method == "louvain"){
    if (verbose){
      print("Finding Neighbors...")
    }
    snn_<- FindNeighbors(object = pc_embedding,
      nn.method = "rann",
      verbose = verbose_more)$snn
    if (verbose){
      print("Finding Clusters...")
    }
    cluster_ <- FindClusters(snn_,
      resolution = resolution,
      verbose = verbose_more)[[1]]
    cluster_ <- as.numeric(as.character(cluster_))
    if (subcluster){
      
      N_cluster<-length(unique(cluster_))
      subcluster_<-rep(NA,length(cluster_))
      for ( i in unique(cluster_)){
        if (verbose){
          print(paste0("Processing cluster ",i))
        }
        index_i<-which(cluster_ == i)
        if (verbose_more){
          print("Finding Neighbors...")
        }
        snn_<- FindNeighbors(object = pc_embedding[index_i,],
          nn.method = "rann",
          verbose = verbose_more)$snn
        if (verbose_more){
          print("Finding Clusters...")
        }
        subcluster_label <- FindClusters(snn_,
          resolution = resolution,
          verbose = verbose_more)[[1]]
        subcluster_[index_i]<-paste0(i,subcluster_label)
      }
      subcluster_ <- as.numeric(as.character(subcluster_))
    }
  }
  if (method == "spectral"){
    num_example<-nrow(pc_embedding)
    N_center<-min(1000,num_example%/%5)
    set.seed(42)
    suppressMessages(
      spectral_eg<-Spectrum(t(
        as.matrix(pc_embedding)),
        method = 1,showres = verbose_more,
        FASP = T,FASPk = N_center))
    cluster_<-spectral_eg$allsample_assignments
    cluster_<-as.numeric(as.character(cluster_))
    if (subcluster){
      N_cluster<-length(unique(cluster_))
      subcluster_<-rep(NA,length(cluster_))
      for ( i in unique(cluster_)){
        #print(paste0("Processing cluster ",i))
        index_i<-which(cluster_ == i)
        if (length(index_i) > min(100,N_center)){
          suppressMessages(
            spectral_eg<-Spectrum(t(
              as.matrix(pc_embedding[index_i,])),
              method = 1,FASP = T,
              FASPk = N_center,
              showres = F))
        }else{
          suppressMessages(
            spectral_eg<-Spectrum(t(
              as.matrix(pc_embedding[index_i,])),
              method = 1,FASP = T,
              FASPk = length(index_i)%/%3,
              showres = F))
        }
        subcluster_[index_i]<-paste0(i,
          spectral_eg$allsample_assignments)
      }
      subcluster_<-as.numeric(as.character(subcluster_))
    }
  }
  if (subcluster){
    newList<-list("cluster" = cluster_,
      "subcluster" = subcluster_)
    return(newList)
  }else{
    newList<-list("cluster" = cluster_)
    return(newList)
  }
}

#' SubPCEmbedding
#'
#' Based on clustering, calculate the subPC embedding for each cluster
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Seurat FindVariableFeatures ScaleData RunPCA
#'
#'
#' @examples
#' a<-1
#'
#'
SubPCEmbedding<-function(
  expr_matrix,
  cluster_label,
  assay_for_var_features = "normalizedcount",
  npcs=20,
  nfeatures =2000,
  return_hvg=FALSE,
  verbose = FALSE
  #  process_min_size = 0
){
  # get unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster
  N_label<-length(label_index)
  
  # get PC for each cluster
  expr_cluster_pc<-list()
  hvg_list<-list()
  hvg_name_list<-list()
  for ( i in c(1:N_label)){
    index_i<-which(cluster_label == label_index[i])
    #if(length(index_i) > process_min_size){
    if (verbose){
      print(paste0("Process cluster ",
        label_index[i]," Yes",
        "(size=",
        length(index_i),")"))
    }
    # expr_matrix is initial expression matrix (gene*cell)
    expr_tmp<-expr_matrix[,index_i]
    if(assay_for_var_features == "rawcount"){
      expr_tmp_hvg <- FindVariableFeatures(
        expr_tmp,verbose = F)$vst.variance.standardized
    }else if(assay_for_var_features == "normalizedcount"){
      expr_tmp_hvg <- ExpFindVariableFeatures(
        expr_tmp,verbose = F)
    }else{
      stop("Please select assay_for_var_features from c('rawcount','normalizedcount')")
    }
    tmp_hvg<-savis_nth(x = expr_tmp_hvg,
      k = nfeatures)
    if(assay_for_var_features == "rawcount"){
      expr_tmp<- NormalizeData(expr_tmp,verbose = F)
    }
    hvg_list[[i]]<-tmp_hvg
    hvg_name_list[[i]]<-rownames(expr_matrix)[tmp_hvg]
    expr_tmp<-expr_tmp[tmp_hvg,]
    expr_tmp <-ScaleData(expr_tmp,verbose = F)
    rm(expr_tmp_hvg,tmp_hvg,index_i)
    
    suppressWarnings(
      expr_tmp <- RunPCA(
        object = expr_tmp,
        features = rownames(expr_matrix),
        npcs = npcs,
        verbose = F)@cell.embeddings)
    
    # Get subgroup PC with scale factor(default to be 1)
    expr_cluster_pc[[i]]<-expr_tmp[,c(1:npcs)]
    
    
    #}
    #else{
    #  if (verbose){
    #    print(paste0("Process cluster ",
    #      label_index[i]," No",
    #      "(size=",
    #      length(index_i),")"))
    #  }
    
    # if the cluster size less than process_min_size,
    # just fill in NA
    #  expr_cluster_pc[[i]]<-NA
    #expr_cluster_pc[[i]]<-matrix(0,
    # nrow = length(index_i),
    # ncol = npcs)
    #}
    #####
  }
  if (return_hvg){
    newList<-list("sub_PC_list" = expr_cluster_pc,
      "hvg_index"=hvg_list,
      "hvg_name"=hvg_name_list)
    return(newList)
  }else{
    return(expr_cluster_pc)
  }
}


#' CombinePC
#'
#' Combine PC and subPCs
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#'
#'
#' @examples
#' a<-1
#'
#'
CombinePC<-function(PC_embedding,
  cluster_label,
  sub_PC_list,
  scale_factor=NULL
){
  label_index<-sort(as.numeric(unique(
    as.character(cluster_label))))
  N_label<-length(label_index)
  if (N_label != length(sub_PC_list)){
    stop("cluster number is not equal to sub PC list length")
  }
  npcs<-ncol(sub_PC_list[[1]])
  PC_embedding_combined<-PC_embedding
  
  sub_PC_list_supp<-matrix(0,nrow = nrow(PC_embedding),
    ncol = npcs)
  sub_PC_list_supp<-data.frame(sub_PC_list_supp)
  rownames(sub_PC_list_supp)<-rownames(PC_embedding)
  
  
  for (i in c(1:N_label)){
    index_i<-which(cluster_label == label_index[i])
    if (is.null(scale_factor[1])){
      if (!is.na(sub_PC_list[[i]][1])){
        sub_PC_list_supp[index_i,]<-sub_PC_list[[i]]
      }else{
        sub_PC_list_supp[index_i,]<-0
        #PC_embedding[index_i,]
      }
    }else{
      if (!is.na(sub_PC_list[[i]][1])){
        sub_PC_list_supp[index_i,]<-sub_PC_list[[i]]*scale_factor[i]
      }else{
        sub_PC_list_supp[index_i,]<-0
        #PC_embedding[index_i,]
      }
    }
    
  }
  
  PC_embedding_combined<-cbind(PC_embedding_combined,
    sub_PC_list_supp)
  
  return(PC_embedding_combined)
}


#' AdaptiveCombine
#'
#' Adaptively Combine PC and subPCs
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#'
#'
#' @examples
#' a<-1
#'
#'
#'
AdaptiveCombine<-function(expr_matrix,
  combined_embedding,
  cluster_label,
  cluster_label_i,
  npcs=20,
  nfeatures=2000,
  scale_factor_separation = 3,
  process_min_size=0
){
  expr_matrix_pca<-combined_embedding[,(1:(npcs))]
  combined_embedding_sub<-combined_embedding[,((npcs+1):(2*npcs))]
  # sorted unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  
  sub_PC_list_i<-lapply(1:N_label, function(i){
    
    index_i<-which(cluster_label == label_index[i])
    if(length(unique(cluster_label_i[[i]]))>1 &
        length(index_i) >= process_min_size){
      sub_PC_list_i<-SubPCEmbedding(
        expr_matrix = expr_matrix[,index_i],
        cluster_label = cluster_label_i[[i]],
        npcs = npcs,
        nfeatures = nfeatures)
      
      combined_embedding_i<-CombinePC(
        PC_embedding = combined_embedding_sub[index_i,],
        cluster_label = cluster_label_i[[i]],
        sub_PC_list = sub_PC_list_i)
      
      combined_embedding_i<-ScaleFactor(
        combined_embedding = combined_embedding_i,
        cluster_label = cluster_label_i[[i]],
        npcs = npcs,
        center_method = "mean",
        scale_factor_separation=scale_factor_separation)$combined_embedding
      combined_embedding_i<-as.matrix(combined_embedding_i)
      colnames(combined_embedding_i)<-c(paste0("subPC",1:npcs),
        paste0("subsubPC",1:npcs))
    }else{
      combined_embedding_i<-matrix(0,
        nrow = length(index_i),
        ncol = npcs)
      combined_embedding_i<-cbind(combined_embedding_sub[index_i,],
        combined_embedding_i)
      combined_embedding_i<-as.matrix(combined_embedding_i)
      colnames(combined_embedding_i)<-c(paste0("subPC",1:npcs),
        paste0("subsubPC",1:npcs))
    }
    combined_embedding_i
  })
  
  combined_embedding<-CombinePC(PC_embedding = expr_matrix_pca,
    cluster_label = cluster_label,
    sub_PC_list = sub_PC_list_i)
  
  
  return(combined_embedding)
}

#' FormCombinedEmbedding
#'
#' FormCombinedEmbedding
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#'
#'
#' @examples
#' a<-1
#'
#'
#'
FormCombinedEmbedding<-function(
  expr_matrix,
  expr_matrix_pca,
  cluster_label,
  assay_for_var_features = "normalizedcount",
  npcs=20,
  nfeatures =2000,
  center_method = "mean",
  scale_factor_separation=3
){
  sub_PC_list<-SubPCEmbedding(
    expr_matrix = expr_matrix,
    cluster_label = cluster_label,
    npcs = npcs,
    nfeatures = nfeatures,
    assay_for_var_features = assay_for_var_features)
  combined_embedding<-CombinePC(
    PC_embedding = expr_matrix_pca,
    cluster_label = cluster_label,
    sub_PC_list = sub_PC_list)
  combined_embedding<-ScaleFactor(
    combined_embedding = combined_embedding,
    cluster_label = cluster_label,
    npcs = npcs,
    center_method = center_method,
    scale_factor_separation = scale_factor_separation)$combined_embedding
  combined_embedding<-as.matrix(combined_embedding)
  return(combined_embedding)
}


#' FormAdaptiveCombineList
#'
#' FormAdaptiveCombineList
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom stats t.test
#'
#' @examples
#' a<-1
#'
#'
#'
#'
FormAdaptiveCombineList<-function(
  expr_matrix,
  expr_matrix_pca,
  max_stratification,
  stratification_count,
  scale_factor_separation,
  resolution,
  cluster_method,
  npcs,
  nfeatures,
  process_min_size,
  assay_for_var_features = "normalizedcount",
  differentail_gene_cutoff = 20,
  do_cluster = TRUE,
  cluster_label = NULL,
  check_differential =TRUE,
  verbose = FALSE
){
  colnames(expr_matrix_pca)<-paste0("Layer",stratification_count,"PC",1:npcs)
  if(nrow(expr_matrix_pca) < process_min_size){
    newList<-list("cluster_label"= -1,
      "combined_embedding"=expr_matrix_pca)
    return(newList)
  }
  if(nrow(expr_matrix_pca)!=ncol(expr_matrix)){
    stop("expr_matrix_pca and expr_matrix do not match")
  }
  N_gene<-nrow(expr_matrix)
  if(is.null(cluster_label[1])){
    do_cluster<-TRUE
  }
  if(do_cluster){
    cluster_label<-DoCluster(
      pc_embedding = expr_matrix_pca,
      method = cluster_method,
      resolution = resolution)$cluster
  }
  # sorted unique cluster label
  label_index<-sort(as.numeric(
    unique(as.character(cluster_label))))
  # number of cluster label
  N_label<-length(label_index)
  
  ## Compare npcs with the cluster size
  size_cluster<-c()
  for ( i in 1:N_label){
    size_cluster<-c(size_cluster,
      sum(cluster_label == label_index[i]))
  }
  if(sum(size_cluster <= npcs) > 0){
    warning("Produce a cluster whose size is less than npcs: Combine it with other cluster, Or please adjust the npcs or resolution")
    cur_index<-which(size_cluster <= npcs)
    pca_center<-t(sapply(1:N_label, function(i){
      index_i<-which(cluster_label == label_index[i])
      colMeans(expr_matrix_pca[index_i,])
    }))
    sample_index_dist<-pdist(pca_center[-cur_index,],pca_center[cur_index,])@dist
    sample_index_dist<-matrix(sample_index_dist,nrow = sum(size_cluster <= npcs))
    comb_list<-sapply( 1:length(cur_index), function(i){
      #c(cur_index[i],nth(x=sample_index_dist[i,],
      #  k = 1,
      #  num.of.nths = 2,
      #  descending = F,
      #  index.return = T)[,1]) 
      c(cur_index[i],which.min(sample_index_dist[i,]))
    })
    for( i in 1:ncol(comb_list)){
      index_1<-which(cluster_label == label_index[comb_list[1,i]])
      cluster_label[index_1]<-label_index[comb_list[2,i]]
    }
    # update label_index
    # sorted unique cluster label
    label_index<-sort(as.numeric(
      unique(as.character(cluster_label))))
    # number of cluster label
    N_label<-length(label_index)
    
    ## Compare npcs with the cluster size
    size_cluster<-c()
    for ( i in 1:N_label){
      size_cluster<-c(size_cluster,
        sum(cluster_label == label_index[i]))
    }
  }
  
  if(check_differential){
    if(N_label == 1){
      newList<-list("cluster_label"= -1,
        "combined_embedding"=expr_matrix_pca)
      return(newList)
    }else{
      ### Begin of else 
      cur_label<- N_label
      while(cur_label>1){
        for (i in 1:(cur_label-1)){
          print(paste0("Current:",label_index[i],"vs",label_index[cur_label]))
          index_1<-which(cluster_label == label_index[i])
          index_2<-which(cluster_label == label_index[cur_label])
          #S<-cbind(expr_matrix[,index_1],expr_matrix[,index_2])
          #S<-CreateSeuratObject(S)
          #label_diff<-c(rep(0,length(index_1)),rep(1,length(index_2)))
          #Idents(S)<-factor(label_diff)
          #.<-capture.output(marker_diff<-FindMarkers(S,
          #  test.use = "t",
          #  ident.1 = unique(S@active.ident)[1],
          # ident.2 =  unique(S@active.ident)[2]))
          #print(sum(marker_diff[,5]<0.05))
          #if (sum(marker_diff[,5]<=0.05) < 5){
          #print(paste0("Don't find differential genes between",label_index[i],"vs",label_index[cur_label]))
          #  cluster_label[index_2]<-label_index[i]
          #  label_index<-sort(as.numeric(
          #    unique(as.character(cluster_label))))
          #  N_label<-length(label_index)
          #  break
          #}
          cur_differential_gene<-0
          cur_gene<-1
          while (cur_differential_gene < differentail_gene_cutoff & 
              cur_gene <= N_gene){
            t_res<-t.test(expr_matrix[cur_gene,index_1],
              expr_matrix[cur_gene,index_2])
            if (is.na(t_res$p.value)){
              t_res$p.value<-1
            }
            if (t_res$p.value <0.05){
              cur_differential_gene <- cur_differential_gene+1
            }
            cur_gene<-cur_gene+1
          }
          if(cur_differential_gene < differentail_gene_cutoff){
            print(paste0("Don't find differential genes between",label_index[i],"vs",label_index[cur_label]))
            cluster_label[index_2]<-label_index[i]
            label_index<-sort(as.numeric(
              unique(as.character(cluster_label))))
            N_label<-length(label_index)
            break
          }
        }
        cur_label<-cur_label - 1
        
      }
      ### End of else 
    }
    
    
  }
  
  if(N_label == 1){
    newList<-list("cluster_label"= -1,
      "combined_embedding"=expr_matrix_pca)
    return(newList)
  }
  
  combined_embedding<-FormCombinedEmbedding(
    expr_matrix=expr_matrix,
    expr_matrix_pca=expr_matrix_pca,
    cluster_label=cluster_label,
    assay_for_var_features = assay_for_var_features,
    npcs=npcs,
    nfeatures=nfeatures,
    scale_factor_separation=scale_factor_separation
  )
  if(max_stratification == stratification_count + 1){
    newList<-list("cluster_label"= cluster_label,
      "combined_embedding"=combined_embedding)
    return(newList)
  }  
  colnames(combined_embedding)<-c(paste0("Layer",stratification_count,"PC",1:npcs),
    paste0("Layer",(stratification_count+1),"PC",1:npcs))
  
  cluster_label_list<-list()
  combined_embedding_list<-list()
  work_adaptive<-FALSE
  max_ncol_sub<-0
  for ( i in 1:N_label){
    index_i<-which(cluster_label == label_index[i])
    tmp<-FormAdaptiveCombineList(
      expr_matrix = expr_matrix[,index_i],
      expr_matrix_pca = combined_embedding[index_i,(npcs+1):(2*npcs)],
      assay_for_var_features = assay_for_var_features,
      max_stratification = max_stratification,
      stratification_count = stratification_count + 1,
      scale_factor_separation = scale_factor_separation,
      resolution = resolution,
      cluster_method = cluster_method,
      npcs = npcs,
      nfeatures = nfeatures,
      process_min_size = process_min_size,
      check_differential = check_differential,
      verbose = verbose
    )
    cluster_label_list[[i]]<-tmp$cluster_label
    combined_embedding_list[[i]]<-tmp$combined_embedding
    max_ncol_sub<-max(max_ncol_sub,ncol(tmp$combined_embedding))
    if (sum(tmp$cluster_label!= -1) != 0){
      work_adaptive<-TRUE
    }else{
      if(ncol(tmp$combined_embedding)!=npcs){
        #print(ncol(tmp$combined_embedding))
        stop("combined_embedding size is wrong")
      }  
    }
    rm(tmp)
  }
  if (work_adaptive){
    
    if(max_ncol_sub < (2*npcs)){
      stop("combined_embedding size is wrong1")
    }else if (max_ncol_sub == (2*npcs)){
      ## cluster label list is vector
      cluster_label_sub<-rep(-1,length(cluster_label))
      combined_embedding_sub<-matrix(0,
        nrow = nrow(expr_matrix_pca),
        ncol = max_ncol_sub)
      colnames(combined_embedding_sub)<-c(paste0("Layer",(stratification_count+1),"PC",1:npcs),
        paste0("Layer",(stratification_count+2),"PC",1:npcs))
      for ( i in 1:N_label){
        index_i<-which(cluster_label == label_index[i])
        cluster_label_sub[index_i]<-cluster_label_list[[i]]
        
        combined_embedding_sub[index_i,
          1:ncol(combined_embedding_list[[i]])]<-
          combined_embedding_list[[i]]
      }
      newList<-list("cluster_label"= cbind(cluster_label,cluster_label_sub),
        "combined_embedding"=cbind(expr_matrix_pca,combined_embedding_sub))
      return(newList)
      
    }else{
      ## cluster label list is matrix
      max_label_count<-0
      for ( i in 1:N_label){
        if(is.null(dim(cluster_label_list[[i]]))){
          max_label_count<-max(max_label_count,1 )
        }else{
          max_label_count<-max(max_label_count,ncol(cluster_label_list[[i]]))
        }
      }
      cluster_label_sub<-matrix(-1,
        nrow = length(cluster_label),
        ncol = max_label_count)
      combined_embedding_sub<-matrix(0,
        nrow = nrow(expr_matrix_pca),
        ncol = max_ncol_sub)
      
      for ( i in 1:N_label){
        index_i<-which(cluster_label == label_index[i])
        if(is.null(dim(cluster_label_list[[i]]))){
          cluster_label_sub[index_i,1]<-cluster_label_list[[i]]
        }else{
          cluster_label_sub[index_i,1:ncol(cluster_label_list[[i]])
            ]<-cluster_label_list[[i]]
        }
        combined_embedding_sub[index_i,
          1:ncol(combined_embedding_list[[i]])]<-
          combined_embedding_list[[i]]
        if(ncol(combined_embedding_list[[i]]) == max_ncol_sub){
          colnames(combined_embedding_sub)<-colnames(combined_embedding_list[[i]])
        }
      }
      newList<-list("cluster_label"= cbind(cluster_label,cluster_label_sub),
        "combined_embedding"=cbind(expr_matrix_pca,combined_embedding_sub))
      return(newList)
    }
  }else{
    
    newList<-list("cluster_label"= cluster_label,
      "combined_embedding"=combined_embedding)
    return(newList)
  }
}




#' adaDimPlot
#'
#' Adaptively Plot the UMAP Embedding
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom dplyr group_by summarise
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats median
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
#'
adaDimPlot<-function(
  umap_embedding,
  label,
  pt.size=0,
  show.legend=TRUE,
  seed.use = 42,
  color.mode = 1
){
  set.seed(seed.use)
  shuffle_index<-sample(1:nrow(umap_embedding))
  umap_embedding<-umap_embedding[shuffle_index,]
  label<-label[shuffle_index]
  
  umap_embedding<-data.frame(umap_embedding)
  if(is.null(colnames(umap_embedding))[1]){
    colnames(umap_embedding)<-paste0("UMAP_",1:ncol(umap_embedding)) 
    xynames<-colnames(umap_embedding)
  }else{
    xynames<-colnames(umap_embedding)
    colnames(umap_embedding)<-paste0("UMAP_",1:ncol(umap_embedding)) 
  }
  umap_embedding$label<-factor(label)
  
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal,
    qual_col_pals$maxcolors,rownames(qual_col_pals)))
  if(color.mode==1){
    col_vector <- unique(col_vector)
    col_vector[4]<-"#ffd000" 
  }
  cpnum<-length(unique(label))%/%length(col_vector)+1
  col_vector<-rep(col_vector,cpnum)
  gg<-ggplot(umap_embedding)+
    geom_point(aes(x = UMAP_1,
      y = UMAP_2,
      color = label),
      size = pt.size)+
    theme(legend.title = element_blank())+
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.key=element_blank())+
    guides(color = guide_legend(override.aes =
        list(size=3)))+
    labs(x = xynames[1],y=xynames[2])+
    scale_colour_manual(values =
        col_vector[c(1:length(unique(label)))])
  if (!show.legend){
    gg<-gg+theme(legend.position="none")
  }
  
  centers<-group_by(umap_embedding,label)
  centers<-summarise(centers,x = median(x = UMAP_1),
    y = median(x = UMAP_2),.groups = 'drop')
  
  gg <- gg +
    geom_text_repel(data = centers,
      mapping = aes(x = x, y = y,
        label = label), size = 4,max.overlaps = 100)
  gg
}


#' adaDimPlot2
#'
#' Adaptively Plot the UMAP Embedding
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom dplyr group_by summarise
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats median
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
#'
adaDimPlot2<-function(
  umap_embedding,
  label,
  pt.size=0
){
  set.seed(42)
  shuffle_index<-sample(1:nrow(umap_embedding))
  umap_embedding<-umap_embedding[shuffle_index,]
  label<-label[shuffle_index]
  umap_embedding<-data.frame(umap_embedding)
  if(is.null(colnames(umap_embedding))[1]){
    colnames(umap_embedding)<-paste0("UMAP_",1:ncol(umap_embedding)) 
    xynames<-colnames(umap_embedding)
  }else{
    xynames<-colnames(umap_embedding)
    colnames(umap_embedding)<-paste0("UMAP_",1:ncol(umap_embedding)) 
  }
  umap_embedding$label<-label
  
  gg<-ggplot(umap_embedding)+
    geom_point(aes(x = UMAP_1,
      y = UMAP_2,
      color = label),
      size = pt.size)+
    theme(legend.title = element_blank())+
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.key=element_blank())+
    #scale_colour_gradientn(colors=rainbow(15)[c(12:1,14,15)])+
    scale_colour_gradientn(colors=c("grey", "red","black"))+
  labs(x = xynames[1],y=xynames[2])
  gg
}



#' ARIEvaluate
#'
#' Adaptively Plot the UMAP Embedding
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom Spectrum Spectrum
#' @importFrom stats kmeans
#' @importFrom mclust adjustedRandIndex
#'
#'
#' @examples
#' a<-1
#'
#'
#'
#'
#'
ARIEvaluate<-function(
  test_data,
  label,
  method = "louvain",
  resolution = 0.5,
  fixk = NULL
){
  if (method == "louvain"){
    snn_<- FindNeighbors(object = test_data,
      nn.method = "rann",
      verbose = F)$snn
    cluster_label <- FindClusters(snn_,
      resolution = resolution,
      verbose = F)[[1]]
  }
  if (method == "spectral"){
    
    set.seed(42)
    N_center<-nrow(test_data)
    if (N_center > 5000){
      N_center <- 1000
    }else{
      N_center<-ceiling(nrow(test_data)/5)
    }
    if (is.null(fixk)){
      suppressMessages(
        spectrual_eg<-Spectrum(t(as.matrix(
          test_data)),
          method = 1,showres = F,
          FASP = T,FASPk = N_center))
    }else{
      suppressMessages(
        spectrual_eg<-Spectrum(t(as.matrix(
          test_data)),
          method = 3,showres = F,
          FASP = T,FASPk = N_center,fixk = fixk))
    }
    
    cluster_label<-spectrual_eg$allsample_assignments
  }
  if (method == "kmeans"){
    res<-kmeans(x = as.matrix(
      test_data),centers = fixk,iter.max = 100)
    cluster_label<-res$cluster
  }
  
  
  c(length(unique(cluster_label)),
    adjustedRandIndex(cluster_label,label))
}

#' SeuratLPCA
#'
#' Get PC matrix of single cell RNAseq matrix using limited memory with Seurat pipline
#'
#' @details Some details
#' @param expr_matrix character
#'
#' @return nothing useful
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#'
#' @export
#'
#' @examples
#' a<-1
#'
#'
#'
SeuratLPCA<-function(expr_matrix,assay_for_var_features = "rawcount",npcs=20,nfeatures=2000){
  if ( assay_for_var_features == "rawcount" ){
    expr_matrix_hvg <- FindVariableFeatures(
      expr_matrix,
      verbose = F)$vst.variance.standardized
    expr_matrix<-NormalizeData(
      expr_matrix,
      verbose = F)
  }else if ( assay_for_var_features == "normalizedcount" ){
    expr_matrix<-NormalizeData(
      expr_matrix,
      verbose = F)
    expr_matrix_hvg <- ExpFindVariableFeatures(
      expr_matrix,
      verbose = F)
  }else{
    stop("Please select assay_for_var_features from c('rawcount','normalizedcount')")
  }
  hvg<-savis_nth(x = expr_matrix_hvg,
    k = nfeatures)
  
  expr_matrix<-expr_matrix[hvg,]
  expr_matrix <- ScaleData(
    expr_matrix,
    verbose = F)
  suppressWarnings(expr_matrix_pca <- RunPCA(
    object = expr_matrix,
    features = rownames(expr_matrix),
    npcs = npcs,
    verbose = F)@cell.embeddings)
  rm(expr_matrix)
  expr_matrix_pca<-data.frame(expr_matrix_pca)
  expr_matrix_pca<-as.matrix(expr_matrix_pca)
  return(expr_matrix_pca)
}
