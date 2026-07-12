#' @title Neighbour composition matrix
#' @name nnMat
#' @description
#' From a metadata file, given a list of cells with the corrisponding annotated cell type, build a neighbor matrix reporting the number of occurrences for each cell type
#'
#' @param data a Seurat object containing Kandinsky data
#' @param label character string specifying the variable name to be used to define cell annotation groups
#' @param method character string specifying the method to be used to identify neighbours.
#' Must be one of the following:
#' 'Q': queen contiguity method,check for contact (not overlap) between any edge or side od two polygons (refers to the queen movement rule in chess). Currently only applicable for Visium/Visium-HD data;
#' 'C': centroid-based method, use maximum centroid distance threshold to identify spot/cell neighbours;
#' 'D': Delaunay triangulation method;
#' 'K': KNN method, define k closest neighbours to each spot/cell;
#' 'M': membrane-based method, check for the occurrence of a physical contact/intersection within a distance threshold between cell boundaries. Not applicable in the case of Visium spots.
#' When argument `method` is not specified or is set to `NULL`, the function will use `nb` slot already stored in Kandinsky object to create neighbour matrix
#' @param snap numeric, extra distance accepted between polygon borders for contiguity relation. Applied when 'method' is set to `Q`. Can't work with point geometries like cell centroids.
#' @param layers numeric, number of concentric contiguous layers to include in spot neighbourhood. Only Applied when 'nb.method = Q'. Default is 1.
#' @param k numeric, number of nearest neighbours to be searched. Applied when `method` is set to `K`.
#' @param d.max numeric, maximum centroid or membrane distance threshold. Applied when `method` is set to `C` or `M`.
#' @param soi boolean, whether or not filter Delaunay network to keep sphere of influence (SOI) graph. Default is FALSE.
#' @param ids_anno optional, character string specifying any extra variable to be added to the final distance matrix
#' @param return.seurat boolean, whether returning neighbour matrix alone (FALSE) or the input Seurat object with a new `nnMat` slot as part of the Kandinsky data
#' @family nnMat
#' @importFrom magrittr %>%
#' @import spatialreg
#' @export
nnMat = function(data,label=NULL,method=c("Q", "C", "D", "K", "M"),layers=1,snap=NULL,k=NULL,d.max=NULL,soi=F,ids_anno=NULL,return.seurat=T){
  data@meta.data$ids = rownames(data@meta.data)
  poly = sf::st_drop_geometry(populate_sf(data,vars=c('ids',label,ids_anno),return.seurat = F))
  if(length(method) >1){
    method = NULL
  }
  if(is.null(method)){
    message('Creating neighbour matrix using pre-defined neighbour network')
    mat = as(KanData(data,'nb'),'CsparseMatrix')
    nb.method = KanData(data,'nb.type')
  }else{
    warning('Creating neighbour matrix using new neighbour criteria')
    if(method == 'Q'){
      if(KanData(data,'platform') %in% c('visium','visium_hd')){
        mat = queen_nb(KanData(data,'sf'),snap=snap,layers=layers)
        nb.method = 'Q'
      }else{
        stop('method "Q" can only be applied to spot-based technologies (Visium/Visium-HD)')
      }
    }
    if(method == 'K'){
      mat = knn_nb(KanData(data,'sf'),k=k)
      nb.method = paste0('K_',k)
    }
    if(method == 'C'){
      mat = centroid_nb(KanData(data,'sf'),d.max=d.max)
      nb.method = paste0('C_',d.max)
    }
    if(method == 'M'){
      mat = membrane_nb(KanData(data,'sf'),d.max=d.max)
      nb.method = paste0('M_',d.max)
    }
    if(method == 'D'){
	mat = tri_nb(KanData(data,'sf'),soi=soi)
    	nb.method = paste0('D_',paste0('soi',soi))
    }
    message('Neighbour method: ',nb.method)
    mat = as(mat,'CsparseMatrix')
  }
  labels = as.factor(as.character(poly[[label]]))
  Z=Matrix::sparse.model.matrix(~labels-1)
  colnames(Z) = gsub('labels','',colnames(Z))
  neigh_mat = mat %*% Z
  neigh_mat = as.data.frame(neigh_mat)
  neigh_mat$tot_nn = rowSums(neigh_mat)
  ids_anno = ids_anno %||% label
  neigh_mat[[ids_anno]]=poly[[ids_anno]]

#  neigh = Matrix::which(mat!=0,arr.ind=TRUE)  %>% as.data.frame()
#  #Add cell ids and neighboring cell types to neigh object
#  neigh$cell_ID = poly[neigh$row,][['ids']]
#  neigh$col.type = as.character(poly[neigh$col,][[colnames(poly)[colnames(poly) == label]]])
#  #Count the occurrences for each cell type among the neighbors of each cell
#  neigh_mat = neigh %>%
#    dplyr::group_by(.data$cell_ID,.data$col.type) %>%
#    dplyr::mutate(n = dplyr::row_number()) %>%
#    dplyr::ungroup() %>%
#    dplyr::select(c(.data$cell_ID,.data$col.type,.data$n))
#  #From tabular to matrix format
#  neigh_mat = reshape2::dcast(neigh_mat,cell_ID~col.type,value.var='n',fun.aggregate = length)
#  rownames(neigh_mat) = neigh_mat$cell_ID
#  neigh_mat$cell_ID = NULL
#  #Be sure that all column in the final matrix are in numeric format
#  neigh_mat = neigh_mat %>% dplyr::mutate_all(as.numeric)
#  #Convert any possible NA value to 0
#  neigh_mat[is.na(neigh_mat)] = 0
#  if(nrow(neigh_mat) < nrow(poly)){
#    missing = setdiff(rownames(poly),rownames(neigh_mat))
#    missing_df = data.frame(matrix(0,nrow=length(missing),ncol=ncol(neigh_mat)))
#    rownames(missing_df) = missing
#    colnames(missing_df)= colnames(neigh_mat)
#    neigh_mat = rbind(neigh_mat,missing_df)
#  }
#  neigh_mat = neigh_mat[rownames(poly),]
#  neigh_mat$tot_nn = rowSums(neigh_mat)
#  ids_anno = ids_anno %||% label
#  neigh_mat[[ids_anno]] = poly[[ids_anno]]
  if(inherits(data,'Seurat')& return.seurat ==T){
    KanData(data,'nnMat') = list(nb.method = nb.method,col.anno = label, nnMat = neigh_mat)
    return(data)
  }else{
    return(neigh_mat)
  }
}

#' @title Interrogate neighbour matrix
#' @name nn_query
#' @description
#' Annotate single cells based on specific features of their neighbor composition
#' @param data a neighbour matrix (nnMat) created with `knnMat`/`rnnMat`/`tnnMat` or a Seurat object with nb matrix stored as a KanData slot
#' @param query boolean string specifying the filtering criteria to annotate single cells. numerical and logical operator can be used.
#' Specify conditions as you were calling dplyr filter() function.
#' Please note that only information already included in the nb matrix can be considered for the query.
#' If you are planning to query your matrix using column with a numeric column name, e.g. 0/1/2/3..., specify the name of that column with backticks to avoid wrong numeric comparisons (e.g., 1 >= 1 will give a positive outcome independently on the values stored in the `1` column)
#' @param anno_name character string indicating the name that will be used to store the output annotation as a column in the neighbour matrix and in Seurat meta.data
#' @param diffexpr whether performing (TRUE) or not (FALSE) differential expression analysis between newly defined cell subgroups.
#' @param label name of variable to use for selecting cell type of interest to perform differential expression
#' @param which cell type label(s) to select from `label` variable
#' @param ... other arguments from FindMarkers() function passed to method when `diffexpr` is set to TRUE
#' @family nnMat
#' @returns indexed character vector with the new annotation for all single cells. Cells not respecting the query conditions will be labelled as 'N', while cells filtered through the query will be labelled as 'Y'
#' @export
nn_query = function(data = NULL,query="tot_nn > 4",anno_name = 'query_class',diffexpr=F,label=NULL,which=NULL,...){
  if(inherits(data,'dgCMatrix')){
    mat = as.data.frame(as.matrix(data))
  }else if(inherits(data,'Seurat')){
    if(is.null(KanData(data)) | is.null(KanData(data,'nnMat'))){
      stop('please provide a neighbourhood matrix created with nnMat or a Seurat object with neighbourhood matrix stored within "nnMat" Kandinsky slot')
    }else{
      mat = KanData(data,'nnMat')$nnMat
      if(inherits(mat,'dgCMatrix')){
        mat = as.data.frame(as.matrix(mat))
      }
    }
  }else{
    mat = data
  }
  query_ids = eval(parse(text=paste0("mat %>% dplyr::filter (",query,")")))
  mat[[anno_name]] = 'N'
  if(nrow(query_ids)>0){
	  mat[rownames(query_ids),][[anno_name]]='Y'
  }
  if(inherits(data,'dgCMatrix') | is.data.frame(data)){
    return(mat)
  }else if(inherits(data,'Seurat')){
    mat = mat[colnames(data),]
    KanData(data,'nnMat')$nnMat = mat
    data@meta.data[[anno_name]] = mat[[anno_name]]
    if(diffexpr==T){
      suppressWarnings({data = data[,data@meta.data[[label]] %in% which]})
      suppressMessages({data = UpdateSeuratObject(data)})
      Idents(data) = anno_name
      ct1 = length(data[[anno_name]][data[[anno_name]] == 'Y'])
      ct2 = length(data[[anno_name]][data[[anno_name]] == 'N'])
      message('Differential expression between ',ct1,' ',anno_name,'-Y and ',ct2,' ',anno_name,'-N cells')
      de = Seurat::FindMarkers(data,ident.1='Y',ident.2='N',...)
      return(de)
    }else{
      return(data)
    }
  }
}


#' @title Identify single cell/spot clusters on the basis of their neighbourhood composition
#' @name nbCluster
#' @description
#' Use model based clustering on single-cell neighbour composition to define cellular clusters
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @param n_clust numeric vector indicating the possible niche numbers to be tested to find the optimal number
#' @param seed numeric, random seed for reproducibility. Default is set to 347548
#' @returns Seurat object with new `nbClust_` variable column stored in the Seurat meta.data
#' @importFrom mclust mclustBIC Mclust
#' @seealso [cellcharter_clust()], [mesa_clust()], [spatialleiden_clust()], [spatopic_clust()], [kmeans_clust()], [giotto_clust()], [hoodscanR_clust()] for alternative clustering methods
#' @export
nbCluster = function(seurat=NULL,n_clust=3:12,seed=347548){
 # if(length(style) >1){
#    style = NULL
#  }
#  if(is.null(style)){
#    message('Setting default clustering method "mclust"')
#  }
#  style = style %||% 'mclust'
  nnmat = KanData(seurat,'nnMat')$nnMat
  method = KanData(seurat,'nnMat')$nb.method
  set.seed(seed)
  tot_nn = nnmat[,'tot_nn']
  tot_nn[tot_nn == 0] = 1
  nnmat[,'tot_nn'] = NULL
  nnmat = nnmat %>% dplyr::select_if(is.numeric)
  scaled_nnmat = Matrix::Diagonal(x=(1/tot_nn),names = T) %*% as(as.matrix(nnmat),'CsparseMatrix')
  rownames(scaled_nnmat) = rownames(nnmat)
  scaled_nnmat = as.data.frame(as.matrix(scaled_nnmat))
 # if(style=='mclust'){
    message('Performing model-based clustering using mclust...')
    clusts = mclust::Mclust(scaled_nnmat,G=n_clust)
    seurat[[paste0("nbClust.",method)]] = as.character(clusts$classification[colnames(seurat)])
#  }else if(style=='greed'){
#    if (!requireNamespace('greed', quietly = TRUE)) {
#      stop("Please install greed to perform greedy clustering")
#    }
#    message('Performing model-based clustering using greed...')
#    clusts = greed::greed(scaled_nnmat,model=greed::DiagGmm(),alg = greed::Seed(),K=max(n_clust))
#    clusts = greed::clustering(clusts)
#    names(clusts) = rownames(na.omit(scaled_nnmat))
#    clusts[setdiff(colnames(seurat),names(clusts))] = NA
#    seurat[[paste0("nbClust.",method)]] = as.character(clusts[colnames(seurat)])
#  }
  return(seurat)
}

##Alternative function for NB clustering. To be commented...
#test_cluster2 = function(seurat=NULL,k_range = 4:15,max_clust=20,seed=347548){
#  nnmat = KanData(seurat,'nnMat')$nnMat
#  tot_nn = nnmat[,'tot_nn']
#  tot_nn[tot_nn == 0] = 1
#  nnmat[,'tot_nn'] = NULL
#  nnmat = nnmat %>% dplyr::select_if(is.numeric)
#  scaled_nnmat = Matrix::Diagonal(x=(1/tot_nn),names = T) %*% as(as.matrix(nnmat),'CsparseMatrix')
#  rownames(scaled_nnmat) = rownames(nnmat)
#  scaled_nnmat = as.data.frame(as.matrix(scaled_nnmat))
#  clusters = lapply(k_range,function(x){
#    set.seed(seed)
#    clust = as.data.frame(kmeans(scaled_nnmat,centers=x,iter.max = 20)$cluster)
#    colnames(clust) = paste0('K_',x)
#    return(clust)
#  })
#  clusters = purrr::reduce(clusters,cbind) %>% mutate_all(as.factor)
#  set.seed(seed)
#  clusters = greed::greed(clusters,model=greed::Lca(),alg = greed::Seed(), K=max_clust)
#  clusters = greed::clustering(clusters)
#  names(clusters) = rownames(na.omit(scaled_nnmat))
#  clusters[setdiff(colnames(seurat),names(clusters))] = NA
#  seurat[["LCA_nbClust"]] = as.character(clusters[colnames(seurat)])
#  return(seurat)
#}


#' @title Calculate Shannon entropy of cell neighbourhood composition
#' @name nb_entropy
#' @description
#' Given cell type composition of single-cell neighbourhood, calculate Shannon entropy associated with each neighbourhood
#' @param seurat seurat object containing Kandinsky data slot with pre-computed neighbourhood composition matrix in `nnMat` slot
#' @returns seurat object with new metadata variable reporting neighbourhood entropy for each cell
#' @export
nb_entropy = function(seurat){
  nnmat = KanData(seurat,'nnMat')$nnMat
  tot_nn = nnmat[,'tot_nn']
  tot_nn[tot_nn == 0] = 1
  nnmat[,'tot_nn'] = NULL
  nnmat = nnmat %>% dplyr::select_if(is.numeric)
  scaled_nnmat = Matrix::Diagonal(x=(1/tot_nn)) %*% as(as.matrix(nnmat),'CsparseMatrix')
  seurat$nb_entropy = apply(scaled_nnmat,1,function(x){
    logv = log2(x)
    logv = ifelse(logv == -Inf,0,logv)
    return(-sum(x*logv))})
  return(seurat)
}



####################################################################################
########################= CELLCHARTER ANALYSIS SECTION =###############################
#################################= START HERE =###########################################
####################################################################################
#' @title Reticulate-based CellCharter clustering
#' @name cellcharter_clust
#' @description
#' CellCharter is called thrugh reticulate and interacts with Seurat/Kandinsky object after anndata conversion with anndataR
#' @details
#' CellCharter must be installed within a python environment together with its dependencies (squidpy,scanpy,scvi)
#' @param data a Seurat object containing Kandinsky data
#' @param sample_key character string specifying a variable stored in the Seurat object to use as sample/batch annotation. If not NULL, neighbour networks will be defined separately for each sample/batch. Default is NULL.
#' @param n_clust numeric, number of expected clusters. Can be a single number or a numeric vector if multiple cluster numbers needs to be tested
#' @param mode molecular data type. Must be one between 'rna' and 'protein'. If set to 'rna', the function will expect to find integer count data. Default is 'rna'.
#' @param preprocess boolean, whether molecular data still needs to be normalized. Default is TRUE.
#' @param seed numeric, random seed.
#' @param python_path path of python environment with CellCharter and its dependencies installed
#' @param use_kandinsky_nb boolean, whether to use neighbour network pre-computed by Kandinsky instead of default CellCharter Delaunay network. Default is TRUE.
#' @param embedding name of dimension reduction embedding to use for CellCharter clustering. If NULL, the function will call SCVI through reticulate. Default is NULL.
#' @param max_runs numeric, maximum number of repetitions for each value of number of clusters.
#' @param convergence_tol numeric, convergence tolerance for the clustering stability. If the Mean Absolute Percentage Error between consecutive iterations is below convergence_tol the algorithm stops without reaching max_runs.
#' @return updated Seurat object with 'cellcharter_clusters' annotation added to metadata slot
#' @export
cellcharter_clust = function(data=NULL,
                             sample_key=NULL,
                             n_clust=6:10,mode=c('rna','protein'),preprocess=T,
                             seed=347548,python_path=NULL,use_kandinsky_nb=T,
                             embedding = NULL,
                             max_runs=10,
                             convergence_tol=0.001){
  message('Initializing CellCharter with reticulate')
  python_path = python_path %||% Sys.getenv('RETICULATE_PYTHON')
  if(is.null(python_path) & Sys.getenv('RETICULATE_PYTHON') == "") stop('Please specificy a path for the python environment containing CellCharter with the parameter "python_path".\nExample: /path/to/my/version/of/python')
  if(Sys.getenv('RETICULATE_PYTHON') == ""){Sys.setenv(RETICULATE_PYTHON = python_path)}
  if(!requireNamespace("reticulate",quietly=T)) stop('You need to install reticulate to run CellCharter.\n Run install.packages("reticulate")')
  if(!requireNamespace("anndataR",quietly=T)) stop('You need to install anndataR to run CellCharter.\n Run BiocManager::install("anndataR") or conda install bioconda::anndataR if you are working within a conda environment')
  message('Loading CellCharter and its dependencies...')
  np = reticulate::import('numpy')
  sq = reticulate::import('squidpy')
  cc = reticulate::import('cellcharter')
  sc = reticulate::import('scanpy')
  if(is.null(embedding)){
    scvi = reticulate::import('scvi')
  }
  message('Converting Seurat object to anndata...')
  coords = sf::st_coordinates(sf::st_centroid(sf::st_geometry(KanData(data,'sf'))))
  colnames(coords) = c('X','Y')
  adata = reticulate::r_to_py(anndataR::as_AnnData(data,assay_name=DefaultAssay(data),output_class='ReticulateAnnData'))
  adata$obsm['spatial']= coords
  if(use_kandinsky_nb){
    adata$obsp['spatial_connectivities'] = reticulate::r_to_py(as(KanData(data,'nb'),'CsparseMatrix'))
    adata$obsp['spatial_distances'] = adata$obsp['spatial_connectivities']
  }
  adata$X = adata$layers['counts']
  adata$layers["counts"] = adata$X$copy()
  if(is.null(sample_key)){
    sample_key = sample_key %||% 'sample'
    adata$obs['sample'] = 'sample1'
  }
  adata$obs[sample_key]=adata$obs[sample_key]$astype('category')
  message('Running CellCharter...')
  if(length(mode) >0){mode = NULL}
  mode = mode %||% 'rna'
  if(!use_kandinsky_nb){
    sq$gr$spatial_neighbors(adata, library_key=sample_key, coord_type='generic', delaunay=T,percentile=99)
  }
  if(mode=='rna'){
    if(preprocess){
      sc$pp$normalize_total(adata, target_sum=1e6)
      sc$pp$log1p(adata)
    }
    if(is.null(embedding)){
      scvi$model$SCVI$setup_anndata(
        adata, 
        layer="counts", 
        batch_key=sample_key,
      )
      scvi$settings$seed = as.integer(seed)
      model = scvi$model$SCVI(adata)
      model$train(early_stopping=T, enable_progress_bar=T)
      adata$obsm['X_scVI'] = model$get_latent_representation(adata)
      adata$obsm['X_scVI'] = adata$obsm['X_scVI']$astype(np$float32)
    }else{
      adata$obsm['X_scVI'] = adata$obsm[embedding]$astype(np$float32)
    }
    cc$gr$aggregate_neighbors(adata, n_layers=as.integer(3), use_rep='X_scVI', out_key='X_cellcharter', sample_key=sample_key)
  }else if(mode=='protein'){
    if(preprocess){
      adata$raw = adata$copy()
      for(sample in adata$obs[sample_key]$cat$categories){
        adata$X[adata$obs[sample_key] == sample] = sc$pp$scale(adata[adata$obs[sample_key] == sample], copy=T)$X
      }
    }
    scvi$settings$seed = as.integer(seed)
    if(is.null(embedding)){
      model = cc$tl$TRVAE(adata)
      model$train(early_stopping=T, enable_progress_bar=T)
      adata$obsm['X_TRVAE'] = model$get_latent(adata$X, adata$obs[sample_key])
      adata$obsm['X_TRVAE'] = adata$obsm['X_TRVAE']$astype(np$float32)
    }else{
      adata$obsm['X_TRVAE'] = adata$obsm[embedding]$astype(np$float32)
    }
    cc$gr$aggregate_neighbors(adata, n_layers=as.integer(3), use_rep='X_TRVAE', out_key='X_cellcharter', sample_key=sample_key)
  }
  if(length(n_clust)>1){
    autok = cc$tl$ClusterAutoK(
      n_clusters=as.integer(n_clust),
      max_runs=as.integer(max_runs),
      convergence_tol=convergence_tol,
    )
    autok$fit(adata, use_rep='X_cellcharter')
    data$cellcharter_clusters = autok$predict(adata, use_rep='X_cellcharter')
    message('Optimal number of clusters: ',length(unique(as.character(data$cellcharter_clusters))))
  }else{
    clust = cc$tl$Cluster(
      n_clusters=as.integer(n_clust),
      random_state = as.integer(seed)
    )
    clust$fit(adata, use_rep='X_cellcharter')
    data$cellcharter_clusters = clust$predict(adata, use_rep='X_cellcharter')
  }
  data
}

####################################################################################
########################= SPATIALLEIDEN ANALYSIS SECTION=###############################
#################################= START HERE =###########################################
####################################################################################
#' @title Reticulate-based SpatialLeiden clustering
#' @name spatialleiden_clust
#' @description
#' SpatialLeiden is called thrugh reticulate and interacts with Seurat/Kandinsky object after anndata conversion with anndataR
#' @details
#' SpatialLeiden must be installed within a python environment together with its dependencies (squidpy,scanpy,leidenalg)
#' @param data a Seurat object containing Kandinsky data
#' @param n_clust numeric, number of expected clusters
#' @param nfeatures numeric, number of variable features to use for dimensionality reduction (PCA)
#' @param seed numeric, random seed
#' @param python_path path of python environment with SpatialLeiden and its dependencies installed
#' @param resolution numeric, resolution for the latent space and spatial layer, respectively.
#' @param embedding name of dimension reduction embedding to use for SpatialLeiden clustering. If NULL, the function will run PCA with scanpy through reticulate. Default is NULL.
#' @param layer_ratio numeric, the ratio of the weighting of the layers; latent space vs spatial. A higher ratio will increase relevance of the spatial neighbors and lead to more spatially homogeneous clusters
#' @return updated Seurat object with 'spatialleiden_clusters' annotation added to metadata slot
#' @export
spatialleiden_clust = function(data=NULL,
                               n_clust=NULL,nfeatures=2000,
                               seed=347548,python_path=NULL,resolution=1,
                               embedding = NULL,
                               layer_ratio=1.5){
  python_path = python_path %||% Sys.getenv('RETICULATE_PYTHON')
  if(is.null(python_path) & Sys.getenv('RETICULATE_PYTHON') == "") stop('Please specificy a path for the python environment containing SpatialLeiden with the parameter "python_path".\nExample: /path/to/my/version/of/python')
  if(Sys.getenv('RETICULATE_PYTHON') == ""){Sys.setenv(RETICULATE_PYTHON = python_path)}
  if(!requireNamespace("reticulate",quietly=T)) stop('You need to install reticulate to run SpatialLeiden.\n Run install.packages("reticulate")')
  if(!requireNamespace("anndataR",quietly = T)) stop('You need to install anndataR to run SpatialLeiden.\n Run BiocManager::install("anndataR") or conda install bioconda::anndataR if you are working within a conda environment')
  message('Loading SpatialLeiden and its dependencies...')
  sq = reticulate::import('squidpy')
  sc = reticulate::import('scanpy')
  sl = reticulate::import('spatialleiden')
  search_resolution = sl$search_resolution
  message('Converting Seurat object to anndata...')
  coords = sf::st_coordinates(sf::st_centroid(sf::st_geometry(KanData(data,'sf'))))
  colnames(coords) = c('X','Y')
  data = NormalizeData(data)
    if(is.null(embedding)){
	    data = FindVariableFeatures(nfeatures = nfeatures) %>%
	    ScaleData() %>%
	    RunPCA(data)
    }
  adata = reticulate::r_to_py(anndataR::as_AnnData(data,assay_name=DefaultAssay(data),output_class='ReticulateAnnData'))
  adata$obsm['spatial']= coords
  adata$obsp['spatial_connectivities'] = reticulate::r_to_py(as(KanData(data,'nb'),'CsparseMatrix'))
  adata$obsp['spatial_distances'] = adata$obsp['spatial_connectivities']
  adata$X = adata$layers['counts']
  message('Running SpatialLeiden...')
  seed = as.integer(seed)
  if(is.null(embedding)){
    sc$pp$neighbors(adata, random_state=seed,use_rep='pca')
  }else{
    sc$pp$neighbors(adata, random_state=seed,use_rep = embedding)
  }
  if(is.null(n_clust)){
    sl$spatialleiden(adata, layer_ratio=layer_ratio, directed=c(F, T), random_state=seed,resolution=resolution)
  }else{
    resolutions = search_resolution(adata,as.integer(n_clust),
                                    latent_kwargs=list(directed=F,random_state=seed),
                                    spatial_kwargs=list(directed=F,random_state=seed,layer_ratio=layer_ratio),
    )
  }
  data$spatialleiden_clusters = reticulate::py_to_r(adata$obs['spatialleiden'])
  data
}


####################################################################################
########################= MESA ANALYSIS SECTION =###############################
#################################= START HERE =###########################################
####################################################################################
#' @title Reticulate-based MESA clustering
#' @name mesa_clust
#' @description
#' MESA is called thrugh reticulate and interacts with Seurat/Kandinsky object after anndata conversion with anndataR
#' @details
#' MESA must be installed within a python environment together with its dependencies (scanpy,sklearn.cluster)
#' @param data a Seurat object containing Kandinsky data 
#' @param k numeric, number of nearest neighbours to consider for MESA gene expression aggregation 
#' @param n_clust numeric, number of expected clusters
#' @param seed numeric, random seed
#' @param python_path path of python environment with MESA and its dependencies installed
#' @return updated Seurat object with 'mesa_clusters' annotation added to metadata slot
#' @export
mesa_clust = function(data=NULL,
                      k=10,
                      n_clust=10,
                      seed=347548,python_path=NULL){
  message('Initializing MESA with reticulate')
  if(is.null(python_path)) stop('Please specificy a path for the python environment containing MESA with the parameter "python_path".\nExample: /path/to/my/version/of/python')
  if(Sys.getenv('RETICULATE_PYTHON') == ""){Sys.setenv(RETICULATE_PYTHON = python_path)}
  if(!requireNamespace("reticulate",quietly=T)) stop('You need to install reticulate to run MESA.\n Run install.packages("reticulate")')
  if(!requireNamespace("anndataR",quietly=T)) stop('You need to install anndataR to run MESA.\n Run BiocManager::install("anndataR") or conda install bioconda::anndataR if you are working within a conda environment')
  np = reticulate::import('numpy')
  pd = reticulate::import('pandas')
  sc = reticulate::import('scanpy')
  mesa = reticulate::import('mesa')
  multiomics = mesa$multiomics
  sklearn = reticulate::import('sklearn.cluster')
  KMeans = sklearn$KMeans
  
  coords = sf::st_coordinates(sf::st_centroid(sf::st_geometry(KanData(data,'sf'))))
  colnames(coords) = c('X','Y')
  coords = reticulate::r_to_py(coords)
  exp_mat = pd$DataFrame(reticulate::r_to_py(LayerData(data,'data'))$toarray())
  colnames(exp_mat) = colnames(data)
  rownames(exp_mat) =rownames(data)
  knn = multiomics$multiomics_spatial$get_spatial_knn_indices(coords,n_neighbors=as.integer(k))
  avg_exp = multiomics$multiomics_spatial$get_avg_expression_neighbors(exp_mat,knn)
  avg_exp = np$stack(avg_exp)
  kmeans_rna = KMeans(n_clusters=as.integer(n_clust), random_state=as.integer(seed))$fit(avg_exp)
  data$mesa_clusters = as.character(kmeans_rna$labels_)
  data
}

####################################################################################
########################= SpatialTopic ANALYSIS SECTION =###############################
#################################= START HERE =###########################################
####################################################################################

#install.packages('SpaTopic')
#' @title SpaTopic clustering
#' @name spatopic_clust
#' @description
#' SpatialTopic clustering performed by calling the R package SpaTopic through Kandinsky
#' @details
#' SpaTopic must be installed within the same environmnent of Kandinsky
#' @param data Seurat object containing Kandinsky data 
#' @param label character string specifying the variable name to be used to defne cell annotation groups
#' @param sample_key character string specifying a variable stored in the Seurat object to use as sample/batch annotation. If not NULL, neighbour networks will be defined separately for each sample/batch. Default is NULL.
#' @param n_clust numeric, number of expected clusters. Can be a single number or a numeric vector if multiple cluster numbers needs to be tested
#' @param sigma Default is 50. The lengthscale of the Nearest-neighbor Exponential Kernel. Sigma controls the strength of decay of correlation with distance in the kernel function. Please check the paper for more information. Need to be adjusted based on the image resolution
#' @param region_radius Default is 400. The radius for each grid square when sampling region centers for each image. Need to be adjusted based on the image resolution and pattern complexity.
#' @param kneigh Default is 5. Only consider the top 5 closest region centers for each cell.
#' @param trace Default is TRUE Compute and save log likelihood, Ndk, Nwk for every posterior samples. Useful when you want to use DIC to select number of topics, but it is time consuming to compute the likelihood for every posterior samples.
#' @param thin Default is 20. Key parameter in Gibbs sampling. Collect a posterior sample for every thin=20 iterations.
#' @param burnin Default is 1000. Key parameter in Gibbs sampling. Start to collect posterior samples after 1000 iterations. You may increase the number of iterations for burn-in for highly complex tissue images.
#' @param niter Default is 200. Key parameter in Gibbs sampling. Number of posterior samples collected for model inference.
#' @param seed Default is 347548 Random seed.
#' @return updated Seurat object with 'spatopic_clusters' annotation added to metadata slot
#' @export
spatopic_clust=function(data=NULL,label=NULL,sample_key=NULL,n_clust=2:9,
                        sigma = 25, region_radius = 500, kneigh=5, trace = TRUE, thin = 20, burnin = 1000, niter = 200,seed=347548){
  if(!requireNamespace("SpaTopic",quietly=T)){
    stop('Please install SpaTopic first by running install.packages("SpaTopic")')
  }
  if(is.null(label)){
    stop('Please specify "label" argument.')
  }
  vars= ifelse(is.null(sample_key),label,c(label,sample_key))
  if(is.list(data)){
    input = lapply(data,function(x){
      polys = populate_sf(x,vars=vars,return.seurat=F)
      polys[,c('X','Y')] = sf::st_coordinates(sf::st_centroid(sf::st_geometry(polys)))
      if(is.null(sample_key)){
        input = polys %>% sf::st_drop_geometry() %>% dplyr::select(X=.data[["X"]],Y=.data[["Y"]],type=.data[[label]]) %>% dplyr::mutate(image='sample')
      }else{
        input = polys %>% sf::st_drop_geometry() %>% dplyr::select(X=.data[["X"]],Y=.data[["Y"]],type=.data[[label]],image= .data[[sample_key]])
      }
      input$cell_ID = rownames(input)
      input})
    input = purrr::reduce(input,rbind)
  }else{
    polys = populate_sf(data,vars=vars,return.seurat=F)
    polys[,c('X','Y')] = sf::st_coordinates(sf::st_centroid(sf::st_geometry(polys)))
    if(is.null(sample_key)){
      input = polys %>% sf::st_drop_geometry() %>% dplyr::select(X=.data[["X"]],Y=.data[["Y"]],type=.data[[label]]) %>% dplyr::mutate(image='sample')
    }else{
      input = polys %>% sf::st_drop_geometry() %>% dplyr::select(X=.data[["X"]],Y=.data[["Y"]],type=.data[[label]],image= .data[[sample_key]])
    }
    input$cell_ID = rownames(input)
    if(length(unique(input$image))>1){
      input = split(input,input$image)
    }
    # add sample ID in a new 'image' column (required by SpaTopic, especially when working with multiple samples)
  }
  #Tuning SpaTopic model following autors recommendation
  if(length(n_clust) <2){
    res <-SpaTopic::SpaTopic_inference(input, ntopics = n_clust,kneigh=kneigh,
                                       sigma = sigma, region_radius = region_radius, 
                                       trace = trace, thin = thin, burnin = burnin, niter = niter,seed=seed)
  }else{
    res = list()
    for(topic in n_clust){
      res[[topic]] <-SpaTopic::SpaTopic_inference(input, ntopics = topic,kneigh=kneigh,
                                                  sigma = sigma, region_radius = region_radius, 
                                                  trace = trace, thin = thin, burnin = burnin, niter = niter,seed=seed)
    }
    Perplexity <- unlist(lapply(res, function(x) x$Perplexity))
    perx_df <- data.frame(
      Topics = n_clust,
      Perplexity = Perplexity
    )
    best_topic = perx_df[perx_df$Perplexity <= min(perx_df$Perplexity),]
    message('Optimal number of clusters: ',best_topic$Topics)
    best_topic
    
    res = res[[best_topic$Topics]]
  }
  if(is.list(data)){
    return(input %>% dplyr::mutate(spatopic_clusters=res$cell_topics))
  }else{
    data@meta.data$spatopic_clusters = as.character(res$cell_topics)
    return(data)
  }
}

####################################################################################
########################= HoodscanR ANALYSIS SECTION =###############################
#################################= START HERE =###########################################
####################################################################################
#Run BiocManager::install("hoodscanR")
#' @title hoodscanR clustering
#' @name hoodscanR_clust
#' @description
#' hoodscanR clustering performed by calling the R package hoodscanR through Kandinsky
#' @details
#' hoodscanR must be installed within the same environmnent of Kandinsky
#' @param data Seurat object containing Kandinsky data 
#' @param k numeric, number of nearest neighbours to consider for hoodscanR clustering
#' @param label character string specifying the variable name to be used to defne cell annotation groups
#' @param n_clust numeric, number of expected clusters. Can be a single number or a numeric vector if multiple cluster numbers needs to be tested
#' @param sample_key character string specifying a variable stored in the Seurat object to use as sample/batch annotation. If not NULL, neighbour networks will be defined separately for each sample/batch. Default is NULL
#' @return updated Seurat object with 'hoodscanR_clusters' annotation added to metadata slot
#' @export
hoodscanR_clust = function(data=NULL,k=100,label=NULL,n_clust=5,sample_key=NULL){
  if((!requireNamespace("hoodscanR",quietly = T))){
    stop('Please install hoodscanR first by running BiocManager::install("hoodscanR")')
  }
  if(is.null(label)){
    stop('Please specify "label" argument.')
  }
  if(!is.null(sample_key)){
    samples = unique(data@meta.data[,sample_key])
  }
  message('Preparing data for hoodscanR...')
  spe = SpatialExperiment::SpatialExperiment(assays=list(counts=LayerData(data,layer='counts')),
                                             rowData=data.frame(row.names=rownames(data)),
                                             colData=data@meta.data,metadata=list(),reducedDims=list(),
                                             altExps=list(),sample_id = 'my_spe',
                                             spatialCoords=sf::st_coordinates(sf::st_centroid(sf::st_geometry(KanData(data,'sf'))))
  )
  rownames(SpatialExperiment::spatialCoords(spe)) <- rownames(KanData(data,'sf'))
  message('Running hoodscanR...')
  if(!is.null(sample_key)){
    nns = list()
    for(s in samples){
      nns[[s]] = hoodscanR::findNearCells(spe[,spe@colData[,sample_key] ==s], k = k,anno_col=label)
    }
    cells = lapply(nns,function(x){x[[1]]})
    dists = lapply(nns,function(x){x[[2]]})
    fnc = list(cells=purrr::reduce(cells,rbind),distance=purrr::reduce(dists,rbind))
  }else{
    fnc <- hoodscanR::findNearCells(spe, k = k,anno_col=label)
  }
  pm <-  hoodscanR::scanHoods(fnc$distance)
  hoods <-  hoodscanR::mergeByGroup(pm, fnc$cells)
  spe <-  hoodscanR::mergeHoodSpe(spe, hoods)
  spe <-  hoodscanR::calcMetrics(spe, pm_cols = colnames(hoods))
  spe <-  hoodscanR::clustByHood(spe, pm_cols = colnames(hoods), k = n_clust)
  data@meta.data[,c('hoodscanR_entropy','hoodscanR_perplexity','hoodscanR_clusters')] = as.data.frame(spe@colData[,c('entropy','perplexity','clusters')])
  data
}

####################################################################################
########################= GIOTTO ANALYSIS SECTION=###############################
#################################= START HERE =###########################################
####################################################################################
#' @title Giotto clustering
#' @name giotto_clust
#' @description
#' Giotto clustering performed by calling the R package Giotto through Kandinsky
#' @details
#' Giotto must be installed within the same environmnent of Kandinsky
#' @param data Seurat object containing Kandinsky data 
#' @param n_clust numeric, number of expected clusters. Can be a single number or a numeric vector if multiple cluster numbers needs to be tested
#' @param ngenes numeric, number of top variable genes to consider for HMRF domain identification Default is 250. If dataset contains less than 250 genes, all genes will be considered.
#' @param betas numeric vector of beta values to use for Giotto HMRF domain identification. Currently only supports two beta values at a time
#' @param seed numeric, random seed.
#' @param cores numeric, number of cores to use for parallelization
#' @param python_path path of python environment with MESA and its dependencies installed
#' @return updated Seurat object with 'giotto_clusters' annotation added to metadata slot for each beta value
#' @export
giotto_clust = function(data=NULL,n_clust=5,ngenes=NULL,betas=c(0,5),seed=347548,cores=1,python_path = NULL){
  if((!requireNamespace("Giotto", quietly = TRUE))){
    stop('Please install Giotto first by running pak::pkg_install("giotto-suite/Giotto")')
  }
warnings('Please install graphcoloring and smfishHmrf packages before using this function:\ndevtools::install_github("saurfang/graphcoloring")\ninstall.packages("smfishHmrf")\n')
Giotto::createGiottoInstructions(python_path = python_path)
genv_exists <- Giotto::checkGiottoEnvironment()
if(!genv_exists){
    # The following command need only be run once to install the Giotto environment
	warning('python env not found. Giotto will create/search its default environment')
	Giotto::installGiottoEnvironment()
}
if(all(Layers(data) != 'data')){
    stop('Please normalize data before running this function by calling NormalizeData()')
  }
  message('Preparing data for Giotto...')
  data$cell_ID = colnames(data)
  
  #Create Giotto object
  spat_coord = sf::st_coordinates(sf::st_centroid(KanData(data,'sf')$geometry))
  spat_coord = as.data.frame(spat_coord)
  rownames(spat_coord) = colnames(data)
  spat_coord$cell_ID <- rownames(spat_coord)
  colnames(spat_coord) <- c("sdimx", "sdimy", "cell_ID")
  adj = as(KanData(data,'nb'),'CsparseMatrix')
  s <- Matrix::summary(adj)
  
  networkDT <- data.table::data.table(
    from = rownames(adj)[s[,1]],
    to   = colnames(adj)[s[,2]],
    weight = s[,3]
  )
  giotto <- Giotto::createGiottoObject(
    expression = list(normalized=LayerData(data,'data')),
    spatial_locs = spat_coord,cell_metadata=data.table::as.data.table(data@meta.data)
  )   
  giotto <- Giotto::processExpression(giotto,
                                      expression_values = "normalized",
                                      Giotto::scaleParam('zscore'),
                                      name = "scaled")
  
  giotto <- Giotto::createSpatialNetwork(giotto,
                                         method = "Delaunay", 
                                         name = "spatial_network"
  )
  snet <- Giotto::getSpatialNetwork(
    giotto,
    name = "spatial_network"
  )
  
  slot(snet, "networkDT") <- networkDT
  
  gobject <- Giotto::setSpatialNetwork(
    giotto,
    x = snet,
    name = "spatial_network"
  ) 
  giotto = Giotto::identifyTMAcores(giotto,spatial_network_name='spatial_network') # useful to split networks from independent samples
  if(nrow(data)<=250){
    ngenes = nrow(data)
  }
  ngenes = ngenes %||% 250
  giotto <- Giotto::binSpect(giotto, 
                             bin_method = "rank",
                             calc_hub = TRUE, 
                             hub_min_int = 5,
                             spatial_network_name = "spatial_network",
                             return_gobject = TRUE,
                             cores=cores
  )
  HMRF_init_obj <- Giotto::initHMRF_V2(gobject = giotto, cl.method = "km",spatial_network_name = "spatial_network",nstart=50,k=n_clust)
  res = Giotto::doHMRF_V2(HMRF_init_obj = HMRF_init_obj, betas = c(betas[1], betas[2], 2))
  res2 = res[[2]]$prob
  res = res[[1]]$prob
  clusters2  = (reshape2::melt(res2)) %>% dplyr::group_by(.data[["Var1"]]) %>% dplyr::filter(.data[["value"]]==max(.data[["value"]])) %>% dplyr::ungroup()
  colnames(clusters2) = c('cell_ID','giotto_clusters2','prob')
  clusters  = (reshape2::melt(res)) %>% dplyr::group_by(.data[["Var1"]]) %>% dplyr::filter(.data[["value"]]==max(.data[["value"]])) %>% dplyr::ungroup()
  colnames(clusters) = c('cell_ID','giotto_clusters','prob')
  clusts = merge(clusters[,1:2],clusters2[,1:2],by='cell_ID')
  meta = merge(data@meta.data,clusts,by='cell_ID')
  rownames(meta) = meta$cell_ID
  data@meta.data[,c('giotto_clusters','giotto_clusters2')] = meta[colnames(data),c('giotto_clusters','giotto_clusters2')]
  data
}

####################################################################################
########################= ISCHIA/K-Means ANALYSIS SECTION =###############################
#################################= START HERE =###########################################
####################################################################################
#' @title K-Means clustering
#' @name kmeans_clust
#' @description
#' ISCHIA-inspired K-Means clustering
#' @details
#' When more than one cluster configurations are provided, the function chooses the one minimizing Within-Cluster Sum of Squares (WSS), as implemented in the R package ISCHIA 
#' @param data Seurat object containing Kandinsky data 
#' @param label character string specifying the variable name to be used to defne cell annotation groups
#' @param n_clust numeric, number of expected clusters. Can be a single number or a numeric vector if multiple cluster numbers needs to be tested
#' @param seed numeric, random seed.
#' @return updated Seurat object with 'kmeans_clusters' annotation added to metadata slot
#' @export
kmeans_clust = function(data=NULL,n_clust=1:10,label=NULL,seed=347548){
  if(is.null(label)){
    stop('Please specify "label" argument.')
  }
  nnmat = nnMat(data,label=label,return.seurat = F)
  set.seed(seed)
  tot_nn = nnmat[, "tot_nn"]
  tot_nn[tot_nn == 0] = 1
  nnmat[, "tot_nn"] = NULL
  nnmat = nnmat %>% dplyr::select_if(is.numeric)
  scaled_nnmat = Matrix::Diagonal(x = (1/tot_nn), names = T) %*% 
    as(as.matrix(nnmat), "CsparseMatrix")
  rownames(scaled_nnmat) = rownames(nnmat)
  scaled_nnmat = as.data.frame(as.matrix(scaled_nnmat))
  if(length(n_clust)>1){
    sample.wss = NULL
    for (i in n_clust){
      sample.fit=kmeans(scaled_nnmat, centers = i)
      sample.wss <- c(sample.wss, sample.fit$tot.withinss)
    }
    sample.wss = diff(sample.wss)*n_clust[-1]
    best_k = n_clust[-1][sample.wss==min(sample.wss)]
    message('Optimal number of clusters: ',best_k)
    clusts = kmeans(scaled_nnmat, centers = best_k)
  }else{
    clusts = kmeans(scaled_nnmat, centers = n_clust)
  }
  data$kmeans_clusters = clusts$cluster
  data
}




