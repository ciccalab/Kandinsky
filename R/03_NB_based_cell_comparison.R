#' @title Neighbour composition matrix
#' @name nnMat
#' @description
#' From a metadata file, given a list of cells with the corrisponding annotated cell type, build a neighbor matrix reporting the number of occurrences for each cell type
#'
#' @param data a Seurat object containing Kandinsky data
#' @param method character string specifying the method to be used to identify neighbours.
#' Must be one of the following:
#' 'Q': queen contiguity method,check for contact (not overlap) between any edge or side od two polygons (refers to the queen movement rule in chess). Currently only applicable for Visium/Visium-HD data
#' 'C': centroid-based method, use maximum centroid distance threshold to identify spot/cell neighbours
#' 'K': KNN method, define k closest neighbours to each spot/cell
#' 'M': membrane-based method, check for the occurrence of a physical contact/intersection within a distance threshold between cell boundaries. Not applicable in the case of Visium spots.
#' When argument `method` is not specified or is set to `NULL`, the function will use `nb` slot already stored in Kandinsky object to create neighbour matrix
#' @param snap numeric, extra distance accepted between polygon borders for contiguity relation. Applied when 'method' is set to `Q`. Can't work with point geometries like cell centroids.
#' @param k numeric, number of nearest neighbours to be searched. Applied when `method` is set to `K`
#' @param d.max numeric, maximum centroid or membrane distance threshold. Applied when `method` is set to `C` or `M`
#' @param label character string specifying the variable name to be used to defne cell annotation groups
#' @param ids_anno optional, character string specifying any extra variable to be added to the final distance matrix
#' @param return.seurat boolean, whether returning neighbour matrix alone (FALSE) or the input Seurat object with a new `nnMat` slot as part of the Kandinsky data
#' @family nnMat
#' @importFrom magrittr %>%
#' @import spatialreg
#' @export
nnMat = function(data,method=c('Q','K','C','M'),snap=NULL,k=NULL,d.max=NULL,label='final_anno',ids_anno=NULL,return.seurat=T){
  data@meta.data$ids = rownames(data@meta.data)
  poly = populate_sf(data,vars=c('ids',label,ids_anno),return.seurat = F)
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
        mat = queen_nb(poly,snap=snap)
        nb.method = 'Q'
      }else{
        stop('method "Q" can only be applied to spot-based technologies (Visium/Visium-HD)')
      }
    }
    if(method == 'K'){
      mat = knn_nb(poly,k=k)
      nb.method = paste0('K_',k)
    }
    if(method == 'C'){
      mat = centroid_nb(poly,d.max=d.max)
      nb.method = paste0('C_',d.max)
    }
    if(method == 'M'){
      mat = membrane_nb(poly,d.max=d.max)
      nb.method = paste0('M_',d.max)
    }
    message('Neighbour method: ',nb.method)
    mat = as(mat,'CsparseMatrix')
  }
  neigh = Matrix::which(mat!=0,arr.ind=TRUE)  %>% as.data.frame()
  #Add cell ids and neighboring cell types to neigh object
  neigh$cell_ID = poly[neigh$row,][['ids']]
  neigh$col.type = as.character(poly[neigh$col,][[colnames(poly)[colnames(poly) == label]]])
  #Count the occurrences for each cell type among the neighbors of each cell
  neigh_mat = neigh %>%
    dplyr::group_by(.data$cell_ID,.data$col.type) %>%
    dplyr::mutate(n = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(c(.data$cell_ID,.data$col.type,.data$n))
  #From tabular to matrix format
  neigh_mat = reshape2::dcast(neigh_mat,cell_ID~col.type,value.var='n',fun.aggregate = length)
  rownames(neigh_mat) = neigh_mat$cell_ID
  neigh_mat$cell_ID = NULL
  #Be sure that all column in the final matrix are in numeric format
  neigh_mat = neigh_mat %>% dplyr::mutate_all(as.numeric)
  #Convert any possible NA value to 0
  neigh_mat[is.na(neigh_mat)] = 0
  if(nrow(neigh_mat) < nrow(poly)){
    missing = setdiff(rownames(poly),rownames(neigh_mat))
    missing_df = data.frame(matrix(0,nrow=length(missing),ncol=ncol(neigh_mat)))
    rownames(missing_df) = missing
    colnames(missing_df)= colnames(neigh_mat)
    neigh_mat = rbind(neigh_mat,missing_df)
  }
  neigh_mat = neigh_mat[rownames(poly),]
  neigh_mat$tot_nn = rowSums(neigh_mat)
  ids_anno = ids_anno %||% label
  neigh_mat[[ids_anno]] = poly[[ids_anno]]
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
  query_ids = eval(parse(text=paste0("mat %>% dplyr::filter (",query,") %>% rownames()")))
  mat[[anno_name]] = 'N'
  mat[query_ids,][[anno_name]]='Y'
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


#' @title Identify single cell/spot clusters on the basis of their neighbourhood compositin
#' @name nbCluster
#' @description
#' Use model based clustering on single-cell neighbour composition to define cellular clusters
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @param n_clust numeric vector indicating the possible niche numbers to be tested to find the optimal number
#' @param seed numeric, random seed for reproducibility. Default is set to 347548
#' @returns Seurat object with new `nbClust_` variable column stored in the Seurat meta.data
#' @importFrom mclust mclustBIC Mclust
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
