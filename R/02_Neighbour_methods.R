#' @title queen contiguity neighbour identification
#' @name queen_nb
#' @description
#' Define neighbour relationships according to Queen criterion for polygon contiguity.
#'
#' @details
#' Based on poly2nb() function from spdep package.
#' Designed to work specifically with Visium / Visium HD data only.
#' Each Visium spot/bin will have a neighbour link with only the first nearest ring of surrounding spots (Visium, n=6) or bins (Visium HD, n=8)
#'
#' @param poly a sf data frame with polygon geometry
#' @param snap numeric, extra distance accepted between polygon borders for contiguity relation
#' @param layers numeric, number of concentric contiguous layers to include in spot neighbourhood. Only Applied when `nb.method = Q`. Default is 1
#' @returns weighted neighbour matrix in a listw object from spdep package
#' @export
#' @family nb_funcs
queen_nb = function(poly,snap=1,layers=1){
  message('Queen neighbour method is designed to work with Visium / Visium HD data only')
  if(!is.data.frame(poly)){
    nbs = lapply(poly,function(p){
      is.poly = substring(class(sf::st_geometry(p))[1], 5) == 'POLYGON'
      if(!is.poly){
        stop('Input geometries should be polygons')
      }
      nb = spdep::poly2nb(p,snap=snap*layers,row.names = rownames(p))
      as(spdep::nb2listw(nb,style = 'B',zero.policy=T),'CsparseMatrix')
    })
    nb_ids = c()
    for(n in 1:length(nbs)){
      nb_ids = c(nb_ids,rownames(nbs[[n]]))
    }
    nbs = Matrix::bdiag(nbs)
    rownames(nbs) = nb_ids
    colnames(nbs) = nb_ids
    (spdep::mat2listw(nbs,row.names=rownames(nbs),style = 'B',zero.policy=T))
  }else{
    is.poly = substring(class(sf::st_geometry(poly))[1], 5) == 'POLYGON'
    if(!is.poly){
      stop('Input geometries should be polygons')
    }
  nb = spdep::poly2nb(poly,snap=snap*layers,row.names = rownames(poly))
  spdep::nb2listw(nb,style = 'B',zero.policy=T)
  }
}

#' @title membrane-based neighbour identification
#' @name membrane_nb
#' @description
#' Define neighbour relationships according to contact distance treshold between polygon borders.
#'
#' @details
#' Based on polygon intersection check done with st_intersects() function from sf package.
#' When minimum distance is higher than 0, a round buffer of size max.dist/2 is applied to all polygons before calling st_intersects()
#'
#' @param poly a sf data frame with polygon geometry
#' @param d.max numeric, maximum distance accepted between polygon borders to define neighbour relation
#' @returns weighted neighbour matrix in a listw object from spdep package
#' @export
#' @family nb_funcs
membrane_nb = function(poly,d.max=0){
  if(!is.data.frame(poly)){
    nbs = lapply(poly,function(p){
      if(d.max>0){
        p = sf::st_buffer(p,dist = d.max/2)
      }
      nb = sf::st_intersects(p)
      nb = lapply(seq_len(length(nb)),function(x) nb[[x]][nb[[x]] != x])
      attr(nb,'class') = c('sgbp','list')
      #Convert neighbors list to a nb object
      attrs = attributes(nb)
      nb = lapply(nb, function(i) { if(length(i) == 0L) 0L else i } )
      attributes(nb) = attrs
      class(nb) = "nb"
      attr(nb,'region.id') = rownames(p)
      as(spdep::nb2listw(nb,style='B',zero.policy=T),'CsparseMatrix')
    })
    nb_ids = c()
    for(n in 1:length(nbs)){
      nb_ids = c(nb_ids,rownames(nbs[[n]]))
    }
    nbs = Matrix::bdiag(nbs)
    rownames(nbs) = nb_ids
    colnames(nbs) = nb_ids
    (spdep::mat2listw(nbs,row.names=rownames(nbs),style = 'B',zero.policy=T))
  }else{
  if(d.max>0){
    poly = sf::st_buffer(poly,dist = d.max/2)
  }
  nb = sf::st_intersects(poly)
  nb = lapply(seq_len(length(nb)),function(x) nb[[x]][nb[[x]] != x])
  attr(nb,'class') = c('sgbp','list')
  #Convert neighbors list to a nb object
  attrs = attributes(nb)
  nb = lapply(nb, function(i) { if(length(i) == 0L) 0L else i } )
  attributes(nb) = attrs
  class(nb) = "nb"
  attr(nb,'region.id') = rownames(poly)
  spdep::nb2listw(nb,style='B',zero.policy=T)
  }
}

#' @title centroid distance neighbour identification
#' @name centroid_nb
#' @description
#' Define neighbour relationships according to radial distance threshold between polygon centroids.
#'
#' @details
#' Based on dnearneigh function from spdep package
#'
#' @param poly a sf data frame with polygon geometry
#' @param d.max numeric, maximum distance accepted between polygon centroids to define neighbour relation
#' @returns weighted neighbour matrix in a listw object from spdep package
#' @export
#' @family nb_funcs
centroid_nb = function(poly,d.max=40){
  if(!is.data.frame(poly)){
    nbs = lapply(poly,function(p){
      is.poly = substring(class(sf::st_geometry(p))[1], 5) == 'POLYGON'
      if(is.poly){
        suppressWarnings({p = sf::st_centroid(p)})
      }
      nb = spdep::dnearneigh(p,d1=0,d2=d.max,row.names = rownames(p))
      as(spdep::nb2listw(nb,style='B',zero.policy=T),'CsparseMatrix')
    })
    nb_ids = c()
    for(n in 1:length(nbs)){
      nb_ids = c(nb_ids,rownames(nbs[[n]]))
    }
    nbs = Matrix::bdiag(nbs)
    rownames(nbs) = nb_ids
    colnames(nbs) = nb_ids
    (spdep::mat2listw(nbs,row.names=rownames(nbs),style = 'B',zero.policy=T))
  }else{
    is.poly = substring(class(sf::st_geometry(poly))[1], 5) == 'POLYGON'
    if(is.poly){
      suppressWarnings({poly = sf::st_centroid(poly)})
    }
    nb = spdep::dnearneigh(poly,d1=0,d2=d.max,row.names = rownames(poly))
    spdep::nb2listw(nb,style='B',zero.policy=T)
  }
}


#' @title K-nearest neighbour neighbour identification
#' @name knn_nb
#' @description
#' Define neighbour relationships according to K-nearest neighbour algorithm.
#'
#' @details
#' based on knearneigh function from spdep package
#'
#' @param poly a sf data frame with polygon geometry
#' @param k numeric, number of nearest neighbour to select with knn algorithm
#' @returns weighted neighbour matrix in a listw object from spdep package
#' @export
#' @family nb_funcs
knn_nb = function(poly=NULL,k=10){
  if(!is.data.frame(poly)){
    nbs = lapply(poly,function(p){
      is.poly = substring(class(sf::st_geometry(p))[1], 5) == 'POLYGON'
      if(is.poly){
        suppressWarnings({p = sf::st_centroid(p)})
      }
      nb = spdep::knearneigh(p,k=k)
      #knb = asplit(FNN::get.knn(data=sf::st_coordinates(p),k=k)$nn.index,1)
      #attr(knb,'class') = c('sgbp','list')
      ##Convert neighbors list to a nb object
      #attrs = attributes(knb)
      #knb = lapply(knb, function(i) { if(length(i) == 0L) 0L else sort(i) } )
      #attributes(knb) = attrs
      #attr(knb,'region.id') = rownames(p)
      #class(knb) = "nb"
      as(spdep::nb2listw(spdep::knn2nb(nb,row.names = rownames(p)),style='B',zero.policy=T),'CsparseMatrix')
    })
    nb_ids = c()
    for(n in 1:length(nbs)){
      nb_ids = c(nb_ids,rownames(nbs[[n]]))
    }
    nbs = Matrix::bdiag(nbs)
    rownames(nbs) = nb_ids
    colnames(nbs) = nb_ids
    (spdep::mat2listw(nbs,row.names=rownames(nbs),style = 'B',zero.policy=T))
  }else{
    is.poly = substring(class(sf::st_geometry(poly))[1], 5) == 'POLYGON'
    if(is.poly){
      suppressWarnings({poly = sf::st_centroid(poly)})
    }
    nb = spdep::knearneigh(poly,k=k)
    #knb = asplit(FNN::get.knn(data=sf::st_coordinates(poly),k=k)$nn.index,1)
    #attr(knb,'class') = c('sgbp','list')
    ##Convert neighbors list to a nb object
    #attrs = attributes(knb)
    #knb = lapply(knb, function(i) { if(length(i) == 0L) 0L else sort(i) } )
    #attributes(knb) = attrs
    #attr(knb,'region.id') = rownames(poly)
    #class(knb) = "nb"
    spdep::nb2listw(spdep::knn2nb(nb,row.names = rownames(poly)),style='B',zero.policy=T)
  }
}


#' @title Delaunay triangulation neighbour identification
#' @name tri_nb
#' @description
#' Define neighbour relationships according to Delaunay triangulation network
#'
#' @details
#' based on tri.mesh function from tripack package.
#'
#' @param poly a sf data frame with polygon geometry
#' @param soi boolean, whether or not filter Delaunay network to keep sphere of influence (SOI) graph. Default is FALSE
#' @returns weighted neighbour matrix in a listw object from spdep package
#' @export
#' @family nb_funcs
tri_nb = function(poly,soi=F){
  if(!is.data.frame(poly)){
    nbs = lapply(poly,function(p){
      is.poly = substring(class(sf::st_geometry(p))[1], 5) == 'POLYGON'
      if(is.poly){
        suppressWarnings({p = sf::st_centroid(p)})
      }
      nb = tripack::tri.mesh(sf::st_coordinates(p)[,1],sf::st_coordinates(p)[,2])
      nb = tripack::neighbours(nb)
      attr(nb,'class') = c('sgbp','list')
      class(nb) = "nb"
      attr(nb,'region.id') = rownames(p)
      if(soi==T){
        nb = spdep::soi.graph(nb,sf::st_coordinates(p))
        nb = spdep::graph2nb(nb)
      }
      attr(nb,'region.id') = rownames(p)
      as(spdep::nb2listw(nb,style='B',zero.policy=T),'CsparseMatrix')
    })
    nb_ids = c()
    for(n in 1:length(nbs)){
      nb_ids = c(nb_ids,rownames(nbs[[n]]))
    }
    nbs = Matrix::bdiag(nbs)
    rownames(nbs) = nb_ids
    colnames(nbs) = nb_ids
    spdep::mat2listw(nbs,row.names=rownames(nbs),style = 'B',zero.policy=T)
  }else{
  #Define neighbours according to sphere of influence graph
    is.poly = substring(class(sf::st_geometry(poly))[1], 5) == 'POLYGON'
    if(is.poly){
      suppressWarnings({poly = sf::st_centroid(poly)})
    }
  nb = tripack::tri.mesh(sf::st_coordinates(poly)[,1],sf::st_coordinates(poly)[,2])
  nb = tripack::neighbours(nb)
  attr(nb,'class') = c('sgbp','list')
  class(nb) = "nb"
  attr(nb,'region.id') = rownames(poly)
  if(soi==T){
  nb = spdep::soi.graph(nb,sf::st_coordinates(poly))
  nb = spdep::graph2nb(nb)
  }
  attr(nb,'region.id') = rownames(poly)
  spdep::nb2listw(nb,style='B',zero.policy=T)
  }
}



#' @title expand neighbour links to N higher orders
#' @name nb_expand
#' @description
#' Given a Kandinsky object containing a neighbour network, this function will identify, for each cell/node, higher order neighbours that are distant n (maxorder) number of links from each other.
#' For more information, see nblag() function in the spdep package.
#'
#' @param seurat a Seurat object containing Kandinsky data slot
#' @param maxorder numeric, link distance considered to define higher order neighbours
#' @param cumul boolean, whether returning updated neighbour network with only higher order relationships (FALSE) or with the cumulative union of first and higher order neighbour links (TRUE)
#' @returns seurat object with updated 'nb' neighbours in Kandinsky data
#' @export
#' @family nb_funcs
nb_expand = function(seurat,maxorder=2,cumul=T){
  if(cumul == F){
    KanData(seurat,'nb') = spdep::nb2listw(spdep::nblag(KanData(seurat,'nb')$neighbours,maxlag=maxorder)[[maxorder]],style='B',zero.policy = T)
  }else{
    KanData(seurat,'nb') = spdep::nb2listw(spdep::nblag_cumul(spdep::nblag(KanData(seurat,'nb')$neighbours,maxlag=maxorder)),style='B',zero.policy = T)
  }
  return(seurat)
}

#' @title Update definition of Kandinsky neighbour network
#' @name nb_update
#' @description
#' Given a Seurat object containing Kandinsky data, this function will modify the neighbour network
#' generated by Kandinsky according to new conditions specified by the user.
#'
#' @param seurat a Seurat object containing Kandinsky data slot
#' @param sample_key character string specifying a variable stored in the Seurat object to use as sample/batch annotation. If not NULL, neighbour networks will be defined separately for each sample/batch. Default is NULL.
#' @param nb.method character string specifying the method to be used to create a `nb` neighbour object.
#' Must be one of the following:
#' 'Q': queen contiguity method,check for contact (not overlap) between any edge or side od two polygons (refers to the queen movement rule in chess). Currently only applicable for Visium/Visium-HD data
#' 'C': centroid-based method, use maximum centroid distance threshold to identify spot/cell neighbours
#' 'D': Delaunay triangulation
#' 'K': KNN method, define k closest neighbours to each spot/cell
#' 'M': membrane-based method, check for the occurrence of a physical contact/intersection within a distance threshold between cell boundaries. Not applicable in the case of Visium spots.
#' @param snap numeric, maximum accepted distance between Visium spots or Visium-HD bins to define contiguity relationships. Only Applied when `nb.method = Q`.
#' @param layers numeric, number of concentric contiguous layers to include in spot neighbourhood. Only Applied when `nb.method = Q`. Default is 1
#' @param d.max numeric, maximum distance accepted between polygon centroids to define neighbour relation
#' @param k numeric, number of nearest neighbour to select with knn algorithm
#' @param soi boolean, whether or not filter Delaunay network to keep sphere of influence (SOI) graph. Default is FALSE
#' @returns seurat object with updated 'nb' neighbours in Kandinsky data
#' @export
#' @family nb_funcs
nb_update = function(seurat=NULL,sample_key=NULL,nb.method=c('K','C','D','M','Q'),snap=NULL,layers=1,d.max=20,k=10,soi=F){
  if(KanData(seurat,'platform') == 'visium'){
    snap = snap %||% ((KanData(seurat,'spot_distance')*sqrt(3))/2)
  }
  if(KanData(seurat,'platform') == 'visium_hd'){
    snap = snap %||% (KanData(seurat,'spot_distance')*0.51*sqrt(2))
  }
  if(!is.null(sample_key)){
    message('Building separate neighbour networks for each "',sample_key,'" id')
    polys = populate_sf(seurat,vars=sample_key,return.seurat = F)
    polys = split(polys,polys[[sample_key]])
  }else{
    polys = KanData(seurat,'sf')
  }
  if(nb.method == 'K'){
    KanData(seurat,'nb') = knn_nb(polys,k = k)
    KanData(seurat,'nb.type') = paste0('K_',k)
  }else if(nb.method=='C'){
    KanData(seurat,'nb') = centroid_nb(polys,d.max=d.max)
    KanData(seurat,'nb.type') = paste0('C_',d.max)
  }else if(nb.method == 'D'){
    KanData(seurat,'nb') = tri_nb(polys,soi=soi)
    KanData(seurat,'nb.type') = paste0('D_',paste0('soi',soi))
  }else if(nb.method == 'M'){
    KanData(seurat,'nb') = membrane_nb(polys,d.max=d.max)
    KanData(seurat,'nb.type') = paste0('M_',d.max)
  }else if(nb.method == 'Q'){
    if(!(KanData(seurat,'platform') %in% c('visium','visium_hd'))){
      warning('neighbour method "Q" is currently only recommended with Visium/Visium-HD or grid-based data. Please consider a different method')
    }
    KanData(seurat,'nb') = queen_nb(polys,snap=snap,layers=layers)
    KanData(seurat,'nb.type') = 'Q'
  }else{
    stop('nb.method parameter must be either "Q", "C", "D", "K", or "M"')
  }
  return(seurat)
}

#' @title Get Kandinsky neighbourhood summary metrics
#' @name nb_summary
#' @description
#' Print a summary of neighbour counts and sizes detected by Kandinsky
#'
#' @param seurat a Seurat object containing Kandinsky data
#' @returns text summary of neighbour metrics
#' @export
nb_summary = function(seurat=NULL){
  if(KanData(seurat,'nb.type')=='Q'){
    method = 'Queen contiguity (Q)'
  }
  if(stringr::str_detect(KanData(seurat,'nb.type'),'C')){
    method = paste0('Centroid distance (C), d.max = ',strsplit(KanData(seurat,'nb.type'),split='_')[[1]][[2]])
  }
  if(stringr::str_detect(KanData(seurat,'nb.type'),'M')){
    method = paste0('Membrane distance (M), d.max = ',strsplit(KanData(seurat,'nb.type'),split='_')[[1]][[2]])
  }
  if(stringr::str_detect(KanData(seurat,'nb.type'),'K')){
    method = paste0('K-nearest neighbour (K), k = ',strsplit(KanData(seurat,'nb.type'),split='_')[[1]][[2]])
  }
  if(stringr::str_detect(KanData(seurat,'nb.type'),'D_')){
    method = paste0('Delaunay triangulation (D), soi = ',gsub('D_soi','',KanData(seurat,'nb.type')))
  }
  ncells = ncol(seurat)
  mat = as(KanData(seurat,'nb'),'CsparseMatrix')
  nlinks = length(mat@x)
  avg_links = round(nlinks/ncells,2)
  pct_links = nlinks/(ncells^2)
  nbs = quantile(Matrix::rowSums(mat))
  missing_nb = length(colnames(seurat)[which(Matrix::rowSums(mat) ==0)])
  message('Neighbour method: ',method)
  message('Total cells/spots: ',ncells)
  message('Total neighbour links: ',nlinks)
  message('Minimum number of neighbours per cell/spot: ',nbs[1])
  message('Mean number of neighbours per cell/spot: ',avg_links)
  message('Median number of neighbours per cell/spot: ',nbs[3])
  message('Max number of neighbours per cell/spot: ',nbs[5])
  rm(mat)
  message('Number of cells/spots with no neighbours: ',missing_nb)
}


#' @title Downsize neighbourhoods
#' @name nb_downsize
#' @description
#' Starting from a list of Kandinsky neighbourhoods, this function will select a random subset of cells/spots within each neighbourhood given a maximum expected neighbourhood size specified by the user
#' @details
#' Neighbourhoods will be used in the form of an adjacency matrix. If neighbourhood relationships are symmetric (as for neighbourhoods created with Queen contiguity and Centroid/Membrane distance methods), the resulting neighbourhood matrix must be symmetric even after the downsizing. Therefore, an additional step to ensure matrix symmetry is applied after the downsizing, and for this reason some neighbourhoods might still include a number of cells/spots higher than the user-provided threshold.
#'
#' @param seurat a Seurat object containing Kandinsky data slot
#' @param exp_links numeric, maximum number of expected cells/spots within each neighbourhood after the downsizing
#' @param seed numeric, random seed for reproducibility
#' @param keep_symmetric boolean, whether to force adjacency matrix symmetry after subsampling. If NULL, it will be set to TRUE for proximity-based neighbour methods (Q/C/M), and FALSE for remaining methods (K/D)
#' @return seurat object with updated neighbourhoods
#' @importFrom utils head
#' @export
nb_downsize = function(seurat=NULL,exp_links=6,seed=347548,keep_symmetric=NULL){
  mat= as(KanData(seurat,'nb'),'CsparseMatrix')
  nb.method = KanData(seurat,'nb.type')
  keep_symmetric = keep_symmetric %||% stringr::str_detect(nb.method, "Q|C|M")
  avg_size= median(Matrix::rowSums(mat))
  message('Current median number of neighbour links per cell is ',avg_size,'. Trying to reduce to ',exp_links,'...')
  #Forced symmetry tends to increase neighbour links. Adjust exp_links to account for that
  exp_links = ifelse(keep_symmetric,exp_links-1,exp_links)
  
  #Build edge list and random subsetting per row
  N = ncol(mat)
  set.seed(seed)
  to = rep(seq_len(N), times = diff(mat@p))
  from = mat@i + 1L
  edges = cbind(from, to)
  
  split_idx = split(seq_len(nrow(edges)), edges[,1])
  keep_idx = unlist(
    lapply(split_idx, function(ix){
      n = length(ix)
      k = min(exp_links, n)
      if (n <= k) ix else sample(ix, k)
    }),
    use.names= F
  )
  
  edges_sub = edges[keep_idx, , drop = FALSE]
  #Build sparse adjacency and symmetrize by addition
  mat = Matrix::sparseMatrix(i = edges_sub[,1], j = edges_sub[,2], x = 1L, dims = c(N, N))
  if(keep_symmetric){
    message("Enforcing symmetric neighbourhoods...")
    mat = as((mat | Matrix::t(mat)), "CsparseMatrix")
    #mat <- Matrix::drop0(mat)
    deg = Matrix::colSums(mat)
    if (any(deg > exp_links)) {
      message('Trying to trim degrees > ', exp_links+1,'...')
      mat <- Matrix::summary(mat)
      j_idx <- mat$j
      set.seed(seed)
      ## Randomly drop excess edges per column (without loops)
      order_j <- order(j_idx)
      j_idx <- j_idx[order_j]
      i_idx <- mat$i[order_j]
      run_len <- rle(j_idx)
      starts <- cumsum(c(1, head(run_len$lengths, -1)))
      ends <- cumsum(run_len$lengths)
      keep <- rep(TRUE, length(i_idx))
      for (k in seq_along(starts)) {
        s <- starts[k]; e <- ends[k]
        n <- e - s + 1
        if (n > exp_links) {
          drop <- sample.int(n, n - exp_links)
          keep[s + drop - 1L] <- FALSE
        }
      }
      # Subset edges
      i_idx <- i_idx[keep]
      j_idx <- j_idx[keep]
      
      # Rebuild symmetric matrix in one pass
      mat <- Matrix::sparseMatrix(i = i_idx, j = j_idx, 
                                  x = 1L, dims = c(N, N))
      # re-enforce symmetry once
      mat <- as((mat | Matrix::t(mat)), "CsparseMatrix")
    }
  }
  mat <- Matrix::drop0(mat)
  new_deg <- Matrix::rowSums(mat)
  message('Returning updated neighour network with a median number of ',median(new_deg),' links per cell')
  KanData(seurat,'nb') = suppressWarnings(spdep::mat2listw(mat,row.names=rownames(mat),style = 'B',zero.policy=T))
  return(seurat)
}

