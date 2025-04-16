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
#' @returns weighted neighbour matrix in a listw object from spdep package
#' @export
#' @family nb_funcs
queen_nb = function(poly,snap=NULL){
  message('Queen neighbour method is designed to work with Visium / Visium HD data only')
  if(levels(droplevels(sf::st_geometry_type(poly))) != 'POLYGON'){
    stop('Input geometries should be polygons')
    #poly = sf::st_buffer(poly,dist=buffer)
  }
  nb = spdep::poly2nb(poly,snap=snap,row.names = poly[['spot_ID']])
  nb = spdep::nb2listw(nb,style = 'B',zero.policy=T)
  return(nb)
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
  if(levels(droplevels(sf::st_geometry_type(poly))) == 'POLYGON'){
    suppressWarnings({poly = sf::st_centroid(poly)})
  }
  nb = spdep::dnearneigh(poly,d1=0,d2=d.max,row.names = rownames(poly))
  nb = spdep::nb2listw(nb,style='B',zero.policy=T)
  return(nb)
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
knn_nb = function(poly,k=20){
  if(levels(droplevels(sf::st_geometry_type(poly))) == 'POLYGON'){
    suppressWarnings({poly = sf::st_centroid(poly)})
  }
  nb = spdep::knearneigh(poly,k=k)
  nb = spdep::knn2nb(nb,row.names = rownames(poly))
  nb = spdep::nb2listw(nb,style='B',zero.policy=T)
  return(nb)
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
  #Define neighbours according to sphere of influence graph
  if(levels(droplevels(sf::st_geometry_type(poly))) == 'POLYGON'){
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
  nb = spdep::nb2listw(nb,style='B',zero.policy=T)
  return(nb)
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
#' @param nb.method character string specifying the method to be used to create a `nb` neighbour object.
#' Must be one of the following:
#' 'Q': queen contiguity method,check for contact (not overlap) between any edge or side od two polygons (refers to the queen movement rule in chess). Currently only applicable for Visium/Visium-HD data
#' 'C': centroid-based method, use maximum centroid distance threshold to identify spot/cell neighbours
#' 'D': Delaunay triangulation
#' 'K': KNN method, define k closest neighbours to each spot/cell
#' 'M': membrane-based method, check for the occurrence of a physical contact/intersection within a distance threshold between cell boundaries. Not applicable in the case of Visium spots.
#' @param snap numeric, maximum accepted distance between Visium spots or Visium-HD bins to define contiguity relationships. Only Applied when `nb.method = Q`.
#' @param d.max numeric, maximum distance accepted between polygon centroids to define neighbour relation
#' @param k numeric, number of nearest neighbour to select with knn algorithm
#' @param soi boolean, whether or not filter Delaunay network to keep sphere of influence (SOI) graph. Default is FALSE
#' @returns seurat object with updated 'nb' neighbours in Kandinsky data
#' @export
#' @family nb_funcs
nb_update = function(seurat=NULL,nb.method=c('K','C','D','M','Q'),snap=NULL,d.max=20,k=10,soi=F){
  if(nb.method == 'K'){
    KanData(seurat,'nb') = knn_nb(KanData(seurat,'sf'),k = k)
    KanData(seurat,'nb.type') = paste0('K_',k)
  }else if(nb.method=='C'){
    KanData(seurat,'nb') = centroid_nb(KanData(seurat,'sf'),d.max=d.max)
    KanData(seurat,'nb.type') = paste0('C_',d.max)
  }else if(nb.method == 'D'){
    KanData(seurat,'nb') = tri_nb(KanData(seurat,'sf'),soi=soi)
    KanData(seurat,'nb.type') = paste0('D_',paste0('soi',soi))
  }else if(nb.method == 'M'){
    KanData(seurat,'nb') = membrane_nb(KanData(seurat,'sf'),d.max=d.max)
    KanData(seurat,'nb.type') = paste0('M_',d.max)
  }else if(nb.method == 'Q'){
    if(KanData(seurat,'platform') == 'visium'){
      snap = snap %||% (KanData(seurat,'spot_distance')/100)*40
    }else if(KanData(seurat,'platform') == 'visium_hd'){
      snap = snap %||% KanData(seurat,'spot_distance')
    }else{
      stop('neighbour method "Q" is currently only compatible with Visium/Visium-HD data. Please choose a different method')
    }
    KanData(seurat,'nb') = queen_nb(KanData(seurat,'sf'),snap=snap)
    KanData(seurat,'nb.type') = 'Q'
  }else{
    stop('nb.method parameter must be either "Q", "C", "K", or "M"')
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
    method = paste0('Delaunay triangulation (D), soi = ',gsub('soi','',KanData(seurat,'nb.type')))
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


