#return Seurat with updated nb
#NB. the subset query should refer to objects related to the nb region.ids, like Seurat metadata or kandinsky sf slots

#subset_kandinsky_nb = function(x,query=""){
# KanData(x,'nb') = subset(KanData(x,'nb'),eval(parse(text=query)))
#return(x)}


#' @title Add variables to sf object stored in Kandinsky data
#' @name populate_sf
#' @description Add information to the Kandinsky 'sf' slot. Based on Seurat::FetchData function
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @param vars character vector containing variable names stored in Seurat object to be copied into the KanData `sf` object
#' @param assay name of the Seurat assay to be used to search for variables of interest
#' @param layer name of the Seurat layer to be used to search for variables of interest
#' @param return.seurat boolean, set to FALSE by default. It indicates whether returning Seurat object with updated KanData `sf` containing the new variables (TRUE), or the updated KanData `sf` object only
#' @returns if return.seurat is set to TRUE, a Seurat object with updated KanData `sf` containing the new variables. if return.seurat is set to FALSE, only the updated KanData `sf` object will be returned
#' @export
populate_sf = function(seurat=NULL,vars=NULL,assay=DefaultAssay(seurat),layer='data',return.seurat=F){
  missing = setdiff(vars,colnames(KanData(seurat,'sf')))
  if(length(missing) >0){
    if(return.seurat==F){
      cbind(KanData(seurat,'sf')[colnames(seurat),],
            Seurat::FetchData(seurat,vars=missing,assay=assay,layer=layer,clean=F))
    }else{
      KanData(seurat,'sf') = cbind(KanData(seurat,'sf')[colnames(seurat),],
                                   Seurat::FetchData(seurat,vars=missing,assay=assay,layer=layer,clean=F))
      return(seurat)
    }
  }else{
    if(return.seurat==F){
      return(KanData(seurat,'sf'))
    }else{
      return(seurat)
    }
  }
}



#read Vizgen cell boundaries parquet files
#read_vizgen_polygons = function(file,zlayer=0){
#  file = normalizePath(file)
#  arrow::open_dataset(file) %>% select(ID,EntityID,ZIndex,Geometry,ZLevel) %>%
#    filter(ZIndex %in% zlayer) %>%
#    collect() %>%
#    sf::st_as_sf() %>%
#    st_cast(.,to='POLYGON',warn=F)
#}


#' @title Rescale Kandinsky polygon xy coordinates
#' @name scale_coords
#' @description
#' Apply a shift or a scaling factor to either one or both x and y polygon coordinates stored in Kandinsky `sf` slot
#'
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @param xfact numeric value, shifting or scaling factor to apply to x coordinates. Default value is 0 for shifting, and 1 for scaling transformations
#' @param yfact numeric value, shifting or scaling factort o apply to y coordinates. Default value is 0 for shifting, and 1 for scaling transformations
#' @param fun character indicating the operator to use for coordinate transformation. Must be one of the following:
#' `"*"`: multiply x/y coordinates by xfact/yfact (scaling)
#' `"+"`: add xfact/yfact to x/y coordinates (shifting)
#' `"/"`: divide x/y coordinates by xfact/yfact (scaling)
#' `"-"`: subtract xfact/yfact to x/y coordinates (shifting)
#' @returns Seurat object with updated polygon coordinates in Kandinsky `sf` slot
#'
#' @export
scale_coords = function(seurat,xfact=NULL,yfact=NULL,fun=c('*','+','/','-')){
  if(fun %in% c('*','/')){
    xfact= xfact %||% 1
  }else{
    xfact=xfact %||% 0
  }
  if(is.null(yfact)){
    if(fun %in% c('*','/')){
      yfact=1
    }else{
      yfact=0
    }
  }
  trans = paste0("sf::st_geometry(KanData(seurat,'sf')) ",fun," matrix(c(xfact,yfact),ncol=2)")
  sf::st_geometry(KanData(seurat,'sf')) = eval(parse(text=trans))
  return(seurat)
}

#' @name stitch_samples
#' @title Merge coordinate reference systems from independent samples
#' @description
#' Combine samples coming from independent experiments by stitching cells/spots coordinates through x/y coordinate shifting into a single grid
#' @param data a data frame containing x and y coordinates into separate columns
#' @param x name of the x coordinate column
#' @param y name of the y coordinate column
#' @param sample.id character string specifying the name of the column storing sample identifiers
#' @param max_rows numeric, maximum number of samples to be disposed on a same column of the final sample grid. If set to `NULL`, maximum number of rows will be automatically calculated to obtain a square grid
#' @param max_cols numeric, maximum number of samples to be disposed on a same row of the final sample grid. If set to `NULL`, maximum number of columns will be automatically calculated to obtain a square grid
#' @returns data frame object with new stitched 'x_global' and 'y_global' coordinate columns
#' @export
stitch_samples = function(data = NULL,x=NULL,y=NULL,sample.id=NULL,max_rows=NULL,max_cols=NULL){
  if(is.factor(data[[sample.id]])){data[[sample.id]] = as.character(data[[sample.id]])}
  ids_list = unique(data[[sample.id]])
  if(length(ids_list) <=2){
    grid_order = matrix(c(ids_list,rep(0,2-length(ids_list))),ncol=2,nrow=1)
  }else{
    if(is.null(max_rows) | is.null(max_cols)){
      max_dim = ceiling(sqrt(length(unique(data[[sample.id]]))))
      n_cells = max_dim^2
      grid_order = matrix(c(ids_list,rep(0,n_cells-length(ids_list))),ncol=max_dim,nrow=max_dim)
    }else{
      n_cells = max_rows*max_cols
      grid_order = matrix(c(ids_list,rep(0,n_cells-length(ids_list))),ncol=max_cols,nrow=max_rows)
    }
  }
  grid_order = reshape2::melt(grid_order)
  grid_order = grid_order[grid_order$value != "0",]
  colnames(grid_order) = c('x','y','id')

  colnames(data)[colnames(data) == x] = 'my_x'
  colnames(data)[colnames(data) == y] = 'my_y'
  data_list = split(data,f=data[[sample.id]])

  coords = t(sapply(data_list,function(x){return(c(x1 = min(x$my_x),y1=min(x$my_y),x2 = max(x$my_x),y2=max(x$my_y)))}))
  coords = as.data.frame(coords)

  x_offset = max(coords$x2 - coords$x1)
  y_offset = max(coords$y2 - coords$y1)

  data_list_shift = lapply(data_list,function(x){
    grid_x = grid_order[which(grid_order$id == unique(x[[sample.id]])),]$x
    grid_y = grid_order[which(grid_order$id == unique(x[[sample.id]])),]$y
    x$x_global = x$my_x + (1.75*x_offset* (grid_x -1))
    x$y_global = x$my_y - (1.75*y_offset* (grid_y -1))
    return(x)
  })
  data_shift = purrr::reduce(data_list_shift,rbind)
  return(data_shift)
}

#' @title resample single cells in Seurat object
#' @name resample_cells
#' @description
#' Apply a random resample of cell identifier in Seurat object. Original cell type proportion in Seurat dataset will be preserved
#' when non-spatial sampling strategy is applied. With spatial resampling, rare/sparse cell types will be less penalised from resampling
#' compared to more abundant and spatially dense cell type, with a sampling proportion that tends to be
#' inversely correlated with the initial cell type relative abundance
#'
#' @param seurat Seurat object
#' @param label character string indicating meta data variable containing cell type annotation
#' @param spatial boolean, whether or not applying a spatially-informed resampling
#' @param maxcells numeric, ideal number of cells to resample.
#' When `spatial` is true, an optimal number of cells as close as possible to `maxcells` argument
#' will be sampled by testing multiple resampling resolution parameters
#' @param seed numeric, random seed for reproducibility
#' @param return.seurat boolean, whether returning resampled IDs (FALSE) or the whole seurat object already filtered with resampled IDs (TRUE)
#' @param update.kandinsky boolean, whether or not update Kandinsky data after resampling. Default is TRUE. Only applied when `return.seurat` is TRUE
#' @returns vector of resampled cell IDs or Seurat object with resampled cell subset
#'
#' @export
resample_cells = function(seurat=NULL,label='cell_types',spatial=T,maxcells=10000,seed=347548,return.seurat=T,update.kandinsky=T){
  seurat@meta.data$ids = rownames(seurat@meta.data)
  if(spatial==F){
    all.prop = maxcells/nrow(seurat@meta.data)
    props = prop.table(table(seurat@meta.data[[label]])) %>% as.data.frame()
    colnames(props) = c(label,'proportion')
    min.size =
      list = seurat@meta.data[,c('ids',label)]
    set.seed(seed)
    new_ids = lapply(props[[label]],function(l){
      sample = sample(list[list[[label]] == l,][['ids']],size=round(maxcells*props[props[[label]]==l,]$proportion))
      return(list[list[['ids']] %in% sample,])
    })
    new_ids=rownames(purrr::reduce(new_ids,rbind))
  }else{
    bbox = sf::st_bbox(KanData(seurat,'sf'))
    opt = expand.grid(seq(50,150,10),seq(1,3,1))
    opt = opt[order(opt$Var1),]
    mindist = Inf
    best = NULL
    gap = 0
    message('Starting grid search for best sampling resolution...')
    for(n in seq_len(nrow(opt))){
      grid = (sf::st_make_grid(bbox,n=opt[n,1],square=T))
      #xr = bbox$xmax-bbox$xmin
      #xbin = round(xr/sp.grid.size)
      #yr = bbox$ymax - bbox$ymin
      #ybin = round(yr/sp.grid.size)
      #grid = (sf::st_make_grid(bbox,n=c(xbin,ybin),square=T,))
      poly.grid = sf::st_as_sf(grid)
      poly.grid$ID = rownames(poly.grid)

      neighbors = sf::st_intersects(poly.grid,suppressWarnings({sf::st_centroid(KanData(seurat,'sf'))}))
      melted_nb = reshape2::melt(neighbors)
      melted_nb$cell_ID = sf::st_drop_geometry(KanData(populate_sf(seurat,vars='ids',return.seurat = T),'sf'))[melted_nb$value,'ids']
      melted_nb$cell_types = sf::st_drop_geometry(KanData(populate_sf(seurat,vars=label,return.seurat = T),'sf'))[melted_nb$value,label]
      set.seed(seed)
      sampled_nb = melted_nb %>%
        dplyr::group_by(.data[["L1"]]) %>%
        dplyr::group_modify(~{
          all_labels = unique(.x[["cell_types"]])
          sampled = lapply(all_labels,function(l){
            ids = .x[.x$cell_types == l,]$cell_ID
            return(suppressWarnings({.x[.x$cell_ID == sample(ids,size=min(opt[n,2],nrow(.x[.x$cell_types == l,,drop=F]))),]}))
          })
          sampled=purrr::reduce(sampled,rbind)
          return(sampled)
        })
      dist = abs(maxcells-nrow(sampled_nb))
      if(dist < mindist){
        best = opt[n,,drop=F]
        mindist=dist
      }else{
        gap = gap+1
      }
      if(gap >=5){
        message('optimal number of resampled cells: ',abs(mindist-maxcells))
        break
      }
    }
    message('Resampling ',best[,2],' cells for each cell type across ',best[,1], ' grid sectors')
    grid = (sf::st_make_grid(bbox,n=best[,1],square=T))
    #xr = bbox$xmax-bbox$xmin
    #xbin = round(xr/sp.grid.size)
    #yr = bbox$ymax - bbox$ymin
    #ybin = round(yr/sp.grid.size)
    #grid = (sf::st_make_grid(bbox,n=c(xbin,ybin),square=T,))
    poly.grid = sf::st_as_sf(grid)
    poly.grid$ID = rownames(poly.grid)

    neighbors = sf::st_intersects(poly.grid,suppressWarnings({sf::st_centroid(KanData(seurat,'sf'))}))
    melted_nb = reshape2::melt(neighbors)
    melted_nb$cell_ID = sf::st_drop_geometry(KanData(populate_sf(seurat,vars='ids',return.seurat = T),'sf'))[melted_nb$value,'ids']
    melted_nb$cell_types = sf::st_drop_geometry(KanData(populate_sf(seurat,vars=label,return.seurat = T),'sf'))[melted_nb$value,label]
    set.seed(seed)
    sampled_nb = melted_nb %>%
      dplyr::group_by(.data[["L1"]]) %>%
      dplyr::group_modify(~{
        all_labels = unique(.x[["cell_types"]])
        sampled = lapply(all_labels,function(l){
          ids = .x[.x$cell_types == l,]$cell_ID
          return(suppressWarnings({.x[.x$cell_ID == sample(ids,size=min(best[,2],nrow(.x[.x$cell_types == l,,drop=F]))),]}))
        })
        sampled=purrr::reduce(sampled,rbind)
        return(sampled)
      })
    #sampled_nb = ddply(melted_nb,.(L1),function(x){
    # all_labels = unique(x$cell_types)
    #sampled = lapply(all_labels,function(l){
    # ids = x[x$cell_types == l,]$cell_ID
    #return(x[x$cell_ID == sample(ids,size=1),])
    #})
    #sampled=purrr::reduce(sampled,rbind)
    #return(sampled)
    #})
    new_ids = sampled_nb$cell_ID
  }
  if(return.seurat==T){
    if(update.kandinsky ==T & any(names(seurat@tools) == 'kandata')){
      return(update_kandinsky(suppressWarnings(seurat[,new_ids])))
    }else{
      return(suppressWarnings(seurat[,new_ids]))
    }
  }else{
    return(new_ids)
  }
}

#' @name global_univ_spatcor
#' @title Compute global Moran's I spatial autocorrelation statistic
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @param var character string or vector specifying the variable(s) to use to compute Moran's I statistics
#' @param sim number of Monte Carlo simulations to be run for estimating Moran coefficients significance
#' @param lag integer value indicating the extent of cell/spot neighbours to be considered to calculate Moran statistics. `lag = 1` indicates that only 1st order neighbours will be considered, while `lag=2` indicates that all neighbours of each 1st order neighbour will be also considered for each spot/cell, and so on.
#' @param alt a character string specifying the alternative hypothesis, Must be one of the following: "two.sided", "greater" (default), or "less".
#' @param seed numeric, random seed for reproducibility. Default is set to 347548
#' @returns a data.frame reporting Moran's I spatial autocorrelation statistics with the associated p value for each selected variable.
#' @export
global_univ_spatcor = function(seurat,var=NULL,sim=49,lag=1,alt='greater',seed=347548){
  set.seed(seed)
  if(lag > 1){
    listw = KanData(nb_expand(seurat,maxorder=lag,cumul=T),'nb')
  }else{
    listw = KanData(seurat,'nb')
  }
  if(length(var) > 1){
    moran=lapply(var,function(v){
      mor = spdep::moran.mc(FetchData(seurat,vars=v,layer='data',clean=F)[,v],listw,nsim=sim,zero.policy=T,alternative=alt)
      return(data.frame(broom::tidy(mor)))
      #return(data.frame(variable = v,Moran_I=mor$statistic,mc_rank=mor$parameter,pval=mor$p.value))
    })
    do.call('rbind',(moran)) %>% dplyr::mutate(variable=var)
  }else{
    mor = spdep::moran.mc(FetchData(seurat,vars=var,layer='data',clean=F)[,var],listw,nsim=sim,zero.policy=T,alternative=alt)
    return(data.frame(broom::tidy(mor)))
    #return(data.frame(variable = var,Moran_I=(mor$statistic),mc_rank=(mor$parameter),pval=(mor$p.value)))
  }
}


#' @name global_biv_spatcor
#' @title Compute global bivariate spatial correlation
#' @description
#' Compute bivariate global spatial correlation. You can use either Lee's L (default) or bivariate Moran's I statistics.
#'
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @param var1 character string or vector specifying the first variable(s) to use to compute bivariate spatial correlation.
#' @param var2 character string or vector specifying the second variable(s) to use to compute bivariate spatial correlation
#' @param sim number of Monte Carlo simulations to be run for estimating Moran/Lee coefficients significance
#' @param lag integer value indicating the extent of cell/spot neighbours to be considered to calculate Moran/Lee statistics. `lag = 1` indicates that only 1st order neighbours will be considered, while `lag=2` indicates that all neighbours of each 1st order neighbour will be also considered for each spot/cell, and so on.
#' @param method character string specifying which method to use to estimate spatial correlation. Must be one of the following: "lee" (default), "moran".
#' @param alt a character string specifying the alternative hypothesis, Must be one of the following: "two.sided" (default), "greater" or "less".
#' @param layer character string indicating which Seurat layer will be considered to search for the variable of interest
#' @param seed numeric, random seed for reproducibility. Default is set to 347548
#' @returns a data.frame reporting the chosen spatial statistic coefficient with the associated p value for each possible combination between variables specified with `var1` `var2` parameters.
#' @export
global_biv_spatcor = function(seurat,var1=NULL,var2=NULL,sim=49,lag=1,method=c('lee','moran'),alt='two.sided',layer='data',seed=347548){
  if(is.null(method) | length(method) > 1){
    method='lee'
  }
  if(lag > 1){
    listw = KanData(nb_expand(seurat,maxorder=lag,cumul=T),'nb')
  }else{
    listw = KanData(seurat,'nb')
  }
  all_combs = expand.grid(var1,var2) %>% dplyr::mutate_all(as.character)
  if(method=='moran'){
    set.seed(seed)
    res=lapply(seq_len(nrow(all_combs)),function(x){
      Var1 = all_combs$Var1[x]
      Var2 = all_combs$Var2[x]
      res=spdep::moran_bv(x=Seurat::FetchData(seurat,vars=Var1,layer=layer,clean=F)[,c(Var1)],y=Seurat::FetchData(seurat,vars=Var2,layer=layer,clean=F)[,c(Var2)],listw,nsim=sim)
      pval = mean(abs(res$t) > abs(res$t0))
      res = data.frame(var1=Var1,var2=Var2,moran_stat=res$t0,moran_pval = pval)
      return(res)
    })
  }else if(method=='lee'){
    set.seed(seed)
    res=lapply(seq_len(nrow(all_combs)),function(x){
      Var1 = all_combs$Var1[x]
      Var2 = all_combs$Var2[x]
      res=spdep::lee.mc(x=Seurat::FetchData(seurat,vars=Var1,layer=layer,clean=F)[,c(Var1)],y=Seurat::FetchData(seurat,vars=Var2,layer=layer,clean=F)[,c(Var2)],listw,nsim=sim,zero.policy=T,alternative=alt)
      res = data.frame(var1=Var1,var2=Var2,lee_stat=res$statistic,lee_pval = res$p.value)
      return(res)
    })
  }
  res = purrr::reduce(res,rbind)
  return(res)
}

#' @title Get neighbourhood aggregated gene expression count matrix
#' @name get_nbcounts
#' @description
#' For each cell, aggregate gene expression counts from all cells belonging to that cell's neighbourhood.
#'
#' @param seurat a Seurat object containing Kandinsky data slot
#' @param label character string indicating meta data variable containing cell type annotation
#' @param which character vector indicating for which cell type creating the neighbourhood expression matrix. Cells belonging to the same class(es) specified through `which` will not be included in the expression count aggregation.
#' @details
#' When `label` or `which` parameters are set to NULL, this function will create an expression count neighbourhood matrix considering all cells included in the dataset.
#' Otherwise, the aggregated count matrix will be created for only cells belonging to the cell class(es) specified through the `which` argument, and cells from that same class will not be considered for the gene expression count aggregation across neighbouring cells.
#' Cell class annotation that will be considered for the analysis needs to be specified through the `label` argument.
#'
#' @returns sparse cell (rows) X gene (columns) count matrix.  Each cell count profile is the result of the aggregation of all gene expression counts found across its neighbouring cells.
#'
#' @export
get_nbcounts = function(seurat=NULL,label=NULL,which=NULL){
  nb = as(KanData(seurat,'nb'),'CsparseMatrix')
  if(!is.null(label) & !is.null(which)){
    nb[,seurat@meta.data[[label]] %in% which]=0
    nb=nb[seurat@meta.data[[label]] %in% which,]
  }else{
    warning('No cell type label of interest has been set. The function will consider all cells included within the dataset')
  }
  #nb = Matrix::Diagonal(x=rep(1, nrow(nb))) %*% nb
  #nb@x[nb@x==0] <- 1
  nb = nb %*% Matrix::t(LayerData(seurat,layer='counts'))
  if(!is.null(label) & !is.null(which)){
    rownames(nb) = colnames(seurat)[seurat@meta.data[[label]] %in% which]
  }else{
    rownames(nb) = colnames(seurat)
  }
  return(nb)
}

#################################
###########COSMX UTILS###########
#################################

#' @title Add transcript file path to Kandinsky data in a Seurat object
#' @name AddTxPath
#' @description
#' Add the full path of CosMx/Xenium/Merscope transcript csv data to Kandinsky data
#' @param seurat  a Seurat object containing Kandinsky data (`KanData()`)
#' @param path character string indicating the path for the transcript coordinates file from CoxMx/Xenium/Merscope platforms
#' @return Seurat object with new `tx` slot added to Kandinsky data
#' @export
AddTxPath = function(seurat,path){
  KanData(seurat,'tx') = normalizePath(path)
  return(seurat)
}


#' @title Query and filter Kandinsky transcript file
#' @name GetTxData
#' @description
#' Load and query CosMx/Xenium/Merscope transcript files as `arrow` datasets
#' @details
#' Transcript files are queried and filtered using `arrow` and `dplyr` functionalities to reduce the memory requirements
#' even in the case of datasets with hundreds of millions of transcripts. In the current version,
#' transcript files can be filtered according to cell identifiers, fov identifiers, and genes/probes of interest.
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @param features character vector indicating genes/probes to select from the list of transcript. If `NULL`, no filter will be applied on genes/probes
#' @param fovs character vector indicating the identifiers of fields of view (FOVs) of interest.  If `NULL`, no filter will be applied on FOVs
#' @param cells character vector indicating the identifiers of the single cells of interest.  If `NULL`, no filter will be applied on single cell identifiers
#' @param tx_path character string indicating the path for the transcript coordinates file from CoxMx/Xenium/Merscope platforms. If `NULL`, the function will expect to find the file path stored within the `tx` Kandinsky data slot
#' @returns data frame containing the rows selected from the input transcript file
#' @importFrom arrow schema int64 int32 string open_dataset
#' @importFrom magrittr %>%
#' @export
GetTxData = function(seurat=NULL,features=NULL,fovs=NULL,cells=NULL,tx_path = NULL){
  if(is.null(KanData(seurat,'tx'))){
    seurat = AddTxPath(seurat,tx_path)
  }
  sch = KanData(seurat,'platform')
  if(sch=='cosmx'){
    schema = schema(fov=int64(),cell_ID=string(),cell=string(),x_local_px=double(),y_local_px=double(),x_global_px=double(),y_global_px=double(),z=int64(),target=string(),CellComp=string())
  }else if(sch=='merscope'){
    schema = schema(row=int64(),barcode_id = string(),global_x = double(),global_y=double(),global_z=double(),x=double(),y=double(),fov=int64(),gene=string(),transcript_id = string(),cell_id=int64())
  }else if(sch=='xenium'){
    schema = schema(transcript_id = string(),cell_id=string(),overlaps_nucleus = int64(),feature_name=string(),x_location=double(),y_location=double(),z_location=double(),qv=double(),fov_name=string(),nucleus_distance=double(),codeword_index=double())
  }
  if(stringr::str_detect(KanData(seurat,'tx'),'parquet')){
    tx = arrow::open_dataset(KanData(seurat,'tx'),format='parquet',schema=schema)
  }else{
    tx = arrow::open_dataset(KanData(seurat,'tx'),format='csv',schema=schema,skip = 1)
  }
  if(!is.null(features)){
    if(sch=='cosmx'){
      tx = tx %>% dplyr::filter(.data$target %in% features) %>% dplyr::compute()
    }else if(sch=='merscope'){
      tx = tx %>% dplyr::filter(.data$gene %in% features) %>% dplyr::compute()
    }else if(sch=='xenium'){
      tx = tx %>% dplyr::filter(.data$feature_name %in% features) %>% dplyr::compute()
    }
  }
  if(!is.null(fovs)){
    if(sch != 'xenium'){
      tx = tx %>% dplyr::filter(.data$fov %in% fovs) %>% dplyr::compute()
    }else{
      tx = tx %>% dplyr::filter(.data$fov_name %in% fovs) %>% dplyr::compute()
    }
  }
  if(!is.null(cells)){
    if(sch=='cosmx'){
      tx = tx %>% dplyr::mutate(cell_ID = paste0(.data$fov,'_',.data$cell_ID)) %>% dplyr::filter(.data$cell_ID %in% cells) %>% dplyr::compute()
    }else if(sch=='merscope'){
      tx = tx %>% dplyr::filter(.data$barcode_id %in% cells) %>% dplyr::compute()
    }else if(sch=='xenium'){
      tx = tx %>% dplyr::filter(.data$cell_id %in% cells) %>% dplyr::compute()
    }
  }
  if(sch=='cosmx'){
    dplyr::mutate(x_global_um=.data[["x_global_px"]]*0.12028,y_global_um=.data[["y_global_px"]]*0.12028) %>%
      dplyr::collect(tx)
  }else{
    dplyr::collect(tx)
  }
}



#' @title Create sf geometry object from CosMx Seurat data
#' @name smi2sf
#' @description
#' Create a sf data frame object (package `sf`) using CosMx cell coordinates
#' @details
#' This function will use CosMx cell boundary coordinates  to create cell polygon geometries using the R package `sf`
#'
#' @param seurat a Seurat object containing CosMx data
#' @param poly_file character string indicating the path for the polygon coordinates file. This parameter is ignored if Seurat object has been created using the function `prepare_cosmx_seurat()`
#' @param id name of the variable to use as single cell identifier
#' @param return.seurat boolean, whether returning the input Seurat file with the sf object stored in the Kandinsky `sf` slot or Seurat `misc` slot (`TRUE`) or the nb object alone (`FALSE`). Default value is `TRUE`
#' @returns Seurat object with the new sf object stored in the `misc` slot, or in the Kandinsky `sf` slot, or the sf object alone, depending on the `return.seurat` parameter
#'
#' @export
smi2sf = function(seurat=NULL,poly_file=NULL,id='cell_ID',return.seurat=F){
  if(is.null(poly_file)){
    poly = sfheaders::sf_point(GetTissueCoordinates(seurat),x='y',y='x',keep=T)
    rownames(poly) = sf::st_drop_geometry(poly)[,'cell']
  }else{
    if(!is.object(poly_file)){
      coords = arrow::read_csv_arrow(poly_file)
    }else{
      coords = poly_file
    }
    if(length(
      intersect(unique(coords[,id]),colnames(seurat))
    ) != length(colnames(seurat))){
      warning('Be sure that the id parameter also corresponds to the cell ids in use in Seurat object')
    }
    poly = sfheaders::sf_polygon(coords,x='x_global_px',y='y_global_px',polygon_id = id,keep=T)
    rownames(poly) = sf::st_drop_geometry(poly)[,id]
  }
  if(return.seurat ==T){
    if(!is.null(KanData(seurat))){
      KanData(seurat,'sf') = poly[colnames(seurat),]
    }else{
      seurat@tools$sf = poly[colnames(seurat),]
      message("New sf data stored as tool 'sf' data (i.e., seurat@tools$sf)")
    }
  }else{
    return(poly[colnames(seurat),])
  }
}


#' @title Read CosMx fov position file
#' @name read_fovfile
#' @description
#' Read fov position csv file generated by NanoString CosMx platform
#' @details
#' The function will create a polygon for each CosMx field of view reported in the input file, with coordinates
#' measured in micrometers. The function assumed a standardized FOV size of 0.51x0.51 millimeters, and it does not account
#' for any possible position offset compared to cell polygon coordinates.
#' The output FOV polygons can be used in principle as a reference to arrange CosMx tif/png images as background for data visualization
#' @param fov_path character string indicating the path for the fov coordinates file from CoxMx platform
#' @param buffer numeric value (micrometers) indicating the radial expansion that will be applied to FOV centroid to reconstruct the FOV full size polygon.
#' It should not be modified, as this is the correct parameter to get the standard CosMx FOV size,
#' @param fov_ids vector of FOV identifier to be considered for the final polygon output. If set to `NULL`, all the FOVs reported in the input file will be used
#' @returns a sf data frame containing a geometry column with all FOV polygon coordinates
#' @export
read_fovfile = function(fov_path=NULL,buffer=255,fov_ids=NULL){
  fov = arrow::read_csv_arrow(fov_path)
  if(any(stringr::str_detect(colnames(fov),'fov'))){
    fov$FOV = fov$fov
  }
  if(!is.null(fov_ids)){
    fov = fov[fov$FOV %in% fov_ids,]
  }
  cols = colnames(fov)
  if(any(stringr::str_detect(cols,'x_global_px|y_global_px'))){
    fov = fov %>% dplyr::mutate(x_global_um = .data[["x_global_px"]]*0.12028,
                                y_global_um = .data[["y_global_px"]]*0.12028)
    fov$x_global_um = fov$x_global_um+255
    fov$y_global_um = fov$y_global_um-255
    fov = sfheaders::sf_point(fov,x='x_global_um',y='y_global_um',keep = T)
  }else{
    #X and Y coordinates reported in this file should correspond to the top left angle of the FOV.
    #Add and subtract half of FOV width to X and Y coords, respectively, to obtain centroid coordinates
    fov$centerX_mm = fov$X_mm + 0.255
    fov$centerY_mm = fov$Y_mm - 0.255
    #Convert mm in um coords
    fov = fov %>% dplyr::mutate(centerX_um = .data[["centerX_mm"]]*1000,
                                centerY_um = .data[["centerY_mm"]]*1000)
    fov = sfheaders::sf_point(fov,x='centerX_um',y='centerY_um',keep = T)
    fov[,c('centerX_um','centerY_um')] = suppressWarnings(sf::st_coordinates(
      sf::st_centroid(fov)))
  }
  #Expand fov centroid up to the real cosmx fov size (0.51x0.51 mm)
  fov = sf::st_buffer(fov,dist = buffer,endCapStyle = "SQUARE")
  fov$fov_ID = stringr::str_pad(fov$FOV,width=3,pad='0')
  fov$fov_ID = paste0('F',fov$fov_ID)
  fov$fov_ID_v2 = gsub('F','F00',fov$fov_ID)
  #  fov = sf::st_geometry(fov) %>%
  #    sf::st_union() %>%
  #    sf::st_cast('POLYGON') %>%
  #    sf::st_as_sf() %>%
  # st_make_valid() %>%
  #  st_union()
  return(fov)
}

#' @title Visualize CosMx IFs staining for a list of FOVs
#' @name CosmxFovPlot
#' @description
#' The function return a plot visualizing the immunofluoresence staining images associated with a list of FOVs provided by the user.
#' When no list is provided, all FOVs represented in the input file will be plotted.
#' @details
#' This function is designed to work with CosMx data only.
#'
#'
#' @param fovs data frame containing FOV geometries (square) created via the read_fovfile() function from Kandinsky
#' @param which vector of FOV identifier to be considered for plotting
#' @param img_dir character stringr specifying the full path of the directory containing Cosmx FOV images (usually named "CellComposite", but other type of images can be selected, such as 'CellOverlay','CellLabels', or 'Morphology2D')
#' @param maxres image resolution threshold. Will be applied to each image separately. When several images are selected for plotting, it is advisable to set a low maxres value (within the 500-800 range) to speed up the processing time
#' @returns plot of concatenated FOV images
#' @export
CosmxFovPlot = function(fovs=NULL,which=NULL,img_dir = NULL,maxres=2000){
  fovs = fovs %>% dplyr::arrange(.data[["FOV"]])
  which = which %||% unique(fovs[["FOV"]])
  fovs = fovs %>% dplyr::filter(.data[["FOV"]] %in% which)
  bbox = split(fovs,f=fovs$FOV) %>% lapply(.,sf::st_bbox)
  imgs = list.files(img_dir,full.names=T) %>% sort()
  if(all(stringr::str_detect(imgs,paste(fovs$fov_ID_v2,collapse="|"),negate=T))){
    imgs = imgs[stringr::str_detect(imgs,paste(fovs$fov_ID,collapse="|"))]
  }else{
    imgs = imgs[stringr::str_detect(imgs,paste(fovs$fov_ID_v2,collapse="|"))]
  }
  suppressWarnings({imgs = lapply(imgs,function(x){
    x = terra::rast(x)
    scalef = round(max(dim(x))/maxres,1)
    if(scalef >= 2){
      x = terra::aggregate(x,fact=scalef,fun='mean')
    }
    x = terra::stretch(x)
    RGB(x) = c(1,2,3)
    return(x)
  })
  })
  imgs = lapply(seq_len(length(imgs)),function(x){terra::ext(imgs[[x]]) = terra::ext(bbox[[x]]);return(imgs[[x]])})
  print(ggplot()+theme_classic()+tidyterra::geom_spatraster_rgb(data=purrr::reduce(imgs,merge))+
          theme(axis.line=element_blank(),
                plot.background = element_rect(fill = "black",colour = "black"),
                axis.text = element_text(colour = "white"),
                axis.ticks=element_blank())
  )

}




#################################
###########VISIUM UTILS###########
#################################

#' @title Create sf geometry object from Visium Seurat data
#' @name visium2sf
#' @description
#' Create a sf data frame object (package `sf`) using Visium spots coordinates
#' @details
#' This function will use Visium spot coordinates using an external 10X spot position file or spot coordinates already stored inside the input Seurat object.
#' The function is compatible with either Visium or Visium-HD data format.
#'
#' @param seurat a Seurat object containing Visium/Visium-HD data
#' @param return.seurat boolean, whether returning the input Seurat file with the nb object stored in the `misc` slot (`TRUE`) or the nb object alone (`FALSE`). Default value is `TRUE`
#' @param is.hd boolean, whether or not the Seurat object contain data generated via Visium-HD technology. Default value is `FALSE`
#' @param binsize numeric value indicating which bin resolution will be considered to create the final sf object when working with Visium-HD data
#' @param res character string specifying which version of H&E Visium image have been used to create Kandinsky data. Must be one of the following: `low`, `high`, `full`, where `full` refers to the original H&E full-resolution tiff image.
#' This parameter will be used to select the right scaling factor to align spot and image coordinates
#' @param img Visium H&E image in array or SpatRaster format. Image resolution should match with the one specified for spot coordinates through the `res` parameter
#' @returns Seurat object with the new sf object stored in the `misc` slot or the sf object alone, depending on the `return.seurat` parameter
#' @export
visium2sf = function(seurat = NULL,is.hd=F,return.seurat=T,binsize=16,res=NULL,img=NULL){
  if(is.hd == T){
    hd.assay = names(seurat@assays)[stringr::str_detect(names(seurat@assays),paste0(binsize,'um'))]
    DefaultAssay(seurat) = hd.assay
    img.id = names(seurat@images)[stringr::str_detect(names(seurat@images),paste0(binsize,'um'))]
    res = seurat@images[[img.id]]@scale.factors$img_res
    ymin = seurat@images[[img.id]]@scale.factors$ymin_hires
    fiducial = seurat@images[[img.id]]@scale.factors$fiducial_diameter_fullres
    scalef_h = seurat@images[[img.id]]@scale.factors$tissue_hires_scalef
  }else{
    res = res %||% seurat@images[[1]]@scale.factors$img_res
  }
  #  if(is.null(spot_file)){
  if(res %in% c('full','high','low')){
    if(res == 'full'){
      coords = Seurat::GetTissueCoordinates(seurat,scale=NULL)
    }else if (res == 'low'){
      coords = Seurat::GetTissueCoordinates(seurat,scale='lowres')
    }else{
      coords = Seurat::GetTissueCoordinates(seurat,scale='hires')
      #coords = coords %>% dplyr::mutate_all(list(~(.*3.333333)))
    }
  }else{
    stop('res parameter must be equal to "full", "high", or "low"')
  }
  #  }else{
  #    if(is.hd==F){
  #      coords = utils::read.delim(spot_file,sep=',') %>%
  #        tibble::column_to_rownames(.data,var='barcode')
  #    }else{
  #      coords = arrow::read_parquet(spot_file) %>%
  #        tibble::column_to_rownames(.data,var='barcode')
  #    }
  #    if(res %in% c('full','high','low')){
  #      if(res == 'full'){
  #        coords = coords[colnames(seurat),c('pxl_row_in_fullres','pxl_col_in_fullres')]
  #      }else if (res == 'low'){
  #        coords = coords[colnames(seurat),c('array_row','array_col')]
  #      }else{
  #        coords = coords[colnames(seurat),c('array_row','array_col')]
  #        if(is.hd==F){
  #          #multiply coords by low-to-high resolution scale factor
  #          coords = coords %>% dplyr::mutate_all(list(~(.data*3.333333)))
  #        }else{
  #          coords = coords %>% dplyr::mutate_all(list(~(.data*10)))
  #        }
  #      }
  #    }else{
  #      stop('res parameter must be equal to "full", "high", or "low"')
  #    }
  #  }
  colnames(coords)[c(1,2)] = c('y','x')
  #Align spot and image coordinates
  coords[,'y'] = (max(coords[,'y']) - coords[,'y']) + min(coords[,'y'])
  img.dim = dim(img)
  shift = c(0, img.dim[[1]]) - range(coords[,'y'])
  coords[,'y'] = coords[,'y']+sum(shift)
  #if(res %in% c('high','low')){
  #  if(is.hd==F){
  #    if(res != 'low'){
  #    coords[,'y'] = coords[,'y'] - seurat@images[[1]]@scale.factors$fiducial
  #    }else{
  #      ##account for scale factor
  #      coords[,'y'] = coords[,'y'] - (seurat@images[[1]]@scale.factors$fiducial/3.333331)
  #    }
  #  }else{
  #    ##Visium HD bins/image alignment
  #    ##NB: need to verify whether y axis shift parameters are valid even in other samples
  #    shift = ymin - fiducial*scalef_h
  #    if(res == 'high'){
  #      #First y-axis shift based on fiducial scale factor from 10X output
  #      coords[,'y'] = coords[,'y']+shift
  #      #Second y-axis shift based on visual inspection - NB: might change for different samples, maybe I should put this extra step as optional or user-based
  #      coords[,'y'] = coords[,'y'] - 85
  #    }else{
  #      #First y-axis shift based on fiducial scale factor from 10X output
  #      coords[,'y'] = coords[,'y']+shift/10
  #      #Second y-axis shift based on visual inspection - NB: might change for different samples, maybe I should put this extra step as optional or user-based
  #      coords[,'y'] = coords[,'y'] - 8.5
  #    }
  #  }
  #}
  if(all(colnames(coords) %in% c(c('x','y')))){
    coords$cell = rownames(coords)
  }
  colnames(coords)[colnames(coords) == 'cell'] = 'spot_ID'
  meta = seurat@meta.data
  meta$spot_ID = rownames(meta)
  meta = merge(meta,coords,by='spot_ID')
  rownames(meta) = meta$spot_ID
  sfobj = sfheaders::sf_point(meta,x='x',y='y',keep=T)
  rownames(sfobj) = sfobj$spot_ID
  if(return.seurat==T){
    if(!is.null(KanData(seurat))){
      KanData(seurat,'sf') = sfobj
    }else{
      seurat@tools$sf = sfobj
    }
    seurat$spot_ID=sfobj[["spot_ID"]]
    return(seurat)
  }else{
    return(sfobj)
  }
}

#' @title Create spatraster image from Visium image file
#' @name visium2rast
#' @description
#' Create a SpatRaster object (package `terra`) using Visium image file
#' @details
#' This function will use Visium H&E image file to create a SpatRaster image for downstream analysis and visualization.
#' The function is compatible with either Visium or Visium-HD data format.
#'
#' @param seurat a Seurat object containing Visium/Visium-HD data
#' @param img_path character string indicating the path of the H&E image to be used for Visium/VisiumHD. If set to `NULL`, no image is loaded and stored in Seurat object
#' @param rm_old_img boolean, whether or not removing any Image object already loaded in the Seurat object (`Images()`)
#' @param return.seurat boolean, whether returning the input Seurat file with the nb object stored in the `misc` slot (`TRUE`) or the nb object alone (`FALSE`). Default value is `TRUE`
#' @param max_dim numeric value indicating the ideal maximum pixel dimension accepted for Visium H&E image. If image size is bigger than this value, a scale factor will be applied to lower image pixel resolution/size up to the defined maximum value
#' @returns Seurat object with the new raster image stored in the `misc` slot or the raster image alone, depending on the `return.seurat` parameter
#' @importFrom terra RGB RGB<-
#' @export
visium2rast = function(seurat = NULL,img_path = NULL,rm_old_img = F,return.seurat=T,max_dim = 2000){
  suppressWarnings({
    img = terra::rast(img_path)
    scalef = round(max(dim(img))/max_dim,1)
    if(scalef >= 2){
      img = terra::aggregate(img,fact=scalef,fun='mean')
    }
    RGB(img) = c(1,2,3)
  })
  if(utils::packageVersion('terra') > '1.7.78'){
    img = terra::flip(img)
  }
  if(return.seurat==T){
    if(!is.null(KanData(seurat))){
      KanData(seurat,'img') = img
    }else{
      seurat@tools$img = img
    }
    if(rm_old_img == T){
      message('Removing Seurat image slot data')
      seurat@images[[1]] = NULL
    }
    message("New SpatRaster data stored as tool 'img' data (i.e., seurat@tools$img)")
    return(seurat)
  }else{
    return(img)
  }
}




#' @title Extract texture features from Visium image
#' @name get_visium_textures
#' @description
#' Gray-level co-occurrence matrix (GLCM) texture metrics extraction from Visium H&E image
#'
#' @details
#' Use the `GLCMTextures` package to extract and summarise GLCM texture features for each tissue spot in the Visium dataset.
#' Since GLCM texture extraction works with single-layer rasters, the function will extract texture features separately for each
#' R/G/B raster layer, plus a fourth created by calculating the cell-wise standard deviation across the three RGB layers.
#' For more information about GLCM texture extraction, see https://github.com/ailich/GLCMTextures
#'
#' @param visium a Visium Seurat object containing Kandinsky data (`KanData()`)
#' @param aggr.factor numeric argument indicating how many raster cells consider to average RGB intensity and lower raster resolution.
#' @param is.hd whether Seurat object contains data from standard Visium (FALSE) or Visium HD (TRUE)
#' The higher the aggr.factor, the lower the final raster resolution.
#' @returns Update Seurat object with a new 'texture' dataframe containing GLCM texture features for each spot stored in the Kandinsky data
#' @importFrom terra RGB RGB<-
#' @export
get_visium_textures = function(visium,aggr.factor=15,is.hd=F){
  if(is.null(KanData(visium,'img'))){
    stop('Visium object must contain a spatraster image. Please use visium2rast() function.')
  }
  if(is.null(KanData(visium,'mask'))){
    visium = he_mask(visium)
  }
  img = KanData(visium,'mask')
  if(is.null(KanData(visium,'sf'))){
    stop('Visium object must contain polygon data in Kandinsky slot "sf"')
  }
  #Perform raster quantization on the aggregated version of each raster RGB layer, plus a fourth layer with standard deviations for each cell across the RGB layers
  message('Raster quantization...')
  if(is.null(RGB(img))){
    RGB(img) = 1:3
  }
  rq_equalprob_r<- GLCMTextures::quantize_raster(r = terra::aggregate((img[[1]]),fact=aggr.factor,fun='mean',na.rm=T), n_levels = 32, method = "equal prob")
  rq_equalprob_g<- GLCMTextures::quantize_raster(r = terra::aggregate((img[[2]]),fact=aggr.factor,fun='mean',na.rm=T), n_levels = 32, method = "equal prob")
  #rq_equalprob_b<- GLCMTextures::quantize_raster(r = aggregate((img[[3]]),fact=aggr.factor,fun=mean,na.rm=T), n_levels = 16, method = "equal prob")
  rq_equalprob_gray<- GLCMTextures::quantize_raster(r = terra::aggregate(terra::stretch(terra::colorize(img,to='col',grays=T)),fact=aggr.factor,fun='mean',na.rm=T), n_levels = 32, method = "equal prob")
  message('...done!')

  message('GLCM texture metrics extraction...')
  texturer = GLCMTextures::glcm_textures(rq_equalprob_r, w = c(5,5), n_levels = 32, quantization = "none",na.rm=T)
  textureg = GLCMTextures::glcm_textures(rq_equalprob_g, w = c(5,5), n_levels = 32, quantization = "none",na.rm=T)
  texturegray = GLCMTextures::glcm_textures(rq_equalprob_gray, w = c(5,5), n_levels = 32, quantization = "none",na.rm=T)

  if(is.hd == F){
    cell_texturer_min = terra::extract(texturer,terra::vect(KanData(visium,'sf')),na.rm=T) %>% dplyr::group_by(.data$ID) %>% dplyr::summarise_if(is.numeric,min,na.rm=T) %>% dplyr::ungroup() %>% as.data.frame()
    cell_texturer_med = terra::extract(texturer,terra::vect(KanData(visium,'sf')),na.rm=T) %>% dplyr::group_by(.data$ID) %>% dplyr::summarise_if(is.numeric,median,na.rm=T) %>% dplyr::ungroup() %>% as.data.frame()
    cell_texturer_max = terra::extract(texturer,terra::vect(KanData(visium,'sf')),na.rm=T) %>% dplyr::group_by(.data$ID) %>% dplyr::summarise_if(is.numeric,max,na.rm=T) %>% dplyr::ungroup() %>% as.data.frame()
    colnames(cell_texturer_min)[-1] = paste0(colnames(cell_texturer_min)[-1],"_q0")
    colnames(cell_texturer_med)[-1] = paste0(colnames(cell_texturer_med)[-1],"_q50")
    colnames(cell_texturer_max)[-1] = paste0(colnames(cell_texturer_max)[-1],"_q100")
    cell_texturer = purrr::reduce(list(cell_texturer_min,cell_texturer_med,cell_texturer_max),merge,by='ID')
    rownames(cell_texturer) = rownames(KanData(visium,'sf'))[cell_texturer$ID]
    cell_texturer$ID = NULL
    colnames(cell_texturer) = paste0(colnames(cell_texturer),"_R")

    cell_textureg_min = terra::extract(textureg,terra::vect(KanData(visium,'sf')),na.rm=T) %>% dplyr::group_by(.data$ID) %>% dplyr::summarise_if(is.numeric,min,na.rm=T) %>% dplyr::ungroup() %>% as.data.frame()
    cell_textureg_med = terra::extract(textureg,terra::vect(KanData(visium,'sf')),na.rm=T) %>% dplyr::group_by(.data$ID) %>% dplyr::summarise_if(is.numeric,median,na.rm=T) %>% dplyr::ungroup() %>% as.data.frame()
    cell_textureg_max = terra::extract(textureg,terra::vect(KanData(visium,'sf')),na.rm=T) %>% dplyr::group_by(.data$ID) %>% dplyr::summarise_if(is.numeric,max,na.rm=T) %>% dplyr::ungroup() %>% as.data.frame()
    colnames(cell_textureg_min)[-1] = paste0(colnames(cell_textureg_min)[-1],"_q0")
    colnames(cell_textureg_med)[-1] = paste0(colnames(cell_textureg_med)[-1],"_q50")
    colnames(cell_textureg_max)[-1] = paste0(colnames(cell_textureg_max)[-1],"_q100")
    cell_textureg = purrr::reduce(list(cell_textureg_min,cell_textureg_med,cell_textureg_max),merge,by='ID')
    rownames(cell_textureg) = rownames(KanData(visium,'sf'))[cell_textureg$ID]
    cell_textureg$ID = NULL
    colnames(cell_textureg) = paste0(colnames(cell_textureg),"_G")

    texturegray = GLCMTextures::glcm_textures(rq_equalprob_gray, w = c(5,5), n_levels = 32, quantization = "none",na.rm=T)
    cell_texturegray_min = terra::extract(texturegray,terra::vect(KanData(visium,'sf')),na.rm=T) %>% dplyr::group_by(.data$ID) %>% dplyr::summarise_if(is.numeric,min,na.rm=T) %>% dplyr::ungroup() %>% as.data.frame()
    cell_texturegray_med = terra::extract(texturegray,terra::vect(KanData(visium,'sf')),na.rm=T) %>% dplyr::group_by(.data$ID) %>% dplyr::summarise_if(is.numeric,median,na.rm=T) %>% dplyr::ungroup() %>% as.data.frame()
    cell_texturegray_max = terra::extract(texturegray,terra::vect(KanData(visium,'sf')),na.rm=T) %>% dplyr::group_by(.data$ID) %>% dplyr::summarise_if(is.numeric,max,na.rm=T) %>% dplyr::ungroup() %>% as.data.frame()
    colnames(cell_texturegray_min)[-1] = paste0(colnames(cell_texturegray_min)[-1],"_q0")
    colnames(cell_texturegray_med)[-1] = paste0(colnames(cell_texturegray_med)[-1],"_q50")
    colnames(cell_texturegray_max)[-1] = paste0(colnames(cell_texturegray_max)[-1],"_q100")
    cell_texturegray = purrr::reduce(list(cell_texturegray_min,cell_texturegray_med,cell_texturegray_max),merge,by='ID')
    rownames(cell_texturegray) = rownames(KanData(visium,'sf'))[cell_texturegray$ID]
    cell_texturegray$ID = NULL
    colnames(cell_texturegray) = paste0(colnames(cell_texturegray),"_Gray")
  }else{
    cell_texturer = terra::extract(texturer,terra::vect(KanData(visium,'sf')),na.rm=T,fun='mean')
    colnames(cell_texturer)[-1] = paste0(colnames(cell_texturer)[-1],"_Red")
    rownames(cell_texturer) = rownames(KanData(visium,'sf'))[cell_texturer$ID]
    cell_texturer$ID = NULL

    cell_textureg = terra::extract(textureg,terra::vect(KanData(visium,'sf')),na.rm=T,fun='mean')
    colnames(cell_textureg)[-1] = paste0(colnames(cell_textureg)[-1],"_Green")
    rownames(cell_textureg) = rownames(KanData(visium,'sf'))[cell_textureg$ID]
    cell_textureg$ID = NULL

    cell_texturegray = terra::extract(textureg,terra::vect(KanData(visium,'sf')),na.rm=T,fun='mean')
    colnames(cell_texturegray)[-1] = paste0(colnames(cell_texturegray)[-1],"_Gray")
    rownames(cell_texturegray) = rownames(KanData(visium,'sf'))[cell_texturegray$ID]
    cell_texturegray$ID = NULL
  }
  message('...done!')

  KanData(visium,'texture') = do.call('cbind',list(cell_texturer,cell_textureg,cell_texturegray))#,cell_textureg,cell_textureb))
  return(visium)
}

#################################
###########XENIUM UTILS###########
#################################

#Load morphology OME tiff images from Xenium output with the help of RBioFormats read.image() function
#You should expect to find 4 channels
#' @title Load Xenium OME tiff image
#' @name LoadXeniumImage
#' @description
#' Load morphology OME tiff images from Xenium output with the help of RBioFormats read.image() function
#'
#' @param path character string indicating the path of the Xenium OME tiff image
#' @param channel_cols vector of color names to be assigned to each of the OME image channel
#' @param res integer specifying resolution levels to read.
#' @returns Xenium image as a multi-layer SpatRaster object
#' @export
LoadXeniumImage = function(path = '/Volumes/colcc-ITH2_scRNA_datasets/data/10X_Xenium_CRC_public/Breast_5K/',channel_cols = c('gray70',"limegreen","cyan","indianred3"),res=5){
  file = paste0(path,"morphology_focus/morphology_focus_0001.ome.tif")
  img = terra::rast(RBioFormats::read.image(file,resolution=res,read.metadata=F,normalize=F))
  names(img) = c('DAPI','ATP1A1.ECad.CD45','rRNA_18S','aSMA.VIM')
  #img = flip(img,'vertical')
  img = terra::stretch(img)
  colfunc <- grDevices::colorRampPalette(c("black", channel_cols[1]))
  coltab(img,layer=1) =  data.frame(value = seq(0,255,1),col=colfunc(256))

  colfunc <- grDevices::colorRampPalette(c("black", channel_cols[2]))
  coltab(img,layer=2) =  data.frame(value = seq(0,255,1),col=colfunc(256))

  colfunc <- grDevices::colorRampPalette(c("black", channel_cols[3]))
  coltab(img,layer=3) =  data.frame(value = seq(0,255,1),col=colfunc(256))

  colfunc <- grDevices::colorRampPalette(c("black", channel_cols[4]))
  coltab(img,layer=4) =  data.frame(value = seq(0,255,1),col=colfunc(256))
  #img = flip(img,'vertical')
  return(img)
}


#' @title align Xenium image and cell coordinates
#' @name AlignXeniumCoords
#' @description
#' Adjust Xenium image coordinates to align with cell polygon coordinates
#' @param img Xenium image as a multi-layer SpatRaster object
#' @param poly Xenium cell polygon data.frame
#' @returns list object containing aligned image and cell polygons objects
#' @export
AlignXeniumCoords = function(img=NULL,poly=NULL){
  #poly = st_as_sf(poly)
  bbox = (sf::st_bbox(poly))
  bbox = sf::st_as_sfc(bbox) - matrix(c(bbox['xmin'],bbox['ymin']),ncol=2)
  #st_geometry(poly) =  st_geometry(poly) - matrix(c(bbox['xmin'],bbox['ymin']),nrow=2)
  #bbox = (st_bbox(poly))
  xmax(img) = sf::st_bbox(bbox)['xmax']
  ymax(img) = sf::st_bbox(bbox)['ymax']
  return(list(img=img,poly=poly))
}




#################################
###########G4X UTILS###########
#################################
#' @title Create spatraster image from G4X jp2 image file
#' @name load_g4x_img
#' @description
#' Create a SpatRaster object (package `terra`) using Visium image file
#' @details
#' This function will use G4X H&E image file to create a SpatRaster image for downstream analysis and visualization.
#'
#' @param seurat a Seurat object containing G4X data
#' @param img_path character string indicating the path of the H&E image to be used. If set to `NULL`, the function will look for any image path 'img' stored in tools slot
#' @param rm_old_img boolean, whether or not removing any Image object already loaded in the Seurat object (`Images()`)
#' @param return.seurat boolean, whether returning the input Seurat file with the raster object stored in the `tools` slot (`TRUE`) or the raster object alone (`FALSE`). Default value is `TRUE`
#' @param max_dim numeric value indicating the ideal maximum pixel dimension accepted for Visium H&E image. If image size is bigger than this value, a scale factor will be applied to lower image pixel resolution/size up to the defined maximum value
#' @returns Seurat object with the new raster image stored in the `misc` slot or the raster image alone, depending on the `return.seurat` parameter
#' @importFrom terra RGB RGB<-
#' @export
load_g4x_img = function(seurat = NULL,img_path = NULL,rm_old_img = F,return.seurat=T,max_dim = 2000){
  if(is.null(img_path)){
    check = names(seurat@tools)
    if(any(check == 'img_path')){
      img_path = seurat@tools$img_path
    }
  }
  suppressWarnings({
    img = terra::rast(img_path)
    scalef = round(max(dim(img))/max_dim,1)
    if(scalef >= 2){
      img = terra::aggregate(img,fact=scalef,fun='mean')
    }
    RGB(img) = c(1,2,3)
  })
  if(utils::packageVersion('terra') <= '1.7.78'){
    img = terra::flip(img)
  }
  if(return.seurat==T){
    if(!is.null(KanData(seurat))){
      KanData(seurat,'img') = img
    }else{
      seurat@tools$img = img
    }
    if(rm_old_img == T){
      message('Removing Seurat image slot data')
      seurat@images[[1]] = NULL
    }
    message("New SpatRaster data stored as tool 'img' data (i.e., seurat@tools$img)")
    return(seurat)
  }else{
    return(img)
  }
}
