#' @title Prepare 10X Visium Seurat data
#' @name prepare_visium_seurat
#'
#' @description
#' Initial formatting of Visium data to work with Seurat and Kandinsky
#' @details
#' This function will use 10X Visium raw input data to build a new Seurat object. All Visium input files must be stored in the same folder that will be specified to call the function
#' @param path character string specifying the path to the Visium input files directory
#' @param dataset.id character string that will be used to name the output Seurat object identity
#' @param img.res character string specifying which version of H&E Visium image will be loaded into the Seurat data. Must be one of the following: `low`, `high`, `full`, where `full` refers to the original H&E full-resolution tiff image. If set to `full`, lowres image will be loaded into Seurat, but the path to the full res tif image will be also stored in the final object.
#' @param fullres_path character string specifying the path to the Visium full resolution H&E tif image. The full image path will be stored into the final Seurat object when `img.res` is set to `full`, and will be used by `Kandinsky` functions for downstream analyses
#' @returns a Seurat object containing Visium count matrix, metadata, and additional files/file paths to be used for downstream analysis with Kandinsky functions
#' @importFrom magrittr %>%
#' @importFrom arrow read_csv_arrow
#' @importFrom sfheaders sf_polygon sf_point
#' @export
#' @family prepare_data
prepare_visium_seurat = function(path=NULL,dataset.id='X',img.res=c('full','high','low'),fullres_path=NULL){
  if(substr(path,nchar(path),nchar(path)) == '/'){path = substr(path,1,nchar(path)-1)}
  path = gsub("//","/",path)
  spatial_path = list.files(path,pattern='spatial',full.names=T)
  if(is.null(img.res)| length(img.res) > 1){
    img.res='low'
  }else if(!(img.res %in% c('full','high','low'))){
    warning('Unknown img.res parameter value,setting img.res to "high"')
    img.res='low'
  }
  if(img.res == 'high'){
    img_path = normalizePath(list.files(spatial_path,pattern='hires_image.png',full.names=T))
    img = 'tissue_hires_image.png'
  }else if(img.res == 'low'){
    img_path = normalizePath(list.files(spatial_path,pattern='lowres_image.png',full.names=T))
    img = 'tissue_lowres_image.png'
  }else{
    img_path = normalizePath(fullres_path)
    img = 'tissue_lowres_image.png'
  }
  mat_file = list.files(path)
  mat_file = mat_file[stringr::str_detect(mat_file,'matrix.mtx')]
  if(length(mat_file) == 0){
    path = paste0(path,'/filtered_feature_bc_matrix')
  }
  mat_file = list.files(path)
  mat_file = mat_file[stringr::str_detect(mat_file,'matrix.mtx')]
  if(length(mat_file) == 0){
    stop('Please provide a path containing either Visium count matrix or a subdirectory "filtered_feature_bc_matrix" containing such file')
  }
  message('Preparing gene expression matrix...')
  counts <- Seurat::Read10X(data.dir = path,
                            gene.column = 2,
                            cell.column = 1,
                            unique.features = TRUE,
                            strip.suffix = FALSE)
  message('Creating Seurat object...')
  sp <- Seurat::CreateSeuratObject(counts = counts,
                                   project = '10XVisium',
                                   assay = 'Visium10x')

  DefaultAssay(sp) <- 'Visium10x' #The default assay is set to `Visium10x`
  img <- Seurat::Read10X_Image(image.dir = spatial_path,
                               image.name = img,
                               assay = 'Visium10x',
                               slice = dataset.id,
                               filter.matrix = TRUE)

  # Getting the cell (spots) names
  sp_cells = Seurat::Cells(sp)
  img <- img[sp_cells]

  # Add the image to the Default object `sp`
  message('Loading H&E image...')
  sp[[dataset.id]] <- img
  sp@tools$img_path = img_path
  sp@images[[1]]@scale.factors$img_res = img.res
  return(sp)
}

#' @title Prepare 10X Visium HD Seurat data
#' @name prepare_visiumHD_seurat
#'
#' @description
#' Initial formatting of Visium data to work with Seurat and Kandinsky
#' @details
#' This function will use 10X Visium raw input data to build a new Seurat object. All Visium input files must be stored in the same folder that will be specified to call the function
#' @param path character string specifying the path to the Visium input files directory
#' @param dataset.id character string that will be used to name the output Seurat object identity
#' @param img.res character string specifying which version of H&E Visium image will be loaded into the Seurat data. Must be one of the following: `low`, `high`, `full`, where `full` refers to the original H&E full-resolution tiff image. If set to `full`, lowres image will be loaded into Seurat, but the path to the full res tif image will be also stored in the final object.
#' @param fullres_path character string specifying the path to the Visium full resolution H&E tif image. The full image path will be stored into the final Seurat object when `img.res` is set to `full`, and will be used by `Kandinsky` functions for downstream analyses
#' @param binsize numeric vector specifying the binsize(s) to be used to build the Seurat object. Can be any combination of `2`, `8`, and `16` sizes
#' @returns a Seurat object containing Visium count matrix, metadata, and additional files/file paths to be used for downstream analysis with Kandinsky functions
#' @importFrom magrittr %>%
#' @importFrom arrow read_csv_arrow
#' @importFrom sfheaders sf_polygon sf_point
#' @export
#' @family prepare_data
prepare_visiumHD_seurat = function(path=NULL,dataset.id='X',img.res=c('full','high','low'),fullres_path=NULL,binsize=c(2,8,16)){
  if(substr(path,nchar(path),nchar(path)) == '/'){path = substr(path,1,nchar(path)-1)}
  path = gsub("//","/",path)
  spatial_path = list.files(path,pattern='spatial',full.names=T)
  if(length(spatial_path)>1){
    spatial_path = spatial_path[stringr::str_detect(spatial_path,'tar|gz',negate=T)]
  }
  bins = stringr::str_pad(binsize,width=3,pad='0')
  if(is.null(img.res)| length(img.res) > 1){
    img.res='low'
  }else if(!(img.res %in% c('full','high','low'))){
    warning('Unknown img.res parameter value,setting img.res to "high"')
    img.res='low'
  }
  if(img.res == 'high'){
    img_path = normalizePath(list.files(spatial_path,pattern='hires_image.png',full.names=T))
    img = 'tissue_hires_image.png'
  }else if(img.res == 'low'){
    img_path = normalizePath(list.files(spatial_path,pattern='lowres_image.png',full.names=T))
    img = 'tissue_lowres_image.png'
  }else{
    img_path = normalizePath(fullres_path)
    img = 'tissue_lowres_image.png'
  }
  bin_path = paste0(path,'/binned_outputs')
  if(length(binsize>1)){
    start_bin = max(bins)
  }else{
    start_bin = bins
  }
  all_seurat = lapply(bins,function(b){
    message('Preparing ',b,'um gene expression matrix...')
    newpath = list.files(bin_path,pattern=b,full.names=T)
    spatial_newpath = list.files(newpath,pattern='spatial',full.names=T)
    ##Read scale factor file and calculate square buffer size for hires and lowres coords
    scalef = arrow::read_json_arrow(paste0(spatial_newpath,'/scalefactors_json.json'))
    buf_h = (scalef$spot_diameter_fullres/2)*scalef$tissue_hires_scalef
    ##Read bin position file and calculate minimym y-axis coordinate for downstream image alignment
    pos = arrow::read_parquet(paste0(spatial_newpath,'/tissue_positions.parquet'))[,c(1,2,5,6)]
    pos$y_hires = pos$pxl_row_in_fullres*scalef$tissue_hires_scalef
    pos$y_hires = (max(pos$y_hires) - pos$y_hires) + min(pos$y_hires)
    ymin = min(pos$y_hires) - buf_h
    rm(pos)
    message('Creating ',b,'um Seurat assay...')
    counts <- Seurat::Read10X_h5(filename = paste0(newpath,'/filtered_feature_bc_matrix.h5'),
                                 use.names = T,unique.features = T)
    sp <- Seurat::CreateSeuratObject(counts = counts,
                                     project = 'VisiumHD',
                                     assay = paste0('VisHD_',b,'um'))
    rm(counts)
    DefaultAssay(sp) <- paste0('VisHD_',b,'um')
    img <- Seurat::Read10X_Image(image.dir = spatial_newpath,
                                 image.name = img,
                                 assay = paste0('VisHD_',b,'um'),
                                 slice = paste0(dataset.id,"_",b,'um'),
                                 filter.matrix = TRUE)
    # Getting the cell (spots) names
    DefaultAssay(sp) <- paste0('VisHD_',b,'um')
    sp_cells = Seurat::Cells(sp)
    img <- img[sp_cells]

    # Add the image to the Default object `sp`
    sp[[paste0(dataset.id,"_",b,'um')]] <- img
    sp@images[[1]]@scale.factors[names(scalef)] = scalef
    sp@images[[1]]@scale.factors$buffer_hires = buf_h
    sp@images[[1]]@scale.factors$ymin_hires = ymin
    sp@images[[1]]@scale.factors$img_res = img.res
    rm(img)
    return(sp)
  })
  message('Loading H&E image...')
  all_seurat = merge(x=all_seurat[[1]],all_seurat[-1])
  DefaultAssay(all_seurat) = paste0('VisHD_',start_bin,'um')
  all_seurat@tools$img_path = img_path
  return(all_seurat)
}

#' @title Prepare NanoString Cosmx Seurat data
#' @name prepare_cosmx_seurat
#'
#' @description
#' Initial formatting of CosMx data to work with Seurat and Kandinsky
#' @details
#' This function will use NanoString CosMx raw input data to build a new Seurat object. All Cosmx input files must be stored in the same folder that will be specified to call the function
#' @param path character string specifying the path to the CosMx input files directory
#' @param dataset.id character string that will be used to name the output Seurat object identity
#' @param pattern character string that will be used as a key to identify CosMx input file names in case when input files from multiple CosMx samples/slides are stored in the same directory
#' @param fovs vector of FOV identifiers to subset for downstream analysis.
#' @returns a Seurat object containing CosMx count matrix, metadata, and additional files/file paths to be used for downstream analysis with Kandinsky functions
#' @importFrom magrittr %>%
#' @importFrom arrow read_csv_arrow
#' @importFrom sfheaders sf_polygon sf_point
#' @export
#' @family prepare_data
prepare_cosmx_seurat = function(path=NULL,dataset.id='X',pattern=NULL,fovs=NULL){
  if(substr(path,nchar(path),nchar(path)) == '/'){path = substr(path,1,nchar(path)-1)}
  path = gsub("//","/",path)
  poly = list.files(path,pattern='poly',full.names = T)
  if(length(poly) > 1){
    poly = poly[stringr::str_detect(poly,pattern)]
  }
  meta = list.files(path,pattern='metadata',full.names = T)
  if(length(meta) > 1){
    meta = meta[stringr::str_detect(meta,pattern)]
  }
  mat = list.files(path,pattern='exprMat',full.names = T)
  if(length(mat) > 1){
    mat = mat[stringr::str_detect(mat,pattern)]
  }
  tx = list.files(path,pattern='tx_',full.names = T)
  if(length(tx) > 1){
    tx = tx[stringr::str_detect(tx,pattern)]
  }
  fov = list.files(path,pattern='fov_',full.names = T)
  if(length(fov) > 1){
    fov = fov[stringr::str_detect(fov,pattern)]
  }

  if(length(mat) > 0){
    message('Preparing gene expression matrix...')
    #Do not use arrow to load CosMx matrix (too many columns)
    mat = data.table::fread(mat,sep=',',data.table=F)
    #Create new cell identifiers
    mat = mat %>%
      dplyr::mutate(cell_ID = paste0(.data$fov,'_',.data$cell_ID))
    if(!is.null(fovs)){
      mat = mat %>% dplyr::filter(.data$fov %in% fovs)
    }
    #Set cell ids as rownames and remove the old id columns from the matrix
    rownames(mat) = mat$cell_ID
    mat[stringr::str_detect(colnames(mat), "^[:upper:]",negate=T)] = NULL
    #Convert count matrix in a sparse format (easier to manage) and transpose it in order to have genes as rows and cells as columns
    mat = as(as.matrix(mat),'CsparseMatrix')
    mat = Matrix::t(mat)

    #Divide transcripts from negative probes
    negmat = (mat[stringr::str_detect(rownames(mat),'Neg'),])
    sysmat = (mat[stringr::str_detect(rownames(mat),'System'),])
    mat = mat[stringr::str_detect(rownames(mat),'Neg|System',negate=T),]
  }else{
    stop('CosMx count matrix not found in the specified folder')
  }

  if(length(meta) > 0){
    message('Loading cell metadata...')
    #Prepare metadata
    meta = arrow::read_csv_arrow(meta)
    #Convert px to um
    meta$CenterX_global_um = meta$CenterX_global_px*0.12028
    meta$CenterY_global_um = meta$CenterY_global_px*0.12028
    #Create new cell identifier
    meta = meta %>%
      dplyr::mutate(cell_ID = paste0(.data$fov,'_',.data$cell_ID)) %>%
      as.data.frame()
    if(!is.null(fovs)){
      meta = meta %>% dplyr::filter(.data$fov %in% fovs)
    }

    common = intersect(colnames(mat),meta$cell_ID)

    #Reorder matrix and metadata
    mat = mat[,common]
    if(dim(negmat)[1] > 0){
      negmat =  negmat[,common]
    }
    if(dim(sysmat)[1] > 0){
      sysmat = sysmat[,common]
    }
    rownames(meta) = meta$cell_ID
    meta = meta[common,]
  }else{
    stop('CosMx metadata file not found in the specified folder')
  }
  if(length(poly) > 0){
    message('Loading cell polygon file...')
    poly = arrow::read_csv_arrow(poly)
    poly = poly %>%
      dplyr::mutate(cell_ID = paste0(.data$fov,'_',.data$cellID))
    if(!is.null(fovs)){
      poly = poly %>% dplyr::filter(.data$fov %in% fovs)
    }
    #Convert pixel to um
    poly$x_global_um = poly$x_global_px*0.12028
    poly$y_global_um = poly$y_global_px*0.12028
    poly = sfheaders::sf_polygon(poly,x='x_global_um',y='y_global_um',polygon_id='cell_ID',keep=T)
    rownames(poly) = poly$cell_ID
    poly = poly[rownames(meta),]
  }else{
    warning('Cosmx polygon file not found in the specified folder.\n
            Polygons will be created using cell centroid coordinates in the metadata table')
    poly = sfheaders::sf_point(meta %>% dplyr::select(.data$CenterX_global_px,.data$CenterY_global_px,.data$cell_ID,.data$fov),
                               x='CenterX_global_px',y='CenterY_global_px',polygon_id='cell_ID',keep=T)
  }
  message('Creating Seurat object...')
  #Create Seurat Object with separate CosMx and NegProbes Assays
  cosmx =  Seurat::CreateSeuratObject(counts = mat, meta.data = meta,assay='CosMx')
  if(dim(negmat)[1] > 0){
    cosmx[["negprobes"]] = Seurat::CreateAssayObject(negmat)
  }
  if(dim(sysmat)[1] > 0){
    cosmx[["SystemControl"]] = Seurat::CreateAssayObject(sysmat)
  }
  #Prepare single-cell spatial coordinates info (i.e., centroid coordinates)
  cents =  CreateCentroids(meta %>% dplyr::select(.data$CenterY_global_px,.data$CenterX_global_px,.data$cell_ID))
  segmentations.data <- list(
    "centroids" = cents
  )
  coords =CreateFOV(
    coords = segmentations.data,
    type = c("centroids"),
    molecules = NULL,
    assay = "CosMx"
  )

  cosmx[[dataset.id]] = coords
  cosmx@tools$poly = poly
  if(length(tx) > 0){
    message('Locating cell transcript file...')
    cosmx@tools$tx = normalizePath(tx)
  }
  if(length(fov)>0){
    message('Loading FOV position file...')
    cosmx@tools$fov = read_fovfile(fov)
    cosmx@tools$fov = cosmx@tools$fov[cosmx@tools$fov$FOV %in% cosmx@meta.data$fov,]
  }
  return(cosmx)
}


#' @title Prepare 10X Xenium Seurat data
#' @name prepare_xenium_seurat
#'
#' @description
#' Initial formatting of Xenium data to work with Seurat and Kandinsky
#' @details
#' This function will use 10X Xenium raw input data to build a new Seurat object. All Xenium input files must be stored in the same folder that will be specified to call the function
#' @param path character string specifying the path to the Xenium input files directory
#' @param dataset.id character string that will be used to name the output Seurat object identity
#' @param h5 boolean, whether reading count matrix from 'cell_feature_matrix.h5' (TRUE) file or zipped 'cell_feature_matrix.tar.gz' folder (FALSE)
#' @returns a Seurat object containing Xenium count matrix, metadata, and additional files/file paths to be used for downstream analysis with Kandinsky functions
#' @importFrom magrittr %>%
#' @importFrom arrow read_csv_arrow
#' @importFrom sfheaders sf_polygon sf_point
#' @export
#' @family prepare_data
prepare_xenium_seurat = function(path=NULL,dataset.id='X',h5=TRUE){
  message('Loading cell polygon file...')
  poly = arrow::read_parquet(paste0(path,'cell_boundaries.parquet'))
  poly = sfheaders::sf_polygon(poly,x='vertex_y',y='vertex_x',polygon_id='cell_id',keep=T)
  rownames(poly) = poly$cell_id
  message('Loading cell metadata...')
  meta = as.data.frame(arrow::read_parquet(paste0(path,'cells.parquet')))
  message('Locating cell transcript file...')
  tx = normalizePath(paste0(path,'transcripts.parquet'))
  rownames(meta) = meta$cell_id
  poly = poly[rownames(meta),]
  message('Preparing fov position data...')
  utils::untar(paste0(path,'aux_outputs.tar.gz'),exdir=path)
  fov_file = list.files(paste0(path,'aux_outputs'),pattern='fov_locations.json',full.names=T)
  if(length(fov_file)>1){
    fov_file = list.files(paste0(path,'aux_outputs'),pattern='morphology_fov_locations.json',full.names=T)
  }
  fovs = as.data.frame(jsonlite::read_json(fov_file))
  colnames(fovs) = gsub('fov_locations.','',colnames(fovs))
  fov_ids = sapply(colnames(fovs)[colnames(fovs)!='units'],function(x){strsplit(x,split='\\.')[[1]][[1]]})
  fov_ids = unique(fov_ids)
  fov_meta = lapply(fov_ids,function(x){
    pattern = paste0(x,'\\.')
    fovs = fovs[,stringr::str_detect(colnames(fovs),pattern)]
    colnames(fovs) = gsub(pattern,'',colnames(fovs))
    fovs$fov = x
    return(fovs)
    })
  fov_mask = lapply(fov_meta,function(x){
    width = x$width
    height= x$height
    centroid = x[,c('x','y')]
    corners <- rbind(
      centroid + c(-width/2, -height/2),
      centroid + c(width/2, -height/2),
      centroid + c(width/2, height/2),
      centroid + c(-width/2, height/2),
      centroid + c(-width/2, -height/2)
    )
    poly = sfheaders::sf_polygon(corners,x='y',y='x')
    poly[,colnames(x)] = x
    return(poly)
  })
  rm(fov_meta)
  fov_mask = purrr::reduce(fov_mask,rbind)
  fov_mask[,c('x','y')] = sf::st_drop_geometry(fov_mask)[,c('y','x')]
  fov_mask$fov_number = 1:nrow(fov_mask)
  message('Preparing gene expression matrix...')
  if(h5==F){
    utils::untar(paste0(path,'cell_feature_matrix.tar.gz'),exdir=path)
    mat = Matrix::readMM(paste0(path,'cell_feature_matrix/matrix.mtx.gz'))
    rownames(mat) = utils::read.delim(paste0(path,'cell_feature_matrix/features.tsv.gz'),header=F)[,2]
    colnames(mat) = utils::read.delim(paste0(path,'cell_feature_matrix/barcodes.tsv.gz'),header=F)[,1]
    unlink(paste0(path,'cell_feature_matrix/'),recursive=T)
  }else{
    probes = rhdf5::h5read(paste0(path,'cell_feature_matrix.h5'),name='/matrix/features/name')
    barcodes = rhdf5::h5read(paste0(path,'cell_feature_matrix.h5'),name='/matrix/barcodes')
    counts= rhdf5::h5read(paste0(path,'cell_feature_matrix.h5'),name='/matrix/data')
    index= rhdf5::h5read(paste0(path,'cell_feature_matrix.h5'),name='/matrix/indices')
    indptr= rhdf5::h5read(paste0(path,'cell_feature_matrix.h5'),name='/matrix/indptr')
    shp= rhdf5::h5read(paste0(path,'cell_feature_matrix.h5'),name='/matrix/shape')
    mat = Matrix::sparseMatrix(i=index+1,p=indptr,x=as.numeric(counts),dims=shp,repr='T')
    features = make.unique(probes)
    rownames(mat) = as.character(features)
    colnames(mat) = as.character(barcodes)
  }
  mat = as(mat,'CsparseMatrix')
  rownames(mat) = gsub("_",".",rownames(mat))

  unsmat = mat[stringr::str_detect(rownames(mat),'Unassigned'),]
  deprmat = mat[stringr::str_detect(rownames(mat),'Deprecated'),]
  negmat = mat[stringr::str_detect(rownames(mat),'NegControl'),]
  genmat = mat[stringr::str_detect(rownames(mat),'Region'),]
  mat = mat[stringr::str_detect(rownames(mat),'Unassigned|NegControl|Deprecated|Region',negate=T),]
  message('Creating Seurat object...')

  xenium = Seurat::CreateSeuratObject(counts=mat,meta.data=meta,assay='Xenium')
  rm(mat)
  if(dim(negmat)[1] > 0){
    xenium[["negprobes"]] = Seurat::CreateAssayObject(negmat)
  }
  if(dim(unsmat)[1] > 0){
    xenium[["unassigned"]] = Seurat::CreateAssayObject(unsmat)
  }
  if(dim(deprmat)[1] > 0){
    xenium[["deprecated"]] = Seurat::CreateAssayObject(deprmat)
  }
  if(dim(genmat)[1] > 0){
    xenium[["intergenic"]] = Seurat::CreateAssayObject(genmat)
  }
  Idents(xenium) = 'cell_id'
  #meta$y_centroid = (max(meta$y_centroid)-meta$y_centroid)+min(meta$y_centroid)
  #Prepare single-cell spatial coordinates info (i.e., centroid coordinates)
  cents =  CreateCentroids(meta %>% dplyr::select(.data$y_centroid,.data$x_centroid,.data$cell_id))
  segmentations.data <- list(
    "centroids" = cents
  )
  coords = CreateFOV(
    coords = segmentations.data,
    type = c("centroids"),
    molecules = NULL,
    assay = "Xenium"
  )
  xenium[[dataset.id]] = coords
  xenium@tools$poly = poly
  xenium@tools$tx = tx
  xenium@tools$fov = fov_mask
  return(xenium)
}




#' @title Prepare Vizgen Merscope Seurat data
#' @name prepare_merscope_seurat
#'
#' @description
#' Initial formatting of Merscope data to work with Seurat and Kandinsky
#' @details
#' This function will use Vizgen Merscope raw input data to build a new Seurat object. All Merscope input files must be stored in the same folder that will be specified to call the function
#' @param path character string specifying the path to the Merscope input files directory
#' @param dataset.id character string that will be used to name the output Seurat object identity
#' @param pattern character string that will be used as a key to identify Merscope input file names in case when input files from multiple Merscope samples/slides are stored in the same directory
#' @param fovs vector of FOV identifiers to subset for downstream analysis.
#' @returns a Seurat object containing Merscope count matrix, metadata, and additional files/file paths to be used for downstream analysis with Kandinsky functions
#' @importFrom magrittr %>%
#' @importFrom arrow read_csv_arrow
#' @importFrom sfheaders sf_polygon sf_point
#' @export
#' @family prepare_data
prepare_merscope_seurat = function(path=NULL,dataset.id='X',pattern=NULL,fovs=NULL){
  if(substr(path,nchar(path),nchar(path)) == '/'){path = substr(path,1,nchar(path)-1)}
  path = gsub("//","/",path)
  poly = list.files(path,pattern='boundaries|micron_space',full.names = T)
  if(length(poly) > 1){
    poly = poly[stringr::str_detect(poly,pattern)]
  }
  meta = list.files(path,pattern='metadata',full.names = T)
  if(length(meta) > 1){
    meta = meta[stringr::str_detect(meta,pattern)]
  }
  mat = list.files(path,pattern='cell_by_gene',full.names = T)
  if(length(mat) > 1){
    mat = mat[stringr::str_detect(mat,pattern)]
  }
  tx = list.files(path,pattern='transcripts',full.names = T)
  if(length(tx) > 1){
    tx = tx[stringr::str_detect(tx,pattern)]
  }
  if(length(mat) > 0){
    message('Preparing gene expression matrix...')
    #Do not use arrow to load Merscope matrix (too many columns)
    mat = data.table::fread(mat,sep=',',data.table=F)
    #Create new cell identifiers
    mat = mat %>%
      dplyr::mutate(cell_ID = paste0(dataset.id,'_',.data$cell))
    #Set cell ids as rownames and remove the old id columns from the matrix
    rownames(mat) = mat$cell_ID
    mat[stringr::str_detect(colnames(mat), "^[:upper:]",negate=T)] = NULL

  }else{
    stop('Merscope count matrix not found in the specified folder')
  }

  if(length(meta) > 0){
    message('Loading cell metadata...')
    #Prepare metadata
    meta = arrow::read_csv_arrow(meta)
    #Create new cell identifier
    meta = meta %>%
      dplyr::mutate(cell_ID = paste0(dataset.id,'_',.data$EntityID)) %>%
      as.data.frame()
    if(!is.null(fovs)){
      meta = meta %>% dplyr::filter(.data$fov %in% fovs)
    }
  }else{
    stop('Merscope metadata file not found in the specified folder')
  }
  common = intersect(rownames(mat),meta$cell_ID)

  #Reorder matrix and metadata
  mat = mat[common,]
  #Convert count matrix in a sparse format (easier to manage) and transpose it in order to have genes as rows and cells as columns
  mat = as(as.matrix(mat),'CsparseMatrix')
  mat = Matrix::t(mat)

  #Divide transcripts from blank probes
  blmat = (mat[stringr::str_detect(rownames(mat),'Blank'),])
  mat = mat[stringr::str_detect(rownames(mat),'Blank',negate=T),]

  rownames(meta) = meta$cell_ID
  meta = meta[common,]

  if(length(poly) > 0){
    message('Loading cell polygon file...')
    poly = arrow::read_parquet(poly)
    poly = poly %>% dplyr::filter(.data$ZIndex == 3) %>%
      dplyr::select(.data$ID,.data$EntityID,.data$ZIndex,.data$Geometry,.data$ZLevel,.data$Type) %>%
      dplyr::mutate(cell_ID = paste0(dataset.id,'_',.data$EntityID)) %>%
      as.data.frame() %>%
      sf::st_as_sf()
    rownames(poly) = poly$cell_ID
    poly = suppressWarnings({sf::st_cast(poly,'POLYGON',do_split = F)})
    poly = poly[rownames(meta),]
  }else{
    warning('Merscope polygon file not found in the specified folder.\n
            Polygons will be created using cell centroid coordinates in the metadata table')
    poly = sfheaders::sf_point(meta %>% dplyr::select(.data$center_x,.data$center_y,.data$cell_ID,.data$fov),
                               x='center_x',y='center_y',polygon_id='cell_ID',keep=T)
  }
  message('Creating Seurat object...')
  #Create Seurat Object with separate CosMx and NegProbes Assays
  mers =  Seurat::CreateSeuratObject(counts = mat, meta.data = meta,assay='Merscope')
  if(dim(blmat)[1] > 0){
    mers[["blank"]] = Seurat::CreateAssayObject(blmat)
  }
  #Prepare single-cell spatial coordinates info (i.e., centroid coordinates)
  cents =  CreateCentroids(meta %>% dplyr::select(.data$center_x,.data$center_y,.data$cell_ID))
  segmentations.data <- list(
    "centroids" = cents
  )
  coords =CreateFOV(
    coords = segmentations.data,
    type = c("centroids"),
    molecules = NULL,
    assay = "Merscope"
  )

  mers[[dataset.id]] = coords
  mers@tools$poly = poly
  if(length(tx) > 0){
    message('Locating cell transcript file...')
    mers@tools$tx = normalizePath(tx)
  }
  return(mers)
}

#' @title prepare seurat object starting from any spatial transcriptomic/proteomic data in tabular format
#' @name prepare_seurat_other
#' @description
#' Given a dataframe containing single-cell marker measurements and x/y centroid coordinates, this function create a new Seurat object
#' with single-cell marker expression and metadata. If multiple independent samples are included in the dataset, there is the possibility to
#' merge sample coordinates into a unique space with parameter 'stitch =T'
#' @param data a dataframe containing single-cell measurements, metadata and spatial coordinates of cell centroids
#' @param markers.ids data column names corresponding to marker expression measurements. These will be used to build Seurat count/data matrix
#' @param xcoord data column name corresponding to cell centroid x coordinates
#' @param ycoord data column name corresponding to cell centroid y coordinates
#' @param stitch boolean, whether or not merging cell coordinates coming from independent samples within the dataset. Default is set to FALSE.
#' @param sample.id data column name reporting sample identifiers. Must be specified when parameter stitch is set to TRUE
#' @param assay.name character string that will be used to name the output Seurat assay
#' @returns a Seurat object containing cell spatial expression measurements and metadata
#' @export
prepare_seurat_other = function(data=NULL,markers.ids=NULL,xcoord=NULL,ycoord=NULL,stitch=F,sample.id=NULL,assay.name='X'){
  if(stitch==T){
    if(is.factor(data[[sample.id]])){data[[sample.id]] = as.character(data[[sample.id]])}
    data=stitch_samples(data,x=xcoord,y=ycoord,sample.id=sample.id)
  }
  mat = data[,markers.ids]
  data = data[,setdiff(colnames(data),markers.ids)]
  mat = Matrix::t(as(as.matrix(mat),'CsparseMatrix'))
  seurat = CreateSeuratObject(CreateAssayObject(counts=mat),meta.data = data, assay=assay.name)
  rm(mat)
  rm(data)
  return(seurat)
}




#' @title Prepare proseg Seurat data
#' @name prepare_proseg_seurat
#'
#' @description
#' Initial formatting of proseg re-segmented data to work with Seurat and Kandinsky
#' @details
#' This function will use proseg raw input data (generated from either CosMx, Xenium, or Merscope datasets) to build a new Seurat object. All proseg input files must be stored in the same folder that will be specified to call the function
#' @param path character string specifying the path to the proseg input files directory
#' @param dataset.id character string that will be used to name the output Seurat object identity
#' @param pattern character string that will be used as a key to identify proseg input file names in case when input files from multiple proseg samples/slides are stored in the same directory
#' @param fovs vector of FOV identifiers to subset for downstream analysis.
#' @param tech character string specifying the technology used to generate the data processed with proseg. Must be one of the following: "cosmx", "merscope", "xenium"
#' @returns a Seurat object containing proseg expected count matrix, metadata, and additional files/file paths to be used for downstream analysis with Kandinsky functions
#' @importFrom magrittr %>%
#' @importFrom arrow read_csv_arrow
#' @importFrom sfheaders sf_polygon sf_point
#' @export
#' @family prepare_data
prepare_proseg_seurat = function(path=NULL,dataset.id=NULL,pattern=NULL,fovs=NULL,tech=c('cosmx','merscope','xenium')){
  if(tech>1){tech=NULL}
  if(substr(path,nchar(path),nchar(path)) == '/'){path = substr(path,1,nchar(path)-1)}
  path = gsub("//","/",path)
  poly = list.files(path,pattern='cell-polygons\\.geojson',full.names = T)
  if(length(poly) > 1){
    poly = poly[stringr::str_detect(poly,pattern)]
  }
  meta = list.files(path,pattern='cell-metadata',full.names = T)
  if(length(meta) > 1){
    meta = meta[stringr::str_detect(meta,pattern)]
  }
  mat = list.files(path,pattern='expected-counts',full.names = T)
  if(length(mat) > 1){
    mat = mat[stringr::str_detect(mat,pattern)]
  }
  tx = list.files(path,pattern='transcript',full.names = T)
  if(length(tx) > 1){
    tx = tx[stringr::str_detect(tx,pattern)]
  }
  if(length(meta)>0){
    meta = arrow::read_csv_arrow(meta)
    meta = meta %>% dplyr::mutate(cell_ID = paste0(.data$fov,'_',.data$cell_ID))
    rownames(meta)=meta$cell_ID
  }else{
    stop('cell metadata file not found')
  }
  if(length(poly) > 0){
    json = readLines(poly)
    json = paste0(json,collapse="")
    poly = geojsonsf::geojson_sf(json)
    poly = sf::st_cast(poly,'POLYGON',do_split=F)
    poly$cell_ID = meta$cell_ID
    rownames(poly) = poly$cell_ID
  }
  if(length(mat)){
    #Load proseg expected count matrix
    mat = data.table::fread(mat,sep=',',data.table=F)
    rownames(mat) = rownames(meta)
    #Convert count matrix in a sparse format (easier to manage) and transpose it in order to have genes as rows and cells as columns
    mat = as(as.matrix(mat),'CsparseMatrix')
    mat = Matrix::t(mat)
  }else{
    stop('expected count matrix not found')
  }
  ##Split count matrix between real probes and negative/positive control probes, depending on the technology
  if(tech == 'cosmx'){
    #Divide transcripts from negative probes
    negmat = (mat[stringr::str_detect(rownames(mat),'Neg'),])
    sysmat = (mat[stringr::str_detect(rownames(mat),'System'),])
    mat = mat[stringr::str_detect(rownames(mat),'Neg|System',negate=T),]
  }else if(tech == 'merscope'){
    #Divide transcripts from blank probes
    blmat = (mat[stringr::str_detect(rownames(mat),'Blank'),])
    mat = mat[stringr::str_detect(rownames(mat),'Blank',negate=T),]
  }else if(tech == 'xenium'){
    unsmat = mat[stringr::str_detect(rownames(mat),'Unassigned'),]
    deprmat = mat[stringr::str_detect(rownames(mat),'Deprecated'),]
    negmat = mat[stringr::str_detect(rownames(mat),'NegControl'),]
    genmat = mat[stringr::str_detect(rownames(mat),'Region'),]
    mat = mat[stringr::str_detect(rownames(mat),'Unassigned|NegControl|Deprecated|Region',negate=T),]
  }else{
    stop('proseg currently works on CosMx/Merscope/Xenium data. Specify an appropriate tech argument.')
  }
  #Create Seurat Object
  proseg =  Seurat::CreateSeuratObject(counts = mat, meta.data = meta,assay=paste0(tech,'_proseg'))

  if(tech == 'cosmx'){
    if(dim(negmat)[1] > 0){
      proseg[["negprobes"]] = Seurat::CreateAssayObject(negmat)
    }
    if(dim(sysmat)[1] > 0){
      proseg[["SystemControl"]] = Seurat::CreateAssayObject(sysmat)
    }
  }
  if(tech == 'xenium'){
    if(dim(negmat)[1] > 0){
      proseg[["negprobes"]] = Seurat::CreateAssayObject(negmat)
    }
    if(dim(unsmat)[1] > 0){
      proseg[["unassigned"]] = Seurat::CreateAssayObject(unsmat)
    }
    if(dim(deprmat)[1] > 0){
      proseg[["deprecated"]] = Seurat::CreateAssayObject(deprmat)
    }
    if(dim(genmat)[1] > 0){
      proseg[["intergenic"]] = Seurat::CreateAssayObject(genmat)
    }
  }
  if(tech == 'merscope'){
    if(dim(blmat)[1] > 0){
      proseg[["blank"]] = Seurat::CreateAssayObject(blmat)
    }
  }

  #Prepare single-cell spatial coordinates info (i.e., centroid coordinates)
  cents =  CreateCentroids(meta %>% dplyr::select(.data$centroid_x,.data$centroid_y,.data$cell_ID))
  segmentations.data <- list(
    "centroids" = cents
  )
  coords =CreateFOV(
    coords = segmentations.data,
    type = c("centroids"),
    molecules = NULL,
    assay = paste0(tech,'_proseg')
  )

  proseg[[dataset.id]] = coords
  if(length(poly>0)){
    proseg@tools$poly = poly
  }
  if(length(tx) > 0){
    proseg@tools$tx = normalizePath(tx)
  }
  return(proseg)
}


#' @title Prepare Slide-Seq Seurat data
#' @name prepare_slideseq_seurat
#'
#' @description
#' Initial formatting of Slide-Seq data to work with Seurat and Kandinsky
#' @details
#' This function will use Slide-Seq input data to build a new Seurat object. All Slide-Seq input files must be stored in the same folder that will be specified to call the function.
#' Kandinsky expects to find Slide-Seq files in the same format proposed in the Slide-Seq V2 original study https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary
#' @param path character string specifying the path to the Slide-Seq input files directory
#' @param dataset.id character string that will be used to name the output Seurat object identity
#' @param pattern character string that will be used as a key to identify Slide-Seq input file names in case when input files from multiple Slide-Seq samples/slides are stored in the same directory
#' @returns a Seurat object containing Slide-Seq count matrix, metadata, and spatial coordinates to be used for downstream analysis with Kandinsky functions
#' @importFrom magrittr %>%
#' @importFrom sfheaders sf_point
#' @export
#' @family prepare_data
prepare_slideseq_seurat = function(path = NULL, dataset.id=NULL,pattern=NULL){
  if(substr(path,nchar(path),nchar(path)) == '/'){path = substr(path,1,nchar(path)-1)}
  path = gsub("//","/",path)
  files = list.files(path,full.names=T)
  mat = files[stringr::str_detect(files,'expression')]
  if(length(mat) > 1){
    mat = mat[stringr::str_detect(mat,pattern)]
  }
  spatial = files[stringr::str_detect(files,'location')]
  if(length(spatial) > 1){
    spatial = spatial[stringr::str_detect(spatial,pattern)]
  }
  mat = data.table::fread(mat,data.table = F)
  colnames(mat)[1] = 'gene'
  rownames(mat) = mat$gene
  mat$gene = NULL
  spatial = data.table::fread(spatial,header=F,data.table = F)
  if(is.character(spatial$V2[1]) | is.character(spatial$V3)[1]){
    colnames(spatial) = spatial[1,]
    spatial = spatial[-1,]
    spatial[,2] = as.numeric(spatial[,2])
    spatial[,3] = as.numeric(spatial[,3])
  }
  colnames(spatial)[c(ncol(spatial)-1, ncol(spatial))] = c('x','y')
  spatial = sfheaders::sf_point(spatial,x='x',y='y',keep=T)
  colnames(spatial)[1] = 'barcode'
  rownames(spatial) = spatial$barcode
  spatial = spatial[colnames(mat),]
  spatial[,c('x','y')] = sf::st_coordinates(spatial)
  message('Creating Seurat object...')
  #Create Seurat Object with SlideSeq data
  slseq =  Seurat::CreateSeuratObject(counts = as(as.matrix(mat),'CsparseMatrix'),assay='Slide-Seq')
  #Prepare single-cell spatial coordinates info (i.e., centroid coordinates)
  cents =  CreateCentroids(spatial[,c('x','y','barcode')])
  segmentations.data <- list(
    "centroids" = cents
  )
  coords =CreateFOV(
    coords = segmentations.data,
    type = c("centroids"),
    molecules = NULL,
    assay = "Slide-Seq"
  )
  slseq[[dataset.id]] = coords
  slseq@tools$poly = spatial
  return(slseq)
}





#' @title Prepare Singular Genomics G4X Seurat data
#' @name prepare_g4x_seurat
#'
#' @description
#' Initial formatting of Singular Genomics G4X data to work with Seurat and Kandinsky
#' @details
#' This function will use Singular Genomics G4X raw input data to build a new Seurat object. All G4X input files must be stored in the same folder that will be specified to call the function
#' @param path character string specifying the path to the G4X input files directory
#' @param dataset.id character string that will be used to name the output Seurat object identity
#' @param pattern character string that will be used as a key to identify G4X input file names in case when input files from multiple G4X samples/slides are stored in the same directory
#' @param h5 boolean, whether reading count matrix and metadata from 'feature_matrix.h5' (TRUE) file or remaining files in the 'single_cell_data/' folder (FALSE)
#' @param segmentation boolean, whether building cell polygons using G4X segmentation masks from .npz files (TRUE) or using cell centroids (FALSE). Reading .npz files requires the use of python numpy library through reticulate
#' @returns a Seurat object containing G4X count matrix, metadata, and additional files/file paths to be used for downstream analysis with Kandinsky functions
#' @importFrom magrittr %>%
#' @importFrom arrow read_csv_arrow
#' @importFrom sfheaders sf_polygon sf_point
#' @export
#' @family prepare_data
prepare_g4x_seurat = function(path = NULL, dataset.id=NULL,pattern=NULL,h5=T,segmentation=F){
  if(substr(path,nchar(path),nchar(path)) == '/'){path = substr(path,1,nchar(path)-1)}
  path = gsub("//","/",path)
  sc_data = list.files(paste0(path,'/single_cell_data'),full.names=T)
  img_data = list.files(paste0(path,'/h_and_e'),full.names=T)
  tx = list.files(paste0(path,'/rna'),full.names=T)
  message('Preparing G4X transcript and protein expression data...')
  if(h5 ==T){
    data = sc_data[stringr::str_detect(sc_data,'feature_matrix.h5')]
    mat = rhdf5::h5read(data,name='X')
    mat = as(mat,'CsparseMatrix')
    obs = as.data.frame(rhdf5::h5read(data,name='obs'))
    rownames(obs) = obs$cell_id
    var = as.data.frame(rhdf5::h5read(data,name='var'))
    rownames(var) = var$gene_id
    rownames(mat) = rownames(var)
    colnames(mat) = rownames(obs)
    message('Creating Seurat object...')
    g4x = Seurat::CreateSeuratObject(counts = mat,meta.data=obs,assay='G4X')
    #Prepare single-cell spatial coordinates info (i.e., centroid coordinates)
    cents =  CreateCentroids(obs[,c('cell_x','cell_y','cell_id')])
    segmentations.data <- list(
      "centroids" = cents
    )
    coords =CreateFOV(
      coords = segmentations.data,
      type = c("centroids"),
      molecules = NULL,
      assay = "G4X"
    )
    g4x[[dataset.id]] = coords
  }else{
    meta = as.data.frame(arrow::read_csv_arrow(sc_data[stringr::str_detect(sc_data,'cell_metadata.csv.gz')]))
    rownames(meta) = meta$label
    prot_mat = as.data.frame(arrow::read_csv_arrow(sc_data[stringr::str_detect(sc_data,'cell_by_protein.csv.gz')]))
    rownames(prot_mat) = prot_mat$label
    prot_mat$label = NULL
    tx_mat = as.data.frame(arrow::read_csv_arrow(sc_data[stringr::str_detect(sc_data,'cell_by_transcript.csv.gz')]))
    rownames(tx_mat) = tx_mat$label
    tx_mat$label = NULL
    tx_mat = as(as.matrix(tx_mat),'CsparseMatrix')
    neg_p = colnames(tx_mat)[stringr::str_detect(colnames(tx_mat),'NCP.')]
    neg_s = colnames(tx_mat)[stringr::str_detect(colnames(tx_mat),'NCS.')]
    negp_mat = tx_mat[,neg_p]
    negs_mat = tx_mat[,neg_s]
    tx_mat = tx_mat[,!(colnames(tx_mat) %in% c(neg_p,neg_s))]
    g4x = Seurat::CreateSeuratObject(counts = Matrix::t(tx_mat),meta.data=meta,assay='G4X')
    g4x@meta.data$ID = sapply(rownames(g4x@meta.data),function(x){strsplit(x,split='-')[[1]][[2]]})
    rm(tx_mat)
    message('Creating Seurat object...')
    g4x[["protein"]] = Seurat::CreateAssayObject(Matrix::t(prot_mat))
    rm(prot_mat)
    g4x[["NCP"]] = Seurat::CreateAssayObject(Matrix::t(negp_mat))
    g4x[["NCS"]] = Seurat::CreateAssayObject(Matrix::t(negs_mat))
    rm(neg_p)
    rm(neg_s)
    cents =  CreateCentroids(meta[,c('cell_x','cell_y','cell_id')])
    segmentations.data <- list(
      "centroids" = cents
    )
    coords =CreateFOV(
      coords = segmentations.data,
      type = c("centroids"),
      molecules = NULL,
      assay = "G4X"
    )
    g4x[[dataset.id]] = coords
  }
  message('Preparing cell polygons...')
  if(segmentation == T){
    seg_data = list.files(paste0(path,'/segmentation'),full.names=T)
    seg_file = seg_data[stringr::str_detect(seg_data,'segmentation_mask.npz')]
    if(rlang::is_installed('reticulate')){
      check_numpy = reticulate::py_list_packages() %>% dplyr::filter(.data[["package"]] == 'numpy')
      if(nrow(check_numpy)==0){
        warning('Please install numpy in order to read G4X segmentation file.
                Using cell centroid coordinates instead of segmentation masks..')
        poly = sfheaders::sf_point(g4x@meta.data,y='cell_x',x='cell_y',keep=T)
      }else{
        np = reticulate::import('numpy')
        seg = np$load(seg_file)
        seg = seg$f[["nuclei"]]
        seg = as(seg,'CsparseMatrix')
        poly = Matrix::which(seg != 0, arr.ind = TRUE) %>% as.data.frame()
        poly$ID = seg@x
        rm(seg)
        poly = poly %>% dplyr::arrange(.data[["ID"]])
        poly = sfheaders::sf_polygon(poly,x='col',y='row',polygon_id ='ID',keep=T)
        poly = merge(g4x@tools$poly,g4x@meta.data,by='ID')
        poly$label = rownames(g4x@meta.data)
        rownames(poly) = poly$label
      }
    }else{
      warning('Please install reticulate and numpy in order to read G4X segmentation file.
                Using cell centroid coordinates instead of segmentation masks..')
      poly = sfheaders::sf_point(g4x@meta.data,y='cell_x',x='cell_y',keep=T)
    }
  }else{
    poly = sfheaders::sf_point(g4x@meta.data,y='cell_x',x='cell_y',keep=T)
  }
  if(length(tx) > 0){
    message('Locating cell transcript file...')
    g4x@tools$tx = normalizePath(tx)
  }
  if(length(img_data)>0){
    message('Locating G4X H&E image...')
    img = img_data[stringr::str_detect(img_data,'h_and_e.jp2')]
    g4x@tools$img_path = normalizePath(img)
  }
  g4x@tools$poly = poly
  return(g4x)
}




