#' @importFrom methods setClass setClassUnion setMethod is slotNames slot slot<- new as show getClass
#' @importFrom stats sd quantile kmeans relevel median dist density na.omit chisq.test p.adjust
#' @importFrom rlang .data
#' @importClassesFrom terra SpatRaster SpatRasterCollection
#' @importClassesFrom Matrix CsparseMatrix ddiMatrix
#' @importFrom terra RGB RGB<- coltab<- as.int xmax<- xmin<- ymax<- ymin<-
#' @import Seurat SeuratObject
NULL

##Add spatstat.geom functions to namespace in case spatstat will be integrated in Kandinsky
#spatstat.geom as.ppp marks ppp marks<-

#' @name KanData
#' @title Show Kandinsky data
#' @description
#' Extract all Kandinsky data or specific slots from Seurat object
#'
#' @param object Seurat or Kandinsky object
#' @param which name of Kandinksy slots to extract. If `NULL`, all Kandinsky data will be extracted
#' @param value object or variable to assign to the selected Kandinksy slot
#' @export
setGeneric(name = "KanData", def = function(object, which=NULL) standardGeneric("KanData"))


#' @title Assign data to Kandinsky
#' @name KanData<-
#' @description
#' Assign new data to Kandinsky data or Kandinsky slots
#'
#' @param object Seurat or Kandinsky object
#' @param which name of Kandinksy slots to extract. If `NULL`, all Kandinsky data will be extracted
#' @param value object or variable to assign to the selected Kandinksy slot
#' @export
setGeneric(name = "KanData<-", def = function(object, which=NULL,value) standardGeneric("KanData<-"))




#' @title Save Kandinsky data
#' @name saveKanData
#' @description
#' saveRDS adaptation for Kandinsky data. it can save Kandinsky data alone or as part of a Seurat object
#'
#' @param object a Seurat object containing Kandinsky data, or a Kandinsky object
#' @param file character string specifying the path to the file where to save Kandinsky data. It must be in .rds format
#' @param with.seurat boolean, whether saving the whole Seurat object together with Kandinsky data (`TRUE`) or Kandinsky data alone (`FALSE`). Default value is set to `TRUE`. Ignored when `object` corresponds to a Kandinsky object
#' @export
setGeneric(name = "saveKanData", def = function(object,file=NULL,with.seurat=T) standardGeneric("saveKanData"))


#' @title Read Kandinsky data
#' @name readKanData
#' @description
#' readRDS adaptation for Kandinsky data. it can save Kandinsky data alone or as part of a Seurat object
#'
#' @param file character string specifying the path to the file where Seurat or Kandinsky object is stored. It must be in .rds format
#' @param with.seurat boolean, whether loading the whole Seurat object together with Kandinsky data (`TRUE`) or Kandinsky data alone (`FALSE`). Default value is set to `TRUE`. Ignored when input file contains only a Kandinsky object
#'
#' @export
setGeneric(name = "readKanData", def = function(file=NULL,with.seurat=T) standardGeneric("readKanData"))



#' @title Add Kandinsky data from .rds file to Seurat object
#' @name addKanData
#' @description
#' read .rds file containing Kandinsky object, and add the obejct to already existing Seurat object
#'
#' @param object a Seurat object
#' @param kanfile character string specifying the path to the file where Kandinsky object is stored. It must be in .rds format
#'
#' @export
setGeneric(name = "addKanData", def = function(object, kanfile=NULL) standardGeneric("addKanData"))


#' @title Mask Visium H&E image
#' @name he_mask
#' @description create a mask image to remove image background from Visium H&E images
#' @details
#' Additional details...
#' This function will apply a variance-based threshold to the RGB channels of the H&E image stored within the Kandinsky data
#' to keep only the area of image covering the tissue sample and to remove the rest of the image background. The masking strategy relies on the fact
#' that RGB values from H&E staining will usually show a higher standard deviation compared to the gray-ish Visium slide background.
#' @param object a Seurat object containing Kandinsky data (`KanData()`) or an image of class `SpatRaster` created with the `terra` package
#' @param maxres numeric value indicating the ideal maximum pixel dimension accepted for Visium H&E image. If image size is bigger than this value, a scale factor will be applied to lower image pixel resolution/size up to the defined maximum value
#' @param stretch boolean value, whether or not to apply a linear stretch to increase masked image contrast. Default is FALSE
#' @param sd_thresh numeric value indicating the minimum RGB standard deviation accepted to keep an image pixel as part of the final tissue mask. Lower `sd_thresh` values will give more permissive filtering results.
#' @param crop boolean value, whether or not to crop the final masked image
#' @param crop_area bounding box coordinates to be used to crop masked image. If `crop = T` and `crop_area` is set to `NULL`, Visium spot coordinates will be used to automatically define the bounding box for cropping the image
#' @return Seurat object with new 'masked_img' slot within the Kandinsky data if object is a Seurat object, or a masked `Spatraster` image
#' @rdname he_mask
#' @export
setGeneric(name = "he_mask", def = function(object,maxres=2000,stretch=F,sd_thresh=9,crop=T,crop_area=NULL) standardGeneric("he_mask"))

if (requireNamespace("terra", quietly = TRUE)) {
  PackedSpatRaster = getClass('PackedSpatRaster',where=asNamespace('terra'))
}


#PackedSpatRaster class from terra is not recognized for import.
#I will explicitly extract it from terra and I include it in the superclass KanImg
#Maybe I could even do it for nb/linstw/sf object classes...

#' @title S4 Class Union KanImg
#' @description
#' Create new class for Kandinsky images
#'
#' @importClassesFrom terra SpatRaster SpatRasterCollection
#' @export
PackedSpatRaster = function() {
  suppressWarnings({
    if (requireNamespace("terra", quietly = TRUE)) {
      PackedSpatRaster = getClass('PackedSpatRaster',where=asNamespace('terra'))
    }
    setClassUnion("KanImg",c("NULL","SpatRaster","SpatRasterCollection",PackedSpatRaster@className))
  })
}

PackedSpatRaster()


#' Kandinsky S4 object class
#'
#' Kandinsky is a list object storing spatial data and information for spatial transcriptomic dataset analysed through Seurat.
#'
#' @aliases Kandinsky Kandinsky-class
#' @slot platform sequencing technology
#' @slot nb weighted neighbour list object created with the package`spdep`
#' @slot sf geometry data frame created with the package `sf` and containing cell/spot coordinates
#' @slot nb.type character indicating the neighbouring strategy used to initialize Kandinsky object. Must be one of the following: `Q`, `K`, `C`, `M`
#' @slot nnMat list object containing nearest neighbour matrix created with one of the `nnMat` functions implemented in Kandinsky
#' @slot spot_distance numeric value indicating centroid-to-centroid distance between Visium spots or Visium HD bins
#' @slot texture data.frame containing GLCM texture features extracted from Visium H&E slide image
#' @slot tx character string indicating the path for the transcript coordinates file from CoxMx/Xenium/Merscope platforms
#' @slot fov_mask geometry data frame created with the package `sf` and containing CosMx FOVs bounding box coordinates
#' @slot img Visium H&E or any other tissue slide image stored as a SpatRaster or SpatrasterCollection object via the package `terra`
#' @slot mask masked H&E image obtained via the Kandinsky function `he_mask` and stored as a SpatRaster or SpatrasterCollection object via the package `terra`
#' @export
#' @exportClass Kandinsky
Kandinsky <- setClass(
  Class = "Kandinsky",
  slots = c(
    "platform" = c("character"),
    "nb" = c('ANY'),
    #    "listw" = c('ANY'),
    "sf" = c("ANY"),
    "nb.type" = c("character"),
    "spot_distance" = c("numeric"),
    "texture" = c("data.frame"),
    "nnMat" = c("list"),
    "tx" = c("character"),
    "fov_mask" = c('ANY'),
    "img" = c("KanImg"),
    "mask" = c("KanImg")
    #    "ppp" = c("ppp")
  )
) -> Kandinsky




#####Maybe to pair with a as.Kandinsky() function for Seurat data
#kan2seurat = function(kan){
#  new = suppressWarnings({Seurat::as.Seurat(Seurat::as.SingleCellExperiment(kan))})
#  slots = slotNames(kan)
#  for(s in slots){
#    if(s != 'kandata'){
#      slot(new,s) = slot(kan,s)
#    }else{
#    new@tools$kandata = slot(kan,s)
#    }
#  }
#  return(new)
#}




#' @rdname KanData
#' @export KanData
setMethod(f='KanData',signature='Kandinsky',
          function(object,which=NULL){if(is.null(which)){object}else{slot(object,which)}}
)

#' @rdname KanData
#' @export KanData
setMethod(f='KanData',signature='Seurat',
          function(object,which=NULL){if(is.null(which)){object@tools$kandata}else{slot(object@tools$kandata,which)}}
)

#' @rdname KanData
#' @export KanData<-
setMethod(f='KanData<-',signature='Kandinsky',
          function(object,which=NULL,value){slot(object,which) = value; object}
)

#' @rdname KanData
#' @export KanData<-
setMethod(f='KanData<-',signature='Seurat',
          function(object,which=NULL,value){slot(object@tools$kandata,which) = value; object}
)



#setMethod(f='saveRDS',signature='Kandinsky',
#         function(object = NULL,file = "", ascii = FALSE, version = NULL,
#                 compress = TRUE, refhook = NULL){
#            if(!is.null(object@img)){
#              object@img = terra::wrap(object@img)
#            }
#            if(!is.null(object@mask)){
#              object@img = terra::wrap(object@img)
#            }
#           base::saveRDS(object,file, ascii = ascii, version = version,
#                          compress = compress, refhook = refhook)
#          }
#)



#setMethod(f='readRDS',signature='character',
#          function(file="",refhook = NULL){
#            data = base::readRDS(file,refhook = refhook)
#            if(inherits(data,'Kandinsky')){
#              if(!is.null(object@img)){
#              object@img = terra::wrap(object@img)
#            }
#            if(!is.null(object@mask)){
#              object@img = terra::wrap(object@img)
#            }
#            }
#            return(data)
#          }
#)

#' @rdname saveKanData
#' @export saveKanData
setMethod(f='saveKanData',signature='Kandinsky',
          function(object,file=NULL,with.seurat=T){
            if(!is.null(object@img)){object@img = terra::wrap(object@img)}
            if(!is.null(object@mask)){object@mask = terra::wrap(object@mask)}
            saveRDS(object,file)
          }
)


#' @rdname saveKanData
#' @export saveKanData
setMethod(f='saveKanData',signature='Seurat',
          function(object,file=NULL,with.seurat=T){
            if(any(Tool(object) == 'kandata')){
              if(!is.null(KanData(object,'img'))){KanData(object,'img') = terra::wrap(KanData(object,'img'))}
              if(!is.null(KanData(object,'mask'))){KanData(object,'mask') = terra::wrap(KanData(object,'mask'))}
              if(with.seurat==T){
                saveRDS(object,file)
              }else{
                saveRDS(KanData(object),file)
              }
            }else{
              saveRDS(object, file)
            }
          }
)

##readRDS adaptation for Kandinsky data. it can load Kandinsky data alone or as part of a Seurat object
#' @rdname readKanData
#' @export readKanData
setMethod(f='readKanData',signature='character',
          function(file=NULL,with.seurat=T){
            obj = base::readRDS(file)
            if(inherits(obj,'Seurat')){
              if(any(Tool(obj) == 'kandata')){
                if(!is.null(KanData(obj,'img'))){KanData(obj,'img') = terra::rast(KanData(obj,'img'))}
                if(!is.null(KanData(obj,'mask'))){KanData(obj,'mask') = terra::rast(KanData(obj,'mask'))}
                if(with.seurat==T){
                  return(obj)
                }else{
                  return(KanData(obj))
                }
              }
            }else if(inherits(obj,'Kandinsky')){
              if(!is.null(KanData(obj,'img'))){KanData(obj,'img') = terra::rast(KanData(obj,'img'))}
              if(!is.null(KanData(obj,'mask'))){KanData(obj,'mask') = terra::rast(KanData(obj,'mask'))}
              return(obj)
            }else{
              stop('Data is not a Seurat or Kandinsky object.')
            }
          }
)

#' @rdname addKanData
#' @export addKanData
setMethod(f='addKanData',signature=c('Seurat'),
          function(object,kanfile=NULL){
            kan = readKanData(file=file,with.seurat=F)
            KanData(object) = kan
            return(object)
          }
)





###HE masking function based on RGB channel standard deviation or custom channel thresholds
####Quite robust and flexible even in the case of weird tissue shapes and staining artifacts
####Only tested for HE staining, no IF or IHC

#setMethod(f='he_mask',signature = c('Kandinsky'),
#          function(object,maxres=2000,stretch=F,sd_thresh=9,crop=T,crop_area=NULL){
#            rast = KanData(object,'img')
#            scalef = round(max(dim(rast))/maxres,1)
#            if(scalef >= 2){
#              rast = terra::aggregate(rast,fact=scalef,fun='mean')
#              }
#            RGB(rast) = c(1,2,3)
#            names(rast) = c('R','G','B')
#          #  if(filter_style == 'rgb'){
#           # if(strong_filter ==T){
#            #  mask =
#             #   (abs(rast$R - rast$G) >=35) |
#              #  (abs(rast$R - rast$G) >=20 & rast$B <= 210) |
#               # #(abs(rast$R - rast$G) >=10 & abs(rast$R - rast$B) <=5) |
#                #(abs(rast$R - rast$G) <=20 & rast$B - rast$R >= 30) |
#                #(abs(rast$R - rast$G) <=8 & abs(rast$R - rast$B) >= 8 & rast$R < 190) |
#                #(abs(rast$R - rast$G) >=10 & rast$R > 230 & rast$B > 210)
#            #}else{
#             # mask = (abs(rast$R - rast$G) >=35) |
#              #  (abs(rast$R - rast$G) >=20 & rast$B <= 210) |
#               # #(abs(rast$R - rast$G) >=10 & abs(rast$R - rast$B) <=5) |
#                #(abs(rast$R - rast$G) <=20 & rast$B - rast$R >= 30) |
#                #(abs(rast$R - rast$G) <=8 & abs(rast$R - rast$B) >= 8 & rast$R < 190) |
#                #(abs(rast$R - rast$G) >=10 & rast$R > 220 & rast$B > 210)
#            #}}
#            #else if(filter_style == 'sd'){
#            sd_rast = rast
#            sd_rast$SD = apply(sd_rast[,,1:3],1,sd)
#            mask = sd_rast$SD > sd_thresh
#            rm(sd_rast)
#            #}
#            mask = as.int(mask)
#            if(stretch==T){
#              KanData(object,'mask') = terra::stretch(terra::mask(rast,mask,maskvalues=0))
#            }else{
#              KanData(object,'mask') = terra::mask(rast,mask,maskvalues=0)
#            }
#            if(crop==T){
#              if(is.null(crop_area)){
#                KanData(object,'mask') = terra::crop(KanData(object,'mask'),sf::st_as_sfc(sf::st_bbox(KanData(object,'sf'))))
#              }else{
#                KanData(object,'mask') = terra::crop(KanData(object,'mask'),crop_area)
#              }
#            }
#            return(object)
#          }
#)


#' @rdname he_mask
#' @aliases he_mask, Seurat_he_mask
#' @export he_mask
setMethod(f='he_mask',signature = c('Seurat'),
          function(object,maxres=2000,stretch=F,sd_thresh=9,crop=T,crop_area=NULL){
            rast = KanData(object,'img')
            scalef = round(max(dim(rast))/maxres,1)
            if(scalef >= 2){
              rast = terra::aggregate(rast,fact=scalef,fun='mean')
            }
            RGB(rast) = c(1,2,3)
            names(rast) = c('R','G','B')
            sd_rast = rast
            sd_rast$SD = apply(sd_rast[,,1:3],1,sd)
            mask = sd_rast$SD > sd_thresh
            rm(sd_rast)
            mask = as.int(mask)
            if(stretch==T){
              KanData(object,'mask') = terra::stretch(terra::mask(rast,mask,maskvalues=0))
            }else{
              KanData(object,'mask') = terra::mask(rast,mask,maskvalues=0)
            }
            if(crop==T){
              if(is.null(crop_area)){
                KanData(object,'mask') = terra::crop(KanData(object,'mask'),sf::st_as_sfc(sf::st_bbox(KanData(object,'sf'))))
              }else{
                KanData(object,'mask') = terra::crop(KanData(object,'mask'),crop_area)
              }
            }
            return(object)
          }
)

#' @rdname he_mask
#' @aliases he_mask, SpatRaster_he_mask
#' @export he_mask
setMethod(f='he_mask',signature = c('SpatRaster'),
          function(object,maxres=2000,stretch=F,sd_thresh=9,crop=T,crop_area=NULL){
            rast = object
            scalef = round(max(dim(rast))/maxres,1)
            if(scalef >= 2){
              rast = terra::aggregate(rast,fact=scalef,fun='mean')
            }
            RGB(rast) = c(1,2,3)
            names(rast) = c('R','G','B')
            sd_rast = rast
            sd_rast$SD = apply(sd_rast[,,1:3],1,sd)
            mask = sd_rast$SD > sd_thresh
            rm(sd_rast)
            mask = as.int(mask)
            if(stretch==T){
              rast = terra::stretch(terra::mask(rast,mask,maskvalues=0))
            }else{
              rast = terra::mask(rast,mask,maskvalues=0)
            }
            if(crop==T){
              rast = terra::crop(rast,crop_area)
            }
            return(rast)
          }
)



###TO OPTIMIZE, I DON'T HOW IT SHOULD WORK EXACTLY....
mergeKanData = function(seurat.list){
  nsamples = length(seurat.list)
  ids = names(seurat.list)
  dim = ceiling(sqrt(nsamples))
  slots = c('nb','sf','listw','img','mask','texture')
  #merge sf
  sf = list()
  for(n in 1:nsamples){
    if(!is.null(KanData(seurat.list[[n]],'sf'))){
      sf[[n]] = KanData(seurat.list[[n]],'sf')
      sf[[n]]$kan_id = paste0(rownames(sf[[n]]),"_",ids[n])
    }
  }
  if(length(sf) > 0){
    sf = purrr::reduce(sf,'rbind')
  }else{sf = NULL}

  #merge listw and create nb objects
  listw = list()
  listw_ids = c()
  for(n in 1:nsamples){
    if(!is.null(KanData(seurat.list[[n]],'listw'))){
      listw[[n]] = as(KanData(seurat.list[[n]],'listw'),'CsparseMatrix')
      listw_ids = c(listw_ids,paste0(rownames(listw[[n]]),"_",ids[n]))
    }
  }
  if(length(listw) > 0){
    merged_listw = Matrix::bdiag(listw)
    rownames(merged_listw) = listw_ids
    colnames(merged_listw) = listw_ids
    merged_listw = spdep::mat2listw(merged_listw,row.names = listw_ids,style='B',zero.policy=T)
    nb = merged_listw[[2]]
  }else{
    listw = NULL
    nb = NULL
  }

  img = list()
  imgnames = c()
  for(n in 1:nsamples){
    if(!is.null(KanData(seurat.list[[n]],'img'))){
      img[[n]] = KanData(seurat.list[[n]],'img')
      imgnames = c(imgnames,ids[n])
    }
  }
  if(length(img) > 0){
    img = terra::sprc(img)
    names(img) = imgnames
  }else{img = NULL}

  mask = list()
  masknames = c()
  for(n in 1:nsamples){
    if(!is.null(KanData(seurat.list[[n]],'mask'))){
      mask[[n]] = KanData(seurat.list[[n]],'mask')
      masknames = c(masknames,ids[n])
    }
  }
  if(length(mask) > 0){
    mask = terra::sprc(mask)
    names(mask) = masknames
  }else{mask = NULL}

  texture = list()
  for(n in 1:nsamples){
    if(!is.null(KanData(seurat.list[[n]],'texture'))){
      texture[[n]] = KanData(seurat.list[[n]],'texture')
      texture[[n]]$kan_id = paste0(rownames(texture[[n]]),"_",ids[n])
    }
  }
  if(length(texture) > 0){
    texture = purrr::reduce(texture,rbind)
  }else{texture = NULL}

  kandinsky = list(
    nb = nb,
    sf = sf,
    listw = merged_listw,
    img = img,
    mask = mask,
    texture=texture
  )
  return(kandinsky)
}

as.kandinsky = function(list){

}


#' @title Initialize Kandinsky data
#' @name kandinsky_init
#' @description initialize Kandinsky data starting from an already existing Seurat object.
#' @details
#' The new data generated with this function will be added to the tool slot of Seurat object
#' as a new 'kandata' list object. `kandinsky_init` will create a fixed set of input data required for downstream
#' spatial analysis, plus additional objects depending on the specific sequencing platform generating the data stored in Seurat.
#' If Seurat object has been created through one of the `prepare_*_seurat` functions, some of the input files required by `kandinsky_init`
#' will be automatically detected even if the corresponding parameters are set to `NULL`.
#' @param seurat a Seurat object
#' @param tech character string specifying the platform used to generate the sequencing data. Must be one of the following: `visium`, `visium_hd`, `cosmx`, `xenium`, `merscope`, `other`.
#' @param img character string indicating the path of the H&E image to be used for Visium/VisiumHD. If set to `NULL`, no image is loaded and stored in the final Kandinsky data
#' @param img_maxdim numeric value indicating the ideal maximum pixel dimension accepted for Visium H&E image. If image size is bigger than this value, a scale factor will be applied to lower image pixel resolution/size up to the defined maximum value
#' @param res character string specifying which version of H&E Visium image will be loaded into the Kandinsky data. Must be one of the following: `low`, `high`, `full`, where `full` refers to the original H&E full-resolution tiff image
#' @param binsize numeric value indicating which bin resolution will be considered to create the final sf object when working with Visium-HD data
#' @param tx_path character string indicating the path for the transcript coordinates file from CoxMx/Xenium/Merscope platforms
#' @param fov_path character string indicating the path for the fov coordinates file from CoxMx platform
#' @param poly_path character string indicating the path for the polygon coordinates file from CoxMx/Xenium platform
#' @param nb.method character string specifying the method to be used to create a `nb` neighbour object.
#' Must be one of the following:
#' 'Q': queen contiguity method,check for contact (not overlap) between any edge or side od two polygons (refers to the queen movement rule in chess). Currently only applicable for Visium/Visium-HD data
#' 'C': centroid-based method, use maximum centroid distance threshold to identify spot/cell neighbours
#' 'K': KNN method, define k closest neighbours to each spot/cell
#' 'M': membrane-based method, check for the occurrence of a physical contact/intersection within a distance threshold between cell boundaries. Not applicable in the case of Visium spots.
#' @param k numeric, number of nearest neighbours to be set when `nb.method = K`
#' @param d.max numeric, maximum centroid distance threshold to be set when `nb.method = C | M`
#' @param hd.snap numeric, scaling factor used on minimal distance between Visium-HD bins to define neighbour relationships. Only Applied when `nb.method = Q`. Higher `hd.snap` values will give lower distance thresholds.
#' @param ids_other character string specifying variable name to be used as cell identifiers when argument tech is set to "other"
#' @param xcoord_other character string specifying variable name to be used as x coordinates when argument tech is set to "other"
#' @param ycoord_other character string specifying variable name to be used as y coordinates when argument tech is set to "other"
#' @export
kandinsky_init = function(seurat=NULL,tech='visium',
                          img=NULL,img_maxdim=2000,
                          res=c('low','high','full'),
                          binsize=16,
                          tx_path=NULL,
                          fov_path = NULL,
                          poly_path = NULL,
                          nb.method = c('Q','C','K','M'),
                          k=20,d.max=40,
                          hd.snap=10,
                          ids_other = NULL,
                          xcoord_other=NULL,
                          ycoord_other=NULL){
  if(!inherits(seurat,'Seurat')){
    stop('data must be a Seurat object')
  }
  tech = intersect(tech,c('visium','visium_hd','cosmx','xenium','merscope','other'))
  if(length(tech) ==0){
    stop('You must specify a technology compatible with Kandinsky: visium/visium_hd/cosmx/xenium/merscope/other')
  }
  kandinsky = list()
  kandinsky$platform = tech
  if(!is.null(img)){
    message('Building Kandinsky slot "img"...')
    kandinsky$img = visium2rast(seurat,img_path = img,rm_old_img = F,return.seurat=F,max_dim=img_maxdim)
  }else if(is.null(img) & !is.null(seurat@tools$img_path)){
    message('Building Kandinsky slot "img"...')
    kandinsky$img = visium2rast(seurat,img_path = seurat@tools$img_path,rm_old_img = F,return.seurat=F,max_dim=img_maxdim)
  }
  message('Building Kandinsky slot "sf"...')
  if(tech == 'visium_hd'){
    kandinsky$sf = visium2sf(seurat,return.seurat=F,is.hd=T,res=res,binsize=binsize,img=kandinsky$img)
  }else if(tech == 'visium'){
    kandinsky$sf = visium2sf(seurat,return.seurat=F,is.hd=F,res=res,img=kandinsky$img)
  }else if(tech %in% c('cosmx','xenium','merscope')){
    if(is.null(seurat@tools$poly)){
      kandinsky$sf = smi2sf(seurat = seurat,poly_file = poly_path,id='cell_ID',return.seurat = F)
    }else{
      kandinsky$sf = seurat@tools$poly[colnames(seurat),]
      seurat@tools$poly = NULL
    }
  }
  if(tech == 'other'){
    kandinsky$sf = sfheaders::sf_point(seurat@meta.data,x=xcoord_other,y=ycoord_other,keep=T)
    rownames(kandinsky$sf) = rownames(seurat@meta.data)
  }
  message('Building Kandinsky slot "nb"...')
  if(tech %in% c('visium','visium_hd')){
    ##Calculate minimal distance between any adjacent spot or bin pairs
    knn = spdep::knearneigh(kandinsky$sf,k = 1)
    knn = spdep::knn2nb(knn,row.names = kandinsky$sf[['spot_ID']])
    mindist = min(unlist(lapply(spdep::nbdists(knn,sf::st_coordinates(kandinsky$sf)),min)))
    #Use minimal distance to expand spot/bin centroids to their real size
    if(tech == 'visium'){
      diameter = (mindist/100)*65
      kandinsky$sf = sf::st_buffer(kandinsky$sf,dist=diameter/2)
      if(nb.method=='Q'){
        ##Add snap term for contiguity check. The real spot/spot distance should be ~35um, but we can set the snap term to 40um
        ## In this way we account for any spot misalignment (i.e., spot-spot distance might be slightly higher than 35)
        snap = (mindist/100)*40
      }
    }else{
      kandinsky$sf = sf::st_buffer(kandinsky$sf,dist=mindist/2,endCapStyle = 'SQUARE')
      if(nb.method=='Q'){
        ##Add snap term for contiguity check. There shouldn't be almost any empty space between adjacent Visium-HD bins.
        ##Therefore, we can set a very small snap term (= 1/10 of minimal bin distance)
        snap = (mindist/hd.snap)
      }
    }
  }

  if(tech %in% c('visium','visium_hd') & !(nb.method %in% c('Q','C','K'))){
    stop('For Visium/Visium-HD data, only the following neighbour methods can be selected: Q(Queen), C(centroid distance), K(K-nearest neighbour)')
  }

  if(length(nb.method) > 1){
    nb.method = NULL
  }
  if(is.null(nb.method)){
    warning('No neighbour method specified by the user. Setting default method to K for single-cell data or Q for spot-based data')
  }
  if(tech %in% c('visium','visium_hd')){
    nb.method = nb.method %||% 'Q'
  }else{
    nb.method = nb.method %||% 'K'
  }

  if(nb.method == 'K'){
    kandinsky$nb = knn_nb(kandinsky$sf,k = k)
    kandinsky$nb.type = paste0('K_',k)
  }else if(nb.method=='C'){
    kandinsky$nb = centroid_nb(kandinsky$sf,d.max=d.max)
    kandinsky$nb.type = paste0('C_',d.max)
  }else if(nb.method == 'M'){
    kandinsky$nb = membrane_nb(kandinsky$sf,d.max=d.max)
    kandinsky$nb.type = paste0('M_',d.max)
  }else if(nb.method == 'Q'){
    kandinsky$nb = queen_nb(kandinsky$sf,snap=snap)
    kandinsky$nb.type = 'Q'
  }else{
    stop('nb.method parameter must be either "Q", "C", "K", or "M"')
  }
  if(is.null(seurat@tools$tx)){
    message('Building Kandinsky slot "tx"...')
    if(is.null(tx_path)){
      kandinsky$tx = NULL
    }else{
      kandinsky$tx = normalizePath(tx_path)
    }
  }else{
    message('Building Kandinsky slot "tx"...')
    kandinsky$tx = seurat@tools$tx
    seurat@tools$tx = NULL
  }
  if(is.null(seurat@tools$fov)){
    message('Building Kandinsky slot "fov_mask"...')
    if(is.null(fov_path)){
      kandinsky$fov_mask = NULL
    }else{
      fov_info = intersect(c('fov','fov_name'),colnames(seurat@meta.data))
      kandinsky$fov_mask = read_fovfile(fov_path,buffer=255,fov_ids = unique(seurat@meta.data[[fov_info]]))
    }
  }else{
    message('Building Kandinsky slot "fov_mask"...')
    kandinsky$fov_mask = seurat@tools$fov
    seurat@tools$fov = NULL
  }

  if(tech %in% c('visium','visium_hd')){
    kandinsky$spot_distance = mindist
  }
  seurat@tools$kandata = Kandinsky()
  for(n in names(kandinsky)){
    KanData(seurat,n) = kandinsky[[n]]
  }
  message('Kandinsky data loaded into Seurat object!')
  return(seurat)
}


###Define a function do
###(i.e., if some cells/spots have been filtered out after initiation of Kandinsky, you will need to update Kandinsky data to take it into account)
#' @title Update Kandinsky data
#' @name update_kandinsky
#' @description data objects included in the KanData slots based on the real Seurat object.
#' By default Kandinsky data will not be subsetted together with Seurat during normal filtering steps using `subset()` or row/column names subset (e.g. `seurat[,ids]`).
#' Use `update_kandinsky()` to check and fix discrepancies between Seurat and KanData cell/spot ids
#' @param seurat a Seurat object containing Kandinsky data (`KanData()`)
#' @returns Seurat object with updated Kandinsky data
#' @export
setGeneric(name = "update_kandinsky", def = function(seurat=NULL) standardGeneric("update_kandinsky"))


#' @rdname update_kandinsky
#' @aliases update_kandinsky, Seurat_update_kandinsky
#' @export update_kandinsky
setMethod(f='update_kandinsky',signature = c('Seurat'),function(seurat=NULL){
  ids = colnames(seurat)
  if(length(setdiff(rownames(KanData(seurat,'sf')),ids)) == 0){
    return(seurat)
  }else{
    message('Updating Kandinsky data')
    KanData(seurat,'sf') = KanData(seurat,'sf')[ids,]
    mat = as(KanData(seurat,'nb'),'CsparseMatrix')
    mat = mat[ids,ids]
    suppressWarnings({
      small_nb = spdep::mat2listw(mat,row.names=rownames(mat),style = 'B',zero.policy=T)
    })
    KanData(seurat,'nb') = small_nb
    rm(mat)
    rm(small_nb)
    if(length(KanData(seurat,'nnMat')) >0){
      KanData(seurat,'nnMat')$nnMat = KanData(seurat,'nnMat')$nnMat[ids,]
    }
    if(nrow(KanData(seurat,'texture')) >0){
      KanData(seurat,'texture') = KanData(seurat,'texture')[ids,]
    }
    return(seurat)
  }
  }
)




##https://guide.macports.org/chunked/installing.macports.html be sure to have MacPorts installed on your mac
#then,from the terminal --> sudo port install zstd

#Sys.setenv( ARROW_WITH_ZSTD = "ON")
#Sys.setenv(ARROW_S3 = "ON")
#Sys.setenv(ARROW_GCS = "ON")
