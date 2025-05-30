% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_Kandinsky_Class_funcs.R
\name{kandinsky_init}
\alias{kandinsky_init}
\title{Initialize Kandinsky data}
\usage{
kandinsky_init(
  seurat = NULL,
  tech = "visium",
  img = NULL,
  img_maxdim = 2000,
  res = c("low", "high", "full"),
  binsize = 16,
  tx_path = NULL,
  fov_path = NULL,
  poly_path = NULL,
  nb.method = c("Q", "C", "D", "K", "M"),
  k = 20,
  d.max = 40,
  soi = F,
  layers = 1,
  ids_other = NULL,
  xcoord_other = NULL,
  ycoord_other = NULL
)
}
\arguments{
\item{seurat}{a Seurat object}

\item{tech}{character string specifying the platform used to generate the sequencing data. Must be one of the following: `visium`, `visium_hd`, `cosmx`, `xenium`, `merscope`, `other`.}

\item{img}{character string indicating the path of the H&E image to be used for Visium/VisiumHD. If set to `NULL`, no image is loaded and stored in the final Kandinsky data}

\item{img_maxdim}{numeric value indicating the ideal maximum pixel dimension accepted for Visium H&E image. If image size is bigger than this value, a scale factor will be applied to lower image pixel resolution/size up to the defined maximum value}

\item{res}{character string specifying which version of H&E Visium image will be loaded into the Kandinsky data. Must be one of the following: `low`, `high`, `full`, where `full` refers to the original H&E full-resolution tiff image}

\item{binsize}{numeric value indicating which bin resolution will be considered to create the final sf object when working with Visium-HD data}

\item{tx_path}{character string indicating the path for the transcript coordinates file from CoxMx/Xenium/Merscope platforms}

\item{fov_path}{character string indicating the path for the fov coordinates file from CoxMx platform}

\item{poly_path}{character string indicating the path for the polygon coordinates file from CoxMx/Xenium platform}

\item{nb.method}{character string specifying the method to be used to create a `nb` neighbour object.
Must be one of the following:
'Q': queen contiguity method,check for contact (not overlap) between any edge or side od two polygons (refers to the queen movement rule in chess). Currently only applicable for Visium/Visium-HD data
'C': centroid-based method, use maximum centroid distance threshold to identify spot/cell neighbours
'D': Delaunay triangulation method
'K': KNN method, define k closest neighbours to each spot/cell
'M': membrane-based method, check for the occurrence of a physical contact/intersection within a distance threshold between cell boundaries. Not applicable in the case of Visium spots.}

\item{k}{numeric, number of nearest neighbours to be set when `nb.method = K`}

\item{d.max}{numeric, maximum centroid distance threshold to be set when `nb.method = C | M`}

\item{soi}{boolean, whether or not filter Delaunay network to keep sphere of influence (SOI) graph. Default is FALSE}

\item{layers}{numeric, number of concentric contiguous layers to include in spot neighbourhood. Only Applied when `nb.method = Q`. Default is 1}

\item{ids_other}{character string specifying variable name to be used as cell identifiers when argument tech is set to "other"}

\item{xcoord_other}{character string specifying variable name to be used as x coordinates when argument tech is set to "other"}

\item{ycoord_other}{character string specifying variable name to be used as y coordinates when argument tech is set to "other"}
}
\description{
initialize Kandinsky data starting from an already existing Seurat object.
}
\details{
The new data generated with this function will be added to the tool slot of Seurat object
as a new 'kandata' list object. `kandinsky_init` will create a fixed set of input data required for downstream
spatial analysis, plus additional objects depending on the specific sequencing platform generating the data stored in Seurat.
If Seurat object has been created through one of the `prepare_*_seurat` functions, some of the input files required by `kandinsky_init`
will be automatically detected even if the corresponding parameters are set to `NULL`.
}
