% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Prepare_seurat_data.R
\name{prepare_seurat_other}
\alias{prepare_seurat_other}
\title{prepare seurat object starting from any spatial transcriptomic/proteomic data in tabular format}
\usage{
prepare_seurat_other(
  data = NULL,
  markers.ids = NULL,
  xcoord = NULL,
  ycoord = NULL,
  stitch = F,
  sample.id = NULL,
  assay.name = "X"
)
}
\arguments{
\item{data}{a dataframe containing single-cell measurements, metadata and spatial coordinates of cell centroids}

\item{markers.ids}{data column names corresponding to marker expression measurements. These will be used to build Seurat count/data matrix}

\item{xcoord}{data column name corresponding to cell centroid x coordinates}

\item{ycoord}{data column name corresponding to cell centroid y coordinates}

\item{stitch}{boolean, whether or not merging cell coordinates coming from independent samples within the dataset. Default is set to FALSE.}

\item{sample.id}{data column name reporting sample identifiers. Must be specified when parameter stitch is set to TRUE}

\item{assay.name}{character string that will be used to name the output Seurat assay}
}
\value{
a Seurat object containing cell spatial expression measurements and metadata
}
\description{
Given a dataframe containing single-cell marker measurements and x/y centroid coordinates, this function create a new Seurat object
with single-cell marker expression and metadata. If multiple independent samples are included in the dataset, there is the possibility to
merge sample coordinates into a unique space with parameter 'stitch =T'
}
