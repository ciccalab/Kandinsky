% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Prepare_seurat_data.R
\name{prepare_slideseq_seurat}
\alias{prepare_slideseq_seurat}
\title{Prepare Slide-Seq Seurat data}
\usage{
prepare_slideseq_seurat(path = NULL, dataset.id = NULL, pattern = NULL)
}
\arguments{
\item{path}{character string specifying the path to the Slide-Seq input files directory}

\item{dataset.id}{character string that will be used to name the output Seurat object identity}

\item{pattern}{character string that will be used as a key to identify Slide-Seq input file names in case when input files from multiple Slide-Seq samples/slides are stored in the same directory}
}
\value{
a Seurat object containing Slide-Seq count matrix, metadata, and spatial coordinates to be used for downstream analysis with Kandinsky functions
}
\description{
Initial formatting of Slide-Seq data to work with Seurat and Kandinsky
}
\details{
This function will use Slide-Seq input data to build a new Seurat object. All Slide-Seq input files must be stored in the same folder that will be specified to call the function.
Kandinsky expects to find Slide-Seq files in the same format proposed in the Slide-Seq V2 original study https://singlecell.broadinstitute.org/single_cell/study/SCP815/sensitive-spatial-genome-wide-expression-profiling-at-cellular-resolution#study-summary
}
\seealso{
Other prepare_data: 
\code{\link{prepare_cosmx_seurat}()},
\code{\link{prepare_g4x_seurat}()},
\code{\link{prepare_merscope_seurat}()},
\code{\link{prepare_proseg_seurat}()},
\code{\link{prepare_visiumHD_seurat}()},
\code{\link{prepare_visium_seurat}()},
\code{\link{prepare_xenium_seurat}()}
}
\concept{prepare_data}
