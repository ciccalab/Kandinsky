% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_Kandinsky_Class_funcs.R
\name{saveKanData}
\alias{saveKanData}
\alias{saveKanData,Kandinsky-method}
\alias{saveKanData,Seurat-method}
\title{Save Kandinsky data}
\usage{
saveKanData(object, file = NULL, with.seurat = T)

\S4method{saveKanData}{Kandinsky}(object, file = NULL, with.seurat = T)

\S4method{saveKanData}{Seurat}(object, file = NULL, with.seurat = T)
}
\arguments{
\item{object}{a Seurat object containing Kandinsky data, or a Kandinsky object}

\item{file}{character string specifying the path to the file where to save Kandinsky data. It must be in .rds format}

\item{with.seurat}{boolean, whether saving the whole Seurat object together with Kandinsky data (`TRUE`) or Kandinsky data alone (`FALSE`). Default value is set to `TRUE`. Ignored when `object` corresponds to a Kandinsky object}
}
\description{
saveRDS adaptation for Kandinsky data. it can save Kandinsky data alone or as part of a Seurat object
}
