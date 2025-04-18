% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_Kandinsky_Class_funcs.R
\name{KanData<-}
\alias{KanData<-}
\title{Assign data to Kandinsky}
\usage{
KanData(object, which = NULL) <- value
}
\arguments{
\item{object}{Seurat or Kandinsky object}

\item{which}{name of Kandinksy slots to extract. If `NULL`, all Kandinsky data will be extracted}

\item{value}{object or variable to assign to the selected Kandinksy slot}
}
\description{
Assign new data to Kandinsky data or Kandinsky slots
}
