% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07_Visualization.R
\name{PropPlot}
\alias{PropPlot}
\title{Proportion barplot}
\usage{
PropPlot(seurat = NULL, var.1 = NULL, var.2 = NULL, cols = NULL)
}
\arguments{
\item{seurat}{a Seurat object containing Kandinsky data}

\item{var.1}{character string, name of the first variable}

\item{var.2}{character string, name of the second variable}

\item{cols}{color palette to use for plotting}
}
\value{
stacked barplot created with ggplot2
}
\description{
Create a stacked barplot reporting the proportion covered by var.2 classes within each var.1 class
}
