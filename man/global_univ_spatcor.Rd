% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_Utils.R
\name{global_univ_spatcor}
\alias{global_univ_spatcor}
\title{Compute global Moran's I spatial autocorrelation statistic}
\usage{
global_univ_spatcor(
  seurat,
  var = NULL,
  sim = 49,
  lag = 1,
  alt = "greater",
  seed = 347548
)
}
\arguments{
\item{seurat}{a Seurat object containing Kandinsky data (`KanData()`)}

\item{var}{character string or vector specifying the variable(s) to use to compute Moran's I statistics}

\item{sim}{number of Monte Carlo simulations to be run for estimating Moran coefficients significance}

\item{lag}{integer value indicating the extent of cell/spot neighbours to be considered to calculate Moran statistics. `lag = 1` indicates that only 1st order neighbours will be considered, while `lag=2` indicates that all neighbours of each 1st order neighbour will be also considered for each spot/cell, and so on.}

\item{alt}{a character string specifying the alternative hypothesis, Must be one of the following: "two.sided", "greater" (default), or "less".}

\item{seed}{numeric, random seed for reproducibility. Default is set to 347548}
}
\value{
a data.frame reporting Moran's I spatial autocorrelation statistics with the associated p value for each selected variable.
}
\description{
Compute global Moran's I spatial autocorrelation statistic
}
