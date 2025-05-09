% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_Hot_cold_spot_analysis_.R
\name{hotspot_analysis}
\alias{hotspot_analysis}
\title{Hotspot analysis using Getis-Ord Gi* statistics}
\usage{
hotspot_analysis(
  seurat = NULL,
  feature = NULL,
  layer = "data",
  sim = 999,
  padj.thresh = 0.05,
  lag = 1,
  seed = 347548,
  prune = F
)
}
\arguments{
\item{seurat}{a Seurat object containing Kandinsky data slot(`KanData()`)}

\item{feature}{character string specifying the name of the feature to use to compute Getis-Ord Gi* statistics}

\item{layer}{character string indicating which Seurat layer will be considered to search for the variable of interest}

\item{sim}{number of Monte Carlo simulations to be run for estimating Local Gi* coefficients significance}

\item{padj.thresh}{numeric value indicating the significance threshold to be applied to the adjusted pvalues resulting from Monte Carlo simulations}

\item{lag}{integer value indicating the extent of cell/spot neighbours to be considered to calculate Getis-Ord Gi* statistics. `lag = 1` indicates that only 1st order neighbours will be considered, while `lag=2` indicates that all neighbours of each 1st order neighbour will be also considered for each spot/cell, and so on.
When the average number or 1st order neighbours is limited as in the case of Visium spots (that is, by default, only spots/bins composing the ring immediately surrounding each spot/bin), increasing the lag might help in reducing the variability of the final results. Default value is set to 1.}

\item{seed}{numeric, random seed for reproducibility. Default is set to 347548}

\item{prune}{boolean, whether ignoring any spot/cell annotated as hot or cold spot with no neighbouring hot or cold spots, respectively (i.e., singlets)}
}
\value{
updated Seurat object with two new columns stored in the meta.data table: `feature_name_GI` and `feature_name_clust`, reporting the local Gi* statistics and the association with a hot- or cold-spot for each spot/cell, respectively
}
\description{
Apply Getis-Ord Gi statistic to any numeric variable stored in Seurat object and calculate spatial local hot- and cold-spots for that variable
}
\details{
The function will use the `localG_perm()` function from the `spdep` package to calculate, for a numeric variable of interest (e.g. expression level of a gene), the local Getis-Ord Gi* statistic for each spot/cell analysed.
The `hotspot` function from the same package wil be then use to define local hot-spots or cold-spots on the basis of the spot/cell Gi statistics. These results should be interpreted as follows:
- Spots / cells classified as hot-spots show significantly high values of the variable of interest, as well as their neighbouring spots/cells
- Spots / cells classified as cold-spots show significantly low values of the variable of interest, as well as their neighbouring spots/cells
For more information about Getis-Ord Gi* statistics implemented in `spdep`, see https://r-spatial.github.io/spdep/reference/localmoran.html
}
