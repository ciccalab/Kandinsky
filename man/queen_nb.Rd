% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_Neighbour_methods.R
\name{queen_nb}
\alias{queen_nb}
\title{queen contiguity neighbour identification}
\usage{
queen_nb(poly, snap = NULL, layers = 1)
}
\arguments{
\item{poly}{a sf data frame with polygon geometry}

\item{snap}{numeric, extra distance accepted between polygon borders for contiguity relation}

\item{layers}{numeric, number of concentric contiguous layers to include in spot neighbourhood. Only Applied when `nb.method = Q`. Default is 1}
}
\value{
weighted neighbour matrix in a listw object from spdep package
}
\description{
Define neighbour relationships according to Queen criterion for polygon contiguity.
}
\details{
Based on poly2nb() function from spdep package.
Designed to work specifically with Visium / Visium HD data only.
Each Visium spot/bin will have a neighbour link with only the first nearest ring of surrounding spots (Visium, n=6) or bins (Visium HD, n=8)
}
\seealso{
Other nb_funcs: 
\code{\link{centroid_nb}()},
\code{\link{knn_nb}()},
\code{\link{membrane_nb}()},
\code{\link{nb_expand}()},
\code{\link{nb_update}()},
\code{\link{tri_nb}()}
}
\concept{nb_funcs}
