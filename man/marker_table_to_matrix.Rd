% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{marker_table_to_matrix}
\alias{marker_table_to_matrix}
\title{Convert marker table to gene x cell type matrix}
\usage{
marker_table_to_matrix(marker_table, known_genes, weighted = TRUE)
}
\arguments{
\item{marker_table}{Marker data.frame generated by MetaMarkers.}

\item{known_genes}{List of genes to include in the matrix, may contain genes
that are absent from the marker lists, or omit genes that are present.}

\item{weighted}{Boolean. If TRUE, all coefficients in the matrix are renormalized
by marker set size, otherwise all coefficients are 1.}
}
\value{
A gene x cell type matrix containing marker sets.
}
\description{
Convert marker table to gene x cell type matrix
}
