% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{summarize_fold_change}
\alias{summarize_fold_change}
\title{Compute signal-to-noise ratio (fold change) of cell types based on marker sets.}
\usage{
summarize_fold_change(scores, true_labels)
}
\arguments{
\item{scores}{A marker score matrix as returned by score_cells or compute_marker_enrichment.}

\item{true_labels}{Ground truth labels for the target cells.}
}
\value{
A marker set x true label fold change matrix. For a given combination, the
marker scores are used as the statistic and the true labels are used to
define the two classes to compare (do cells belong to the given cell type?).
}
\description{
Compute signal-to-noise ratio (fold change) of cell types based on marker sets.
}
