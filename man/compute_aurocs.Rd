% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aurocs.R
\name{compute_aurocs}
\alias{compute_aurocs}
\title{Efficient computation of AUROCs (vectorized for predictors and categories).}
\usage{
compute_aurocs(predictors, label_matrix, compute_tie_correction = FALSE)
}
\arguments{
\item{predictors}{Matrix where each column is a predictor and each row is a sample.}

\item{label_matrix}{One-hot encoded matrix where columns are categories and each row is a sample.
The number of rows must be identical to the number of rows in predictors.
1 indicates that the sample on this row belongs to the category on this column.}

\item{compute_tie_correction}{Boolean. If TRUE, for each AUROC, compute classical
tie correction (only useful for p-value computation).}
}
\value{
An AUROC matrix of size #predictors x #categories, containing
all (predictor, category) combinations.
}
\description{
Efficient computation of AUROCs (vectorized for predictors and categories).
}
