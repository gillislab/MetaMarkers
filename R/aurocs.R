
#' Convert character vector to one-hot encoded matrix, with one column per
#' unique string.
#'
#' @param cell_type Character vector.
#' @param scale_columns Boolean, should columns be normalized to sum to 1?
#'
#' @export
design_matrix = function(cell_type, scale_columns=FALSE) {
  factors = levels(as.factor(cell_type))
  if (length(factors) > 1) {
    result = stats::model.matrix(~as.factor(cell_type)-1)
  } else {
    result = matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) = factors
  if (scale_columns) {
    result = scale(result, center = FALSE, scale = colSums(result))
  }
  return(result)
}

#' Efficient computation of AUROCs (vectorized for predictors and categories).
#'
#' @param predictors Matrix where each column is a predictor and each row is a sample.
#' @param label_matrix One-hot encoded matrix where columns are categories and each row is a sample.
#' The number of rows must be identical to the number of rows in predictors.
#' 1 indicates that the sample on this row belongs to the category on this column.
#' @param compute_tie_correction Boolean. If TRUE, for each AUROC, compute classical
#' tie correction (only useful for p-value computation).
#' @return An AUROC matrix of size #predictors x #categories, containing 
#' all (predictor, category) combinations.
compute_aurocs = function(predictors, label_matrix, compute_tie_correction = FALSE) {
    label_matrix = as.matrix(label_matrix)
    n_positives = colSums(label_matrix)
    n_negatives = nrow(label_matrix) - n_positives
    if (methods::is(predictors, "dgCMatrix")) {
        # we shift all ranks after the matrix multiplication to keep
        # the predictor matrix sparse
        ranks = rank_sparse(predictors)
        sum_of_positive_ranks = as.matrix(Matrix::crossprod(label_matrix, ranks)) +
            outer(n_positives, rank_zero(predictors))
    } else {
        predictors = as.matrix(predictors)
        ranks = matrixStats::colRanks(predictors, ties.method = "average", preserveShape=TRUE)
        sum_of_positive_ranks = crossprod(label_matrix, ranks)
        colnames(sum_of_positive_ranks) = colnames(predictors)
    }
    result = (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives

    if (compute_tie_correction) {
        tie_correction = compute_tie_correction(ranks)
    } else {
        tie_correction = 0
    }
    return(list(aurocs = result, tie_corrections = tie_correction))
}

# these 2 functions only rank non-zeros, implicitly shifting the matrix of ranks
# to keep the matrix sparse according to the formula:
#   true_ranks(M) = rank_zero(M) + rank_sparse(M), where:
#     + rank_zero(M) = (#zeros + 1)/2 = ((nrow(M) - diff(M@p)) + 1) / 2
#     + rank_sparse(M) = rank(nonzero element) + (#zeros - 1)/2
# faster than solution using base::tapply, probably data.table would be faster
rank_sparse = function(M) {
    nnz = diff(M@p)
    ranks = tibble::tibble(x = M@x, j = rep.int(1:ncol(M), nnz)) %>%
        dplyr::group_by(.data$j) %>%
        dplyr::mutate(rank =
            c(matrixStats::colRanks(as.matrix(.data$x), ties.method = "average"))
        )
    R = M
    R@x = ranks$rank + rep.int((nrow(M)-nnz-1)/2, nnz)
    return(R)
}

rank_zero = function(M) {
    return(((nrow(M) - diff(M@p)) + 1) / 2)
}

# For the following functions, see
#
#   https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Normal_approximation_and_tie_correction
#
# The tie correction effectively computes lost variance because of ties (compared to discrete uniform).
# Computing the Wikipedia formula naively is slow, this method is equivalent and fast.
compute_tie_correction = function(ranks) {
    if (methods::is(ranks, "dgCMatrix")) {
        observed_var = colVars_sparse(ranks)
    } else {
        observed_var = matrixStats::colVars(as.matrix(ranks))
    }
    max_var = stats::var(seq_len(nrow(ranks)))
    return((max_var-observed_var) * 12 / nrow(ranks))
}

colVars_sparse = function(M) {
    result = (Matrix::colMeans(M**2) - Matrix::colMeans(M)**2)*nrow(M)/(nrow(M)-1)
    return(result)
}

auroc_p_value = function(aurocs, label_matrix, two_tailed = TRUE, tie_correction = 0, log.p = FALSE) {
    p = colSums(label_matrix)
    n = nrow(label_matrix) - p
  
    # Careful: NAs may arise from tie_correction (predictor with 0 variance)
    if (length(tie_correction) > 1) {
        Z = (aurocs - 0.5) * sqrt(12*n*p)
        Z = t(t(Z) / sqrt(nrow(label_matrix)+1-tie_correction))
    } else {
        Z = (aurocs - 0.5) / sqrt((nrow(label_matrix)+1-tie_correction)/(12*n*p))
    }
  
    result = Z
    if (two_tailed) {
        is_not_na = !is.na(Z)
        result[Z<=0 & is_not_na] = stats::pnorm(Z[Z<=0 & is_not_na], log.p = log.p) * 2
        result[Z>0 & is_not_na] = stats::pnorm(Z[Z>0 & is_not_na], lower.tail=FALSE, log.p = log.p) * 2
    } else {
        result = stats::pnorm(Z, lower.tail=FALSE, log.p = log.p)
    }
    return(result)
}
