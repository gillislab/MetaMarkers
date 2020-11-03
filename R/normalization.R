
#' Convert expression matrix to CPM (library size normalization).
#'
#' @param M Expression matrix (with genes as rows) in sparse or dense format.
#' @return Matrix in same format as input.
#'
#' @export
convert_to_cpm = function(M) {
    normalization_factor = Matrix::colSums(M) / 1000000
    if (methods::is(M, "dgCMatrix")) {
        M@x = M@x / rep.int(normalization_factor, diff(M@p))
        return(M)
    } else {
        return(scale(M, center = FALSE, scale = normalization_factor))
    }
}