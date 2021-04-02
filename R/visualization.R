
#' Plot expression of markers in cell type as violin plots.
#'
#' @param expression_matrix Expression matrix (with genes as rows) in sparse or
#' dense format.
#' @param genes Vector of gene names.
#' @param cell_type_label Vector of cell type names (must be same length as
#' number of columns in expression_matrix).
#' @return A ggplot2 object.
#'
#' @export
plot_marker_expression = function(expression_matrix, genes, cell_type_label) {
    expr = Matrix::t(expression_matrix[genes,]) %>%
        as.matrix() %>%
        tibble::as_tibble() %>%
        cbind(label = cell_type_label) %>%
        tidyr::pivot_longer(cols = dplyr::all_of(genes), names_to = "gene", values_to = "cpm") %>%
        dplyr::mutate(gene = factor(.data$gene, levels = genes))

    result = expr %>%
        ggplot2::ggplot(ggplot2::aes(x = .data$label, y = .data$cpm+1, fill = .data$label)) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::scale_y_log10() +
        ggplot2::facet_grid(gene ~ .) +
        ggplot2::guides(fill = FALSE) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(y = "Gene expression", x = "Cell type")
    return(result)
}

#' Plot markers scores in UMAP space (one facet per marker set).
#'
#' @param scores A marker score matrix as returned by score_cells or compute_marker_enrichment.
#' @param umap_coordinates A (cell x 2) data.frame containing UMAP coordinates, where the ordering
#'  of cells is the same as in the score matrix
#' @param normalize_scores Boolean. If TRUE, marker scores are renormalized between
#'  0 and 1 for all marker sets to obtain similar color ranges across facets.
#' @param point_size Numeric value indicating the dot size.
#' @param rasterize_umap Boolean. If TRUE, the UMAP part of the figure is rasterized.
#'  This is particularly useful for large datasets.
#' @return A ggplot2 object.
#'
#' @export
plot_marker_scores = function(scores, umap_coordinates, normalize_scores=FALSE, point_size=1, rasterize_umap=FALSE) {
    colnames(umap_coordinates) = c("umap_1", "umap_2")
    to_plot = tibble::as_tibble(t(as.matrix(scores)), .name_repair = "minimal") %>%
        dplyr::bind_cols(umap_coordinates) %>%
        tidyr::pivot_longer(c(-.data$umap_1, -.data$umap_2),
                            names_to = "cell_type", values_to = "score")
    
    if (normalize_scores) {
        to_plot = to_plot %>%
            dplyr::group_by(.data$cell_type) %>%
            dplyr::mutate(score = (.data$score-min(.data$score))/(max(.data$score)-min(.data$score)))
    }
        
    result = ggplot2::ggplot(to_plot, ggplot2::aes(x = .data$umap_1, y = .data$umap_2, col = .data$score))
    if (rasterize_umap) {
        result = result + ggrastr::geom_point_rast(size=point_size)
    } else {
        result = result + ggplot2::geom_point(size=point_size)
    }
    result = result +
        ggplot2::facet_wrap(~ .data$cell_type) +
        ggplot2::scale_color_gradient(low = "gray100", high = "darkgreen") +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::labs(x = "UMAP1", y = "UMAP2", col = "Marker score")
    return(result)
}

#' Plot cell type assignments in UMAP space.
#'
#' @param assignments A cell type assignment table as returned by assign_cells.
#' @param umap_coordinates A (cell x 2) data.frame containing UMAP coordinates, where the ordering
#'  of cells is the same as in the score matrix
#' @param enrichment_threshold Numeric value defining a threshold below which a
#'  cell should be considered unassigned (QC metric). 1 by default.
#' @param point_size Numeric value indicating the dot size.
#' @param rasterize_umap Boolean. If TRUE, the UMAP part of the figure is rasterized.
#'  This is particularly useful for large datasets.
#' @return A ggplot2 object.
#'
#' @export
plot_assignments = function(assignments, umap_coordinates, enrichment_threshold = 1, point_size=1, rasterize_umap=FALSE) {
    colnames(umap_coordinates) = c("umap_1", "umap_2")
    to_plot = assignments %>%
        dplyr::bind_cols(umap_coordinates) %>%
        dplyr::mutate(predicted = as.character(.data$predicted)) %>%
        dplyr::mutate(predicted = ifelse(.data$enrichment < enrichment_threshold, NA, .data$predicted)) %>%
        dplyr::mutate(predicted = ifelse(.data$predicted == "unassigned", NA, .data$predicted)) %>%
        dplyr::mutate(predicted = factor(.data$predicted))
    
    result = ggplot2::ggplot(to_plot, ggplot2::aes(x = .data$umap_1, y = .data$umap_2, col = .data$predicted))
    if (rasterize_umap) {
        result = result + ggrastr::geom_point_rast(size=point_size)
    } else {
        result = result + ggplot2::geom_point(size=point_size)
    }
    result = result +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::scale_color_hue(na.value = "gray80") +
        ggplot2::labs(x = "UMAP1", y = "UMAP2", col = "Cell type")
    return(result)    
}

#' Compute UMAP coordinates for a target dataset.
#'
#' @param expr Expression matrix (with genes as rows) of the target dataset in
#'  sparse or dense format. We recommend providing log-normalized counts for
#'  carefully selected highly variable genes.
#' @return A tibble with one row per cell and 3 columns containing the sample
#'  name, UMAP coordinate 1 and UMAP coordinate 2.
#'
#' @export
compute_umap = function(expr) {
    result = uwot::umap(t(as.matrix(expr))) %>%
        tibble::as_tibble(.name_repair = "minimal") %>%
        rlang::set_names(c("umap_1", "umap_2")) %>%
        tibble::add_column(sample = colnames(expr), .before = "umap_1")
    return(result)
}
