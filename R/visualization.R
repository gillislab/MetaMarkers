
#' @export
plot_marker_scores = function(scores, umap_coordinates, normalize_scores=FALSE) {
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
        
    result = to_plot %>%
        ggplot2::ggplot(ggplot2::aes(x = .data$umap_1, y = .data$umap_2, col = .data$score)) +
        ggplot2::geom_point(size=0.1) +
        ggplot2::facet_wrap(~ .data$cell_type) +
        ggplot2::scale_color_gradient(low = "gray100", high = "darkgreen") +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::labs(x = "UMAP1", y = "UMAP2", col = "Marker score")
    return(result)
}

plot_assignments = function(assignments, umap_coordinates, enrichment_threshold = 1) {
    colnames(umap_coordinates) = c("umap_1", "umap_2")
    to_plot = assignments %>%
        dplyr::bind_cols(umap_coordinates) %>%
        dplyr::mutate(predicted = ifelse(.data$enrichment < enrichment_threshold, NA, .data$predicted)) %>%
        dplyr::mutate(predicted = ifelse(.data$predicted == "unassigned", NA, .data$predicted))
    
    result = to_plot %>%
        ggplot2::ggplot(ggplot2::aes(x = .data$umap_1, y = .data$umap_2, col = .data$predicted)) +
        ggplot2::geom_point(size=0.1) +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::scale_color_hue(na.value = "gray80") +
        ggplot2::labs(x = "UMAP1", y = "UMAP2", col = "Cell type")
    return(result)    
}

compute_umap = function(expr) {
    result = uwot::umap(t(as.matrix(expr))) %>%
        tibble::as_tibble(.name_repair = "minimal") %>%
        rlang::set_names(c("umap_1", "umap_2")) %>%
        tibble::add_column(sample = colnames(expr), .before = "umap_1")
    return(result)
}