
#' @export
plot_marker_scores = function(scores, umap_coordinates, normalize_scores=FALSE) {
    colnames(umap_coordinates) = c("umap_1", "umap_2")
    to_plot = as_tibble(t(as.matrix(scores))) %>%
        bind_cols(umap_coordinates) %>%
        pivot_longer(c(-umap_1, -umap_2), names_to = "cell_type", values_to = "score")
    
    if (normalize_scores) {
        to_plot = to_plot %>%
            group_by(cell_type) %>%
            mutate(score = (score-min(score))/(max(score)-min(score)))
    }
        
    result = to_plot %>%
        ggplot(aes(x = umap_1, y = umap_2, col = score)) +
        geom_point(size=0.1) +
        facet_wrap(~ cell_type) +
        scale_color_gradient(low = "gray100", high = "darkgreen") +
        theme_bw(base_size = 20) +
        labs(x = "UMAP1", y = "UMAP2", col = "Marker score")
    return(result)
}

plot_assignments = function(assignments, umap_coordinates, enrichment_threshold = 1) {
    colnames(umap_coordinates) = c("umap_1", "umap_2")
    to_plot = assignments %>%
        bind_cols(umap_coordinates) %>%
        mutate(predicted = ifelse(enrichment < enrichment_threshold, NA, predicted)) %>%
        mutate(predicted = ifelse(predicted == "unassigned", NA, predicted))
    
    result = to_plot %>%
        ggplot(aes(x = umap_1, y = umap_2, col = predicted)) +
        geom_point(size=0.1) +
        theme_bw(base_size = 20) +
        scale_color_hue(na.value = "gray80") +
        labs(x = "UMAP1", y = "UMAP2", col = "Cell type")
    return(result)    
}

compute_umap = function(expr) {
    result = uwot::umap(t(as.matrix(expr))) %>%
        as_tibble(.name_repair = "minimal") %>%
        set_names(c("umap_1", "umap_2")) %>%
        add_column(sample = colnames(expr), .before = "umap_1")
    return(result)
}