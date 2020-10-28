
plot_pareto_markers = function(meta_markers, cell_type_name, min_recurrence = 1,
                               min_auroc = 0.5, min_fc = 0) {
    max_recurrence = max(meta_markers$recurrence)
    my_blue = RColorBrewer::brewer.pal(n = 9, "Blues")[c(3,9)]
    my_palette = c(grDevices::colorRampPalette(my_blue)(max_recurrence-min_recurrence), "black")
    names(my_palette) = min_recurrence:max_recurrence
    
    filtered_markers = meta_markers %>%
        dplyr::select(.data$cell_type, .data$gene, .data$recurrence,
                      .data$auroc, .data$fold_change_detection) %>%
        dplyr::filter(.data$cell_type == cell_type_name & .data$recurrence >= min_recurrence) %>%
        dplyr::filter(.data$auroc > min_auroc & .data$fold_change_detection > min_fc) %>%
        dplyr::mutate(recurrence = factor(.data$recurrence)) %>%
        dplyr::group_by(.data$cell_type) %>%
        dplyr::mutate(is_pareto = is_pareto_front(.data$auroc, .data$fold_change_detection))

    filtered_markers %>%
        ggplot2::ggplot(ggplot2::aes(x = log2(.data$fold_change_detection),
                                     y = .data$auroc,
                                     label = .data$gene)) +
        ggplot2::geom_point(ggplot2::aes(col = .data$recurrence)) +
        ggplot2::geom_line(data = dplyr::filter(filtered_markers, .data$is_pareto)) +
        ggplot2::geom_label(data = dplyr::filter(filtered_markers, .data$is_pareto),
                            ggplot2::aes(col = .data$recurrence),
                            show.legend = FALSE, size = 8) +
        ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = 3, linetype = "dashed") +
        ggplot2::labs(x = "log2(Fold change) of detection rate", y = "AUROC",
                      col = "# Datasets") +
        ggplot2::scale_color_manual(values = my_palette)
}

is_pareto_front = function(x, y, tolerance = 0) {
    ids = seq_along(x)
    result = data.frame(x, y, ids) %>%
        dplyr::arrange(dplyr::desc(.data$x), dplyr::desc(.data$y)) %>%
        dplyr::mutate(max_y = cummax(.data$y), is_pareto = .data$y >= .data$max_y - tolerance) %>%
        dplyr::select(.data$ids, .data$is_pareto) %>%
        dplyr::arrange(.data$ids) %>%
        dplyr::pull(.data$is_pareto)
    return(result)
}

plot_pareto_summary = function(meta_markers) {
    meta_markers %>%
        dplyr::group_by(.data$cell_type) %>%
        dplyr::filter(is_pareto_front(.data$auroc, .data$fold_change_detection)) %>%
        ggplot2::ggplot(ggplot2::aes(x = log2(.data$fold_change_detection),
                                     y = .data$auroc,
                                     col = .data$cell_type,
                                     label = .data$gene)) +
        ggplot2::geom_line() +
        ggplot2::geom_label(show.legend=FALSE, size=6) +
        ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = 3, linetype = "dashed") +
        ggplot2::labs(x = "log2(Fold change) of detection rate", y = "AUROC") +
        ggplot2::theme(legend.title = ggplot2::element_blank())
}

get_pareto_markers = function(meta_markers, cell_type_name) {
    meta_markers %>%
        dplyr::filter(.data$cell_type == cell_type_name) %>%
        dplyr::filter(is_pareto_front(.data$auroc, .data$fold_change_detection)) %>%
        dplyr::arrange(.data$fold_change_detection) %>%
        dplyr::pull(.data$gene)
}

