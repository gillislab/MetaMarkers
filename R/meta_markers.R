
#' Extract meta-analytic markers from multiple datasets.
#'
#' @param marker_lists Named list, where elements are marker statistics obtained
#' with `compute_markers` and names are names of the dataset from which markers
#' have been extracted.
#' @param order_by Secondary statistic by which meta-markers should be ranked.
#' @param fdr_threshold FDR threshold for a gene to be considered Differentially
#' Expressed (DE).
#' @param fc_threshold Fold change threshold for a gene to be DE.
#' @param detection_threshold Detection rate threshold for a gene to be DE.
#' @param detailed_stats Boolean. By default, only output a list of best markers,
#' alternatively output additional statistics, such as average AUROC, FC, etc.
#' @param common_genes_only Boolean. Keep only genes that are common to all
#' datasets?
#' @param check_duplicates Boolean. Check and remove duplicated gene names? In
#' theory, this step was already performed in `compute_markers` and does not
#' need to be performed again (time consuming).
#'
#' @return A tibble containing ranked meta-markers and (if desired) average DE
#' statistics for each cell type. Within each cell type, meta-markers are
#' ranked by recurrence (# datasets in which gene is DE), then by a secondary
#' statistics, as specified in "order_by" (AUROC by default).
#'
#' @export
make_meta_markers = function(marker_lists, order_by = "auroc",
                             fdr_threshold = 0.05, fc_threshold = 4,
                             detection_threshold = 0,
                             detailed_stats=FALSE,
                             common_genes_only=TRUE,
                             check_duplicates=FALSE) {
    if (check_duplicates) {
        marker_lists = lapply(marker_lists, remove_duplicated_genes)
    }
    
    all_markers = dplyr::bind_rows(marker_lists, .id = "dataset")
    if (common_genes_only) {
        all_markers = dplyr::filter(
            all_markers, .data$gene %in% find_common_genes(marker_lists)
        )
    }
    
    all_markers = all_markers %>%
        dplyr::mutate(is_de = .data$log_fdr <= log(fdr_threshold) &
                      .data$fold_change >= fc_threshold,
                      .data$detection_rate >= detection_threshold)
    
    result = all_markers %>%
        dplyr::group_by(.data$group, .data$cell_type, .data$gene) %>%
        dplyr::summarize(recurrence = sum(.data$is_de),
                         auroc = mean(.data$auroc),
                         fold_change = gmean(.data$fold_change),
                         fold_change_detection = gmean(.data$fold_change_detection),
                         expression = gmean(.data$average_expression),
                         precision = mean(.data$precision),
                         recall = mean(.data$detection_rate),
                         population_size = mean(.data$population_size),
                         n_datasets = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$group, .data$cell_type) %>%
        dplyr::arrange_at(c("recurrence", order_by), dplyr::desc, .by_group=TRUE) %>%
        dplyr::mutate(rank = 1:dplyr::n()) %>%
        dplyr::select(.data$group, .data$cell_type, .data$rank, dplyr::everything()) %>%
        dplyr::ungroup()

    if (detailed_stats) {
        recurrence_table = all_markers %>%
            dplyr::select(.data$group, .data$cell_type, .data$gene,
                          .data$dataset, .data$is_de) %>%
            tidyr::pivot_wider(id_cols = dplyr::everything(),
                               names_from = .data$dataset,
                               values_from = .data$is_de)
        result = dplyr::inner_join(result, recurrence_table,
                                   by = c("group", "cell_type", "gene"))
    } else {
        result = dplyr::select(result, .data$group, .data$cell_type,
                               .data$gene, .data$rank)
    }
    return(result)
}

# Identify common genes across multiple marker tables
find_common_genes = function(marker_lists) {
    genes = lapply(marker_lists, function(marker_stats) {
        unique(marker_stats$gene)
    })
    return(Reduce(intersect, genes))
}

# Geometric mean
gmean = function(x) {
    exp(mean(log(x)))
}

#' Export meta-analytic markers to file (including header with time stamp).
#'
#' @param meta_markers Meta-marker table obtained with `make_meta_marker`.
#' @param filename File where to write meta-marker table.
#' @param dataset_names Names of datasets used to compute meta-markers.
#' @param order_by Secondary statistic by which meta-markers have been ranked.
#' @param fdr_threshold FDR threshold used to identify Differentially
#' Expressed (DE) genes.
#' @param fc_threshold Fold change threshold used to identify DE genes.
#' @param detection_threshold Detection rate threshold used to identify DE genes.
#' @param gzip Boolean. Should output file be zipped to spare memory?
#'
#' @export
export_meta_markers = function(meta_markers, filename, dataset_names,
                               order_by = "auroc", fdr_threshold = 0.05,
                               fc_threshold = 4, detection_threshold = 0,
                               gzip=TRUE) {
    write_meta_header(filename, order_by, fdr_threshold, fc_threshold,
                      detection_threshold, dataset_names)
    meta_markers %>%
#        dplyr::mutate_if(is.double, round, 2) %>%
        data.table::fwrite(filename, append=TRUE, row.names=FALSE, col.names=TRUE)
    if (gzip) {
        R.utils::gzip(filename, overwrite=TRUE)
    }
}

# Header for meta-marker file
# (including time stamp, dataset names and recurrence criteria)
write_meta_header = function(filename, order_by = "auroc", fdr_threshold = 0.05,
                             fc_threshold = 4, detection_threshold = 0,
                             dataset_names) {
    lines = paste0(
        "# Meta marker list generated on ", Sys.Date(), ". ",
        "Ordered by recurrence among DE genes (FC >= ", fc_threshold,
        ", FDR <= ", fdr_threshold, ", detection_rate >= ",
        detection_threshold, "), then ", order_by, ". ",
        "Based on: ", paste(dataset_names, collapse = ", "), ". ",
        "All stats are averaged across datasets, last columns indicate in which dataset the gene is DE."
    )
    writeLines(lines, filename)
}

#' Export meta-analytic markers to file (one file per cell type).
#'
#' @param meta_markers Meta-marker table obtained with `make_meta_marker`.
#' @param output_prefix File prefix, typically a path or a path + prefix.
#' Each cell type will be exported to a standardized file name, consisting of
#' output_prefix, a standardized version of the cell type name 
#' using `make.names`) and a .csv extension.
#' @param dataset_names Names of datasets used to compute meta-markers.
#' @param order_by Secondary statistic by which meta-markers have been ranked.
#' @param fdr_threshold FDR threshold used to identify Differentially
#' Expressed (DE) genes.
#' @param fc_threshold Fold change threshold used to identify DE genes.
#' @param detection_threshold Detection rate threshold used to identify DE genes.
#' @param gzip Boolean. Should output file be zipped to spare memory?
#'
#' @export
export_meta_markers_by_cell_type = function(
        meta_markers, output_prefix, dataset_names, order_by = "auroc",
        fdr_threshold = 0.05, fc_threshold = 4, detection_threshold = 0,
        gzip=TRUE) {
    for (ct in unique(meta_markers$cell_type)) {
        filename = paste0(output_prefix, make.names(ct), ".csv")
        meta_markers %>%
            dplyr::filter(.data$cell_type == ct) %>%
            export_meta_markers(filename, dataset_names, order_by, fdr_threshold,
                                fc_threshold, detection_threshold, gzip)
    }
}

#' Import meta-marker statistics.
#'
#' @param filename File name.
#' @param header_size Size of header (# rows) containing meta-information such
#' as time stamp and list of dataset names (for compatibility purposes only).
#'
#' @export
read_meta_markers = function(filename, header_size=1) {
    result = data.table::fread(filename, header=TRUE, skip=header_size)
    return(tibble::as_tibble(result))
}
