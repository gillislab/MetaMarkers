
#' Compute differential expression statistics for given dataset and cell types,
#' stratified by groups.
#'
#' @param expression Expression matrix (may be sparse).
#' @param cell_type_labels Character vector providing
#' cell type names for each sample in the expression matrix.
#' @param group_labels Character vector providing hierarchical grouping for
#' cell types (one group name for each sample in the expression matrix).
#' @param two_tailed Boolean. If FALSE, only upregulated genes are considered
#' significant (ROC test).
#' @param tie_correction Boolean. For the ROC test, should tie correction be
#' applied? Note that skipping tie correction is slightly conservative.
#' @param genes_are_rows Boolean. In the expression matrix, were genes provided
#' as rows (as in SingleCellExperiment or Seurat objects)?
#'
#' @return A tibble containing basic differential expression statistics for
#' all cell types and genes. All statistics are 1-vs-rest within groups. Note
#' genes with duplicate names will be removed.
#'
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr %>%
compute_markers = function(expression, cell_type_labels,
                           group_labels = rep("all", length(cell_type_labels)),
                           two_tailed = FALSE, tie_correction = FALSE,
                           genes_are_rows=TRUE) {
    result = lapply(purrr::set_names(unique(group_labels)), function(group) {
        keep_cell = group_labels == group
        if (genes_are_rows) {
            expr = expression[, keep_cell]
        } else {
            expr = expression[keep_cell,]
        }
        label_matrix = design_matrix(cell_type_labels[keep_cell])
        if (ncol(label_matrix) > 1) {
            compute_markers_(expr, label_matrix, two_tailed,
                             tie_correction, genes_are_rows)
        } else {
            warning(paste0("Skipping group ", group,
                    ": only 1 cluster (", colnames(label_matrix) ,")."))
            return(NULL)
        }
    })
    result = result[!sapply(result, is.null)]
    result = dplyr::bind_rows(result, .id = "group")
    return(result)
}

# Actual implementation of marker computation (for a given group)
compute_markers_ = function(expression, cell_type_labels, two_tailed = FALSE,
                            tie_correction = FALSE, genes_are_rows=TRUE) {
    if (is.vector(cell_type_labels)) {
        cell_type_labels = design_matrix(as.character(cell_type_labels))
    }
    if (genes_are_rows) {
        expression = Matrix::t(expression)
    }
    
    aurocs = compute_aurocs(expression, cell_type_labels,
                            compute_tie_correction = tie_correction)
    p_values = auroc_p_value(
        aurocs$aurocs, cell_type_labels, two_tailed,
        tie_correction = aurocs$tie_correction, log.p = TRUE
    )
    log_fdr = matrix(my_fdr(p_values, log.p = TRUE), nrow = nrow(p_values),
                     dimnames = dimnames(p_values))
    
    population_size = colSums(cell_type_labels)
    population_fraction = population_size / sum(population_size)
    
    average = average_expression(expression, cell_type_labels)
    fold_change = (average$positives+1) / (average$negatives+1)
    uncentered_var = average_expression(expression**2, cell_type_labels)
    standard_error = sqrt((uncentered_var$positives - average$positives**2) / population_size)
          
    binary_expression = expression > 0
    m_binary = average_expression(binary_expression, cell_type_labels)
    fc_binary = (m_binary$positives + 0.001) / (m_binary$negatives + 0.001)
    n_expressing_cells = as.matrix(Matrix::crossprod(cell_type_labels, binary_expression))
    binary_precision = t(t(n_expressing_cells) / c(my_col_sums(binary_expression)))
    binary_recall = n_expressing_cells / population_size

    result = tidy_stat(fold_change, "fold_change")
    # WARNING: the indexing technique is efficient, but only works if
    # gene names are unique!!!
    full_indices = as.matrix(result[, c("cell_type", "gene")])
    result = result %>% dplyr::mutate(
        auroc = aurocs$aurocs[full_indices],
        log_fdr = log_fdr[full_indices],
        population_size = population_size[result$cell_type],
        population_fraction = population_fraction[result$cell_type],
        average_expression = average$positives[full_indices],
        se_expression = standard_error[full_indices],
        detection_rate = m_binary$positives[full_indices],
        fold_change_detection = fc_binary[full_indices],
        precision = binary_precision[full_indices],
        recall = binary_recall[full_indices]
    ) %>%
        dplyr::arrange(.data$cell_type, dplyr::desc(.data$auroc))
    
    # remove duplicate genes (see WARNING above)
    # NOTE: more efficient to do it now than in the original matrix
    genes = colnames(expression)
    duplicated_genes = unique(genes[duplicated(genes)])
    result = dplyr::filter(result, !(.data$gene %in% duplicated_genes))

    return(result)
}

# Compute average expression for each cell type + average background expression
average_expression = function(expression, design_matrix) {
    scaled_positives = scale(design_matrix, center=FALSE, scale=colSums(design_matrix))
    scaled_negatives = scale(1-design_matrix, center=FALSE, scale=colSums(1-design_matrix))
    return(list(
        positives = as.matrix(Matrix::crossprod(scaled_positives, expression)),
        negatives = as.matrix(Matrix::crossprod(scaled_negatives, expression))
    ))
}

# More efficient version of colSums (at least for sparse matrices)
my_col_sums = function(M) {
    as.matrix(rep(1, nrow(M)) %*% M)
}

# Adpated from p.adjust (with the possibility to provide log p-values)
my_fdr = function(p_values, log.p = FALSE) {
    i = length(p_values):1L
    o = order(p_values, decreasing = TRUE)
    n = length(p_values)
    result = rep(0, n)
    if (log.p) {
        result[o] = pmin(0, cummin(log(n) - log(i) + p_values[o]))
    } else {
        result[o] = pmin(1, cummin(n/i * p_values[o]))
    }
    return(result)
}

# Convert stat matrix to tidy format
tidy_stat = function(stat, stat_name, row_name="cell_type", col_name="gene") {
    stat %>%
        tibble::as_tibble(rownames=row_name, .name_repair = "minimal") %>%
        tidyr::pivot_longer(cols=-{{row_name}}, names_to=col_name, values_to=stat_name)
}

# Compute markers blockwise (dense blocks of the original sparse matrix)
compute_markers_blockwise = function(expression, cell_type_labels,
                                     group_labels = rep("all", length(cell_type_labels)),
                                     two_tailed = FALSE, tie_correction = FALSE,
                                     n_genes = NULL, genes_are_rows = TRUE) {
    if (genes_are_rows) {
        expression = Matrix::t(expression)
    }
    if (is.null(n_genes)) {
        n_genes = floor(2e9 / nrow(expression))
    }
    
    blocks = unique(c(seq(1, ncol(expression), n_genes), ncol(expression)))
    result = lapply(seq_len(length(blocks) - 1), function(i) {
        compute_markers(
            as.matrix(expression[, blocks[i]:blocks[i+1]]),
            cell_type_labels, group_labels, two_tailed,
            tie_correction, genes_are_rows=FALSE
        )
    })
    result = dplyr::bind_rows(result) %>%
        dplyr::arrange(.data$cell_type, .data$gene)

    # remove duplicate genes (see WARNING in compute_markers_)
    # NOTE: more efficient to do it now than in the original matrix
    genes = colnames(expression)
    duplicated_genes = unique(genes[duplicated(genes)])
    result = dplyr::filter(result, !(.data$gene %in% duplicated_genes))
    return(result)
}

#' Remove duplicate gene from a marker table (individual dataset)
#'
#' @param marker_stats Tibble obtained by calling `compute_marker_stats`.
#' 
#' @return Tibble with marker stats, all duplicated genes being removed.
#'
#' @export
remove_duplicated_genes = function(marker_stats) {
    duplicated_genes = marker_stats %>%
        dplyr::group_by(.data$group, .data$cell_type, .data$gene) %>%
        dplyr::tally() %>%
        dplyr::filter(.data$n > 1) %>%
        dplyr::pull(.data$gene) %>%
        unique()
    return(dplyr::filter(marker_stats, !(.data$gene %in% duplicated_genes)))
}


#' Export marker statistics to file.
#'
#' This function will export marker statistics in a simple format, including a short
#' summary header with cell types analyzed and date of export.
#'
#' @param marker_stats Tibble obtained by calling `compute_marker_stats`.
#' @param filename File name.
#' @param gzip Boolean. Should output file be zipped to spare memory?
#'
#' @export
export_markers = function(marker_stats, filename, gzip=TRUE) {
    write_header(unique(marker_stats$cell_type), filename)
    data.table::fwrite(marker_stats, filename, append=TRUE,
                       row.names=FALSE, col.names=TRUE)
    if (gzip) {
        R.utils::gzip(filename, overwrite=TRUE)
    }
}

# Generic header (time stamp, cell types)
write_header = function(cell_types, filename) {
    lines = paste0(
        "# Differential expression statistics generated on ", Sys.Date(), ". ",
        "Cell types: ", paste(cell_types, collapse = ", "), "."
    )
    writeLines(lines, filename)
}

#' Export markers statistics to file, one file per group.
#'
#' This function will export marker statistics in a simple format, including a short
#' summary header with cell types analyzed and date of export.
#'
#' @param de_stats Tibble obtained by calling `compute_markers`.
#' @param output_prefix File prefix, typically a path or a path + prefix.
#' Each group will be exported to a standardized file name, consisting of
#' output_prefix, a standardized version of the group name 
#' using `make.names`) and a .csv extension.
#' @param gzip Boolean. Should output file be zipped to spare memory?
#'
#' @export
export_markers_by_group = function(de_stats, output_prefix, gzip=TRUE) {
    for (group in unique(de_stats$group)) {
        filename = paste0(output_prefix, make.names(group), ".csv")
        export_markers(dplyr::filter(de_stats, group == group),
                       filename, gzip)
    }
}

#' Export marker statistics to file, one file per cell type.
#'
#' This function will export marker statistics in a simple format, including a short
#' summary header with cell types analyzed and date of export.
#'
#' @param de_stats Tibble obtained by calling `compute_markers`.
#' @param output_prefix File prefix, typically a path or a path + prefix.
#' Each cell type will be exported to a standardized file name, consisting of
#' output_prefix, a standardized version of the cell type name 
#' using `make.names`) and a .csv extension.
#' @param gzip Boolean. Should output file be zipped to spare memory?
#'
#' @export
export_markers_by_cell_type = function(de_stats, output_prefix, gzip=TRUE) {
    for (ct in unique(de_stats$cell_type)) {
        filename = paste0(output_prefix, make.names(ct), ".csv")
        export_markers_for_cell_type(de_stats, ct, filename, gzip)
    }
}

#' Export marker statistics of given cell type to file.
#'
#' This function will export marker statistics in a simple format, including a short
#' summary header with cell types analyzed and date of export.
#'
#' @param de_stats Tibble obtained by calling `compute_markers`.
#' @param cell_type Name of cell type of interest (contained in de_stats).
#' @param filename File name.
#' @param gzip Boolean. Should output file be zipped to spare memory?
#'
#' @export
export_markers_for_cell_type = function(de_stats, cell_type, filename, gzip=TRUE) {
    all_cell_types = unique(de_stats$cell_type)
    outgroups = all_cell_types[all_cell_types != cell_type]
    write_cell_type_header(cell_type, outgroups, filename)
    data.table::fwrite(dplyr::filter(de_stats, cell_type == cell_type),
                       filename, append=TRUE, row.names=FALSE, col.names=TRUE)
    if (gzip) {
        R.utils::gzip(filename, overwrite=TRUE)
    }
}
 
# Cell type header (time stamp, reference cell type, outgroups)
write_cell_type_header = function(reference, outgroups, filename) {
    lines = paste0(
        "# Differential expression statistics generated on ", Sys.Date(), ". ",
        "Reference: ", reference, ", outgroups: ", paste(outgroups, collapse = ", "), "."
    )
    writeLines(lines, filename)
}

#' Import marker statistics.
#'
#' @param filename File name.
#' @param header_size Size of header (# rows) containing meta-information such
#' as time stamp and list of cell types (for compatibility purposes only).
#'
#' @export
read_markers = function(filename, header_size=1) {
    result = data.table::fread(filename, header=TRUE, skip=header_size)
    return(tibble::as_tibble(result))
}
