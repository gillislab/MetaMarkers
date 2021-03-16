#' @param marker_set A marker table containing the following columns: group, cell_type and gene.
#' @export
score_cells = function(expr, marker_set) {
    if (is.data.frame(marker_set)) {
        marker_matrix = marker_table_to_matrix(marker_set, rownames(expr))
    } else if (is.matrix(marker_set)) {
        marker_matrix = marker_set
    } else {
        marker_matrix = marker_list_to_matrix(marker_set, rownames(expr))        
    }
    as.matrix(Matrix::crossprod(marker_matrix, expr))
}

#' @export
marker_list_to_matrix = function(marker_list, known_genes, weighted=TRUE) {
    marker_list %>%
        tibble::enframe("cell_type", "gene") %>%
        tidyr::unnest(c("gene")) %>%
        marker_table_to_matrix(known_genes, weighted)
}

#' @export
marker_table_to_matrix = function(marker_table, known_genes, weighted=TRUE) {
    marker_table = marker_table %>%
        dplyr::mutate(cell_type = paste(.data$group, .data$cell_type, sep = "|")) %>%
        dplyr::select(.data$cell_type, .data$gene) %>%
        dplyr::filter(.data$gene %in% known_genes)
    if (weighted) {
        x = marker_table %>%
            dplyr::group_by(cell_type) %>%
            dplyr::mutate(w = 1/dplyr::n()) %>%
            dplyr::pull(.data$w)
    } else {
        x = 1
    }
    cell_type = as.factor(marker_table$cell_type)
    result = Matrix::sparseMatrix(
        i = match(marker_table$gene, known_genes),
        j = as.numeric(cell_type),
        x = x,
        dims = c(length(known_genes), length(levels(cell_type))),
        dimnames = list(known_genes, levels(cell_type))
    )
    return(result)
}

#' @export
compute_marker_enrichment = function(scores, by_group=TRUE) {
    if (by_group) {
        group = get_group(rownames(scores))
        if(any(group == "")) {
            warning("No group information: assuming that all cell types are at the same level.")
            by_group = FALSE
        }
    }
    if (by_group) {
        tscores = t(scores+0.0001)
        result = lapply(unique(group), function(g) {
            scores_g = tscores[, group == g, drop=FALSE]
            subresult = scores_g / rowMeans(scores_g)
            dimnames(subresult) = dimnames(scores_g)
            subresult
        })
        result = t(do.call(cbind, result))
    } else {
        result = scores+0.0001
        result = matrixStats::t_tx_OP_y(result, colMeans(result), "/")
        dimnames(result) = dimnames(scores)
    }
    return(result)
}

#' @export
get_group = function(cell_type_label) {
    if (all(grepl("\\|", cell_type_label))) {
        sapply(strsplit(cell_type_label, "|", TRUE), "[", 1)
    } else {
        return(rep("", length(cell_type_label)))
    }
}

#' @export
get_cell_type = function(cell_type_label) {
    if (all(grepl("\\|", cell_type_label))) {
        sapply(lapply(strsplit(cell_type_label, "|", TRUE), "[", -1), paste, collapse="|")
    } else {
        return(cell_type_label)
    }
}

#' @export
assign_cells = function(scores, group_assignment=NULL) {
    group = get_group(rownames(scores))
    unique_groups = unique(group)
    if (is.null(group_assignment) & length(unique_groups)>1) {
        warning("Detected group information in the score matrix (",
                paste(unique_groups, collapse = ", "),
                ") but no group assignments were provided. ",
                "Usually this means that marker sets were intended to use ",
                "hierarchically and groups needed to be assigned first. ",
                "Proceeding without group information (any cell can be any cell type).")
    }
    if (!is.null(group_assignment)) {
        unique_assignments = unique(group_assignment)
        unknown_assignments = unique_assignments[!(unique_assignments %in% unique_groups)]
        if (any(unique_groups == "")) {
            stop("Group assignments provided but no group information in the ",
                 "score matrix. Did you provide group information to score_cells?")
        }
        if (length(unknown_assignments) > 0) {
            warning("Some group assignments (",
                    paste(unknown_assignments, collapse = ", "),
                    ") do not match groups in the score matrix (",
                    paste(unique_groups, collapse = ", "),
                    ") and will result in NA predictions.")
        }
    }
    
    if (is.null(group_assignment)) {
        result = assign_cells_(scores)
    } else {
        result = lapply(unique_assignments, function(group_name) {
            keep_cell = group_assignment == group_name
            keep_ct = group == group_name
            assign_cells_(scores[keep_ct, keep_cell, drop=FALSE]) %>%
                tibble::add_column(cell_id = seq_along(keep_cell)[keep_cell])
        }) %>%
            dplyr::bind_rows() %>%
            dplyr::arrange(.data$cell_id) %>%
            dplyr::select(-.data$cell_id)
    }
    return(result)
}

assign_cells_ = function(scores) {
    if (nrow(scores) == 0) {
        return(data.frame(group = rep(NA, ncol(scores)), predicted = NA,
                          score = NA, enrichment = NA))
    }
    tscores = t(scores)
    first_index = max.col(tscores, ties.method = "first")
    first_value = scores[cbind(first_index, seq_len(ncol(scores)))]
    tscores[cbind(seq_len(ncol(scores)), first_index)] = 0
    second_index = max.col(tscores, ties.method = "first")
    second_value = scores[cbind(second_index, seq_len(ncol(scores)))]

    label = rownames(scores)[first_index]
    group_label = get_group(label)
    label = get_cell_type(label)
    enrichment = (first_value+0.0001) / colMeans(scores+0.0001)
    is_tie = first_value == second_value
    label[is_tie] = "unassigned"
    #enrichment[is_tie] = 0
    
    return(data.frame(group = group_label, predicted = label,
                      score = first_value, enrichment = enrichment))
}

#' @export
summarize_auroc = function(scores, true_labels) {
    if (!is.matrix(true_labels)) {
        true_labels = design_matrix(true_labels)
    }
    result = compute_aurocs(t(scores), true_labels)$aurocs
    return(result)
}

#' @export
summarize_fold_change = function(scores, true_labels) {
    if (!is.matrix(true_labels)) {
        true_labels = design_matrix(true_labels)
    }
    scaled_positives = scale(true_labels, center=FALSE, scale=colSums(true_labels))
    scaled_negatives = 1-true_labels
    scaled_negatives = scale(scaled_negatives, center=FALSE, scale=colSums(scaled_negatives))
    positive_expression = as.matrix(scores %*% scaled_positives)
    negative_expression = as.matrix(scores %*% scaled_negatives)
    return(positive_expression / negative_expression)
}

#' @export
summarize_precision_recall = function(scores, true_labels, score_threshold) {
    if (!is.matrix(true_labels)) {
        true_labels = design_matrix(true_labels)
    }
    
    if (length(score_threshold) > 1) {
        result = lapply(score_threshold, function(t) {
            summarize_precision_recall(scores, true_labels, t)
        })
        result = dplyr::bind_rows(result)
        return(result)
    }
    
    binary_scores = scores > score_threshold
    pp = rowSums(binary_scores)
    tp = binary_scores %*% true_labels
    fp = pp - tp
    p = colSums(true_labels)
    n = nrow(true_labels) - p
    
    result = list(
        precision = tp / pp,
        recall = t(t(tp) / p),
        fpr = t(t(fp) / n)
    )
    result$f1 = 2/(1/result$precision + 1/result$recall)
    result$balanced_accuracy = (1 + result$recall - result$fpr) / 2
    result = tibble::as_tibble(sapply(result, as.vector)) %>%
        tibble::add_column(marker_set = rep(rownames(scores), ncol(true_labels)), .before=1) %>%
        tibble::add_column(true_label = rep(colnames(true_labels), each = nrow(scores)), .after=1) %>%
        tibble::add_column(score_threshold = score_threshold, .before=1)
    return(result)
}
