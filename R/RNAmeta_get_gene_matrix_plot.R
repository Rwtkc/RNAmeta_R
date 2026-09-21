.rnameta_gene_matrix_feature_ranges <- function(feature) {
  if (inherits(feature, "GRangesList") || inherits(feature, "CompressedGRangesList")) {
    ranges <- BiocGenerics::unlist(feature)
    ranges$tx_name <- names(ranges)
  } else {
    ranges <- feature
    columns <- names(S4Vectors::mcols(ranges))
    if (!"tx_name" %in% columns) {
      if (!is.null(names(ranges)) && any(nzchar(names(ranges)))) ranges$tx_name <- names(ranges)
      else if ("tx_id" %in% columns) ranges$tx_name <- as.character(ranges$tx_id)
      else ranges$tx_name <- NA_character_
    }
  }
  keep <- !is.na(ranges$tx_name) & nzchar(as.character(ranges$tx_name))
  ranges <- ranges[keep]
  if (length(ranges)) ranges$tx_name <- as.character(ranges$tx_name)
  ranges
}

.rnameta_gene_matrix_gene_sets_online <- function(preliminary_analysis_data) {
  peaks <- .rnameta_preliminary_value(preliminary_analysis_data, "peak_gr_list", 2L)
  features <- .rnameta_preliminary_value(preliminary_analysis_data, "txdb_features", 7L)
  gff <- .rnameta_preliminary_value(preliminary_analysis_data, "gff", 3L)
  feature_sets <- list(
    promoters = features$promoters,
    fiveUTRs = features$five_utrs,
    cds = features$cds,
    threeUTRs = features$three_utrs,
    introns = features$introns,
    stopCodon = features$stop_codon,
    transcripts = features$transcripts
  )
  feature_sets <- lapply(feature_sets, .rnameta_gene_matrix_feature_ranges)
  transcript_meta <- data.table::as.data.table(as.data.frame(gff))
  transcript_meta <- transcript_meta[as.character(type) == "transcript" & !is.na(transcript_id) & nzchar(as.character(transcript_id)) & !is.na(gene_name) & nzchar(as.character(gene_name)), .(transcript_id = as.character(transcript_id), gene_name = as.character(gene_name))]
  transcript_meta <- transcript_meta[, .(gene_name = gene_name[[1L]]), by = transcript_id]
  gene_lookup <- stats::setNames(transcript_meta$gene_name, transcript_meta$transcript_id)
  lapply(peaks, function(sample_peaks) {
    query <- GenomicRanges::resize(sample_peaks, width = 1L, fix = "center")
    transcript_ids <- unique(unlist(lapply(feature_sets, function(feature_ranges) {
      if (!length(feature_ranges)) return(character())
      overlaps <- GenomicRanges::findOverlaps(query, feature_ranges, ignore.strand = TRUE)
      if (!length(overlaps)) return(character())
      first <- !duplicated(S4Vectors::queryHits(overlaps))
      unique(as.character(feature_ranges$tx_name[S4Vectors::subjectHits(overlaps)[first]]))
    }), use.names = FALSE))
    gene_names <- unname(gene_lookup[transcript_ids])
    sort(unique(as.character(stats::na.omit(gene_names[nzchar(gene_names)]))))
  })
}

.rnameta_gene_matrix_sets <- function(preliminary_analysis_data, set_group_name = NULL) {
  can_use_online <- is.list(preliminary_analysis_data) &&
    !is.null(.rnameta_preliminary_value(preliminary_analysis_data, "peak_gr_list", 2L)) &&
    !is.null(.rnameta_preliminary_value(preliminary_analysis_data, "txdb_features", 7L)) &&
    !is.null(.rnameta_preliminary_value(preliminary_analysis_data, "gff", 3L))
  sets <- if (can_use_online) {
    .rnameta_gene_matrix_gene_sets_online(preliminary_analysis_data)
  } else {
    shared <- attr(preliminary_analysis_data, "peak_annotation", exact = TRUE)
    entries <- if (!is.null(shared)) shared$components$gene else {
      if (!is.list(preliminary_analysis_data) || !length(.rnameta_preliminary_value(preliminary_analysis_data, "cct_list", 1L))) stop("Provide preliminary analysis data.")
      lapply(.rnameta_preliminary_value(preliminary_analysis_data, "cct_list", 1L), function(entry) entry[[3]])
    }
    lapply(entries, function(entry) {
      if (!"gene_name" %in% names(entry)) stop("Gene Matrix annotation lacks gene_name.")
      sort(unique(as.character(stats::na.omit(entry$gene_name))))
    })
  }
  count <- length(sets)
  if (count < 2L || count > 5L) stop("Gene Matrix requires between 2 and 5 samples; received ", count, ".")
  labels <- if (is.null(set_group_name)) paste0("Group", seq_len(count)) else as.character(set_group_name)
  if (length(labels) != count || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) stop("set_group_name must contain one unique nonempty name per sample.")
  names(sets) <- labels
  sets
}
#' Prepare Gene Matrix intersections
#' @param preliminary_analysis_data Preliminary Analysis output.
#' @param set_group_name Optional sample labels.
#' @return Gene sets and all nonempty set combinations.
#' @export
RNAmeta_gene_matrix_data <- function(preliminary_analysis_data = NULL, set_group_name = NULL) {
  sets <- .rnameta_gene_matrix_sets(preliminary_analysis_data, set_group_name)
  labels <- names(sets)
  intersections <- unlist(lapply(seq_along(labels), function(size) {
    lapply(utils::combn(labels, size, simplify = FALSE), function(combo) {
      genes <- if (length(combo) == 1L) sets[[combo]] else Reduce(intersect, sets[combo])
      list(sets = unname(as.list(combo)), size = length(genes), genes = unname(as.list(genes)))
    })
  }), recursive = FALSE)
  list(sampleLabels = labels, geneSets = sets, intersections = intersections)
}

.rnameta_gene_matrix_payload <- function(data) {
  labels <- data$sampleLabels; colors <- c("#d94841", "#2563eb", "#f59e0b", "#0f9d78", "#7c3aed")[seq_along(labels)]
  list(type = "gene_matrix", title = "Gene Matrix", sampleLabels = unname(labels),
       colors = stats::setNames(as.list(colors), labels), intersections = data$intersections)
}

#' Plot Gene Matrix intersections
#' @param preliminary_analysis_data Preliminary Analysis output or RNAmeta_gene_matrix_data() result.
#' @param set_group_name Optional sample labels.
#' @param width,height Optional widget dimensions.
#' @return An htmlwidget.
#' @export
RNAmeta_get_gene_matrix_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL, width = NULL, height = NULL) {
  data <- if (is.list(preliminary_analysis_data) && all(c("sampleLabels", "intersections") %in% names(preliminary_analysis_data))) preliminary_analysis_data else RNAmeta_gene_matrix_data(preliminary_analysis_data, set_group_name)
  .rnameta_chart_widget(.rnameta_gene_matrix_payload(data), width, height)
}

#' Export a Gene Matrix plot
#' @param data Preliminary Analysis output, Gene Matrix data, or its widget.
#' @param file PNG or PDF output path.
#' @param set_group_name Optional sample labels.
#' @param width,height,aspect_ratio,scale Export dimensions and PNG scale.
#' @return Output path invisibly.
#' @export
RNAmeta_gene_matrix_export <- function(data, file, set_group_name = NULL, width = 1180, height = NULL, aspect_ratio = c(1120, 700), scale = 2) {
  config <- if (inherits(data, "htmlwidget") && identical(data$x$config$type, "gene_matrix")) data$x$config else {
    prepared <- if (is.list(data) && all(c("sampleLabels", "intersections") %in% names(data))) data else RNAmeta_gene_matrix_data(data, set_group_name)
    .rnameta_gene_matrix_payload(prepared)
  }
  .rnameta_chart_export(config, file, width, height, aspect_ratio, scale)
}
