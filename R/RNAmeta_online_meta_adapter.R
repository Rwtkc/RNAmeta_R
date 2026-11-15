.rnameta_online_meta <- function(bed_file, txdb_file, plot_tx_type, sample_model, head_or_tail,
    set_group_name, set_ambiguity, set_sample_number, tx_five_utr_min_length, tx_cds_min_length,
    tx_three_utr_min_length, tx_long_ncrna_min_length, tx_lncrna_overlap_mrna,
    tx_promoter_length, tx_tail_length, tx_ambiguity, tx_primary_only, map_filter_transcript,
    enable_confidence_interval, overlap_index, site_length_index, adjustment_factor,
    confidence_interval_resampling_time, confidence_interval_range, seed) {
  seed <- .rnameta_validate_seed(seed)
  allowed <- c(tx = "tx", mRNA = "mrna", ncRNA = "ncrna")
  if (!length(plot_tx_type) || anyNA(plot_tx_type) || any(!plot_tx_type %in% names(allowed))) {
    stop("plot_tx_type must contain one or more of: tx, mRNA, ncRNA.")
  }
  group_count <- if (is.null(set_group_name)) length(bed_file) else length(set_group_name)
  prepared <- .bed_file_test_for_gene_meta(bed_file, .txdb_file_test(txdb_file), group_count)
  on.exit(unlink(prepared, force = TRUE), add = TRUE)
  txdb <- AnnotationDbi::loadDb(txdb_file)
  sites <- .getStGroup(prepared, set_group_name = set_group_name)
  controls <- list(stAmblguity = set_ambiguity, stSampleNum = set_sample_number,
    stSampleModle = sample_model, txfiveutrMinLength = tx_five_utr_min_length,
    txcdsMinLength = tx_cds_min_length, txthreeutrMinLength = tx_three_utr_min_length,
    txlongNcrnaMinLength = tx_long_ncrna_min_length, txlncrnaOverlapmrna = as.character(tx_lncrna_overlap_mrna),
    txpromoterLength = if (head_or_tail) tx_promoter_length else 0L,
    txtailLength = if (head_or_tail) tx_tail_length else 0L,
    txPrimaryOnly = as.character(tx_primary_only), mapFilterTranscript = as.character(map_filter_transcript),
    overlapIndex = overlap_index, siteLengthIndex = site_length_index)
  features <- .rnameta_meta_cached(paste("online-meta-features", .rnameta_file_fingerprint(txdb_file)),
    function() .rnameta_online_meta_plot_build_annotation_features(txdb))
  widgets <- lapply(unique(plot_tx_type), function(type) {
    controls$pltTxType <- unname(allowed[[type]])
    key <- paste("online-guitar-20261029", .rnameta_file_fingerprint(txdb_file),
      .rnameta_object_fingerprint(controls), sep = ":")
    guitar <- .rnameta_meta_cached(key, function() .rnameta_online_meta_plot_make_guitar_txdb(txdb, controls, features))
    if (!controls$pltTxType %in% guitar$txTypes) stop("No eligible transcripts for ", type, ".")
    normalized <- lapply(seq_along(sites), function(i) {
      map_key <- paste("online-meta-map", key, .rnameta_granges_fingerprint(sites[[i]]), seed)
      .rnameta_meta_cached(map_key, function() .rnameta_with_seed(.rnameta_seed_offset(seed, i - 1L), {
        mapped <- .rnameta_online_meta_plot_granges_list_map_to_transcripts(sites[[i]], map_filter_transcript,
          guitar[[controls$pltTxType]]$tx)
        points <- .rnameta_online_meta_plot_build_sample_point_ranges(mapped, controls)
        .rnameta_online_meta_plot_normalize(points, guitar, controls)
      }))
    })
    positions <- stats::setNames(lapply(normalized, `[[`, "position"), names(sites))
    weights <- stats::setNames(lapply(normalized, `[[`, "weight"), names(sites))
    density <- .rnameta_with_seed(seed, .rnameta_online_meta_plot_generate_density_ci(positions, weights,
      adjustment_factor, enable_confidence_interval, confidence_interval_resampling_time, confidence_interval_range))
    if (!nrow(density)) stop("Fewer than two valid sampled positions for density estimation.")
    if ("confidenceDown" %in% names(density)) {
      names(density)[names(density) == "confidenceDown"] <- "confidence_down"
      names(density)[names(density) == "confidenceUp"] <- "confidence_up"
    }
    payload <- .rnameta_meta_payload(density, guitar[[controls$pltTxType]]$componentWidthAverage_pct,
      paste("Distribution on", switch(type, tx = "Transcript", mRNA = "mRNA", ncRNA = "ncRNA")),
      show_ci = enable_confidence_interval, head_or_tail = head_or_tail,
      tx_promoter_length = tx_promoter_length, tx_tail_length = tx_tail_length)
    widget <- .rnameta_chart_widget(payload)
    attr(widget, "density_data") <- density
    widget
  })
  names(widgets) <- unique(plot_tx_type)
  if (length(widgets) == 1L) widgets[[1L]] else structure(widgets, class = c("RNAmeta_meta_plots", "list"))
}
