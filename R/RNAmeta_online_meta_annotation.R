# Numerical kernels from deployed RNAmeta Online; see inst/online-algorithms.json.

.rnameta_online_meta_plot_generate_checking_ranges <- function(tx_information, component_types, checking_ranges_number = 500) {
  tx_number <- length(tx_information$names)
  component_type_number <- length(component_types)
  component_width <- matrix(0, tx_number, component_type_number)
  rownames(component_width) <- tx_information$names
  colnames(component_width) <- component_types
  end_point <- component_width
  start_point <- component_width

  for (component_type in component_types) {
    component_width[, component_type] <- sum(IRanges::width(tx_information[[component_type]]))
  }

  component_width_ratio <- component_width / rowSums(component_width)
  component_width_ratio_avg <- colSums(component_width_ratio)
  component_width_average <- floor(component_width_ratio_avg / sum(component_width_ratio_avg) * checking_ranges_number + 0.5)
  end_point <- t(apply(component_width, 1, cumsum))
  start_point[, 1] <- 0

  if (component_type_number > 1) {
    start_point[, seq(2, component_type_number)] <- end_point[, seq_len(component_type_number - 1)]
  }

  start_point <- start_point + 1

  list(
    componentWidth = component_width,
    startPoint = start_point,
    endPoint = end_point,
    componentWidthAverage = component_width_average
  )
}

.rnameta_online_meta_plot_primary_tx_names <- function(tx_lengths) {
  keep <- !is.na(tx_lengths$gene_id) & !is.na(tx_lengths$tx_len)
  split_index <- GenomicRanges::split(seq_len(nrow(tx_lengths))[keep], tx_lengths$gene_id[keep])
  primary_index <- BiocGenerics::unlist(
    lapply(split_index, function(index) index[[which.max(tx_lengths$tx_len[index])]]),
    use.names = FALSE
  )
  tx_lengths$tx_name[primary_index]
}

.rnameta_online_meta_plot_build_annotation_features <- function(txdb) {
  tx_lengths <- as.data.frame(GenomicFeatures::transcriptLengths(txdb))
  tx <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  tx_range <- base::range(tx)
  cds <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  utr5 <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  utr3 <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)

  list(
    tx_lengths = tx_lengths,
    tx = tx,
    tx_range = tx_range,
    tx_overlap_count = GenomicRanges::countOverlaps(tx, tx),
    cds = cds,
    utr5 = utr5,
    utr3 = utr3,
    all_mRNA_names = unique(c(names(utr5), names(utr3), names(cds)))
  )
}

.rnameta_online_meta_plot_extract_component <- function(
  txdb,
  annotation_features = NULL,
  txfiveutrMinLength = 100,
  txcdsMinLength = 100,
  txthreeutrMinLength = 100,
  txlongNcrnaMinLength = 100,
  txlncrnaOverlapmrna = TRUE,
  txpromoterLength = 1000,
  txtailLength = 1000,
  txAmblguity = 5,
  txPrimaryOnly = FALSE,
  pltTxType = c("tx", "mrna", "ncrna")
) {
  if (is.null(annotation_features)) {
    annotation_features <- .rnameta_online_meta_plot_build_annotation_features(txdb)
  }

  component <- list(txTypes = character(0))
  tx_lengths <- annotation_features$tx_lengths
  name_filter_tx <- if (txPrimaryOnly) .rnameta_online_meta_plot_primary_tx_names(tx_lengths) else tx_lengths$tx_name
  tx <- annotation_features$tx
  overlap_count <- annotation_features$tx_overlap_count
  name_filter_tx <- intersect(name_filter_tx, names(tx[overlap_count < (txAmblguity + 2)]))
  tx <- tx[name_filter_tx]
  tx_range <- annotation_features$tx_range[name_filter_tx]
  name_filter_tx <- names(tx_range)[S4Vectors::elementNROWS(tx_range) == 1L]
  tx <- tx[name_filter_tx]
  tx_range <- tx_range[name_filter_tx]
  cds <- annotation_features$cds
  utr5 <- annotation_features$utr5
  utr3 <- annotation_features$utr3
  mRNA_name <- Reduce(
    intersect,
    list(names(utr5)[sum(IRanges::width(utr5)) > txfiveutrMinLength], names(utr3)[sum(IRanges::width(utr3)) > txthreeutrMinLength], names(cds)[sum(IRanges::width(cds)) > txcdsMinLength])
  )
  name_filter_mRNA <- intersect(mRNA_name, name_filter_tx)
  ncRNA <- tx[setdiff(name_filter_tx, annotation_features$all_mRNA_names)]
  ncRNA_name <- names(ncRNA)[sum(IRanges::width(ncRNA)) > txlongNcrnaMinLength]

  if (!isTRUE(txlncrnaOverlapmrna)) {
    overlap_with_mRNA <- GenomicRanges::countOverlaps(ncRNA, tx[name_filter_mRNA])
    ncRNA_name <- intersect(ncRNA_name, names(ncRNA)[overlap_with_mRNA < 1])
  }

  promoter <- GenomicRanges::flank(tx_range, txpromoterLength, start = TRUE)
  tail <- GenomicRanges::flank(tx_range, txtailLength, start = FALSE)
  tx_gr <- BiocGenerics::unlist(tx)
  S4Vectors::mcols(tx_gr) <- NULL
  tx_with_flank <- GenomicRanges::reduce(GenomicRanges::split(c(BiocGenerics::unlist(promoter), tx_gr, BiocGenerics::unlist(tail)), names(c(BiocGenerics::unlist(promoter), tx_gr, BiocGenerics::unlist(tail)))))

  if ("tx" %in% pltTxType) {
    component$txTypes <- c(component$txTypes, "tx")
    component$tx <- list(
      componentTypes = c("promoter", "rna", "tail"),
      names = name_filter_tx,
      txWithFlank = tx_with_flank,
      tx = tx,
      promoter = promoter,
      rna = tx,
      tail = tail,
      txRange = tx_range
    )
  }

  if ("mrna" %in% pltTxType && length(name_filter_mRNA) > 0) {
    component$txTypes <- c(component$txTypes, "mrna")
    component$mrna <- list(
      componentTypes = c("promoter", "utr5", "cds", "utr3", "tail"),
      names = name_filter_mRNA,
      txWithFlank = tx_with_flank[name_filter_mRNA],
      tx = tx[name_filter_mRNA],
      promoter = promoter[name_filter_mRNA],
      utr5 = utr5[name_filter_mRNA],
      cds = cds[name_filter_mRNA],
      utr3 = utr3[name_filter_mRNA],
      tail = tail[name_filter_mRNA]
    )
  }

  if ("ncrna" %in% pltTxType && length(ncRNA_name) > 0) {
    component$txTypes <- c(component$txTypes, "ncrna")
    component$ncrna <- list(
      componentTypes = c("promoter", "ncrna", "tail"),
      names = ncRNA_name,
      txWithFlank = tx_with_flank[ncRNA_name],
      tx = tx[ncRNA_name],
      promoter = promoter[ncRNA_name],
      ncrna = tx[ncRNA_name],
      tail = tail[ncRNA_name]
    )
  }

  component
}

.rnameta_online_meta_plot_make_guitar_txdb <- function(txdb, controls, annotation_features = NULL) {
  tx_component <- .rnameta_online_meta_plot_extract_component(
    txdb = txdb,
    annotation_features = annotation_features,
    txfiveutrMinLength = controls$txfiveutrMinLength,
    txcdsMinLength = controls$txcdsMinLength,
    txthreeutrMinLength = controls$txthreeutrMinLength,
    txlongNcrnaMinLength = controls$txlongNcrnaMinLength,
    txlncrnaOverlapmrna = identical(controls$txlncrnaOverlapmrna, "TRUE"),
    txpromoterLength = controls$txpromoterLength,
    txtailLength = controls$txtailLength,
    txAmblguity = controls$stAmblguity,
    txPrimaryOnly = identical(controls$txPrimaryOnly, "TRUE"),
    pltTxType = controls$pltTxType
  )
  guitar_txdb <- list(txTypes = tx_component$txTypes)

  for (tx_type in tx_component$txTypes) {
    component_types <- tx_component[[tx_type]]$componentTypes
    checking_ranges <- .rnameta_online_meta_plot_generate_checking_ranges(tx_component[[tx_type]], component_types, 500)
    guitar_txdb[[tx_type]] <- list(
      tx = tx_component[[tx_type]]$txWithFlank,
      txLength = sum(IRanges::width(tx_component[[tx_type]]$txWithFlank)),
      componentWidth = checking_ranges$componentWidth,
      componentWidthPtc = checking_ranges$componentWidth / rowSums(checking_ranges$componentWidth),
      startPoint = checking_ranges$startPoint,
      endPoint = checking_ranges$endPoint,
      componentWidthAverage = checking_ranges$componentWidthAverage,
      componentWidthAverage_pct = checking_ranges$componentWidthAverage / sum(checking_ranges$componentWidthAverage)
    )
    cumulative_start <- cumsum(guitar_txdb[[tx_type]]$componentWidthAverage_pct)
    guitar_txdb[[tx_type]]$componentStartAverage_pct <- c(0, head(cumulative_start, -1))
    names(guitar_txdb[[tx_type]]$componentStartAverage_pct) <- component_types
  }

  guitar_txdb
}
