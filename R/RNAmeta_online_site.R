# Numerical kernels from deployed RNAmeta Online; see inst/online-algorithms.json.

.rnameta_online_site_profile_center_peaks <- function(peaks) {
  GenomicRanges::resize(peaks, width = 1, fix = "center")
}

.rnameta_online_site_profile_heatmap_peaks <- function(peaks) {
  GenomicRanges::resize(peaks, width = 15, fix = "center")
}

.rnameta_online_site_profile_query_table <- function(query_ranges, sample_name) {
  peak_table <- data.table::as.data.table(query_ranges)
  peak_table[, index := .I]
  peak_table[, Sample := sample_name]
  peak_table
}

.rnameta_online_site_profile_empty_hits <- function(sample_name, feature_label, x_min, x_max, flank_size) {
  data.table::data.table(
    seqnames = character(),
    start = integer(),
    end = integer(),
    width = integer(),
    strand = character(),
    Sample = sample_name,
    feature = feature_label,
    rel_pos = numeric(),
    flankSize = flank_size,
    xmin = x_min,
    xmax = x_max
  )
}

.rnameta_online_site_profile_collect_boundary_hits <- function(query_ranges, windows, sample_name, feature_label, rel_offset, x_min, x_max, flank_size) {
  hits <- GenomicFeatures::mapToTranscripts(query_ranges, windows, ignore.strand = FALSE)

  if (!length(hits)) {
    return(.rnameta_online_site_profile_empty_hits(sample_name, feature_label, x_min, x_max, flank_size))
  }

  peak_table <- .rnameta_online_site_profile_query_table(query_ranges, sample_name)
  hits_table <- data.table::as.data.table(hits)
  data.table::setnames(hits_table, c("seqnames", "xHits"), c("transcriptID", "index"))

  peak_hits <- peak_table[hits_table, on = "index"]
  peak_hits[, c("index", "transcriptsHits") := NULL]
  peak_hits[, feature := feature_label]
  peak_hits[, rel_pos := i.start - rel_offset]
  peak_hits[, flankSize := flank_size]
  peak_hits[, xmin := x_min]
  peak_hits[, xmax := x_max]
  peak_hits
}

.rnameta_online_site_profile_prepare_transcription_windows <- function(species, txdb_features, txlens, flank_size = 1000) {
  .rnameta_online_site_profile_cached_boundary_windows("transcription", species, flank_size, function() {
    # One true transcript span per isoform, including its UTRs.
    transcripts <- BiocGenerics::unlist(base::range(txdb_features$exons))
    boundary_window <- function(fix) {
      anchor <- GenomicRanges::resize(transcripts, width = 1, fix = fix)
      symmetric <- GenomicRanges::resize(anchor, width = 2 * flank_size + 1, fix = "center")
      GenomicRanges::resize(symmetric, width = 2 * flank_size, fix = "start")
    }
    coverage_tss <- boundary_window("start")
    coverage_tes <- boundary_window("end")
    tx_names <- unique(as.character(data.table::as.data.table(txlens)$tx_name))
    heatmap_tss <- coverage_tss[names(coverage_tss) %in% tx_names]
    heatmap_tes <- coverage_tes[names(coverage_tes) %in% tx_names]

    list(
      coverage_tss = coverage_tss,
      coverage_tes = coverage_tes,
      heatmap_tss = heatmap_tss,
      heatmap_tes = heatmap_tes,
      x_positions = seq(-flank_size, flank_size - 1),
      x_domain = c(-flank_size, flank_size - 1),
      rel_offset = flank_size + 1L,
      flank_size = flank_size
    )
  })
}

.rnameta_online_site_profile_prepare_translation_windows <- function(species, txlens, exons, flank_size = 100) {
  .rnameta_online_site_profile_cached_boundary_windows("translation", species, flank_size, function() {
    txlens_dt <- data.table::as.data.table(txlens)

    use_cds_only <- identical(species, "sce_R64") && !any(
      txlens_dt$utr5_len > flank_size & txlens_dt$utr3_len > flank_size &
        txlens_dt$cds_len > flank_size, na.rm = TRUE
    )
    if (use_cds_only) {
      usable <- txlens_dt[cds_len > flank_size, .(tx_name, utr5_len, cds_len, utr3_len)]
      tss_ranges <- usable[, .(
        seqnames = tx_name,
        start = utr5_len + 1,
        end = utr5_len + flank_size,
        strand = "*"
      )]
      tes_ranges <- usable[, .(
        seqnames = tx_name,
        start = utr5_len + cds_len - flank_size + 1,
        end = utr5_len + cds_len,
        strand = "*"
      )]
      x_domain <- c(1, flank_size)
      x_positions <- seq_len(flank_size)
    } else {
      usable <- txlens_dt[utr5_len > flank_size & utr3_len > flank_size & cds_len > flank_size, .(tx_name, utr5_len, cds_len, utr3_len)]
      tss_ranges <- usable[, .(
        seqnames = tx_name,
        start = utr5_len + 1 - flank_size,
        end = utr5_len + 1 + flank_size,
        strand = "*"
      )]
      tes_ranges <- usable[, .(
        seqnames = tx_name,
        start = utr5_len + cds_len - flank_size,
        end = utr5_len + cds_len + flank_size,
        strand = "*"
      )]
      x_domain <- c(-flank_size, flank_size)
      x_positions <- seq(-flank_size, flank_size)
    }

    tss_windows <- GenomicRanges::makeGRangesFromDataFrame(tss_ranges, keep.extra.columns = TRUE)
    tes_windows <- GenomicRanges::makeGRangesFromDataFrame(tes_ranges, keep.extra.columns = TRUE)
    names(tss_windows) <- as.character(tss_ranges$seqnames)
    names(tes_windows) <- as.character(tes_ranges$seqnames)
    GenomeInfoDb::seqlevels(tss_windows) <- unique(c(GenomeInfoDb::seqlevels(tss_windows), names(exons)))
    GenomeInfoDb::seqlevels(tes_windows) <- unique(c(GenomeInfoDb::seqlevels(tes_windows), names(exons)))

    list(
      tss_windows = tss_windows,
      tes_windows = tes_windows,
      x_positions = x_positions,
      x_domain = x_domain,
      rel_offset = if (use_cds_only) 0L else flank_size + 1L,
      flank_size = flank_size
    )
  })
}

.rnameta_online_site_profile_prepare_splicesite_windows <- function(species, exons_by_tx, flank_size = 100) {
  .rnameta_online_site_profile_cached_boundary_windows("splicesite", species, flank_size, function() {
    exons_flat <- BiocGenerics::unlist(exons_by_tx)
    exons_flat <- exons_flat[IRanges::width(exons_flat) > flank_size]
    transcript_ranges <- BiocGenerics::unlist(base::range(exons_by_tx))

    ss5p <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(exons_flat),
      ranges = IRanges::IRanges(
        ifelse(as.character(GenomicRanges::strand(exons_flat)) == "+", BiocGenerics::end(exons_flat), BiocGenerics::start(exons_flat)),
        ifelse(as.character(GenomicRanges::strand(exons_flat)) == "+", BiocGenerics::end(exons_flat), BiocGenerics::start(exons_flat))
      ),
      strand = GenomicRanges::strand(exons_flat)
    )
    names(ss5p) <- names(exons_flat)

    ss3p <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(exons_flat),
      ranges = IRanges::IRanges(
        ifelse(as.character(GenomicRanges::strand(exons_flat)) == "+", BiocGenerics::start(exons_flat), BiocGenerics::end(exons_flat)),
        ifelse(as.character(GenomicRanges::strand(exons_flat)) == "+", BiocGenerics::start(exons_flat), BiocGenerics::end(exons_flat))
      ),
      strand = GenomicRanges::strand(exons_flat)
    )
    names(ss3p) <- names(exons_flat)

    transcript5p <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(transcript_ranges),
      ranges = IRanges::IRanges(
        ifelse(as.character(GenomicRanges::strand(transcript_ranges)) == "+", BiocGenerics::end(transcript_ranges), BiocGenerics::start(transcript_ranges)),
        ifelse(as.character(GenomicRanges::strand(transcript_ranges)) == "+", BiocGenerics::end(transcript_ranges), BiocGenerics::start(transcript_ranges))
      ),
      strand = GenomicRanges::strand(transcript_ranges)
    )

    transcript3p <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(transcript_ranges),
      ranges = IRanges::IRanges(
        ifelse(as.character(GenomicRanges::strand(transcript_ranges)) == "+", BiocGenerics::start(transcript_ranges), BiocGenerics::end(transcript_ranges)),
        ifelse(as.character(GenomicRanges::strand(transcript_ranges)) == "+", BiocGenerics::start(transcript_ranges), BiocGenerics::end(transcript_ranges))
      ),
      strand = GenomicRanges::strand(transcript_ranges)
    )

    ss5p <- ss5p[setdiff(seq_along(ss5p), unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(ss5p, transcript5p))))]
    ss3p <- ss3p[setdiff(seq_along(ss3p), unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(ss3p, transcript3p))))]

    windows_5p <- GenomicRanges::resize(
      GenomicRanges::resize(ss5p, width = flank_size * 2 + 1L, fix = "center"),
      width = flank_size * 3 + 1L,
      fix = "start"
    )
    windows_3p <- GenomicRanges::resize(
      GenomicRanges::resize(ss3p, width = flank_size * 4 + 1L, fix = "center"),
      width = flank_size * 3 + 1L,
      fix = "start"
    )

    list(
      windows_5p = windows_5p,
      windows_3p = windows_3p,
      offset_5p = flank_size + 1L,
      offset_3p = flank_size * 2 + 1L
    )
  })
}

.rnameta_online_site_profile_map_to_transcript_space <- function(source_file, species, peaks, exons, cache_tag = "default") {
  .rnameta_online_site_profile_cached_transcript_space_peaks(source_file, species, cache_tag = cache_tag, function() {
    transcript_peaks <- GenomicFeatures::mapToTranscripts(peaks, exons)
    if (length(transcript_peaks)) {
      GenomicRanges::strand(transcript_peaks) <- "*"
    }
    transcript_peaks
  })
}

.rnameta_online_site_profile_evenly_sample_windows <- function(windows, sample_n = 10000) {
  if (is.null(windows) || !length(windows) || length(windows) <= sample_n) {
    return(windows)
  }

  sampled_index <- unique(round(seq(1, length(windows), length.out = sample_n)))
  windows[sampled_index]
}

.rnameta_online_site_profile_score_matrix_rna <- function(target, windows) {
  widths <- unique(IRanges::width(windows))
  if (length(widths) > 1L) stop("Heatmap windows must have equal widths.")
  ncols <- if (length(widths)) widths[[1]] else 0L
  if (!length(windows)) return(matrix(numeric(), 0L, ncols))
  chromosomes <- as.character(GenomeInfoDb::seqnames(windows))
  limits <- GenomeInfoDb::seqlengths(windows)[chromosomes]
  valid <- IRanges::start(windows) >= 1L & (is.na(limits) | IRanges::end(windows) <= limits)
  windows <- windows[valid]
  chromosomes <- chromosomes[valid]
  result <- matrix(0, nrow = length(windows), ncol = ncols)
  if (!length(windows) || !length(target)) return(result)
  # Pad coverage to the full annotation window; never shrink its coordinate axis.
  coverage_widths <- tapply(IRanges::end(windows), chromosomes, max)
  target_levels <- GenomeInfoDb::seqlevels(target)
  coverage_widths <- coverage_widths[intersect(names(coverage_widths), target_levels)]
  target <- GenomeInfoDb::keepSeqlevels(target, names(coverage_widths), pruning.mode = "coarse")
  target_rle <- GenomicRanges::coverage(target)
  for (chromosome in names(coverage_widths)) {
    indices <- which(chromosomes == chromosome)
    signal <- target_rle[[chromosome]]
    required <- coverage_widths[[chromosome]]
    if (length(signal) < required) signal <- c(signal, S4Vectors::Rle(0, required - length(signal)))
    views <- IRanges::Views(signal, start = IRanges::start(windows)[indices], end = IRanges::end(windows)[indices])
    result[indices, ] <- t(IRanges::viewApply(views, as.vector))
  }
  minus <- as.character(GenomicRanges::strand(windows)) == "-"
  if (any(minus)) result[minus, ] <- result[minus, ncols:1L, drop = FALSE]
  result
}

.rnameta_online_site_profile_density_points <- function(hit_table, x_min, x_max, feature_order = NULL, point_count = 512) {
  if (is.null(hit_table) || !nrow(hit_table)) {
    return(data.table::data.table(
      Panel = character(),
      Sample = character(),
      Position = numeric(),
      Density = numeric()
    ))
  }

  density_rows <- lapply(GenomicRanges::split(hit_table, by = c("feature", "Sample"), keep.by = TRUE), function(sample_hits) {
    density_values <- sample_hits$rel_pos[is.finite(sample_hits$rel_pos)]
    if (!length(density_values)) {
      density_values <- c(0, 0)
    } else if (length(density_values) < 2) {
      density_values <- rep(density_values[[1]], 2)
    }

    density_curve <- stats::density(
      density_values,
      from = x_min,
      to = x_max,
      n = point_count,
      bw = if (length(unique(density_values)) < 2) 1 else "nrd0",
      na.rm = TRUE
    )

    data.table::data.table(
      Panel = unique(as.character(sample_hits$feature))[[1]],
      Sample = unique(as.character(sample_hits$Sample))[[1]],
      Position = as.numeric(density_curve$x),
      Density = as.numeric(density_curve$y)
    )
  })

  density_data <- data.table::rbindlist(density_rows, use.names = TRUE, fill = TRUE)

  if (!is.null(feature_order) && length(feature_order)) {
    density_data[, Panel := factor(Panel, levels = feature_order)]
    data.table::setorder(density_data, Panel, Sample, Position)
    density_data[, Panel := as.character(Panel)]
  }

  density_data
}

.rnameta_online_site_profile_prepare_heatmap_matrix <- function(tag_matrix) {
  if (is.null(tag_matrix)) return(matrix(numeric(), nrow = 0, ncol = 0))
  if (!nrow(tag_matrix)) return(tag_matrix)

  row_order <- order(rowSums(tag_matrix))
  tag_matrix <- tag_matrix[row_order, , drop = FALSE]
  non_zero_rows <- rowSums(tag_matrix) > 0

  tag_matrix <- tag_matrix[non_zero_rows, , drop = FALSE]

  tag_matrix
}
