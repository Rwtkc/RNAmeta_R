# Local adapters for the deployed Online numerical kernels.
.rnameta_online_site_profile_cached_boundary_windows <- function(kind, species, flank_size, builder) {
  builder()
}

.rnameta_online_site_profile_cached_transcript_space_peaks <- function(source_file, species, cache_tag, builder) {
  builder()
}

.rnameta_online_peak_distribution_first_hit_index <- function(index_vector) which(!duplicated(index_vector))

.rnameta_online_raw_peaks <- function(preliminary) {
  shared <- attr(preliminary, "peak_annotation", exact = TRUE)
  if (!is.null(shared$raw_peaks)) return(shared$raw_peaks)
  .rnameta_preliminary_value(preliminary, "peak_gr_list", 2L)
}

.rnameta_online_profile_data <- function(preliminary, kind) {
  features <- .rnameta_preliminary_value(preliminary, "txdb_features", 7L)
  txlens <- .rnameta_preliminary_value(preliminary, "tx_lengths", 10L)
  peaks <- .rnameta_online_raw_peaks(preliminary)
  flank <- .rnameta_preliminary_value(preliminary,
    if (kind == "transcription") "flanking_size" else "short_flanking_size",
    if (kind == "transcription") 6L else 9L)
  shared <- attr(preliminary, "peak_annotation", exact = TRUE)
  species <- if (!is.null(shared$annotation_bundle$species_id)) shared$annotation_bundle$species_id else "local"
  annotation_key <- paste(species, .rnameta_granges_fingerprint(features$exons),
    .rnameta_object_fingerprint(txlens), sep = ":")
  # Keep the special yeast policy identifiable while including the annotation
  # fingerprint in the outer cache key for every local reference bundle.
  windows <- .rnameta_site_profile_cached(paste("online-windows", kind, annotation_key, flank, sep = ":"), function() {
    if (kind == "transcription") {
      .rnameta_online_site_profile_prepare_transcription_windows(annotation_key, features, txlens, flank)
    } else if (kind == "translation") {
      .rnameta_online_site_profile_prepare_translation_windows(species, txlens, features$exons, flank)
    } else {
      .rnameta_online_site_profile_prepare_splicesite_windows(annotation_key, features$exons, flank)
    }
  })
  sample_names <- names(peaks)
  if (is.null(sample_names)) sample_names <- paste0("Group", seq_along(peaks))
  samples <- lapply(seq_along(peaks), function(i) {
    key <- paste("online-profile-20261029", kind, annotation_key, flank, sample_names[[i]],
      .rnameta_granges_fingerprint(peaks[[i]]), sep = ":")
    .rnameta_site_profile_cached(key, function() {
      query <- peaks[[i]]
      if (kind == "translation") {
        query <- GenomicFeatures::mapToTranscripts(query, features$exons)
        GenomicRanges::strand(query) <- "*"
      }
      if (kind == "splicesite") {
        starts <- windows$windows_5p; ends <- windows$windows_3p
        labels <- c("5PSS", "3PSS"); offsets <- c(windows$offset_5p, windows$offset_3p)
        limits <- list(c(-flank, flank * 2), c(-flank * 2, flank))
      } else {
        starts <- if (kind == "transcription") windows$coverage_tss else windows$tss_windows
        ends <- if (kind == "transcription") windows$coverage_tes else windows$tes_windows
        labels <- paste0(if (kind == "transcription") "Transcription" else "Translation", c("_start_site", "_end_site"))
        offsets <- rep(windows$rel_offset, 2L); limits <- rep(list(windows$x_domain), 2L)
      }
      hits <- lapply(seq_len(2L), function(j) .rnameta_online_site_profile_collect_boundary_hits(
        query, if (j == 1L) starts else ends, sample_names[[i]], labels[[j]], offsets[[j]],
        limits[[j]][[1]], limits[[j]][[2]], flank))
      heat <- NULL
      if (kind != "splicesite") {
        target <- if (kind == "transcription") .rnameta_online_site_profile_heatmap_peaks(peaks[[i]]) else {
          value <- GenomicFeatures::mapToTranscripts(.rnameta_online_site_profile_center_peaks(peaks[[i]]), features$exons)
          GenomicRanges::strand(value) <- "*"
          value
        }
        heat_starts <- if (kind == "transcription") windows$heatmap_tss else starts
        heat_ends <- if (kind == "transcription") windows$heatmap_tes else ends
        matrices <- lapply(list(heat_ends, heat_starts), function(ranges) {
          ranges <- .rnameta_online_site_profile_evenly_sample_windows(ranges, 4000L)
          .rnameta_online_site_profile_prepare_heatmap_matrix(.rnameta_online_site_profile_score_matrix_rna(target, ranges))
        })
        heat <- list(windows$x_positions, windows$x_positions, matrices[[1L]], matrices[[2L]])
      }
      list(hits = data.table::rbindlist(hits, use.names = TRUE, fill = TRUE), heatmap = heat)
    })
  })
  hits <- data.table::rbindlist(lapply(samples, `[[`, "hits"), use.names = TRUE, fill = TRUE)
  if (kind == "splicesite") return(hits)
  list(hits, stats::setNames(lapply(samples, `[[`, "heatmap"), sample_names))
}

.rnameta_online_exon_components <- function(peaks, features, txlens, sample_name) {
  exons <- .rnameta_online_peak_exon_prepare_exon_ranges(features$exons)
  transcripts <- .rnameta_online_peak_exon_prepare_transcript_ranges(features$transcripts)
  lookup <- .rnameta_online_peak_exon_prepare_transcript_lookup(txlens, features$exons)
  list(exon = .rnameta_online_peak_exon_compute_exon_size_data(peaks, exons, sample_name),
       transcript = .rnameta_online_peak_exon_compute_transcript_stats(peaks, transcripts, lookup, sample_name))
}
