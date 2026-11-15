# Numerical kernels from deployed RNAmeta Online; see inst/online-algorithms.json.

.rnameta_online_peak_exon_prepare_transcript_lookup <- function(txlens, exons = NULL) {
  txlens_dt <- data.table::as.data.table(txlens)

  if (!"tx_name" %in% names(txlens_dt)) {
    row_names <- rownames(as.data.frame(txlens))
    if (length(row_names) == nrow(txlens_dt) && any(nzchar(row_names))) {
      txlens_dt[, tx_name := row_names]
    }
  }

  if (!"nexon" %in% names(txlens_dt) && !is.null(exons)) {
    exon_counts <- data.table::data.table(
      tx_name = names(exons),
      nexon = S4Vectors::elementNROWS(exons)
    )
    txlens_dt <- merge(txlens_dt, exon_counts, by = "tx_name", all.x = TRUE)
  }

  txlens_dt[, tx_len := suppressWarnings(as.numeric(tx_len))]
  txlens_dt[, nexon := suppressWarnings(as.numeric(nexon))]
  unique(txlens_dt[, .(tx_name, tx_len, nexon)], by = "tx_name")
}

.rnameta_online_peak_exon_prepare_exon_ranges <- function(exons) {
  exon_ranges <- BiocGenerics::unlist(exons)
  exon_ranges$tx_name <- names(exon_ranges)
  exon_ranges$length <- IRanges::width(exon_ranges)

  exon_counts <- data.table::data.table(name = exon_ranges$tx_name)[, .N, by = name]
  exon_count_values <- exon_counts$N
  names(exon_count_values) <- exon_counts$name

  # exonsBy(..., by = "tx") supplies exon_rank in transcript order on both
  # strands. Reversing negative-strand ranks would swap First and Last.
  exon_ranges$rank <- as.integer(exon_ranges$exon_rank)

  exon_total <- unname(exon_count_values[as.character(exon_ranges$tx_name)])
  exon_loc <- rep("Middle", length(exon_ranges))
  exon_loc[exon_ranges$rank == 1] <- "First"
  exon_loc[!is.na(exon_total) & exon_ranges$rank == exon_total] <- "Last"

  exon_ranges$exonLoc <- factor(exon_loc, levels = .rnameta_online_peak_exon_location_levels())
  exon_ranges$exonID <- paste0(exon_ranges$tx_name, ":", exon_ranges$rank, recycle0 = TRUE)
  exon_ranges
}

.rnameta_online_peak_exon_prepare_transcript_ranges <- function(transcripts) {
  transcript_ranges <- transcripts

  if (!"tx_name" %in% names(S4Vectors::mcols(transcript_ranges))) {
    if ("tx_name" %in% names(S4Vectors::mcols(transcript_ranges))) {
      transcript_ranges$tx_name <- as.character(transcript_ranges$tx_name)
    } else if (!is.null(names(transcript_ranges)) && any(nzchar(names(transcript_ranges)))) {
      transcript_ranges$tx_name <- names(transcript_ranges)
    } else if ("tx_id" %in% names(S4Vectors::mcols(transcript_ranges))) {
      transcript_ranges$tx_name <- as.character(transcript_ranges$tx_id)
    } else {
      transcript_ranges$tx_name <- NA_character_
    }
  }

  transcript_ranges[!is.na(transcript_ranges$tx_name) & nzchar(transcript_ranges$tx_name)]
}

.rnameta_online_peak_exon_first_hit_subject_index <- function(peaks, feature_ranges) {
  overlaps <- GenomicRanges::findOverlaps(peaks, feature_ranges, ignore.strand = TRUE)

  if (!length(overlaps)) {
    return(integer())
  }

  subject_hits <- S4Vectors::subjectHits(overlaps[.rnameta_online_peak_distribution_first_hit_index(S4Vectors::queryHits(overlaps))])
  unique(subject_hits)
}

.rnameta_online_peak_exon_location_levels <- function() {
  c("First", "Middle", "Last")
}

.rnameta_online_peak_exon_compute_exon_size_data <- function(peaks, exon_ranges, sample_name) {
  subject_index <- .rnameta_online_peak_exon_first_hit_subject_index(peaks, exon_ranges)

  if (!length(subject_index)) {
    return(data.table::data.table(
      exonID = character(),
      tx_name = character(),
      rank = numeric(),
      length = numeric(),
      exonLoc = character(),
      Sample = character()
    ))
  }

  exon_info <- data.table::as.data.table(as.data.frame(S4Vectors::mcols(exon_ranges[subject_index])))
  output_columns <- intersect(
    c("exonID", "exon_id", "exon_name", "exon_rank", "tx_name", "length", "rank", "exonLoc"),
    names(exon_info)
  )
  exon_info <- unique(exon_info[, ..output_columns], by = "exonID")
  exon_info[, Sample := sample_name]
  exon_info[, exonLoc := factor(exonLoc, levels = .rnameta_online_peak_exon_location_levels())]
  exon_info[order(Sample, exonLoc, rank)]
}

.rnameta_online_peak_exon_compute_transcript_stats <- function(peaks, transcript_ranges, transcript_lookup, sample_name) {
  subject_index <- .rnameta_online_peak_exon_first_hit_subject_index(peaks, transcript_ranges)

  if (!length(subject_index)) {
    return(data.table::data.table(
      tx_name = character(),
      tx_len = numeric(),
      nexon = numeric(),
      Sample = character()
    ))
  }

  transcript_ids <- unique(as.character(transcript_ranges$tx_name[subject_index]))
  transcript_stats <- transcript_lookup[tx_name %in% transcript_ids, .(tx_name, tx_len, nexon)]
  transcript_stats <- transcript_stats[!is.na(tx_name)]
  transcript_stats[, Sample := sample_name]
  transcript_stats[order(Sample, tx_name)]
}

.rnameta_online_peak_gene_size_compute_lengths <- function(peaks, transcript_ranges, txlens, sample_name) {
  transcript_ranges_unstranded <- transcript_ranges
  GenomicRanges::strand(transcript_ranges_unstranded) <- "*"
  overlaps <- GenomicRanges::findOverlaps(peaks, transcript_ranges_unstranded)

  if (!length(overlaps)) {
    return(data.table::data.table(
      tx_name = character(),
      tx_len = numeric(),
      Sample = character()
    ))
  }

  hit_index <- which(!duplicated(S4Vectors::queryHits(overlaps)))
  subject_index <- S4Vectors::subjectHits(overlaps[hit_index])
  transcript_hits <- unique(as.character(transcript_ranges$tx_name[subject_index]))

  txlens_dt <- data.table::as.data.table(txlens)
  lengths_dt <- txlens_dt[tx_name %in% transcript_hits, .(tx_name, tx_len)]
  lengths_dt <- unique(lengths_dt, by = "tx_name")
  lengths_dt <- lengths_dt[!is.na(tx_len) & tx_len > 0]
  lengths_dt[, Sample := sample_name]
  lengths_dt[order(Sample, tx_len)]
}
