.rnameta_peak_annotation_cache <- new.env(hash = TRUE, parent = emptyenv())

.rnameta_peak_annotation_sample <- function(bed, group_name, txdb_file, txdb_path, txdb_features, tx_length_file, gff_file, genomic_annotation_priority, annotation_key, use_cache = TRUE) {
  key <- paste("peak-annotation-online-20261029-v2", group_name, .rnameta_input_fingerprint(bed), annotation_key, sep = ":")
  if (isTRUE(use_cache) && .rnameta_cache_has(.rnameta_peak_annotation_cache, key)) return(.rnameta_cache_get(.rnameta_peak_annotation_cache, key))
  original <- .batch_peak_reader(group_name = group_name, txdb_file = txdb_file, bed_file = bed, txdb_file_path = txdb_path)
  motif <- if (length(original) && max(GenomicRanges::width(original)) < 15L) GenomicRanges::resize(original, width = 15L, fix = "center") else original
  peak <- GenomicRanges::resize(original, width = 1L, fix = "center")
  peak_list <- stats::setNames(list(peak), group_name)
  peak_dt <- stats::setNames(list(data.table::as.data.table(peak)), group_name)
  components <- .batch_genomic_annotation(1L, peak_gr = peak_list, txdb_features = txdb_features, genomic_annotation_priority = genomic_annotation_priority, tx_length_file = tx_length_file, gff_file = gff_file, peak_dt = peak_dt, group_name = group_name)
  gene_type_counts <- attr(components, "gene_type_counts", exact = TRUE)
  result <- list(components = components, raw_peak = original, peak = peak, motif_peak = motif, gene_type_counts = gene_type_counts)
  if (isTRUE(use_cache)) .rnameta_cache_set(.rnameta_peak_annotation_cache, key, result)
  result
}

#' Compute and cache shared peak annotation results
#' @param bed_file BED paths or GRanges objects.
#' @param txdb_file,tx_length_file,gff_file Annotation files.
#' @param use_cache Reuse annotation bundles and per-sample peak annotations.
#' @return A named shared annotation result used by Preliminary Analysis and plots.
#' @export
RNAmeta_peak_annotation_data <- function(bed_file = NULL, txdb_file = NULL, tx_length_file = NULL, gff_file = NULL, use_cache = TRUE) {
  if (!length(bed_file)) stop("bed_file must contain at least one BED path or GRanges object.")
  number_of_bins <- 100L
  priority <- c("stop_codon", "start_codon", "CDS", "5UTR", "3UTR", "Intron", "Transcript", "Promoter", "Intergenic")
  bundle <- .rnameta_annotation_bundle(txdb_file, tx_length_file, gff_file, number_of_bins, use_cache)
  annotation_key <- paste("annotation-online-20261029", .rnameta_file_fingerprint(txdb_file), .rnameta_file_fingerprint(tx_length_file), .rnameta_file_fingerprint(gff_file), number_of_bins, sep = ":")
  samples <- lapply(seq_along(bed_file), function(i) .rnameta_peak_annotation_sample(bed_file[[i]], paste0("Group", i), bundle$txdb, txdb_file, bundle$txdb_features, data.table::copy(bundle$txlens), bundle$gff, priority, annotation_key, use_cache))
  sample_names <- paste0("Group", seq_along(samples)); names(samples) <- sample_names
  component_names <- c("distribution", "exon", "gene", "detail", "transcript")
  components <- stats::setNames(lapply(seq_along(component_names), function(index) stats::setNames(lapply(samples, function(sample) sample$components[[index]]), sample_names)), component_names)
  structure(list(samples = sample_names, sample_results = lapply(samples, `[[`, "components"), components = components, peaks = stats::setNames(lapply(samples, `[[`, "peak"), sample_names), raw_peaks = stats::setNames(lapply(samples, `[[`, "raw_peak"), sample_names), motif_peaks = stats::setNames(lapply(samples, `[[`, "motif_peak"), sample_names), gene_type_counts = stats::setNames(lapply(samples, `[[`, "gene_type_counts"), sample_names), annotation_bundle = bundle, annotation_key = annotation_key), class = "RNAmeta_peak_annotation")
}
