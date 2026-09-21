# Gene Type is intentionally independent of the expensive all-module pre-analysis.
.rnameta_gene_type_features <- function(gff, txdb = NULL) {
  metadata <- S4Vectors::mcols(gff)
  field <- if ("transcript_type" %in% names(metadata)) "transcript_type" else "transcript_biotype"
  if (!field %in% names(metadata)) stop("GFF annotation has no transcript biotype column.")
  values <- as.character(metadata[[field]])
  keep <- !is.na(values) & nzchar(values)
  if (!is.null(txdb)) {
    keep <- keep & as.character(GenomeInfoDb::seqnames(gff)) %in% GenomeInfoDb::seqlevels(txdb)
  }
  features <- gff[keep]
  S4Vectors::mcols(features) <- S4Vectors::DataFrame(transcript_type = values[keep])
  features
}

.rnameta_gene_type_bundle <- function(txdb_file, gff_file, use_cache = TRUE) {
  if (!is.logical(use_cache) || length(use_cache) != 1L || is.na(use_cache)) {
    stop("use_cache must be TRUE or FALSE.")
  }
  paths <- lapply(list(txdb_file, gff_file), .rnameta_annotation_file_state)
  key <- paste(c("gene-type", vapply(paths, `[[`, character(1), "key")), collapse = "\n")
  state <- lapply(paths, function(entry) entry[c("size", "mtime")])
  if (use_cache && .rnameta_cache_has(.rnameta_annotation_cache, key)) {
    cached <- .rnameta_cache_get(.rnameta_annotation_cache, key)
    if (identical(cached$state, state)) return(cached$bundle)
  }
  txdb <- AnnotationDbi::loadDb(paths[[1L]]$path)
  gff <- .rnameta_load_rda_object(paths[[2L]]$path)
  bundle <- list(txdb = txdb, features = .rnameta_gene_type_features(gff, txdb))
  if (use_cache) {
    .rnameta_cache_set(.rnameta_annotation_cache, key, list(state = state, bundle = bundle))
  }
  bundle
}

.rnameta_gene_type_read_bed <- function(input, txdb) {
  if (inherits(input, "GRanges")) {
    peaks <- input
  } else {
    .rnameta_annotation_file_state(input)
    bed <- data.table::fread(input, header = FALSE, sep = "\t", showProgress = FALSE)
    if (ncol(bed) < 3L) stop("BED input must have at least three columns.")
    starts <- suppressWarnings(as.numeric(bed[[2L]]))
    ends <- suppressWarnings(as.numeric(bed[[3L]]))
    if (anyNA(starts) || anyNA(ends) || any(starts < 0 | ends <= starts) ||
        any(starts != floor(starts) | ends != floor(ends)) ||
        any(ends > .Machine$integer.max)) {
      stop("BED start/end must be non-negative integers with end > start.")
    }
    strands <- if (ncol(bed) >= 6L) as.character(bed[[6L]]) else rep("*", nrow(bed))
    strands[is.na(strands) | !strands %in% c("+", "-", "*")] <- "*"
    peaks <- GenomicRanges::GRanges(
      seqnames = as.character(bed[[1L]]),
      ranges = IRanges::IRanges(start = as.integer(starts + 1), end = as.integer(ends)),
      strand = strands
    )
  }
  valid <- GenomeInfoDb::seqlevels(txdb)
  seq_names <- as.character(GenomeInfoDb::seqnames(peaks))
  alternate <- ifelse(startsWith(seq_names, "chr"), sub("^chr", "", seq_names), paste0("chr", seq_names))
  seq_names[!seq_names %in% valid] <- alternate[!seq_names %in% valid]
  keep <- seq_names %in% valid
  peaks <- peaks[keep]
  peaks <- GenomicRanges::GRanges(
    seqnames = seq_names[keep],
    ranges = IRanges::ranges(peaks),
    strand = as.character(GenomicRanges::strand(peaks))
  )
  GenomicRanges::resize(peaks, width = 1L, fix = "center")
}

.rnameta_gene_type_counts <- function(peaks, features, sample_name) {
  empty <- data.table::data.table(V1 = character(), N = integer(), Sample = character())
  if (!length(peaks) || !length(features)) return(empty)
  hits <- GenomicRanges::findOverlaps(peaks, features)
  if (!length(hits)) return(empty)
  values <- S4Vectors::mcols(features)$transcript_type[S4Vectors::subjectHits(hits)]
  counts <- data.table::as.data.table(table(as.character(values)))
  data.table::setnames(counts, c("V1", "N"))
  counts[, Sample := sample_name]
  counts
}

#' Fast Gene Type counts from BED and annotations
#'
#' Compute Gene Type counts without running the full preliminary analysis.
#' BED inputs use zero-based, half-open coordinates; each peak is counted at
#' its center. Multiple overlapping GFF features contribute multiple counts.
#' @param bed_file A BED path, a GRanges object, or a list of either, one per sample.
#' @param txdb_file Path to a TxDb SQLite annotation file.
#' @param gff_file Path to an RDA file containing one annotated GRanges object.
#' @param set_group_name Optional labels in input order; defaults to Group1, Group2, etc.
#' @param use_cache Reuse the in-memory Gene Type annotation bundle when unchanged.
#' @return A data.table with columns V1 (biotype), N (count), Sample, and Source_File;
#'   it can be passed directly to RNAmeta_gene_type_plot().
#' @export
RNAmeta_gene_type_data <- function(bed_file, txdb_file, gff_file,
                                   set_group_name = NULL, use_cache = TRUE) {
  files <- if (inherits(bed_file, "GRanges")) list(bed_file) else if (is.character(bed_file)) as.list(bed_file) else bed_file
  if (!is.list(files) || !length(files)) stop("bed_file must contain at least one BED path or GRanges object.")
  if (is.null(set_group_name)) set_group_name <- paste0("Group", seq_along(files))
  if (!is.character(set_group_name) || length(set_group_name) != length(files) ||
      anyNA(set_group_name) || any(!nzchar(set_group_name)) || anyDuplicated(set_group_name)) {
    stop("set_group_name must contain one distinct, nonempty label per BED input.")
  }
  bundle <- .rnameta_gene_type_bundle(txdb_file, gff_file, use_cache)
  results <- lapply(seq_along(files), function(i) {
    peaks <- .rnameta_gene_type_read_bed(files[[i]], bundle$txdb)
    counts <- .rnameta_gene_type_counts(peaks, bundle$features, set_group_name[[i]])
    source_name <- if (is.character(files[[i]])) basename(files[[i]]) else set_group_name[[i]]
    counts[, Source_File := source_name]
    counts
  })
  data.table::rbindlist(results, use.names = TRUE)
}
