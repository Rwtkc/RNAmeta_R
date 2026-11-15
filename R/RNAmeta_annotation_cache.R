# Annotation bundles are cached in memory only; no files are written.
.rnameta_annotation_cache <- new.env(hash = TRUE, parent = emptyenv())

.rnameta_annotation_file_state <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Annotation inputs must be existing file paths.")
  }
  # Keep the original path for file I/O: R on Windows may reject Unicode absolute paths
  # under a C locale even when a relative path to the same file is readable.
  info <- file.info(path)
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (isTRUE(info$isdir) || is.na(info$size) || is.na(info$mtime)) {
    stop("Annotation input is not a readable file: ", path)
  }
  list(path = path, key = resolved, size = info$size, mtime = as.numeric(info$mtime))
}

.rnameta_load_rda_object <- function(path) {
  destination <- new.env(parent = emptyenv())
  object_names <- load(path, envir = destination)
  if (length(object_names) != 1L) {
    stop("Expected exactly one object in annotation file: ", path)
  }
  get(object_names[[1L]], envir = destination, inherits = FALSE)
}

.rnameta_build_annotation_bundle <- function(txdb_file, tx_length_file, gff_file, number_of_bins) {
  txdb <- AnnotationDbi::loadDb(txdb_file)
  txlens <- data.table::as.data.table(.rnameta_load_rda_object(tx_length_file))
  gene_max_tx_lengths <- txlens[, lapply(.SD, max), by = gene_id]
  features <- .get_txdb_features_from_GRanges(
    txdb_file = txdb, tx_length_file = txlens, number_of_bins = number_of_bins
  )
  gff <- .rnameta_load_rda_object(gff_file)
  columns <- names(S4Vectors::mcols(gff))
  columns[columns == "gene_biotype"] <- "gene_type"
  columns[columns == "transcript_biotype"] <- "transcript_type"
  names(S4Vectors::mcols(gff)) <- columns
  gene_type_features <- .rnameta_gene_type_features(gff, txdb)
  features$annotation_index <- .rnameta_build_annotation_index(features, gene_type_features)

  list(species_id = sub("\\.txdb\\.sqlite$", "", basename(txdb_file)), txdb = txdb, txlens = txlens, gene_max_tx_lengths = gene_max_tx_lengths,
       txdb_features = features, gff = gff)
}

.rnameta_annotation_bundle <- function(txdb_file, tx_length_file, gff_file,
                                       number_of_bins = 100L, use_cache = TRUE) {
  paths <- lapply(list(txdb_file, tx_length_file, gff_file), .rnameta_annotation_file_state)
  number_of_bins <- as.integer(number_of_bins)
  if (length(number_of_bins) != 1L || is.na(number_of_bins) || number_of_bins < 1L) {
    stop("number_of_bins must be a positive integer.")
  }
  if (!is.logical(use_cache) || length(use_cache) != 1L || is.na(use_cache)) {
    stop("use_cache must be TRUE or FALSE.")
  }
  key <- paste(c(vapply(paths, `[[`, character(1), "key"), number_of_bins), collapse = "\n")
  state <- lapply(paths, function(entry) entry[c("size", "mtime")])
  if (use_cache && .rnameta_cache_has(.rnameta_annotation_cache, key)) {
    cached <- .rnameta_cache_get(.rnameta_annotation_cache, key)
    if (identical(cached$state, state)) return(cached$bundle)
  }
  bundle <- .rnameta_build_annotation_bundle(
    paths[[1L]]$path, paths[[2L]]$path, paths[[3L]]$path, number_of_bins
  )
  if (use_cache) {
    .rnameta_cache_set(.rnameta_annotation_cache, key, list(state = state, bundle = bundle))
  }
  bundle
}

#' Clear cached annotation bundles
#'
#' Clears only the in-memory annotation cache in the current R process. Source
#' annotation files and analysis results are not deleted.
#' @return Invisibly returns the number of cached bundles removed.
#' @export
RNAmeta_clear_annotation_cache <- function() {
  counts <- RNAmeta_clear_cache("annotation")
  invisible(sum(counts))
}
