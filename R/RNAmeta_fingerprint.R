.rnameta_file_fingerprint <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !file.exists(path)) stop("Cannot fingerprint missing file: ", paste(path, collapse = ""))
  unname(tools::md5sum(path)[[1L]])
}

.rnameta_object_fingerprint <- function(object) {
  path <- tempfile(pattern = "rnameta_fingerprint_", fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  saveRDS(object, path, version = 2, compress = FALSE)
  .rnameta_file_fingerprint(path)
}

.rnameta_granges_fingerprint <- function(granges) {
  if (inherits(granges, "GenomicRangesList")) granges <- BiocGenerics::unlist(granges, use.names = FALSE)
  if (is.list(granges) && !inherits(granges, "GenomicRanges")) return(paste(vapply(granges, .rnameta_granges_fingerprint, character(1)), collapse = ":"))
  .rnameta_object_fingerprint(list(
    seqnames = as.character(GenomeInfoDb::seqnames(granges)),
    start = as.integer(IRanges::start(granges)), end = as.integer(IRanges::end(granges)),
    strand = as.character(BiocGenerics::strand(granges))
  ))
}

.rnameta_input_fingerprint <- function(input) {
  if (is.character(input) && length(input) == 1L && !is.na(input) && file.exists(input)) return(.rnameta_file_fingerprint(input))
  if (inherits(input, "GenomicRanges")) return(.rnameta_granges_fingerprint(input))
  .rnameta_object_fingerprint(input)
}
