#' RNAmeta package
#'
#' Genomic interval annotation, profiling and D3 visualization.
#'
#' @keywords internal
#' @importFrom data.table := setnames
#' @importFrom dplyr %>% filter group_by
#' @importFrom stringr str_match
#' @importFrom BiocGenerics width
#' @importFrom methods as is
#' @importFrom stats ave end start
#' @importFrom utils head read.table tail
"_PACKAGE"

utils::globalVariables(c(
  ".", ".I", ".N", ".SD", ":=", "anno", "annotation", "cds_len",
  "chr", "Count", "end", "exonID", "exonLoc", "feature", "Feature",
  "flank_size", "Frequency", "gene", "gene_id", "gene_name", "i.start", "Index",
  "location", "N", "name", "nexon", "peak", "peak_name", "Percentage",
  "rel_pos", "Sample", "seqnames", "Source_File", "start", "strand",
  "Transcript", "transcript_id", "transcriptID", "tx_len", "tx_name",
  "type", "utr3_len", "utr5_len", "V1", "V2", "V3", "V4", "V5",
  "V6", "width", "xmax", "xmin"
))