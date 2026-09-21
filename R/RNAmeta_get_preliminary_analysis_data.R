#' @title Get Preliminary Analysis Data
#'
#' @param bed_file A list of BED file paths or GRanges objects. Each entry in the list
#' represents a set of genomic regions for analysis.
#' @param txdb_file A file path to a TxDb (transcript database) file. The TxDb file
#' contains genomic annotations related to transcripts, exons, and other genomic features.
#' @param tx_length_file A file path to a data file containing transcript length data.
#' This file should include information about the lengths of different transcripts.
#' @param gff_file A file path to a GFF (General Feature Format) file that provides
#' genomic annotations, including gene and transcript information. The file should
#' contain columns such as gene ID and feature types (e.g., "gene", "CDS", "exon").
#' @param use_cache Reuse the in-memory annotation bundle when inputs are unchanged.
#'   Set to FALSE to force a fresh load.
#' @returns A list containing the following elements:
#' - `cct_list`: A list of results from batch genomic annotation processes for each
#'   BED file.
#' - `peak_gr_list`: A list of `GRanges` objects representing peaks for each group.
#' - gff: The processed GFF GRanges object.
#' - `peak_gr_for_motifscan`: A list of `GRanges` objects resized for motif scanning.
#' - `motif_parameter`: A vector containing motif scanning parameters extracted from input.
#' - `flanking_size`: The flanking size (default: 1000 bp) used for genomic analyses.
#' - `txdb_features`: A list of transcript-related genomic features (e.g., exons, promoters, CDS) extracted from the TxDb file.
#' - `gene_max_tx_lengths`: A data table containing the maximum transcript lengths
#'   for each gene, aiding in gene length normalization in downstream analysis.
#' - `short_flanking_size`: A smaller flanking size (default: 100 bp) used in certain genomic analysis steps.
#' - tx_lengths: A processed transcript-length table.
#'
#' The named fields are preferred. Legacy positional access at positions 1 through 10 remains supported.
#'
#' @export
#'
RNAmeta_get_preliminary_analysis_data <- function(bed_file = NULL, txdb_file = NULL, tx_length_file = NULL, gff_file = NULL, use_cache = TRUE){
  motif_parameter <- BiocGenerics::unlist(strsplit("6,7,8_1000", "_"))
  shared <- RNAmeta_peak_annotation_data(bed_file, txdb_file, tx_length_file, gff_file, use_cache)
  bundle <- shared$annotation_bundle
  preliminary_analysis_data <- list(
    cct_list = shared$sample_results,
    peak_gr_list = shared$peaks,
    gff = bundle$gff,
    peak_gr_for_motifscan = shared$motif_peaks,
    motif_parameter = motif_parameter,
    flanking_size = 1000,
    txdb_features = bundle$txdb_features,
    gene_max_tx_lengths = data.table::copy(bundle$gene_max_tx_lengths),
    short_flanking_size = 100,
    tx_lengths = data.table::copy(bundle$txlens)
  )
  attr(preliminary_analysis_data, "peak_annotation") <- shared
  class(preliminary_analysis_data) <- c("RNAmeta_preliminary", "list")
  preliminary_analysis_data
}
