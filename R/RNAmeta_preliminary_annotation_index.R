.rnameta_annotation_feature_part <- function(region, annotation_type) {
  if (inherits(region, "GenomicRangesList")) {
    counts <- S4Vectors::elementNROWS(region)
    gr <- BiocGenerics::unlist(region, use.names = FALSE)
    tx_name <- rep(names(region), counts)
  } else {
    gr <- region
    metadata <- S4Vectors::mcols(gr)
    tx_name <- if ("tx_name" %in% names(metadata)) as.character(metadata$tx_name) else names(gr)
  }
  if (!length(gr)) return(GenomicRanges::GRanges())
  if (is.null(tx_name) || length(tx_name) != length(gr)) tx_name <- rep("", length(gr))
  tx_name[is.na(tx_name)] <- ""
  rank <- rep(NA_integer_, length(gr)); total <- rep(NA_integer_, length(gr))
  if (annotation_type %in% c("CDS", "Intron", "ExonCDS")) {
    if (inherits(region, "GenomicRangesList")) {
      group <- rep.int(seq_along(counts), counts)
      rank <- sequence(counts)
      total <- counts[group]
      nonempty <- which(counts > 0L)
      group_starts <- cumsum(c(1L, head(counts, -1L)))[nonempty]
      negative_group <- rep(FALSE, length(counts))
      negative_group[nonempty] <- as.character(GenomicRanges::strand(gr[group_starts])) == "-"
      negative <- negative_group[group]
      rank[negative] <- total[negative] - rank[negative] + 1L
      if (annotation_type == "ExonCDS" && "exon_rank" %in% names(S4Vectors::mcols(gr))) rank <- as.integer(gr$exon_rank)
    } else {
      group <- match(tx_name, unique(tx_name))
      rank <- ave(seq_along(gr), group, FUN = seq_along)
      total_by_group <- tabulate(group)
      total <- total_by_group[group]
      negative <- as.character(GenomicRanges::strand(gr)) == "-"
      rank[negative] <- total[negative] - rank[negative] + 1L
      if (annotation_type == "ExonCDS" && "exon_rank" %in% names(S4Vectors::mcols(gr))) rank <- as.integer(gr$exon_rank)
    }
  }
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    annotation_type = rep(annotation_type, length(gr)),
    tx_name = tx_name,
    feature_rank = rank,
    feature_total = total,
    feature_length = as.integer(GenomicRanges::width(gr)),
    annotation_value = rep(NA_character_, length(gr))
  )
  gr
}

.rnameta_gene_type_index_part <- function(features) {
  if (!length(features)) return(GenomicRanges::GRanges())
  values <- as.character(S4Vectors::mcols(features)$transcript_type)
  gr <- features
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    annotation_type = rep("GeneType", length(gr)), tx_name = rep("", length(gr)),
    feature_rank = rep(NA_integer_, length(gr)), feature_total = rep(NA_integer_, length(gr)),
    feature_length = as.integer(GenomicRanges::width(gr)), annotation_value = values
  )
  gr
}
.rnameta_build_annotation_index <- function(txdb_features, gene_type_features = NULL) {
  parts <- list(
    .rnameta_annotation_feature_part(txdb_features$exons, "ExonCDS"),
    .rnameta_annotation_feature_part(txdb_features$introns, "Intron"),
    .rnameta_annotation_feature_part(txdb_features$cds, "CDS"),
    .rnameta_annotation_feature_part(txdb_features$three_utrs, "three_utr"),
    .rnameta_annotation_feature_part(txdb_features$five_utrs, "five_utr"),
    .rnameta_annotation_feature_part(txdb_features$stop_codon, "stop_codon"),
    .rnameta_annotation_feature_part(txdb_features$promoters, "Promoter"),
    .rnameta_annotation_feature_part(txdb_features$start_codon, "start_codon"),
    .rnameta_annotation_feature_part(txdb_features$transcripts, "Transcript"),
    .rnameta_gene_type_index_part(gene_type_features)
  )
  index <- Reduce(c, parts, init = GenomicRanges::GRanges())
  unstranded <- index
  if (length(unstranded)) {
    ignore_strand <- as.character(S4Vectors::mcols(unstranded)$annotation_type) != "GeneType"
    GenomicRanges::strand(unstranded)[ignore_strand] <- "*"
  }
  list(stranded = index, unstranded = unstranded)
}

.rnameta_annotation_label <- function(type, tx_name, rank, total, feature_length) {
  if (type == "Intron") return(paste0("Intron (", tx_name, ", intron ", rank, " of ", total, ")"))
  if (type %in% c("CDS", "ExonCDS")) return(paste0("CDS (", tx_name, ", CDS ", rank, " of ", total, ", ", feature_length, "bp)"))
  if (type == "five_utr") return(paste0("UTR5 (", tx_name, ")"))
  if (type == "three_utr") return(paste0("UTR3 (", tx_name, ")"))
  if (type == "stop_codon") return(paste0("Stop Codon (", tx_name, ")"))
  if (type == "start_codon") return(paste0("Start Codon (", tx_name, ")"))
  paste0(type, " (", tx_name, ")")
}

.rnameta_annotation_vectors <- function(peaks, annotation_index, same_strand = FALSE) {
  index <- if (isTRUE(same_strand)) annotation_index$stranded else annotation_index$unstranded
  types <- c("ExonCDS", "Intron", "CDS", "three_utr", "five_utr", "stop_codon", "Promoter", "start_codon", "Transcript")
  output <- stats::setNames(lapply(types, function(type) rep(NA_character_, length(peaks))), types)
  hit_queries <- stats::setNames(lapply(types, function(type) integer()), types)
  if (!length(peaks) || !length(index)) return(list(values = output, queries = hit_queries, gene_types = character()))

  overlaps <- GenomicRanges::findOverlaps(peaks, index)
  if (!length(overlaps)) return(list(values = output, queries = hit_queries, gene_types = character()))
  query <- S4Vectors::queryHits(overlaps); subject <- S4Vectors::subjectHits(overlaps)
  metadata <- S4Vectors::mcols(index)
  hit_type <- as.character(metadata$annotation_type[subject])
  gene_types <- as.character(metadata$annotation_value[subject[hit_type == "GeneType"]])
  gene_types <- gene_types[!is.na(gene_types) & nzchar(gene_types)]
  for (type in types) {
    positions <- which(hit_type == type)
    if (!length(positions)) next
    positions <- positions[!duplicated(query[positions])]
    q <- query[positions]; s <- subject[positions]
    output[[type]][q] <- vapply(seq_along(s), function(i) {
      .rnameta_annotation_label(
        type, as.character(metadata$tx_name[s[[i]]]),
        as.integer(metadata$feature_rank[s[[i]]]),
        as.integer(metadata$feature_total[s[[i]]]),
        as.integer(metadata$feature_length[s[[i]]])
      )
    }, character(1))
    hit_queries[[type]] <- q
  }
  list(values = output, queries = hit_queries, gene_types = gene_types)
}

.get_genomic_annotation <- function(peaks = NULL, txdb_features = NULL, genomic_annotation_priority = NULL, same_strand = FALSE) {
  annotation_index <- txdb_features$annotation_index
  if (is.null(annotation_index)) annotation_index <- .rnameta_build_annotation_index(txdb_features)
  collected <- .rnameta_annotation_vectors(peaks, annotation_index, same_strand)
  values <- collected$values; queries <- collected$queries
  annotation <- rep("Intergenic", length(peaks))
  priority_map <- c(CDS = "CDS", `5UTR` = "five_utr", `3UTR` = "three_utr", Intron = "Intron", Promoter = "Promoter")
  for (priority in rev(genomic_annotation_priority)) {
    type <- unname(priority_map[priority])
    if (!length(type) || is.na(type)) next
    keep <- !is.na(values[[type]])
    annotation[keep] <- values[[type]][keep]
  }
  flag <- rep(FALSE, length(peaks))
  detail <- data.frame(genic = flag, Intergenic = flag, Promoter = flag, five_utr = flag,
                       three_utr = flag, CDS = flag, Intron = flag, stop_codon = flag)
  for (type in c("Promoter", "five_utr", "three_utr", "CDS", "Intron")) {
    if (length(queries[[type]])) detail[queries[[type]], type] <- TRUE
  }
  genic <- detail$CDS | detail$Intron | detail$five_utr | detail$three_utr
  detail$genic <- genic
  detail$Intergenic <- !genic
  gene_type_counts <- if (length(collected$gene_types)) data.table::as.data.table(table(collected$gene_types)) else data.table::data.table(V1 = character(), N = integer())
  if (nrow(gene_type_counts)) data.table::setnames(gene_type_counts, c("V1", "N"))
  list(
    gene_type_counts = gene_type_counts,
    annotation = annotation,
    annotation_exon = values$ExonCDS,
    annotation_stop_codon = values$stop_codon,
    annotation_start_codon = values$start_codon,
    annotation_transcript = values$Transcript,
    detail_genomic_annotation = detail
  )
}
