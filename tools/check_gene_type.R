# Run from the package root: Rscript tools/check_gene_type.R
source("R/RNAmeta_cache_core.R")
source("R/RNAmeta_annotation_cache.R")
source("R/RNAmeta_gene_type_data.R")
source("R/RNAmeta_gene_type_plot_dependency_functions.R")
source("R/RNAmeta_gene_type_widget.R")
source("R/RNAmeta_gene_type_export.R")
source("R/RNAmeta_gene_type_plot.R")
root <- "RNAmeta_test_data"
db <- file.path(root, "osa_IRGSP_1.txdb.sqlite")
gff <- file.path(root, "osa_IRGSP_1.gff.rda")
beds <- file.path(root, c("rice_up.bed", "rice_Down.bed"))
labels <- c("up", "down")
t0 <- proc.time()[["elapsed"]]
result <- RNAmeta_gene_type_data(beds, db, gff, labels)
cat("cold_gene_type_seconds:", round(proc.time()[["elapsed"]] - t0, 3), "\n")
stopifnot(all(c("V1", "N", "Sample") %in% names(result)), identical(unique(result$Sample), labels))
bundle <- .rnameta_gene_type_bundle(db, gff)
raw_gff <- .rnameta_load_rda_object(gff)
for (i in seq_along(beds)) {
  peaks <- .rnameta_gene_type_read_bed(beds[[i]], bundle$txdb)
  start <- proc.time()[["elapsed"]]
  reference <- .peak_type(peaks, raw_gff)
  old_counts <- table(as.character(reference$transcript_biotype))
  old_time <- proc.time()[["elapsed"]] - start
  start <- proc.time()[["elapsed"]]
  new_counts <- .rnameta_gene_type_counts(peaks, bundle$features, labels[[i]])
  new_time <- proc.time()[["elapsed"]] - start
  actual <- stats::setNames(as.integer(new_counts$N), as.character(new_counts$V1))
  stopifnot(identical(as.integer(actual[names(old_counts)]), as.integer(old_counts)))
  cat(labels[[i]], "old_seconds:",round(old_time,3), "new_seconds:",round(new_time,3),
      "count:",sum(new_counts$N), "\n")
}
again <- RNAmeta_gene_type_data(beds, db, gff, labels)
stopifnot(identical(result, again), length(ls(.rnameta_annotation_cache)) == 1L)
cat("cached_result_equal: TRUE\n")
empty <- GenomicRanges::GRanges("chr1",IRanges::IRanges(1,1))
none <- RNAmeta_gene_type_data(empty, db, gff)
stopifnot(nrow(none) == 0L)
cat("empty_overlaps: TRUE\n")
widget <- RNAmeta_gene_type_plot(result, animate = FALSE)
stopifnot(inherits(widget, "htmlwidget"), identical(widget$x$config$type, "grouped_bar"))
cat("interactive_plot_from_counts: TRUE\n")
legacy_peaks <- lapply(beds, function(bed) .rnameta_gene_type_read_bed(bed, bundle$txdb))
names(legacy_peaks) <- c("Group1", "Group2")
legacy <- list(vector("list", length(legacy_peaks)), legacy_peaks, raw_gff)
legacy_widget <- RNAmeta_gene_type_plot(legacy, animate = FALSE)
stopifnot(inherits(legacy_widget, "htmlwidget"))
cat("legacy_plot_compatible: TRUE\n")
multiple_peaks <- c(legacy_peaks[[1]][1], legacy_peaks[[1]][2])
gr_result <- RNAmeta_gene_type_data(multiple_peaks, db, gff)
stopifnot(identical(sum(gr_result$N), sum(.rnameta_gene_type_counts(multiple_peaks, bundle$features, "Group1")$N)))
cat("single_granges_input: TRUE\n")
invalid_path <- tempfile(fileext = ".bed")
writeLines("chr1\tbad\t10", invalid_path)
invalid_result <- try(RNAmeta_gene_type_data(invalid_path, db, gff), silent = TRUE)
unlink(invalid_path)
stopifnot(inherits(invalid_result, "try-error"), grepl("BED start/end",as.character(invalid_result), fixed=TRUE))
cat("invalid_bed_rejected: TRUE\n")
png_file <- tempfile(fileext = ".png")
pdf_file <- tempfile(fileext = ".pdf")
wide_file <- tempfile(fileext = ".png")
exact_file <- tempfile(fileext = ".png")
RNAmeta_gene_type_export(result, png_file)
RNAmeta_gene_type_export(result, pdf_file)
RNAmeta_gene_type_export(result, wide_file, width = 1600, aspect_ratio = c(16, 9), scale = 1)
RNAmeta_gene_type_export(result, exact_file, width = 1400, height = 700, scale = 1)
png_dimensions <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  bytes <- as.integer(readBin(con, "raw", 24))
  c(width = sum(bytes[17:20] * 256^(3:0)), height = sum(bytes[21:24] * 256^(3:0)))
}
stopifnot(file.info(png_file)$size > 1000, file.info(pdf_file)$size > 1000,
          identical(readBin(png_file, "raw", 8), as.raw(c(137,80,78,71,13,10,26,10))),
          identical(rawToChar(readBin(pdf_file, "raw", 4)), "%PDF"),
          identical(unname(png_dimensions(wide_file)), c(1600, 900)),
          identical(unname(png_dimensions(exact_file)), c(1400, 700)))
unlink(c(png_file, pdf_file, wide_file, exact_file))
cat("headless_png_pdf_export: TRUE\n")
cat("custom_aspect_and_dimensions: TRUE\n")
