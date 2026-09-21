library(data.table)
source("R/RNAmeta_gene_type_export.R")
source("R/RNAmeta_chart_core.R")
source("R/RNAmeta_chart_export.R")
source("R/RNAmeta_peak_distribution_plot.R")
source("R/RNAmeta_peak_exon_type_plot.R")
source("R/RNAmeta_peak_exon_number_plot.R")
source("R/RNAmeta_peak_exon_size_plot.R")
source("R/RNAmeta_peak_gene_size_plot.R")

make_sample <- function(name, offset = 0L) {
  list(
    data.table(Feature = c("CDS", "Intron", "UTR3"), Frequency = c(20, 35, 8) + offset, Sample = name),
    data.table(exonID = paste0(name, 1:6), length = c(100, 180, 240, 400, 800, 1600) + offset,
               exonLoc = rep(c("First", "Middle", "Last"), 2), Sample = name),
    data.table(), data.table(),
    data.table(tx_name = paste0(name, 1:6), nexon = c(2, 3, 4, 6, 8, 12),
               tx_len = c(500, 800, 1200, 2000, 3500, 6000) + offset, Sample = name)
  )
}
data <- list(list(make_sample("Group1"), make_sample("Group2", 5L)))
labels <- c("Up", "Down")
plots <- list(
  RNAmeta_peak_distribution_plot(data, labels, animate = FALSE),
  RNAmeta_peak_exon_type_plot(data, labels, animate = FALSE),
  RNAmeta_peak_exon_number_plot(data, labels),
  RNAmeta_peak_exon_size_plot(data, labels),
  RNAmeta_peak_gene_size_plot(data, labels)
)
stopifnot(all(vapply(plots, inherits, logical(1), "htmlwidget")))
out <- tempfile("rnameta-peak-charts-")
dir.create(out)
on.exit(unlink(out, recursive = TRUE, force = TRUE), add = TRUE)
files <- file.path(out, c("distribution.png", "exon-type.pdf", "exon-number.png", "exon-size.pdf", "gene-size.png"))
RNAmeta_peak_distribution_export(data, files[[1]], labels, width = 1600, aspect_ratio = c(16, 9), scale = 1)
RNAmeta_peak_exon_type_export(data, files[[2]], labels)
RNAmeta_peak_exon_number_export(data, files[[3]], labels, scale = 1)
RNAmeta_peak_exon_size_export(data, files[[4]], labels)
RNAmeta_peak_gene_size_export(data, files[[5]], labels, scale = 1)
stopifnot(all(file.exists(files)), all(file.info(files)$size > 1000))
con <- file(files[[1]], "rb"); bytes <- as.integer(readBin(con, "raw", 24)); close(con)
dimensions <- c(sum(bytes[17:20] * 256^(3:0)), sum(bytes[21:24] * 256^(3:0)))
stopifnot(identical(dimensions, c(1600, 900)))
cat("peak_chart_widgets: TRUE\npeak_chart_exports: TRUE\ncustom_dimensions: TRUE\n")
