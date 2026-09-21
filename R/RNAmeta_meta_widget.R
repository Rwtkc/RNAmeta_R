.rnameta_meta_component_payload <- function(component_width, peak, head_or_tail = TRUE,
                                            tx_promoter_length = 1000,
                                            tx_tail_length = 1000) {
  widths <- as.numeric(component_width)
  names(widths) <- names(component_width)
  if (!isTRUE(head_or_tail) && length(names(widths))) widths[names(widths) %in% c("promoter", "tail")] <- 0
  if (sum(widths[is.finite(widths) & widths > 0]) > 0) widths <- widths / sum(widths[is.finite(widths) & widths > 0])
  keep <- is.finite(widths) & widths > 0
  widths <- widths[keep]
  names(widths) <- tolower(names(widths))
  if (!length(widths)) return(list(trackBase = -0.03 * peak, labelY = -0.08 * peak, segments = list(), separators = list()))
  ends <- cumsum(widths)
  starts <- c(0, head(ends, -1)) + 0.001
  labels <- names(widths)
  labels[labels == "promoter"] <- if (isTRUE(head_or_tail)) sprintf("Promoter(%skb)", format(tx_promoter_length / 1000, trim = TRUE)) else ""
  labels[labels == "tail"] <- if (isTRUE(head_or_tail)) sprintf("Tail(%skb)", format(tx_tail_length / 1000, trim = TRUE)) else ""
  labels[labels == "utr5"] <- "5'UTR"
  labels[labels == "utr3"] <- "3'UTR"
  labels[labels == "cds"] <- "CDS"
  labels[labels == "ncrna"] <- "ncRNA"
  labels[labels == "rna"] <- "RNA"
  segments <- lapply(seq_along(widths), function(i) {
    component <- names(widths)[[i]]
    list(
      start = starts[[i]], end = ends[[i]], mid = (starts[[i]] + ends[[i]]) / 2,
      label = labels[[i]], component = component,
      alpha = if (component == "cds") 0.272 else if (component %in% c("ncrna", "rna")) 0.2 else 0.99,
      height = if (component == "cds") 0.01 * peak else if (component %in% c("ncrna", "rna")) 0.01 * peak else if (component %in% c("utr5", "utr3")) 0.005 * peak else 0.002 * peak
    )
  })
  list(
    trackBase = -0.03 * peak,
    labelY = -0.085 * peak,
    segments = segments,
    separators = lapply(head(ends, -1), function(value) list(x = unname(value), y1 = 1.05 * peak, y2 = -0.03 * peak))
  )
}

.rnameta_meta_payload <- function(density_data, component_width, title = "Meta Plot",
                                  show_ci = FALSE, head_or_tail = TRUE,
                                  tx_promoter_length = 1000,
                                  tx_tail_length = 1000) {
  density_data <- data.table::as.data.table(density_data)
  group_col <- if ("group" %in% names(density_data)) "group" else if ("samples" %in% names(density_data)) "samples" else stop("Meta Plot density data lacks a group column.")
  groups <- unique(as.character(density_data[[group_col]])); colors <- .rnameta_chart_palette()
  peak_values <- density_data$density
  if (isTRUE(show_ci) && "confidence_up" %in% names(density_data)) peak_values <- c(peak_values, density_data$confidence_up)
  finite_peak <- peak_values[is.finite(peak_values)]
  peak <- if (length(finite_peak)) max(finite_peak, 1e-12) else 1
  list(
    type = "meta_plot", title = title, yLabel = "Frequency", xDomain = c(0, 1),
    yDomain = c(-0.12 * peak, 1.05 * peak), showCI = isTRUE(show_ci),
    components = .rnameta_meta_component_payload(component_width, peak, head_or_tail, tx_promoter_length, tx_tail_length),
    series = lapply(seq_along(groups), function(i) {
      current <- density_data[as.character(get(group_col)) == groups[[i]]]
      list(name = groups[[i]], color = colors[((i - 1L) %% length(colors)) + 1L], values = lapply(seq_len(nrow(current)), function(j) list(x = current$x[[j]], density = current$density[[j]], confidenceDown = if ("confidence_down" %in% names(current)) current$confidence_down[[j]] else NULL, confidenceUp = if ("confidence_up" %in% names(current)) current$confidence_up[[j]] else NULL)))
    })
  )
}

#' Export a Meta Plot widget
#' @param plot A widget returned by RNAmeta_meta_plot().
#' @param file PNG or PDF path.
#' @param tx_type Transcript type when exporting an RNAmeta_meta_plots collection.
#' @param width Output width.
#' @param height Optional exact height.
#' @param aspect_ratio Width/height ratio when height is NULL.
#' @param scale PNG resolution multiplier.
#' @return Output path invisibly.
#' @export
RNAmeta_meta_plot_export <- function(plot, file, tx_type = NULL, width = 1180, height = NULL, aspect_ratio = c(1180, 620), scale = 2) {
  if (inherits(plot, "RNAmeta_meta_plots")) {
    if (is.null(tx_type) || length(tx_type) != 1L || !tx_type %in% names(plot)) {
      stop("tx_type must select one of: ", paste(names(plot), collapse = ", "))
    }
    plot <- plot[[tx_type]]
  }
  if (!inherits(plot, "htmlwidget") || is.null(plot$x$config) || !identical(plot$x$config$type, "meta_plot")) stop("plot must be returned by RNAmeta_meta_plot().")
  .rnameta_chart_export(plot$x$config, file, width, height, aspect_ratio, scale)
}

.bed_file_test_for_gene_meta <- function(bed_file = NULL, chrom_info = NULL, set_group_name_length = NULL) {
  if (length(bed_file) != set_group_name_length) stop("The number of BED files must match the number of group names.")
  paths <- as.character(unlist(bed_file, use.names = FALSE))
  if (!length(paths) || any(!file.exists(paths))) stop("One or more BED files do not exist.")
  if (any(tolower(tools::file_ext(paths)) != "bed")) stop("All Meta Plot inputs must be BED files.")
  copies <- vapply(seq_along(paths), function(i) {
    target <- tempfile(pattern = paste0("rnameta_meta_", i, "_"), fileext = ".bed")
    if (!file.copy(paths[[i]], target, overwrite = TRUE)) stop("Failed to prepare BED input: ", paths[[i]])
    target
  }, character(1))
  copies
}

.rnameta_meta_cache <- new.env(parent = emptyenv())
.rnameta_meta_cached <- function(key, builder) {
  if (.rnameta_cache_has(.rnameta_meta_cache, key)) return(.rnameta_cache_get(.rnameta_meta_cache, key))
  value <- builder(); .rnameta_cache_set(.rnameta_meta_cache, key, value); value
}
