.rnameta_gene_type_plot_table <- function(data, set_group_name = NULL) {
  if (is.data.frame(data) && all(c("V1", "N", "Sample") %in% names(data))) {
    counts <- data.table::as.data.table(data.table::copy(data))
  } else {
    shared <- attr(data, "peak_annotation", exact = TRUE)
    if (!is.null(shared) && length(shared$gene_type_counts)) {
      counts <- data.table::rbindlist(shared$gene_type_counts, use.names = TRUE, fill = TRUE)
    } else {
      if (!is.list(data) || length(data) < 3L) {
        stop("Provide preliminary analysis data or RNAmeta_gene_type_data() output.")
      }
      peak_gr <- data[[2L]]
      feature_ranges <- .rnameta_gene_type_features(data[[3L]])
      rows <- lapply(seq_along(peak_gr), function(i) {
        sample_name <- names(peak_gr)[[i]]
        if (is.null(sample_name) || is.na(sample_name) || !nzchar(sample_name)) sample_name <- paste0("Group", i)
        .rnameta_gene_type_counts(peak_gr[[i]], feature_ranges, sample_name)
      })
      counts <- data.table::rbindlist(rows, use.names = TRUE)
    }
  }
  counts[, V1 := as.character(V1)]
  counts[, Sample := as.character(Sample)]
  counts[, N := as.numeric(N)]
  counts <- counts[!is.na(V1) & nzchar(V1) & is.finite(N)]
  if (!nrow(counts)) stop("No overlapping transcript biotypes were found for these samples.")
  samples <- unique(counts$Sample)
  if (!is.null(set_group_name)) {
    if (!is.character(set_group_name) || length(set_group_name) != length(samples) ||
        anyNA(set_group_name) || any(!nzchar(set_group_name)) || anyDuplicated(set_group_name)) {
      stop("set_group_name must contain one distinct, nonempty label per sample.")
    }
    counts[, Sample := set_group_name[match(Sample, samples)]]
  }
  counts
}

.rnameta_gene_type_payload <- function(counts, show_labels = FALSE, animate = TRUE) {
  categories <- unique(as.character(counts$V1))
  samples <- unique(as.character(counts$Sample))
  palette <- c("#de7a61", "#3f3f5b", "#81b29d", "#d69779", "#b7c88b")
  max_value <- max(counts$N, na.rm = TRUE)
  y_max <- if (is.finite(max_value) && max_value > 0) max_value * 1.12 else 1
  y_ticks <- pretty(c(0, y_max), n = 5)
  y_ticks <- y_ticks[y_ticks >= 0 & y_ticks <= max(y_ticks[y_ticks <= y_max], y_max)]

  list(
    type = "grouped_bar", title = "Gene Type", xLabel = "", yLabel = "Frequency",
    categories = unname(categories), hoverMinHeight = 18, showLabels = isTRUE(show_labels),
    animate = isTRUE(animate), yTicks = unname(as.numeric(y_ticks)), yDomain = c(0, y_max),
    series = lapply(seq_along(samples), function(i) {
      sample <- samples[[i]]
      sample_data <- counts[Sample == sample]
      values <- sample_data$N[match(categories, sample_data$V1)]
      values[is.na(values)] <- 0
      source_name <- if ("Source_File" %in% names(sample_data) && nrow(sample_data)) {
        as.character(sample_data$Source_File[[1L]])
      } else sample
      list(name = sample, originalName = source_name,
           color = palette[((i - 1L) %% length(palette)) + 1L], values = unname(as.numeric(values)))
    })
  )
}

.rnameta_gene_type_widget <- function(data, set_group_name = NULL, show_labels = TRUE,
                                      animate = TRUE, width = NULL, height = 560,
                                      elementId = NULL) {
  counts <- .rnameta_gene_type_plot_table(data, set_group_name)
  htmlwidgets::createWidget(
    name = "rnametaD3",
    x = list(config = .rnameta_gene_type_payload(counts, show_labels, animate)),
    width = width, height = height, package = "RNAmeta", elementId = elementId,
    sizingPolicy = htmlwidgets::sizingPolicy(
      viewer.defaultWidth = 1180, viewer.defaultHeight = 600,
      browser.fill = FALSE, padding = 0
    )
  )
}
