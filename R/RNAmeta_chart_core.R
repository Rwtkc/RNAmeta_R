.rnameta_chart_palette <- function() {
  c("#de7a5b", "#3f3f5b", "#859b7a", "#c9886d", "#b9c78f", "#6676a0")
}

.rnameta_preliminary_component <- function(data, index, required = character()) {
  if (is.data.frame(data)) {
    result <- data.table::as.data.table(data.table::copy(data))
  } else {
    shared <- attr(data, "peak_annotation", exact = TRUE)
    component_name <- c("distribution", "exon", "gene", "detail", "transcript")[[index]]
    if (!is.null(shared) && !is.null(shared$components[[component_name]])) {
      result <- data.table::rbindlist(shared$components[[component_name]], use.names = TRUE, fill = TRUE)
    } else {
      if (!is.list(data) || !length(data) || !is.list(data[[1]]) || !length(data[[1]])) {
        stop("Provide preliminary analysis data or a prepared data frame.")
      }
      result <- data.table::rbindlist(
        lapply(data[[1]], function(entry) {
          if (!is.list(entry) || length(entry) < index) return(NULL)
          entry[[index]]
        }),
        use.names = TRUE, fill = TRUE
      )
    }
  }
  missing <- setdiff(required, names(result))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
  result
}

.rnameta_clean_sample_names <- function(data) {
  if (!"Sample" %in% names(data)) stop("Missing required column: Sample")
  data[, Sample := sub(".*_", "", as.character(Sample))]
  data
}

.rnameta_rename_samples <- function(data, set_group_name = NULL) {
  samples <- unique(as.character(data$Sample))
  if (is.null(set_group_name)) {
    data[, Sample := factor(as.character(Sample), levels = samples)]
    return(data)
  }
  if (length(set_group_name) != length(samples)) {
    stop(
      "The number of group names does not match the number of samples. ",
      "Expected ", length(samples), ", received ", length(set_group_name), "."
    )
  }
  labels <- stats::setNames(as.character(set_group_name), samples)
  renamed <- unname(labels[as.character(data$Sample)])
  data[, Sample := factor(renamed, levels = as.character(set_group_name))]
  data
}

.rnameta_chart_widget <- function(config, width = NULL, height = NULL,
                                  element_id = NULL) {
  htmlwidgets::createWidget(
    name = "rnametaD3",
    x = list(config = config),
    width = width,
    height = height,
    package = "RNAmeta",
    elementId = element_id
  )
}

.rnameta_boxplot_summary <- function(values, scale_transform = "linear") {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values) & values > if (identical(scale_transform, "log2")) 0 else -Inf]
  if (!length(values)) return(NULL)
  display <- if (identical(scale_transform, "log2")) log2(values) else values
  quartiles <- as.numeric(stats::quantile(display, c(0.25, 0.5, 0.75), names = FALSE, type = 7))
  iqr <- quartiles[[3]] - quartiles[[1]]
  lower_limit <- quartiles[[1]] - 1.5 * iqr
  upper_limit <- quartiles[[3]] + 1.5 * iqr
  inside <- display[display >= lower_limit & display <= upper_limit]
  list(
    count = length(display),
    min = min(display), q1 = quartiles[[1]], median = quartiles[[2]], q3 = quartiles[[3]], max = max(display),
    lowerWhisker = if (length(inside)) min(inside) else quartiles[[1]],
    upperWhisker = if (length(inside)) max(inside) else quartiles[[3]],
    outliers = display[display < lower_limit | display > upper_limit]
  )
}

.rnameta_chart_export <- function(config, file, width = 1180, height = NULL,
                                  aspect_ratio = c(1180, 560), scale = 2) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("file must be a nonempty path.")
  }
  extension <- tolower(tools::file_ext(file))
  if (!extension %in% c("png", "pdf")) stop("file must end in .png or .pdf.")
  if (!dir.exists(dirname(file))) stop("Output directory does not exist: ", dirname(file))
  size <- .rnameta_resolve_export_size(width, height, aspect_ratio)
  svg <- charToRaw(.rnameta_chart_svg(config, size[["width"]], size[["height"]]))
  if (identical(extension, "png")) {
    scale <- as.numeric(scale)
    if (length(scale) != 1L || !is.finite(scale) || scale <= 0) stop("scale must be a positive number.")
    rsvg::rsvg_png(svg, file = file, width = round(size[["width"]] * scale), height = round(size[["height"]] * scale))
  } else {
    rsvg::rsvg_pdf(svg, file = file)
  }
  invisible(normalizePath(file, winslash = "/", mustWork = TRUE))
}
