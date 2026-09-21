.rnameta_xml_escape <- function(value) {
  value <- as.character(value)
  value <- gsub("&", "&amp;", value, fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub('"', "&quot;", value, fixed = TRUE)
  gsub("'", "&apos;", value, fixed = TRUE)
}

.rnameta_svg_band <- function(count, width, padding = 0.18) {
  if (count < 1L) return(list(start = numeric(), width = 0, step = 0))
  step <- width / max(1, count - padding + 2 * padding)
  list(start = step * padding + (seq_len(count) - 1) * step,
       width = step * (1 - padding), step = step)
}

.rnameta_resolve_export_size <- function(width, height = NULL,
                                         aspect_ratio = c(1180, 560)) {
  width <- as.numeric(width)
  if (length(width) != 1L || !is.finite(width) || width < 400) {
    stop("width must be a single number of at least 400 pixels.")
  }
  if (is.null(height)) {
    ratio <- as.numeric(aspect_ratio)
    if (length(ratio) == 2L) ratio <- ratio[[1L]] / ratio[[2L]]
    if (length(ratio) != 1L || !is.finite(ratio) || ratio <= 0) {
      stop("aspect_ratio must be a positive number or c(width, height).")
    }
    height <- width / ratio
  } else {
    height <- as.numeric(height)
  }
  if (length(height) != 1L || !is.finite(height) || height < 300) {
    stop("height must be a single number of at least 300 pixels.")
  }
  c(width = round(width), height = round(height))
}

.rnameta_gene_type_svg <- function(counts, show_labels = TRUE,
                                    width = 1180, height = 560) {
  config <- .rnameta_gene_type_payload(counts, show_labels, animate = FALSE)
  width <- as.numeric(width)
  height <- as.numeric(height)
  scale_x <- width / 1180
  scale_y <- height / 560
  visual_scale <- max(0.65, min(scale_x, scale_y))
  fs <- function(value) round(value * visual_scale, 3)
  sw <- function(value) round(value * visual_scale, 3)
  margin <- c(top = 92 * scale_y, right = 36 * scale_x,
              bottom = 108 * scale_y, left = 82 * scale_x)
  plot_width <- width - margin[["left"]] - margin[["right"]]
  plot_height <- height - margin[["top"]] - margin[["bottom"]]
  categories <- unlist(config$categories, use.names = FALSE)
  series <- config$series
  outer <- .rnameta_svg_band(length(categories), plot_width, 0.18)
  inner <- .rnameta_svg_band(length(series), outer$width, 0.18)
  y_ticks <- as.numeric(unlist(config$yTicks, use.names = FALSE))
  y_max <- max(c(as.numeric(config$yDomain[[2L]]), y_ticks), na.rm = TRUE)
  y_pos <- function(value) margin[["top"]] + plot_height * (1 - value / y_max)
  baseline_y <- margin[["top"]] + plot_height
  axis_color <- "#41503c"
  esc <- .rnameta_xml_escape
  pieces <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s" role="img" aria-label="Gene Type chart">', width, height, width, height),
    sprintf('<rect width="%s" height="%s" fill="#ffffff"/>', width, height),
    sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="#ffffff"/>', margin[["left"]], margin[["top"]], plot_width, plot_height, fs(18))
  )

  for (tick in y_ticks) {
    y <- y_pos(tick)
    label <- if (tick >= 1000) sprintf("%.1fk", tick / 1000) else sprintf("%.0f", tick)
    pieces <- c(
      pieces,
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#dfe4da" stroke-width="%s"/>', margin[["left"]], margin[["left"]] + plot_width, y, y, sw(1)),
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="%s" stroke-width="%s"/>', margin[["left"]] - 6 * visual_scale, margin[["left"]], y, y, axis_color, sw(1.1)),
      sprintf('<text x="%s" y="%s" text-anchor="end" dominant-baseline="middle" fill="%s" font-family="sans-serif" font-size="%s" font-weight="600">%s</text>', margin[["left"]] - 10 * visual_scale, y, axis_color, fs(13), esc(label))
    )
  }

  pieces <- c(
    pieces,
    sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="%s" stroke-width="%s"/>', margin[["left"]], margin[["left"]], margin[["top"]], baseline_y, axis_color, sw(1.2)),
    sprintf('<text x="%s" y="%s" fill="#22301f" font-family="sans-serif" font-size="%s" font-weight="800">Gene Type</text>', margin[["left"]], 80 * scale_y, fs(15)),
    sprintf('<text transform="translate(%s %s) rotate(-90)" fill="%s" font-family="sans-serif" font-size="%s" font-weight="700">Frequency</text>', 26 * scale_x, margin[["top"]] + plot_height / 2, axis_color, fs(15))
  )

  for (category_index in seq_along(categories)) {
    category_x <- margin[["left"]] + outer$start[[category_index]]
    pieces <- c(
      pieces,
      sprintf('<text x="%s" y="%s" text-anchor="middle" fill="%s" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', category_x + outer$width / 2, baseline_y + 28 * visual_scale, axis_color, fs(13), esc(categories[[category_index]]))
    )
    for (series_index in seq_along(series)) {
      value <- as.numeric(unlist(series[[series_index]]$values, use.names = FALSE)[[category_index]])
      bar_x <- category_x + inner$start[[series_index]]
      bar_y <- y_pos(value)
      bar_height <- baseline_y - bar_y
      pieces <- c(
        pieces,
        sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="%s" stroke="#ffffff" stroke-width="%s"/>', bar_x, bar_y, inner$width, max(0, bar_height), esc(series[[series_index]]$color), sw(1.2))
      )
      if (isTRUE(show_labels) && value > 0) {
        pieces <- c(
          pieces,
          sprintf('<text x="%s" y="%s" text-anchor="middle" fill="%s" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', bar_x + inner$width / 2, bar_y - 8 * visual_scale, axis_color, fs(12), esc(format(value, scientific = FALSE, trim = TRUE)))
        )
      }
    }
  }

  # Draw the baseline after the bars, matching the Viewer/D3 layer order.
  pieces <- c(
    pieces,
    sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="%s" stroke-width="%s"/>', margin[["left"]], margin[["left"]] + plot_width, baseline_y, baseline_y, axis_color, sw(1.2))
  )

  names <- vapply(series, `[[`, character(1), "name")
  item_widths <- pmax(110 * visual_scale,
                      nchar(names) * 8.3 * visual_scale + 44 * visual_scale)
  gap <- 24 * visual_scale
  row_gap <- 28 * visual_scale
  rows <- list()
  current <- integer()
  current_width <- 0
  for (i in seq_along(series)) {
    projected <- if (length(current)) current_width + gap + item_widths[[i]] else item_widths[[i]]
    if (length(current) && projected > plot_width) {
      rows[[length(rows) + 1L]] <- current
      current <- integer()
      current_width <- 0
    }
    current <- c(current, i)
    current_width <- if (length(current) == 1L) item_widths[[i]] else current_width + gap + item_widths[[i]]
  }
  if (length(current)) rows[[length(rows) + 1L]] <- current

  legend_y <- height - 72 * scale_y
  for (row_index in seq_along(rows)) {
    indexes <- rows[[row_index]]
    row_width <- sum(item_widths[indexes]) + gap * max(0, length(indexes) - 1L)
    cursor <- margin[["left"]] + (plot_width - row_width) / 2
    for (i in indexes) {
      y <- legend_y + (row_index - 1L) * row_gap
      pieces <- c(
        pieces,
        sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s"/>', cursor, y, fs(16), fs(16), fs(4), esc(series[[i]]$color)),
        sprintf('<text x="%s" y="%s" dominant-baseline="middle" fill="%s" font-family="sans-serif" font-size="%s" font-weight="600">%s</text>', cursor + 24 * visual_scale, y + 8 * visual_scale, axis_color, fs(13), esc(series[[i]]$name))
      )
      cursor <- cursor + item_widths[[i]] + gap
    }
  }

  paste0(paste0(pieces, collapse = ""), "</svg>")
}

#' Export a Gene Type plot to PNG or PDF
#'
#' This exporter does not require RStudio Viewer, a browser, X11 or an active
#' display, and is suitable for headless Linux servers.
#' @param data A counts table returned by `RNAmeta_gene_type_data()`, or a legacy
#'   preliminary analysis list.
#' @param file Output filename ending in `.png` or `.pdf`.
#' @param set_group_name Optional replacement labels in current sample order.
#' @param show_labels Display count labels above bars.
#' @param width Output width in pixels before PNG scaling.
#' @param height Optional explicit output height. When `NULL`, height is derived
#'   from `width` and `aspect_ratio`.
#' @param aspect_ratio Positive width/height ratio, or `c(width, height)`.
#'   Ignored when `height` is supplied.
#' @param scale PNG resolution multiplier; ignored for PDF.
#' @return Invisibly returns the normalized output path.
#' @export
RNAmeta_gene_type_export <- function(data, file, set_group_name = NULL,
                                     show_labels = TRUE, width = 1180,
                                     height = NULL,
                                     aspect_ratio = c(1180, 560),
                                     scale = 2) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("file must be a nonempty path.")
  }
  extension <- tolower(tools::file_ext(file))
  if (!extension %in% c("png", "pdf")) stop("file must end in .png or .pdf.")
  parent <- dirname(file)
  if (!dir.exists(parent)) stop("Output directory does not exist: ", parent)
  size <- .rnameta_resolve_export_size(width, height, aspect_ratio)
  counts <- .rnameta_gene_type_plot_table(data, set_group_name)
  svg <- charToRaw(.rnameta_gene_type_svg(
    counts, show_labels, width = size[["width"]], height = size[["height"]]
  ))
  if (identical(extension, "png")) {
    scale <- as.numeric(scale)
    if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
      stop("scale must be a positive number.")
    }
    rsvg::rsvg_png(
      svg, file = file,
      width = round(size[["width"]] * scale),
      height = round(size[["height"]] * scale)
    )
  } else {
    rsvg::rsvg_pdf(svg, file = file)
  }
  invisible(normalizePath(file, winslash = "/", mustWork = TRUE))
}
