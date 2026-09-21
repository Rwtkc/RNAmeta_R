.rnameta_chart_ticks <- function(minimum, maximum, count = 5L) {
  if (!is.finite(minimum) || !is.finite(maximum)) return(c(0, 1))
  if (minimum == maximum) {
    padding <- if (minimum == 0) 1 else abs(minimum) * 0.1
    minimum <- minimum - padding
    maximum <- maximum + padding
  }
  pretty(c(minimum, maximum), n = count)
}

.rnameta_chart_legend_svg <- function(series, left, top, width, visual_scale) {
  if (!length(series)) return(character())
  esc <- .rnameta_xml_escape
  names <- vapply(series, function(item) as.character(item$name), character(1))
  colors <- vapply(series, function(item) as.character(item$color), character(1))
  item_widths <- pmax(100 * visual_scale, nchar(names) * 8 * visual_scale + 42 * visual_scale)
  gap <- 22 * visual_scale
  total <- sum(item_widths) + gap * max(0, length(series) - 1L)
  cursor <- left + max(0, (width - total) / 2)
  pieces <- character()
  for (i in seq_along(series)) {
    pieces <- c(
      pieces,
      sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s"/>', cursor, top, 15 * visual_scale, 15 * visual_scale, 4 * visual_scale, esc(colors[[i]])),
      sprintf('<text x="%s" y="%s" dominant-baseline="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="600">%s</text>', cursor + 23 * visual_scale, top + 7.5 * visual_scale, 13 * visual_scale, esc(names[[i]]))
    )
    cursor <- cursor + item_widths[[i]] + gap
  }
  pieces
}

.rnameta_bar_chart_svg <- function(config, width, height) {
  esc <- .rnameta_xml_escape
  sx <- width / 1180
  sy <- height / 560
  vs <- sqrt(sx * sy)
  margin <- c(top = 92 * sy, right = 36 * sx, bottom = 108 * sy, left = 82 * sx)
  plot_width <- width - margin[["left"]] - margin[["right"]]
  plot_height <- height - margin[["top"]] - margin[["bottom"]]
  baseline <- margin[["top"]] + plot_height
  categories <- unlist(config$categories, use.names = FALSE)
  series <- config$series
  stacked <- identical(config$type, "stacked_bar")
  all_values <- unlist(lapply(series, function(item) as.numeric(unlist(item$values, use.names = FALSE))))
  y_max <- if (length(config$yDomain) >= 2L) as.numeric(config$yDomain[[2]]) else max(c(1, all_values), na.rm = TRUE)
  if (stacked && is.null(config$yDomain)) y_max <- max(100, y_max)
  ticks <- if (length(config$yTicks)) as.numeric(unlist(config$yTicks)) else .rnameta_chart_ticks(0, y_max)
  y_max <- max(c(y_max, ticks), na.rm = TRUE)
  y <- function(value) margin[["top"]] + plot_height - (as.numeric(value) / y_max) * plot_height
  outer <- .rnameta_svg_band(length(categories), plot_width, 0.18)
  inner <- .rnameta_svg_band(max(1, length(series)), outer$width, 0.18)
  pieces <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s">', width, height, width, height),
    sprintf('<rect width="%s" height="%s" fill="#ffffff"/>', width, height),
    sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="#ffffff"/>', margin[["left"]], margin[["top"]], plot_width, plot_height, 18 * vs)
  )
  for (tick in ticks) {
    yy <- y(tick)
    label <- if (stacked) paste0(round(tick), "%") else if (tick >= 1000) paste0(format(round(tick / 1000, 1), trim = TRUE), "k") else format(round(tick), trim = TRUE)
    pieces <- c(pieces,
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#dfe4da" stroke-width="%s"/>', margin[["left"]], margin[["left"]] + plot_width, yy, yy, vs),
      sprintf('<text x="%s" y="%s" text-anchor="end" dominant-baseline="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="600">%s</text>', margin[["left"]] - 10 * vs, yy, 13 * vs, esc(label))
    )
  }
  pieces <- c(pieces,
    sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', margin[["left"]], margin[["left"]], margin[["top"]], baseline, 1.2 * vs),
    sprintf('<text x="%s" y="%s" fill="#22301f" font-family="sans-serif" font-size="%s" font-weight="800">%s</text>', margin[["left"]], 72 * sy, 15 * vs, esc(config$title %||% "")),
    sprintf('<text transform="translate(%s %s) rotate(-90)" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', 26 * sx, margin[["top"]] + plot_height / 2, 15 * vs, esc(config$yLabel %||% ""))
  )
  for (ci in seq_along(categories)) {
    cx <- margin[["left"]] + outer$start[[ci]]
    pieces <- c(pieces, sprintf('<text x="%s" y="%s" text-anchor="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', cx + outer$width / 2, baseline + 28 * vs, 13 * vs, esc(categories[[ci]])))
    cumulative <- 0
    for (si in seq_along(series)) {
      value <- as.numeric(unlist(series[[si]]$values, use.names = FALSE)[[ci]])
      if (!is.finite(value)) value <- 0
      if (stacked) {
        x0 <- cx + outer$width * 0.22
        bar_width <- outer$width * 0.56
        y0 <- y(cumulative + value)
        bar_height <- y(cumulative) - y(cumulative + value)
        cumulative <- cumulative + value
      } else {
        x0 <- cx + inner$start[[si]]
        bar_width <- inner$width
        y0 <- y(value)
        bar_height <- baseline - y0
      }
      pieces <- c(pieces, sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="%s" stroke="#ffffff" stroke-width="%s"/>', x0, y0, bar_width, max(0, bar_height), esc(series[[si]]$color), 1.2 * vs))
      if (!stacked && isTRUE(config$showLabels) && value > 0) pieces <- c(pieces, sprintf('<text x="%s" y="%s" text-anchor="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', x0 + bar_width / 2, y0 - 8 * vs, 12 * vs, esc(format(value, scientific = FALSE, trim = TRUE))))
    }
  }
  pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', margin[["left"]], margin[["left"]] + plot_width, baseline, baseline, 1.2 * vs))
  pieces <- c(pieces, .rnameta_chart_legend_svg(series, margin[["left"]], height - 72 * sy, plot_width, vs), "</svg>")
  paste0(pieces, collapse = "")
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

.rnameta_box_panel_svg <- function(groups, x, y, width, height, scale_transform, title = NULL, show_y = TRUE, visual_scale = 1) {
  esc <- .rnameta_xml_escape
  summaries <- lapply(groups, function(group) .rnameta_boxplot_summary(unlist(group$values, use.names = FALSE), scale_transform))
  keep <- !vapply(summaries, is.null, logical(1))
  groups <- groups[keep]; summaries <- summaries[keep]
  if (!length(groups)) return(character())
  values <- unlist(lapply(summaries, function(s) c(s$lowerWhisker, s$upperWhisker, s$outliers)))
  ticks <- .rnameta_chart_ticks(min(values, na.rm = TRUE), max(values, na.rm = TRUE))
  ymin <- min(ticks); ymax <- max(ticks)
  ypos <- function(value) y + height - (value - ymin) / (ymax - ymin) * height
  band <- .rnameta_svg_band(length(groups), width, 0.24)
  pieces <- sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="#ffffff"/>', x, y, width, height, 16 * visual_scale)
  for (tick in ticks) {
    yy <- ypos(tick)
    label <- if (identical(scale_transform, "log2")) format(round(2^tick), big.mark = ",", scientific = FALSE, trim = TRUE) else format(round(tick, 1), trim = TRUE)
    pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#dfe4da" stroke-width="%s"/>', x, x + width, yy, yy, visual_scale))
    if (show_y) pieces <- c(
      pieces,
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', x - 6 * visual_scale, x, yy, yy, 1.1 * visual_scale),
      sprintf('<text x="%s" y="%s" text-anchor="end" dominant-baseline="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="600">%s</text>', x - 10 * visual_scale, yy, 13 * visual_scale, esc(label))
    )
  }
  if (show_y) pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', x, x, y, y + height, 1.2 * visual_scale))
  pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', x, x + width, y + height, y + height, 1.2 * visual_scale))
  for (i in seq_along(groups)) {
    group <- groups[[i]]; s <- summaries[[i]]
    cx <- x + band$start[[i]] + band$width / 2
    bw <- band$width * 0.62
    pieces <- c(pieces,
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', cx, cx, ypos(s$lowerWhisker), ypos(s$upperWhisker), 1.5 * visual_scale),
      sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="%s" fill-opacity="0.86" stroke="#ffffff" stroke-width="%s"/>', cx - bw / 2, ypos(s$q3), bw, max(1, ypos(s$q1) - ypos(s$q3)), esc(group$color), 1.2 * visual_scale),
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#ffffff" stroke-width="%s"/>', cx - bw / 2, cx + bw / 2, ypos(s$median), ypos(s$median), 2 * visual_scale),
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', cx - bw * .3, cx + bw * .3, ypos(s$lowerWhisker), ypos(s$lowerWhisker), 1.5 * visual_scale),
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', cx - bw * .3, cx + bw * .3, ypos(s$upperWhisker), ypos(s$upperWhisker), 1.5 * visual_scale),
      sprintf('<text x="%s" y="%s" text-anchor="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', cx, y + height + 27 * visual_scale, 12 * visual_scale, esc(group$name))
    )
    if (length(s$outliers)) for (outlier in s$outliers) pieces <- c(pieces, sprintf('<circle cx="%s" cy="%s" r="%s" fill="%s" stroke="#111111" stroke-width="%s"/>', cx, ypos(outlier), 3 * visual_scale, esc(group$color), visual_scale))
  }
  if (!is.null(title)) pieces <- c(pieces, sprintf('<text x="%s" y="%s" text-anchor="middle" fill="#22301f" font-family="sans-serif" font-size="%s" font-weight="800">%s</text>', x + width / 2, y - 16 * visual_scale, 13 * visual_scale, esc(title)))
  pieces
}

.rnameta_box_chart_svg <- function(config, width, height) {
  esc <- .rnameta_xml_escape
  sx <- width / 1180; sy <- height / 560; vs <- sqrt(sx * sy)
  margin <- c(top = 125 * sy, right = 36 * sx, bottom = 96 * sy, left = 82 * sx)
  plot_width <- width - margin[["left"]] - margin[["right"]]
  plot_height <- height - margin[["top"]] - margin[["bottom"]]
  pieces <- c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s">', width, height, width, height), sprintf('<rect width="%s" height="%s" fill="#ffffff"/>', width, height), sprintf('<text x="%s" y="%s" fill="#22301f" font-family="sans-serif" font-size="%s" font-weight="800">%s</text>', margin[["left"]], 72 * sy, 15 * vs, esc(config$title %||% "")), sprintf('<text transform="translate(%s %s) rotate(-90)" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', 26 * sx, margin[["top"]] + plot_height / 2, 15 * vs, esc(config$yLabel %||% "")))
  if (identical(config$type, "faceted_boxplot")) {
    facets <- config$facets
    gap <- 24 * sx
    panel_width <- (plot_width - gap * (length(facets) - 1L)) / length(facets)
    for (i in seq_along(facets)) pieces <- c(pieces, .rnameta_box_panel_svg(facets[[i]]$groups, margin[["left"]] + (i - 1L) * (panel_width + gap), margin[["top"]], panel_width, plot_height, config$scaleTransform %||% "linear", facets[[i]]$name, i == 1L, vs))
  } else {
    pieces <- c(pieces, .rnameta_box_panel_svg(config$groups, margin[["left"]], margin[["top"]], plot_width, plot_height, config$scaleTransform %||% "linear", NULL, TRUE, vs))
  }
  paste0(c(pieces, "</svg>"), collapse = "")
}

.rnameta_gene_matrix_svg <- function(config, width, height) {
  esc <- .rnameta_xml_escape
  labels <- as.character(unlist(config$sampleLabels, use.names = FALSE))
  rows <- config$intersections
  if (!length(labels) || !length(rows)) stop("Gene Matrix has no intersections to export.")
  sizes <- vapply(rows, function(row) as.numeric(row$size), numeric(1))
  set_counts <- vapply(rows, function(row) length(row$sets), integer(1))
  ids <- vapply(rows, function(row) paste(unlist(row$sets, use.names = FALSE), collapse = "__"), character(1))
  rows <- rows[order(-sizes, set_counts, ids)]
  blocks <- split(rows, ceiling(seq_along(rows) / 8L))
  base_width <- 1120; base_height <- 700
  margin <- c(top = 52, right = 40, bottom = 54, left = 40)
  top_bar_height <- 250; row_height <- 40; column_width <- 132; block_gap <- 68
  block_height <- top_bar_height + 52 + length(labels) * row_height
  legend_columns <- if (length(labels) >= 4L) 2L else min(3L, length(labels))
  legend_item_width <- 320; legend_swatch_width <- 58; legend_offset_x <- 72
  legend_total_width <- legend_columns * legend_item_width
  max_block_columns <- max(vapply(blocks, length, integer(1)))
  max_block_width <- max(420, max_block_columns * column_width)
  viewport_width <- max(base_width, max_block_width + margin[["left"]] + 76 + margin[["right"]], legend_total_width + margin[["left"]] + margin[["right"]])
  legend_rows <- ceiling(length(labels) / legend_columns); legend_height <- legend_rows * 68
  legend_top <- margin[["top"]] + length(blocks) * block_height + max(0, length(blocks) - 1L) * block_gap + 52
  viewport_height <- max(base_height, legend_top + legend_height + margin[["bottom"]])
  colors <- config$colors %||% stats::setNames(as.list(c("#d94841", "#2563eb", "#f59e0b", "#0f9d78", "#7c3aed")[seq_along(labels)]), labels)
  sample_sizes <- stats::setNames(rep(0, length(labels)), labels)
  for (row in rows) { sets <- unlist(row$sets, use.names = FALSE); if (length(sets) == 1L && sets[[1]] %in% labels) sample_sizes[[sets[[1]]]] <- as.numeric(row$size) }
  global_max <- max(c(1, vapply(rows, function(row) as.numeric(row$size), numeric(1))), na.rm = TRUE)
  y_ticks <- pretty(c(0, global_max), n = 5); y_ticks <- y_ticks[y_ticks >= 0]; y_max <- max(y_ticks, global_max)
  pieces <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s" preserveAspectRatio="xMidYMid meet">', width, height, viewport_width, viewport_height),
    sprintf('<rect width="%s" height="%s" fill="#ffffff"/>', viewport_width, viewport_height)
  )
  for (bi in seq_along(blocks)) {
    block <- blocks[[bi]]; block_top <- margin[["top"]] + (bi - 1L) * (block_height + block_gap)
    bar_bottom <- block_top + top_bar_height; matrix_top <- bar_bottom + 52
    block_width <- max(420, length(block) * column_width); left <- (viewport_width - block_width) / 2
    padding_inner <- .28; padding_outer <- .08
    step <- block_width / max(1, length(block) - padding_inner + 2 * padding_outer)
    band_width <- step * (1 - padding_inner); start <- left + step * padding_outer
    ypos <- function(value) bar_bottom - as.numeric(value) / y_max * top_bar_height
    pieces <- c(pieces,
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#74806f" stroke-width="1.2"/>', left, left, block_top, bar_bottom),
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#74806f" stroke-width="1.2"/>', left, left + block_width, bar_bottom, bar_bottom))
    for (tick in y_ticks[y_ticks <= y_max]) pieces <- c(pieces,
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#74806f" stroke-width="1.2"/>', left - 7, left, ypos(tick), ypos(tick)),
      sprintf('<text x="%s" y="%s" text-anchor="end" dominant-baseline="middle" fill="#4f5b4a" font-family="sans-serif" font-size="14" font-weight="700">%s</text>', left - 12, ypos(tick), esc(format(tick, trim = TRUE, scientific = FALSE, big.mark = ","))))
    if (bi == 1L) pieces <- c(pieces, sprintf('<text x="%s" y="%s" fill="#253320" font-family="sans-serif" font-size="15" font-weight="800">Intersection Size</text>', left + 6, block_top - 18))
    for (ri in seq_along(block)) {
      row <- block[[ri]]; sets <- as.character(unlist(row$sets, use.names = FALSE)); value <- as.numeric(row$size)
      x0 <- start + (ri - 1L) * step; xc <- x0 + band_width / 2
      pieces <- c(pieces,
        sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="#859b7a" fill-opacity="0.88"/>', x0, ypos(value), band_width, max(0, bar_bottom - ypos(value))),
        sprintf('<text x="%s" y="%s" text-anchor="middle" fill="#5b4a38" font-family="sans-serif" font-size="%s" font-weight="800">%s</text>', xc, ypos(value) - 10, if (band_width >= 22) 13 else 11, esc(format(value, trim = TRUE, scientific = FALSE, big.mark = ","))))
      active <- which(labels %in% sets)
      active_y <- matrix_top + (active - 1L) * row_height + row_height / 2
      if (length(active_y) > 1L) pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#7a6655" stroke-width="3" stroke-linecap="round"/>', xc, xc, min(active_y), max(active_y)))
      for (li in seq_along(labels)) {
        cy <- matrix_top + (li - 1L) * row_height + row_height / 2; on <- li %in% active
        fill <- if (on) as.character(colors[[labels[[li]]]] %||% "#859b7a") else "#d9ded3"
        pieces <- c(pieces, sprintf('<circle cx="%s" cy="%s" r="7.5" fill="%s" stroke="%s" stroke-width="%s"/>', xc, cy, esc(fill), if (on) "#ffffff" else "#b8c2b1", if (on) 2.2 else 1.1))
      }
    }
  }
  legend_left <- (viewport_width - legend_total_width) / 2 + legend_offset_x
  for (i in seq_along(labels)) {
    column <- (i - 1L) %% legend_columns; row <- (i - 1L) %/% legend_columns
    item_left <- legend_left + column * legend_item_width; item_top <- legend_top + row * 68
    color <- as.character(colors[[labels[[i]]]] %||% "#859b7a")
    pieces <- c(pieces,
      sprintf('<rect x="%s" y="%s" width="%s" height="16" rx="5" fill="%s" fill-opacity="0.9"/>', item_left, item_top + 2, legend_swatch_width, esc(color)),
      sprintf('<text x="%s" y="%s" dominant-baseline="hanging" fill="#253320" font-family="sans-serif" font-size="13" font-weight="700">%s</text>', item_left + legend_swatch_width + 12, item_top + 4, esc(labels[[i]])),
      sprintf('<text x="%s" y="%s" dominant-baseline="hanging" fill="#5b4a38" font-family="sans-serif" font-size="13" font-weight="800">%s</text>', item_left + legend_swatch_width + 12, item_top + 30, esc(format(sample_sizes[[labels[[i]]]], trim = TRUE, scientific = FALSE, big.mark = ","))))
  }
  paste0(c(pieces, "</svg>"), collapse = "")
}
.rnameta_chart_svg <- function(config, width, height) {
  if (config$type %in% c("grouped_bar", "stacked_bar")) return(.rnameta_bar_chart_svg(config, width, height))
  if (config$type %in% c("boxplot", "faceted_boxplot")) return(.rnameta_box_chart_svg(config, width, height))
  if (identical(config$type, "gene_matrix")) return(.rnameta_gene_matrix_svg(config, width, height))
  if (identical(config$type, "site_profile")) return(.rnameta_site_profile_svg(config, width, height))
  if (identical(config$type, "meta_plot")) return(.rnameta_meta_plot_svg(config, width, height))

  stop("Unsupported chart type: ", config$type)
}


.rnameta_site_profile_svg <- function(config, width, height) {
  esc <- .rnameta_xml_escape
  panel <- config$panels[[1L]]
  sx <- width / 1180; sy <- height / 560; vs <- sqrt(sx * sy)
  left <- 112 * sx; right <- 128 * sx; top <- 64 * sy
  pieces <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s">', width, height, width, height),
    sprintf('<rect width="%s" height="%s" fill="#ffffff"/>', width, height),
    sprintf('<text x="%s" y="%s" fill="#22301f" font-family="sans-serif" font-size="%s" font-weight="800">%s</text>', left, 34 * sy, 15 * vs, esc(panel$title %||% ""))
  )
  if (identical(panel$type, "density")) {
    bottom <- 130 * sy; plot_width <- width - left - right; plot_height <- height - top - bottom
    xs <- unlist(lapply(panel$series, function(s) vapply(s$values, function(p) as.numeric(p$x), numeric(1))))
    ys <- unlist(lapply(panel$series, function(s) vapply(s$values, function(p) as.numeric(p$density), numeric(1))))
    xr <- as.numeric(unlist(panel$xDomain %||% range(xs, finite = TRUE)))
    yr <- as.numeric(unlist(panel$yDomain %||% c(0, max(c(ys, 1e-12), na.rm = TRUE) * 1.08)))
    if (diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
    if (diff(yr) == 0) yr <- yr + c(0, max(1, abs(yr[[1]]) * .1))
    xp <- function(value) left + (as.numeric(value) - xr[[1]]) / diff(xr) * plot_width
    yp <- function(value) top + plot_height - (as.numeric(value) - yr[[1]]) / diff(yr) * plot_height
    y_ticks <- as.numeric(unlist(panel$yTicks %||% pretty(yr, n = 5)))
    y_ticks <- y_ticks[y_ticks >= yr[[1]] & y_ticks <= yr[[2]]]
    x_ticks <- pretty(xr, n = 6); x_ticks <- x_ticks[x_ticks >= xr[[1]] & x_ticks <= xr[[2]]]
    format_tick <- function(value) format(signif(value, 4), scientific = abs(value) >= 1e5, trim = TRUE, big.mark = if (abs(value) >= 1000) "," else "")
    for (tick in y_ticks) pieces <- c(pieces,
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#dfe4da" stroke-width="%s"/>', left, left + plot_width, yp(tick), yp(tick), vs),
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', left - 6 * vs, left, yp(tick), yp(tick), 1.2 * vs),
      sprintf('<text x="%s" y="%s" text-anchor="end" dominant-baseline="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', left - 10 * vs, yp(tick), 12 * vs, esc(format_tick(tick))))
    pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', left, left, top, top + plot_height, 1.2 * vs))
    for (tick in x_ticks) pieces <- c(pieces, sprintf('<text x="%s" y="%s" text-anchor="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', xp(tick), top + plot_height + 25 * sy, 12 * vs, esc(format_tick(tick))))
    pieces <- c(pieces,
      sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', left, left + plot_width, top + plot_height, top + plot_height, 1.2 * vs),
      sprintf('<text x="%s" y="%s" text-anchor="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', left + plot_width / 2, top + plot_height + 52 * sy, 15 * vs, esc(panel$xLabel %||% "")),
      sprintf('<text transform="translate(%s %s) rotate(-90)" text-anchor="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', 34 * sx, top + plot_height / 2, 15 * vs, esc(panel$yLabel %||% "")))
    for (guide in unlist(panel$guideLines %||% list())) if (is.finite(guide) && guide >= xr[[1]] && guide <= xr[[2]]) pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s" stroke-dasharray="3 6"/>', xp(guide), xp(guide), top, top + plot_height, 1.2 * vs))
    for (series in panel$series) {
      points <- vapply(series$values, function(point) sprintf('%s,%s', xp(point$x), yp(point$density)), character(1))
      pieces <- c(pieces, sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="%s" stroke-linejoin="round" stroke-linecap="round"/>', paste(points, collapse = " "), esc(series$color), 3.6 * vs))
    }
    if (length(panel$series) > 1L) pieces <- c(pieces, .rnameta_chart_legend_svg(panel$series, left, height - 30 * sy, plot_width, vs))
  } else if (identical(panel$type, "heatmap")) {
    bottom <- 74 * sy; plot_width <- width - left - right; plot_height <- height - top - bottom
    matrix <- do.call(rbind, lapply(panel$values, unlist))
    if (length(matrix)) {
      nr <- nrow(matrix); nc <- ncol(matrix)
      maxv <- as.numeric(stats::quantile(matrix, panel$colorMaxQuantile %||% .92, na.rm = TRUE))
      if (!is.finite(maxv) || maxv <= 0) maxv <- max(matrix, 1, na.rm = TRUE)
      palette <- grDevices::colorRampPalette(c("#ffffff", "#f7cf7a", "#ee9b4b", "#da6a33", "#b23a2c", "#7f1d1d"))(128)
      rw <- plot_width / nc; rh <- plot_height / nr
      pieces <- c(pieces, sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="%s"/>', left, top, plot_width, plot_height, esc(panel$backgroundColor %||% "#ffffff")))
      for (r in seq_len(nr)) for (cc in seq_len(nc)) {
        idx <- min(128L, max(1L, floor(matrix[r, cc] / maxv * 127) + 1L))
        pieces <- c(pieces, sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="%s"/>', left + (cc - 1) * rw, top + (r - 1) * rh, rw + .2, rh + .2, palette[[idx]]))
      }
      y_ticks <- unique(c(1L, max(1L, round(nr / 2)), nr))
      y_pos <- function(value) top + (as.numeric(value) - .5) / nr * plot_height
      for (tick in y_ticks) pieces <- c(pieces,
        sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', left - 8 * vs, left, y_pos(tick), y_pos(tick), 1.2 * vs),
        sprintf('<text x="%s" y="%s" text-anchor="end" dominant-baseline="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', left - 12 * vs, y_pos(tick), 12 * vs, tick))
      pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', left, left, top, top + plot_height, 1.2 * vs))
      xr <- as.numeric(unlist(panel$xDomain %||% c(panel$xMin, panel$xMax)))
      base_ticks <- pretty(xr, n = 5); base_ticks <- base_ticks[base_ticks >= xr[[1]] & base_ticks <= xr[[2]]]
      base_step <- if (length(base_ticks) >= 2L) min(diff(base_ticks)) else diff(xr)
      endpoint_threshold <- max(1, abs(base_step) * .35)
      add_min <- !length(base_ticks) || min(abs(base_ticks - xr[[1]])) > endpoint_threshold
      add_max <- !length(base_ticks) || min(abs(base_ticks - xr[[2]])) > endpoint_threshold
      x_ticks <- unique(c(if (add_min) xr[[1]], base_ticks, if (add_max) xr[[2]]))
      xp <- function(value) left + (as.numeric(value) - xr[[1]]) / diff(xr) * plot_width
      for (tick in x_ticks) pieces <- c(pieces,
        sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', xp(tick), xp(tick), top + plot_height, top + plot_height + 8 * vs, 1.2 * vs),
        sprintf('<text x="%s" y="%s" text-anchor="middle" fill="#41503c" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', xp(tick), top + plot_height + 25 * sy, 12 * vs, esc(format(signif(tick, 5), trim = TRUE, big.mark = ","))))
      pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#41503c" stroke-width="%s"/>', left, left + plot_width, top + plot_height, top + plot_height, 1.2 * vs))
    }
  }
  paste0(c(pieces, "</svg>"), collapse = "")
}

.rnameta_meta_plot_svg <- function(config, width, height) {
  esc <- .rnameta_xml_escape
  sx <- width / 1180; sy <- height / 520; vs <- sqrt(sx * sy)
  margin <- c(top = 58 * sy, right = 28 * sx, bottom = 112 * sy, left = 76 * sx)
  plot_width <- width - margin[["left"]] - margin[["right"]]
  plot_height <- height - margin[["top"]] - margin[["bottom"]]
  x_domain <- as.numeric(unlist(config$xDomain)); y_domain <- as.numeric(unlist(config$yDomain))
  xp <- function(value) margin[["left"]] + (as.numeric(value) - x_domain[[1]]) / diff(x_domain) * plot_width
  yp <- function(value) margin[["top"]] + plot_height - (as.numeric(value) - y_domain[[1]]) / diff(y_domain) * plot_height
  pieces <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s">', width, height, width, height),
    sprintf('<rect width="%s" height="%s" fill="#ffffff"/>', width, height),
    sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="#ffffff"/>', margin[["left"]], margin[["top"]], plot_width, plot_height),
    sprintf('<text x="%s" y="%s" fill="#111111" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', margin[["left"]], 50 * sy, 15 * vs, esc(config$title %||% "")),
    sprintf('<text transform="translate(%s %s) rotate(-90)" fill="#111111" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', 30 * sx, 250 * sy, 15 * vs, esc(config$yLabel %||% ""))
  )
  for (tick in seq(0, 1, .25)) pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#dfe4da" stroke-width="%s"/>', xp(tick), xp(tick), margin[["top"]], margin[["top"]] + plot_height, 1.1 * vs))
  for (tick in pretty(y_domain, n = 5)) pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#dfe4da" stroke-width="%s"/>', margin[["left"]], margin[["left"]] + plot_width, yp(tick), yp(tick), 1.1 * vs), sprintf('<text x="%s" y="%s" text-anchor="end" dominant-baseline="middle" fill="#4c4c4c" font-family="sans-serif" font-size="%s">%s</text>', margin[["left"]] - 9 * vs, yp(tick), 13 * vs, esc(format(signif(tick, 3), trim = TRUE))))
  baseline <- yp(0)
  for (series in config$series) {
    line_points <- vapply(series$values, function(value) sprintf('%s,%s', xp(value$x), yp(value$density)), character(1))
    area_points <- paste(c(sprintf('%s,%s', xp(series$values[[1]]$x), baseline), line_points, sprintf('%s,%s', xp(series$values[[length(series$values)]]$x), baseline)), collapse = " ")
    pieces <- c(pieces, sprintf('<polygon points="%s" fill="%s" fill-opacity="0.18"/>', area_points, esc(series$color)), sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="%s" stroke-linejoin="round"/>', paste(line_points, collapse = " "), esc(series$color), 4 * vs))
    if (isTRUE(config$showCI)) {
      valid <- vapply(series$values, function(value) is.finite(as.numeric(value$confidenceDown)) && is.finite(as.numeric(value$confidenceUp)), logical(1))
      if (all(valid)) {
        upper <- vapply(series$values, function(value) sprintf('%s,%s', xp(value$x), yp(value$confidenceUp)), character(1))
        lower <- rev(vapply(series$values, function(value) sprintf('%s,%s', xp(value$x), yp(value$confidenceDown)), character(1)))
        pieces <- c(pieces, sprintf('<polygon points="%s" fill="%s" fill-opacity="0.12"/>', paste(c(upper, lower), collapse = " "), esc(series$color)))
      }
    }
  }
  components <- config$components %||% list()
  for (separator in components$separators %||% list()) pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#111111" stroke-width="%s" stroke-dasharray="2 6"/>', xp(separator$x), xp(separator$x), yp(separator$y1), yp(separator$y2), 1.6 * vs))
  minimum_height <- c(promoter = 2, tail = 2, utr5 = 4, utr3 = 4, cds = 6, ncrna = 6, rna = 6) * vs
  component_base <- as.numeric(components$trackBase %||% 0)
  label_y <- as.numeric(components$labelY %||% component_base)
  for (segment in components$segments %||% list()) {
    center <- yp(component_base)
    segment_height <- max(minimum_height[[segment$component]] %||% (3 * vs), abs(yp(component_base - segment$height) - yp(component_base + segment$height)))
    x <- xp(segment$start); segment_width <- max(1, xp(segment$end) - x)
    pieces <- c(pieces, sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="#ffffff" fill-opacity="%s" stroke="#111111" stroke-width="%s"/>', x, center - segment_height / 2, segment_width, segment_height, segment$alpha, 1.1 * vs))
    if (nzchar(segment$label %||% "")) pieces <- c(pieces, sprintf('<text x="%s" y="%s" text-anchor="middle" fill="#222222" font-family="sans-serif" font-size="%s" font-weight="700">%s</text>', xp(segment$mid), yp(label_y) + 12 * vs, 13 * vs, esc(segment$label)))
  }
  pieces <- c(pieces, sprintf('<line x1="%s" x2="%s" y1="%s" y2="%s" stroke="#4c4c4c" stroke-width="%s"/>', margin[["left"]], margin[["left"]], margin[["top"]], margin[["top"]] + plot_height, 1.5 * vs), .rnameta_chart_legend_svg(config$series, margin[["left"]], height - 85 * sy, plot_width, vs), '</svg>')
  paste0(pieces, collapse = "")
}
