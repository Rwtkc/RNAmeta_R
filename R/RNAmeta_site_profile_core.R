.rnameta_site_profile_cache <- new.env(parent = emptyenv())

.rnameta_site_profile_cached <- function(key, builder) {
  if (.rnameta_cache_has(.rnameta_site_profile_cache, key)) return(.rnameta_cache_get(.rnameta_site_profile_cache, key))
  value <- builder()
  .rnameta_cache_set(.rnameta_site_profile_cache, key, value)
  value
}

.rnameta_density_points <- function(data, set_group_name = NULL, n = 512L) {
  data <- data.table::as.data.table(data.table::copy(data))
  if (!all(c("feature", "Sample", "rel_pos") %in% names(data))) stop("Site-profile data must contain feature, Sample and rel_pos columns.")
  data <- .rnameta_rename_samples(data, set_group_name)
  groups <- unique(data[, .(feature = as.character(feature), Sample = as.character(Sample))])
  result <- lapply(seq_len(nrow(groups)), function(i) {
    current <- data[as.character(feature) == groups$feature[[i]] & as.character(Sample) == groups$Sample[[i]], as.numeric(rel_pos)]
    current <- current[is.finite(current)]
    if (!length(current)) return(NULL)
    if (all(c("xmin", "xmax") %in% names(data))) {
      xmin <- min(data[as.character(feature) == groups$feature[[i]], as.numeric(xmin)], na.rm = TRUE)
      xmax <- max(data[as.character(feature) == groups$feature[[i]], as.numeric(xmax)], na.rm = TRUE)
    } else if ("flank_size" %in% names(data)) {
      flank <- unique(data[as.character(feature) == groups$feature[[i]], as.numeric(flank_size)])[[1]]
      if (identical(groups$feature[[i]], "5PSS")) { xmin <- -flank; xmax <- 2 * flank } else if (identical(groups$feature[[i]], "3PSS")) { xmin <- -2 * flank; xmax <- flank } else { xmin <- -flank; xmax <- flank }
    } else {
      xmin <- min(current); xmax <- max(current)
    }
    if (!is.finite(xmin)) xmin <- min(current); if (!is.finite(xmax)) xmax <- max(current)
    if (length(unique(current)) < 2L) {
      fit <- stats::density(rep(current[[1]], 2L), from = xmin, to = xmax, n = n, bw = 1)
      x <- fit$x; y <- fit$y
    } else {
      density <- stats::density(current, from = xmin, to = xmax, n = n, adjust = 1)
      x <- density$x; y <- density$y
    }
    data.table::data.table(feature = groups$feature[[i]], Sample = groups$Sample[[i]], x = x, density = y)
  })
  data.table::rbindlist(result, use.names = TRUE, fill = TRUE)
}

.rnameta_downsample_matrix <- function(matrix, max_rows = 160L, max_cols = 360L) {
  matrix <- as.matrix(matrix)
  if (!length(matrix) || !nrow(matrix) || !ncol(matrix)) return(matrix(numeric(), 0, 0))
  row_index <- unique(round(seq(1, nrow(matrix), length.out = min(max_rows, nrow(matrix)))))
  col_index <- unique(round(seq(1, ncol(matrix), length.out = min(max_cols, ncol(matrix)))))
  matrix[row_index, col_index, drop = FALSE]
}

.rnameta_heatmap_panel <- function(key, title, matrix, x_positions = NULL, sample = NULL) {
  matrix <- as.matrix(matrix)
  original_columns <- ncol(matrix)
  row_index <- unique(round(seq(1, nrow(matrix), length.out = min(160L, nrow(matrix)))))
  column_index <- unique(round(seq(1, original_columns, length.out = min(360L, original_columns))))
  matrix <- matrix[row_index, column_index, drop = FALSE]
  if (!is.null(x_positions) && length(x_positions) == original_columns) {
    x_positions <- as.numeric(x_positions[column_index])
  } else {
    x_positions <- seq_len(ncol(matrix))
  }
  list(type = "heatmap", key = key, title = title, sample = sample, sampleName = sample, rows = nrow(matrix), columns = ncol(matrix),
       empty = nrow(matrix) == 0L,
       xMin = min(x_positions), xMax = max(x_positions), xDomain = c(min(x_positions), max(x_positions)),
       backgroundColor = "#ffffff", colorMaxQuantile = 0.92,
       values = lapply(seq_len(nrow(matrix)), function(i) unname(as.numeric(matrix[i, ]))))
}

.rnameta_site_profile_heatmaps <- function(data, samples, set_group_name = NULL, heatmap_sample = "all") {
  if (length(data) < 2L || !is.list(data[[2]]) || !length(data[[2]])) return(list())
  raw <- data[[2]]
  is_legacy <- length(raw) >= 4L && !is.list(raw[[1]])
  heatmaps <- if (is_legacy) stats::setNames(list(raw), samples[[1]]) else raw
  if (is.null(names(heatmaps)) || any(!nzchar(names(heatmaps)))) names(heatmaps) <- paste0("Group", seq_along(heatmaps))
  display_names <- names(heatmaps)
  if (!is.null(set_group_name)) {
    if (is_legacy) {
      display_names <- as.character(set_group_name[[1L]])
    } else {
      if (length(set_group_name) != length(heatmaps)) stop("The number of group names does not match the number of heatmaps.")
      display_names <- as.character(set_group_name)
    }
  }
  if (!is.null(heatmap_sample) && length(heatmap_sample) == 1L && !identical(as.character(heatmap_sample), "all")) {
    index <- if (is.numeric(heatmap_sample)) as.integer(heatmap_sample) else {
      value <- as.character(heatmap_sample); matched <- match(value, display_names)
      if (is.na(matched)) matched <- match(value, names(heatmaps))
      matched
    }
    if (is.na(index) || index < 1L || index > length(heatmaps)) stop("Unknown heatmap_sample: ", heatmap_sample)
    heatmaps <- heatmaps[index]; display_names <- display_names[index]
  }
  stats::setNames(heatmaps, display_names)
}

.rnameta_site_profile_payload <- function(data, set_group_name = NULL, module = "Site Profile", heatmap_sample = "all") {
  hits <- data[[1]]
  density <- .rnameta_density_points(hits, set_group_name)
  samples <- unique(as.character(density$Sample)); features <- unique(as.character(density$feature)); colors <- .rnameta_chart_palette()
  titles <- switch(module,
    Transcription = c(Transcription_start_site = "Distribution on transcription start site", Transcription_end_site = "Distribution on transcription end site"),
    Translation = c(Translation_start_site = "Distribution on Translation start site", Translation_end_site = "Distribution on Translation end site"),
    Splicesite = c(`5PSS` = "Distribution on 5' Splicing Sites (5PSS)", `3PSS` = "Distribution on 3' Splicing Sites (3PSS)"),
    stats::setNames(features, features))
  panels <- lapply(features, function(feature_key) {
    current <- density[as.character(feature) == feature_key]
    max_density <- max(current$density, na.rm = TRUE)
    if (!is.finite(max_density) || max_density <= 0) max_density <- 1
    y_domain <- c(0, max_density * 1.08)
    list(type = "density", key = feature_key, title = unname(titles[[feature_key]] %||% feature_key), xLabel = "Region around the site (0 and 5' -> 3' direction)", yLabel = "Density",
         xDomain = range(current$x, finite = TRUE), yDomain = y_domain, yTicks = pretty(y_domain, n = 5), guideLines = list(0),
         series = lapply(seq_along(samples), function(i) {
           sample <- samples[[i]]; values <- current[Sample == sample]
           list(name = sample, color = colors[((i - 1L) %% length(colors)) + 1L], values = lapply(seq_len(nrow(values)), function(j) list(x = values$x[[j]], density = values$density[[j]])))
         }))
  })
  if (module %in% c("Transcription", "Translation")) {
    heatmaps <- .rnameta_site_profile_heatmaps(data, samples, set_group_name, heatmap_sample)
    for (sample in names(heatmaps)) {
      heat <- heatmaps[[sample]]
      if (!is.list(heat) || length(heat) < 4L) next
      sample_key <- gsub("[^A-Za-z0-9]+", "_", sample)
      panels <- c(panels, list(
        .rnameta_heatmap_panel(paste0("Tss_heatmap_", sample_key), paste0(module, ".Tss.Heatmap - ", sample), heat[[4]], heat[[2]], sample),
        .rnameta_heatmap_panel(paste0("Tes_heatmap_", sample_key), paste0(module, ".Tes.Heatmap - ", sample), heat[[3]], heat[[1]], sample)
      ))
    }
  }
  list(type = "site_profile", title = module, panels = panels, heatmapSamples = names(.rnameta_site_profile_heatmaps(data, samples, set_group_name, "all")),
       legend = lapply(seq_along(samples), function(i) list(name = samples[[i]], color = colors[((i - 1L) %% length(colors)) + 1L])))
}

.rnameta_site_profile_payloads <- function(data, set_group_name = NULL, module = "Site Profile", heatmap_sample = "all") {
  config <- .rnameta_site_profile_payload(data, set_group_name, module, heatmap_sample)
  panels <- config$panels
  if (!length(panels)) stop("No Site Profile panels are available.")
  keys <- make.unique(vapply(panels, function(panel) as.character(panel$key %||% panel$title), character(1)), sep = "_")
  payloads <- lapply(seq_along(panels), function(i) {
    single <- config
    single$panels <- list(panels[[i]])
    if (!identical(panels[[i]]$type, "density")) single$legend <- list()
    single
  })
  stats::setNames(payloads, keys)
}

.rnameta_select_site_profile_payloads <- function(payloads, panel = "all") {
  if (is.null(panel) || (length(panel) == 1L && identical(as.character(panel), "all"))) return(payloads)
  if (is.numeric(panel)) {
    index <- as.integer(panel)
    if (anyNA(index) || any(index < 1L | index > length(payloads))) stop("panel index is out of range.")
  } else {
    panel <- as.character(panel)
    index <- match(panel, names(payloads))
    if (anyNA(index)) stop("Unknown panel: ", paste(panel[is.na(index)], collapse = ", "), ". Available panels: ", paste(names(payloads), collapse = ", "))
  }
  payloads[index]
}

.rnameta_site_profile_plot <- function(data, set_group_name, module, width = NULL, height = NULL, heatmap_sample = "all", panel = "all") {
  payloads <- .rnameta_select_site_profile_payloads(
    .rnameta_site_profile_payloads(data, set_group_name, module, heatmap_sample), panel
  )
  widgets <- lapply(payloads, .rnameta_chart_widget, width = width, height = height)
  if (length(widgets) == 1L) return(widgets[[1L]])
  structure(widgets, class = c("RNAmeta_site_profile_plots", "list"))
}

.rnameta_site_profile_file_names <- function(file, panel_names) {
  extension <- tolower(tools::file_ext(file))
  stem <- tools::file_path_sans_ext(file)
  safe <- gsub("[^A-Za-z0-9._-]+", "_", panel_names)
  paste0(stem, "_", safe, ".", extension)
}

.rnameta_site_profile_export <- function(data, file, set_group_name, module, width = 1180, height = NULL, aspect_ratio = c(1180, 560), scale = 2, heatmap_sample = "all", panel = "all") {
  payloads <- .rnameta_select_site_profile_payloads(
    .rnameta_site_profile_payloads(data, set_group_name, module, heatmap_sample), panel
  )
  files <- if (length(payloads) == 1L) file else .rnameta_site_profile_file_names(file, names(payloads))
  paths <- vapply(seq_along(payloads), function(i) {
    .rnameta_chart_export(payloads[[i]], files[[i]], width, height, aspect_ratio, scale)
  }, character(1))
  names(paths) <- names(payloads)
  invisible(paths)
}
