# Numerical kernels from deployed RNAmeta Online; see inst/online-algorithms.json.

.rnameta_online_meta_plot_granges_list_map_to_transcripts <- function(site, map_filter_transcript = FALSE, transcripts) {
  names(site) <- seq_along(site)
  x_widths <- sum(IRanges::width(site))
  names(x_widths) <- names(site)
  tx_coord <- GenomicFeatures::mapToTranscripts(BiocGenerics::unlist(site), transcripts, ignore.strand = FALSE)
  grouped <- GenomicRanges::split(tx_coord, paste(names(tx_coord), tx_coord$transcriptsHits, sep = "-"))
  reduced_grouped <- GenomicRanges::reduce(grouped)
  mapping_filtered <- BiocGenerics::unlist(reduced_grouped[S4Vectors::elementNROWS(reduced_grouped) == 1L], use.names = TRUE)
  x_hit_tx_hit <- strsplit(names(mapping_filtered), "-", fixed = TRUE)
  S4Vectors::mcols(mapping_filtered) <- S4Vectors::DataFrame(
    xHits = as.numeric(vapply(x_hit_tx_hit, `[`, character(1), 1)),
    txHits = as.numeric(vapply(x_hit_tx_hit, `[`, character(1), 2))
  )
  if (map_filter_transcript) {
    mapping_filtered <- mapping_filtered[IRanges::width(mapping_filtered) == x_widths[mapping_filtered$xHits]]
  }
  out_of_bound <- GenomicRanges:::get_out_of_bound_index(mapping_filtered)
  if (length(out_of_bound) > 0) {
    mapping_filtered <- mapping_filtered[-out_of_bound]
  }
  mapping_filtered
}

.rnameta_online_meta_plot_build_sample_point_ranges <- function(mapped, controls) {
  st_sample_num <- 2 * controls$stSampleNum - 1
  site_width <- IRanges::width(mapped)
  sample_fn <- if (identical(controls$stSampleModle, "random")) {
    function(x, i) if (i == 1) round(x / 2) else sort(sample(x, i, replace = FALSE))
  } else {
    function(x, i) if (i == 1) round(x / 2) else round(seq(1, x - 1, length.out = i))
  }
  site_points <- t(vapply(site_width, sample_fn, i = st_sample_num, numeric(st_sample_num)))
  site_point_vector <- as.numeric(t(site_points))
  point_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    data.frame(
      chr = rep(as.character(GenomicRanges::seqnames(mapped)), each = st_sample_num),
      start = site_point_vector + BiocGenerics::start(rep(mapped, each = st_sample_num)) - 1,
      end = site_point_vector + BiocGenerics::start(rep(mapped, each = st_sample_num))
    )
  )
  overlap_count <- ave(seq_along(mapped), mapped$xHits, FUN = length)
  S4Vectors::mcols(point_ranges) <- S4Vectors::DataFrame(
    sitesLength = rep(site_width, each = st_sample_num),
    xHits = rep(mapped$xHits, each = st_sample_num),
    pointsOverlapTx = rep(overlap_count, each = st_sample_num)
  )
  point_ranges
}

.rnameta_online_meta_plot_normalize <- function(points, guitar_txdb, controls) {
  if (!length(points)) return(list(position = numeric(), weight = numeric()))
  tx_type <- controls$pltTxType
  points_pos_tx <- BiocGenerics::end(points)
  names(points_pos_tx) <- as.character(GenomicRanges::seqnames(points))
  start_point_mat <- guitar_txdb[[tx_type]]$startPoint[names(points_pos_tx), , drop = FALSE]
  start_point_differ <- start_point_mat - points_pos_tx
  point_component <- apply(start_point_differ, 1, function(x) max(which(x < 0)))
  point_position_component <- -start_point_differ[cbind(seq_along(point_component), point_component)]
  component_width_avg <- guitar_txdb[[tx_type]]$componentWidthAverage_pct[point_component]
  component_start_pct <- guitar_txdb[[tx_type]]$componentStartAverage_pct[point_component]
  component_width_mat <- guitar_txdb[[tx_type]]$componentWidth[names(points_pos_tx), , drop = FALSE]
  point_component_width <- component_width_mat[cbind(seq_along(point_component), point_component)]
  point_position_normalized <- point_position_component / point_component_width * component_width_avg + component_start_pct
  names(point_position_normalized) <- points$xHits
  component_pct <- guitar_txdb[[tx_type]]$componentWidthPtc[names(point_component), , drop = FALSE]
  point_component_pct <- BiocGenerics::unlist(lapply(seq_along(point_component), function(i) component_pct[i, ][point_component[[i]]]))
  point_weight <- component_width_avg / (points$pointsOverlapTx ^ controls$overlapIndex) / point_component_pct * (points$sitesLength ^ controls$siteLengthIndex)
  names(point_weight) <- points$xHits
  list(position = point_position_normalized, weight = point_weight)
}

.rnameta_online_meta_plot_generate_density_ci <- function(relative_points, point_weight, adjust, enableCI, CI_ResamplingTime = 1000, CI_interval = c(0.025, 0.975)) {
  density_dataframe <- data.frame()
  for (group_name in names(relative_points)) {
    valid <- is.finite(point_weight[[group_name]]) & point_weight[[group_name]] > 0 & is.finite(relative_points[[group_name]])
    point_weight[[group_name]] <- point_weight[[group_name]][valid]
    relative_points[[group_name]] <- relative_points[[group_name]][valid]
    if (length(relative_points[[group_name]]) < 2L) next
    site_weight <- point_weight[[group_name]] / sum(point_weight[[group_name]])
    site_id <- relative_points[[group_name]]
    fit <- suppressWarnings(stats::density(site_id, adjust = adjust, from = 0, to = 1, n = 256, weight = site_weight))
    tmp <- data.frame(x = fit$x, density = fit$y / (sum(diff(fit$x) * (head(fit$y, -1) + tail(fit$y, -1))) / 2), group = rep(group_name, length(fit$y)))
    if (enableCI) {
      point_index <- seq_along(site_id)
      names(point_index) <- names(site_id)
      point_grouped <- GenomicRanges::split(point_index, names(point_index))
      fit_resampled <- replicate(CI_ResamplingTime, {
        resampled_groups <- sample(names(point_grouped), replace = TRUE)
        resampled_index <- BiocGenerics::unlist(point_grouped[resampled_groups])
        suppressWarnings(stats::density(site_id[resampled_index], adjust = adjust, from = 0, to = 1, n = 256, weight = site_weight[resampled_index])$y)
      })
      fit_quantile <- apply(fit_resampled, 1, stats::quantile, CI_interval)
      tmp$confidenceDown <- fit_quantile[1, ]
      tmp$confidenceUp <- fit_quantile[2, ]
    }
    density_dataframe <- rbind(density_dataframe, tmp)
  }
  density_dataframe
}
