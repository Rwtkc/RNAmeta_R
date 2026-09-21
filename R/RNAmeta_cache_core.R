.rnameta_cache_limit <- function() {
  value <- getOption("RNAmeta.cache.max_entries", 128L)
  value <- suppressWarnings(as.integer(value))
  if (length(value) != 1L || !is.finite(value) || value < 1L) 128L else value
}

.rnameta_cache_order <- function(cache) {
  order <- attr(cache, "rnameta_order", exact = TRUE)
  if (is.null(order)) character() else as.character(order)
}

.rnameta_cache_touch <- function(cache, key) {
  order <- .rnameta_cache_order(cache)
  attr(cache, "rnameta_order") <- c(order[order != key], key)
  invisible(cache)
}

.rnameta_cache_has <- function(cache, key) {
  exists(key, envir = cache, inherits = FALSE)
}

.rnameta_cache_get <- function(cache, key) {
  value <- get(key, envir = cache, inherits = FALSE)
  .rnameta_cache_touch(cache, key)
  value
}

.rnameta_cache_set <- function(cache, key, value) {
  limit <- .rnameta_cache_limit()
  if (!exists(key, envir = cache, inherits = FALSE)) {
    order <- .rnameta_cache_order(cache)
    while (length(order) >= limit) {
      oldest <- order[[1L]]
      if (exists(oldest, envir = cache, inherits = FALSE)) rm(list = oldest, envir = cache)
      order <- order[-1L]
    }
    attr(cache, "rnameta_order") <- order
  }
  assign(key, value, envir = cache)
  .rnameta_cache_touch(cache, key)
  invisible(value)
}

.rnameta_cache_clear_env <- function(cache) {
  keys <- ls(envir = cache, all.names = TRUE)
  if (length(keys)) rm(list = keys, envir = cache)
  attr(cache, "rnameta_order") <- character()
  length(keys)
}

.rnameta_validate_seed <- function(seed) {
  seed <- suppressWarnings(as.integer(seed))
  if (length(seed) != 1L || is.na(seed) || seed < 0L) stop("seed must be one non-negative integer.")
  seed
}

.rnameta_seed_offset <- function(seed, offset = 0L) {
  seed <- .rnameta_validate_seed(seed)
  value <- (as.double(seed) + as.double(offset)) %% .Machine$integer.max
  as.integer(value)
}

.rnameta_with_seed <- function(seed, expression) {
  seed <- .rnameta_validate_seed(seed)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

#' Clear RNAmeta in-memory caches
#'
#' @param scope Cache group to clear: all caches, annotation/peak annotation,
#'   site-profile calculations, or Meta Plot calculations.
#' @return Invisibly returns a named integer vector with removed entry counts.
#' @export
RNAmeta_clear_cache <- function(scope = c("all", "annotation", "site_profile", "meta")) {
  scope <- match.arg(scope)
  caches <- list(
    annotation = list(.rnameta_annotation_cache, .rnameta_peak_annotation_cache),
    site_profile = list(.rnameta_site_profile_cache),
    meta = list(.rnameta_meta_cache)
  )
  selected <- if (identical(scope, "all")) names(caches) else scope
  counts <- vapply(selected, function(name) {
    sum(vapply(caches[[name]], .rnameta_cache_clear_env, integer(1)))
  }, integer(1))
  invisible(counts)
}