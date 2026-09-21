source("R/RNAmeta_cache_core.R")
source("R/RNAmeta_meta_plot_dependency_functions.R")

.rnameta_annotation_cache <- new.env(parent = emptyenv())
.rnameta_peak_annotation_cache <- new.env(parent = emptyenv())
.rnameta_site_profile_cache <- new.env(parent = emptyenv())
.rnameta_meta_cache <- new.env(parent = emptyenv())

old_option <- getOption("RNAmeta.cache.max_entries")
on.exit(options(RNAmeta.cache.max_entries = old_option), add = TRUE)
options(RNAmeta.cache.max_entries = 2L)
.rnameta_cache_set(.rnameta_site_profile_cache, "a", 1L)
.rnameta_cache_set(.rnameta_site_profile_cache, "b", 2L)
stopifnot(identical(.rnameta_cache_get(.rnameta_site_profile_cache, "a"), 1L))
.rnameta_cache_set(.rnameta_site_profile_cache, "c", 3L)
stopifnot(.rnameta_cache_has(.rnameta_site_profile_cache, "a"),
          .rnameta_cache_has(.rnameta_site_profile_cache, "c"),
          !.rnameta_cache_has(.rnameta_site_profile_cache, "b"))

set.seed(99)
before <- .Random.seed
first <- .rnameta_with_seed(7L, sample.int(100L, 20L))
stopifnot(identical(before, .Random.seed), identical(first, .rnameta_with_seed(7L, sample.int(100L, 20L))))

positions <- list(Group1 = stats::setNames(seq(0.05, 0.95, length.out = 30), rep(paste0("site", 1:10), each = 3)))
weights <- list(Group1 = rep(1, 30))
density_one <- .rnameta_with_seed(12L, .generate_density_confidence_interval(
  positions, weights, confidence_interval_resampling_time = 25L,
  adjustment_factor = 1, enable_confidence_interval = TRUE
))
density_two <- .rnameta_with_seed(12L, .generate_density_confidence_interval(
  positions, weights, confidence_interval_resampling_time = 25L,
  adjustment_factor = 1, enable_confidence_interval = TRUE
))
stopifnot(identical(density_one, density_two))

counts <- RNAmeta_clear_cache("all")
stopifnot(all(counts >= 0L), !length(ls(.rnameta_site_profile_cache)), !length(ls(.rnameta_meta_cache)))
cat("cache_lru: TRUE\nseed_reproducible: TRUE\nrng_state_restored: TRUE\ncache_clear_all: TRUE\n")