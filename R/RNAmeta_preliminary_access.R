.rnameta_preliminary_value <- function(data, name, legacy_index) {
  if (!is.list(data)) stop("Provide preliminary analysis data.")
  if (!is.null(names(data)) && name %in% names(data)) return(data[[name]])
  if (length(data) < legacy_index) stop("Preliminary analysis data is missing: ", name)
  data[[legacy_index]]
}