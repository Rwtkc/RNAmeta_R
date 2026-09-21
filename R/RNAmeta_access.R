.onAttach <- function(libname, pkgname) {
  logo <- c(
    " ____  _   _    _    __  __      _        ",
    "|  _ \\| \\ | |  / \\  |  \\/  | ___| |_ __ _ ",
    "| |_) |  \\| | / _ \\ | |\\/| |/ _ \\ __/ _` |",
    "|  _ <| |\\  |/ ___ \\| |  | |  __/ || (_| |",
    "|_| \\_\\_| \\_/_/   \\_\\_|  |_|\\___|\\__\\__,_|"
  )
  online_url <- "https://rnainformatics.cn/RNAmeta-online/"
  github_url <- "https://github.com/Rwtkc/RNAmeta_R"
  email <- "rwtkc@foxmail.com"

  message <- c(
    crayon::green(logo),
    "",
    crayon::bold$green(
      "Welcome to RNAmeta! This package is currently under active development."
    ),
    crayon::silver("Genomic region annotation, profiling and visualization"),
    "",
    paste0("  ", crayon::bold("Online "), "  ", cli::style_hyperlink(online_url, online_url)),
    paste0("  ", crayon::bold("GitHub "), "  ", cli::style_hyperlink(github_url, github_url)),
    paste0("  ", crayon::bold("Contact"), "  ", cli::style_hyperlink(email, paste0("mailto:", email)))
  )
  packageStartupMessage(paste(message, collapse = "\n"))
}