#' Package Metadata Helper
#'
#' @return A named list describing the current package scaffold state.
#' @export
fm_package_info <- function() {
  list(
    package = "fmapsR",
    version = utils::packageVersion("fmapsR"),
    status = "scaffold"
  )
}
