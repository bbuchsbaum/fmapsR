#' Package Metadata Helper
#'
#' @return A named list describing the current package development state.
#' @export
fm_package_info <- function() {
  list(
    package = "fmapsR",
    version = utils::packageVersion("fmapsR"),
    status = "active-development",
    milestones = c("M1", "M2", "M3")
  )
}
