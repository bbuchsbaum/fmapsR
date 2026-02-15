load_network_tools <- function() {
  env <- new.env(parent = globalenv())

  root_hints <- unique(Filter(nzchar, c(
    Sys.getenv("FMAPSR_SOURCE_ROOT", unset = ""),
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )))

  candidate_paths <- unlist(lapply(root_hints, function(root) {
    c(
      file.path(root, "tools", "benchmark_network_sync.R"),
      file.path(root, "fmapsR", "tools", "benchmark_network_sync.R")
    )
  }), use.names = FALSE)

  script <- NULL
  for (p in candidate_paths) {
    pp <- normalizePath(p, winslash = "/", mustWork = FALSE)
    if (file.exists(pp)) {
      script <- pp
      break
    }
  }

  if (is.null(script)) {
    testthat::skip("tools/benchmark_network_sync.R unavailable outside source checkout")
  }

  source(script, local = env)
  env
}

test_that("network benchmark parser preserves frozen defaults and overrides", {
  env <- load_network_tools()

  opts_default <- env$parse_args(character())
  expect_equal(opts_default$scales, c(10L, 25L, 50L))
  expect_equal(opts_default$seed, 2026L)
  expect_identical(opts_default$mode, "robust")

  opts <- env$parse_args(c("--scales", "8,16", "--seed", "99", "--mode", "cycle"))
  expect_equal(opts$scales, c(8L, 16L))
  expect_equal(opts$seed, 99L)
  expect_identical(opts$mode, "cycle")
})

test_that("network benchmark reference targets include baseline scales", {
  env <- load_network_tools()
  targets <- env$reference_targets()

  expect_equal(as.numeric(targets[["10"]]), 0.60)
  expect_equal(as.numeric(targets[["25"]]), 0.55)
  expect_equal(as.numeric(targets[["50"]]), 0.50)
})
