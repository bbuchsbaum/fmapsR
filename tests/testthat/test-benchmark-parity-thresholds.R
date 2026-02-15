load_parity_threshold_tools <- function() {
  env <- new.env(parent = globalenv())

  root_hints <- unique(Filter(nzchar, c(
    Sys.getenv("FMAPSR_SOURCE_ROOT", unset = ""),
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )))

  candidate_paths <- unlist(lapply(root_hints, function(root) {
    c(
      file.path(root, "tools", "recommend_parity_thresholds.R"),
      file.path(root, "fmapsR", "tools", "recommend_parity_thresholds.R")
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
    testthat::skip("tools/recommend_parity_thresholds.R unavailable outside source checkout")
  }

  source(script, local = env)
  env
}

test_that("parity threshold recommender parser accepts tuning flags", {
  env <- load_parity_threshold_tools()

  opts <- env$parse_args(c(
    "--dir", "benchmarks/parity",
    "--quantile", "0.90",
    "--safety", "1.5",
    "--min-n", "5"
  ))

  expect_identical(opts$dir, "benchmarks/parity")
  expect_equal(opts$quantile, 0.90)
  expect_equal(opts$safety, 1.5)
  expect_equal(opts$min_n, 5L)
})

test_that("recommend_threshold falls back on insufficient sample size", {
  env <- load_parity_threshold_tools()

  rec <- env$recommend_threshold(
    values = c(1.0, 1.2),
    current = 2.5,
    q = 0.95,
    safety = 1.25,
    min_n = 3L,
    floor_value = 0.1
  )

  expect_equal(rec$value, 2.5)
  expect_equal(rec$n, 2L)
  expect_identical(rec$reason, "insufficient_data")
})

test_that("collect_parity_rows harmonizes legacy and current report schemas", {
  env <- load_parity_threshold_tools()
  tmp <- tempfile("parity-reco-")
  dir.create(tmp, recursive = TRUE)

  report_legacy <- list(
    scenario_results = data.frame(
      scenario = "easy",
      baseline_available = TRUE,
      r_objective = 2.0,
      py_objective = 1.5,
      stringsAsFactors = FALSE
    )
  )
  report_current <- list(
    scenario_results = data.frame(
      scenario = "noisy",
      baseline_available = TRUE,
      objective_abs_gap = 1.2,
      map_fro_norm_delta = 1e-3,
      map_orth_resid_delta = 2e-3,
      stringsAsFactors = FALSE
    )
  )

  saveRDS(report_legacy, file.path(tmp, "parity-benchmark-legacy.rds"))
  saveRDS(report_current, file.path(tmp, "parity-benchmark-current.rds"))

  rows <- env$collect_parity_rows(tmp)
  expect_equal(nrow(rows), 2L)
  expect_true(all(c("objective_abs_gap", "map_fro_norm_delta", "map_orth_resid_delta") %in% names(rows)))
  easy <- rows[rows$scenario == "easy", , drop = FALSE]
  expect_equal(as.numeric(easy$objective_abs_gap), 0.5)
})
