load_parity_gate_tools <- function() {
  env <- new.env(parent = globalenv())

  root_hints <- unique(Filter(nzchar, c(
    Sys.getenv("FMAPSR_SOURCE_ROOT", unset = ""),
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )))

  candidate_paths <- unlist(lapply(root_hints, function(root) {
    c(
      file.path(root, "tools", "check_parity_gate.R"),
      file.path(root, "fmapsR", "tools", "check_parity_gate.R")
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
    testthat::skip("tools/check_parity_gate.R unavailable outside source checkout")
  }

  source(script, local = env)
  env
}

test_that("parity gate parser supports required scenario override", {
  env <- load_parity_gate_tools()

  opts <- env$parse_args(c("--report", "x.rds", "--required", "easy,noisy"))
  expect_identical(opts$report, "x.rds")
  expect_equal(opts$required, c("easy", "noisy"))
})

test_that("parity gate scenario evaluator detects threshold failures", {
  env <- load_parity_gate_tools()
  thr <- env$benchmark_parity_gate_thresholds()

  row_ok <- data.frame(
    scenario = "easy",
    baseline_available = TRUE,
    accuracy_delta = -0.05,
    geodesic_norm_delta = 0.02
  )
  chk_ok <- env$scenario_gate_row(row_ok, thresholds = thr)
  expect_true(isTRUE(chk_ok$ok))

  row_bad <- data.frame(
    scenario = "noisy",
    baseline_available = TRUE,
    accuracy_delta = -0.30,
    geodesic_norm_delta = 0.02
  )
  chk_bad <- env$scenario_gate_row(row_bad, thresholds = thr)
  expect_false(isTRUE(chk_bad$ok))
  expect_true(grepl("accuracy_delta", chk_bad$reason, fixed = TRUE))
})
