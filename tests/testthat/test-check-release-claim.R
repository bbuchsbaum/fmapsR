load_release_gate_tools <- function() {
  env <- new.env(parent = globalenv())

  root_hints <- unique(Filter(nzchar, c(
    Sys.getenv("FMAPSR_SOURCE_ROOT", unset = ""),
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )))

  candidate_paths <- unlist(lapply(root_hints, function(root) {
    c(
      file.path(root, "tools", "check_release_claim.R"),
      file.path(root, "fmapsR", "tools", "check_release_claim.R")
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
    testthat::skip("tools/check_release_claim.R unavailable outside source checkout")
  }

  source(script, local = env)
  env
}

test_that("release gate parser supports objective and speed overrides", {
  env <- load_release_gate_tools()

  opts <- env$parse_args(c(
    "--report", "x.rds",
    "--required-quality", "easy,noisy",
    "--speed-target-ratio", "0.25",
    "--speed-stability-margin", "0.03",
    "--require-objective",
    "--no-stable-speed"
  ))

  expect_identical(opts$report, "x.rds")
  expect_equal(opts$required_quality, c("easy", "noisy"))
  expect_equal(opts$speed_target_ratio, 0.25)
  expect_equal(opts$speed_stability_margin, 0.03)
  expect_true(isTRUE(opts$require_objective_parity))
  expect_false(isTRUE(opts$require_stable_speed))
})

test_that("release gate thresholds inherit overrides", {
  env <- load_release_gate_tools()

  thr <- env$thresholds_from_args(list(
    required_quality = c("easy"),
    speed_target_ratio = 0.22,
    speed_stability_margin = 0.01,
    require_stable_speed = FALSE,
    require_objective_parity = TRUE
  ))

  expect_equal(thr$required_quality_scenarios, "easy")
  expect_equal(thr$speed_target_ratio, 0.22)
  expect_equal(thr$speed_stability_margin, 0.01)
  expect_false(isTRUE(thr$require_stable_speed))
  expect_true(isTRUE(thr$require_objective_parity))
})
