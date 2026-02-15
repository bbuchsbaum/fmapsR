load_pairwise_tools <- function() {
  env <- new.env(parent = globalenv())

  root_hints <- unique(Filter(nzchar, c(
    Sys.getenv("FMAPSR_SOURCE_ROOT", unset = ""),
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )))

  candidate_paths <- unlist(lapply(root_hints, function(root) {
    c(
      file.path(root, "tools", "benchmark_pairwise.R"),
      file.path(root, "fmapsR", "tools", "benchmark_pairwise.R")
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
    testthat::skip("tools/benchmark_pairwise.R unavailable outside source checkout")
  }

  source(script, local = env)
  env
}

test_that("pairwise benchmark arg parser accepts stability options", {
  env <- load_pairwise_tools()

  defaults <- env$parse_args(character())
  expect_equal(defaults$stability_seed, 2026L)

  opts <- env$parse_args(c(
    "--r-runs", "5",
    "--py-runs", "4",
    "--target-ratio", "0.25",
    "--stability-margin", "0.03",
    "--stability-boot", "250",
    "--stability-seed", "99",
    "--output-dir", "tmp/out"
  ))

  expect_equal(opts$r_runs, 5L)
  expect_equal(opts$py_runs, 4L)
  expect_equal(opts$target_ratio, 0.25)
  expect_equal(opts$stability_margin, 0.03)
  expect_equal(opts$stability_boot, 250L)
  expect_equal(opts$stability_seed, 99L)
  expect_identical(opts$output_dir, "tmp/out")
})

test_that("improvement samples use bootstrap over independent runtime draws", {
  env <- load_pairwise_tools()

  s <- env$compute_improvement_samples(
    r_runtimes = c(0.7, 0.8, 0.9),
    py_runtimes = c(1.0, 1.0, 1.0),
    n_boot = 200L,
    seed = 123L
  )

  expect_equal(length(s), 200L)
  expect_true(all(is.finite(s)))
  expect_true(all(s >= 0.09 & s <= 0.31))
  expect_gt(mean(s), 0.16)
  expect_lt(mean(s), 0.24)

  s2 <- env$compute_improvement_samples(
    r_runtimes = c(0.7, 0.8),
    py_runtimes = c(1.0, 0.0, 2.0),
    n_boot = 120L,
    seed = 7L
  )
  expect_equal(length(s2), 120L)
  expect_true(all(is.finite(s2)))
})



test_that("stability sample preprocessing drops one high outlier when n>=3", {
  env <- load_pairwise_tools()

  meta <- env$prepare_stability_samples(c(0.009, 0.010, 0.090), drop_max = TRUE)
  expect_equal(meta$n_raw, 3L)
  expect_equal(meta$n_used, 2L)
  expect_equal(meta$n_dropped_high, 1L)
  expect_true(all(meta$values < 0.02))

  meta2 <- env$prepare_stability_samples(c(0.009, 0.010), drop_max = TRUE)
  expect_equal(meta2$n_raw, 2L)
  expect_equal(meta2$n_used, 2L)
  expect_equal(meta2$n_dropped_high, 0L)
})

test_that("runtime target assessment distinguishes strong and borderline zones", {
  env <- load_pairwise_tools()

  strong <- env$assess_runtime_target(
    improvement = 0.36,
    improvement_samples = c(0.35, 0.37, 0.34),
    target_ratio = 0.30,
    stability_margin = 0.02
  )
  expect_true(isTRUE(strong$target_met_stable))
  expect_identical(strong$classification, "pass_strong")

  borderline_stable <- env$assess_runtime_target(
    improvement = 0.31,
    improvement_samples = c(0.29, 0.31, 0.33),
    target_ratio = 0.30,
    stability_margin = 0.02
  )
  expect_true(isTRUE(borderline_stable$target_met_stable))
  expect_identical(borderline_stable$classification, "pass_borderline_stable")

  borderline_unstable <- env$assess_runtime_target(
    improvement = 0.31,
    improvement_samples = c(0.24, 0.31, 0.34),
    target_ratio = 0.30,
    stability_margin = 0.02
  )
  expect_false(isTRUE(borderline_unstable$target_met_stable))
  expect_identical(borderline_unstable$classification, "fail_borderline_unstable")

  strong_fail <- env$assess_runtime_target(
    improvement = 0.20,
    improvement_samples = c(0.19, 0.20, 0.22),
    target_ratio = 0.30,
    stability_margin = 0.02
  )
  expect_false(isTRUE(strong_fail$target_met_stable))
  expect_identical(strong_fail$classification, "fail_strong")
})

test_that("pairwise report writer includes stability metadata", {
  env <- load_pairwise_tools()

  report <- list(
    metadata = list(
      timestamp = "2026-02-12 00:00:00 UTC",
      r_version = "4.4.0",
      platform = "x86_64",
      target_ratio = 0.30,
      stability_margin = 0.02,
      stability_boot = 500L
    ),
    r_pipeline = list(
      median_runtime_sec = 0.5,
      mean_runtime_sec = 0.55,
      optimizer = "cg",
      kernel_backend = "cpp",
      median_objective = 1.0,
      mem_bytes = NA_real_
    ),
    backend_comparison = list(
      baseline_backend = "r",
      active_backend = "cpp",
      baseline_runtime_sec = 0.8,
      active_runtime_sec = 0.5,
      active_vs_baseline_speedup = 0.375,
      objective_gap = 0.01
    ),
    pyfm_baseline = list(status = "ok", n_runs = 3L),
    comparison = list(
      pyfm_runtime_sec = 0.7,
      pyfm_objective = 1.1,
      runtime_improvement_ratio = 0.285,
      runtime_improvement_sample_median = 0.29,
      runtime_improvement_sample_q25 = 0.27,
      runtime_improvement_sample_q75 = 0.31,
      runtime_improvement_boot_n = 500L,
      stability_r_runtime_n_raw = 3L,
      stability_r_runtime_n_used = 2L,
      stability_py_runtime_n_raw = 3L,
      stability_py_runtime_n_used = 2L,
      target_ratio = 0.30,
      target_band_low = 0.28,
      target_band_high = 0.32,
      target_met_raw = FALSE,
      target_met_stable = TRUE,
      stability_decision = "pass_borderline_stable",
      baseline_available = TRUE,
      notes = "Borderline but stable in bootstrap spread"
    )
  )

  out_dir <- tempfile("pairwise-report-")
  outputs <- env$write_report(report, output_dir = out_dir)

  expect_true(file.exists(outputs$rds))
  expect_true(file.exists(outputs$markdown))

  md <- readLines(outputs$markdown, warn = FALSE)
  expect_true(any(grepl("Stable target met", md, fixed = TRUE)))
  expect_true(any(grepl("Stability decision", md, fixed = TRUE)))
  expect_true(any(grepl("Runtime improvement bootstrap samples", md, fixed = TRUE)))
  expect_true(any(grepl("Stability sample usage", md, fixed = TRUE)))
})
