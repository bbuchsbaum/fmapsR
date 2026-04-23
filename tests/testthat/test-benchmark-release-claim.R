load_release_claim_tools <- function() {
  env <- new.env(parent = globalenv())

  root_hints <- unique(Filter(nzchar, c(
    Sys.getenv("FMAPSR_SOURCE_ROOT", unset = ""),
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )))

  candidate_paths <- unlist(lapply(root_hints, function(root) {
    c(
      file.path(root, "tools", "benchmark_release_claim.R"),
      file.path(root, "fmapsR", "tools", "benchmark_release_claim.R")
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
    testthat::skip("tools/benchmark_release_claim.R unavailable outside source checkout")
  }

  source(script, local = env)
  env
}

mock_pairwise_report <- function(target_ratio = 0.30, stable = TRUE, runtime_improvement = 0.39) {
  list(
    comparison = list(
      baseline_available = TRUE,
      target_ratio = target_ratio,
      runtime_improvement_ratio = runtime_improvement,
      target_met_raw = runtime_improvement >= target_ratio,
      target_met_stable = stable,
      stability_decision = if (stable) "pass_strong" else "fail_borderline_unstable"
    )
  )
}

mock_parity_report <- function(partial_objective = 4.1) {
  list(
    scenario_results = data.frame(
      scenario = c("easy", "noisy", "partial"),
      baseline_available = c(TRUE, TRUE, TRUE),
      runtime_comparable = c(TRUE, TRUE, TRUE),
      runtime_improvement_ratio = c(0.55, 0.41, 0.44),
      accuracy_delta = c(0.00, -0.03, 0.03),
      geodesic_norm_delta = c(0.00, -0.02, 0.02),
      map_fro_norm_delta = c(0, 0, 0),
      map_orth_resid_delta = c(0, 0, 0),
      objective_abs_gap = c(1.2, 1.4, partial_objective),
      stringsAsFactors = FALSE
    )
  )
}

test_that("release benchmark parser accepts benchmark overrides", {
  env <- load_release_claim_tools()

  opts <- env$parse_args(c(
    "--output-dir", "tmp/release",
    "--parity-scenarios", "easy,noisy,partial",
    "--required-quality", "easy,partial",
    "--pairwise-r-runs", "4",
    "--pairwise-py-runs", "5",
    "--parity-r-runs", "2",
    "--parity-py-runs", "3",
    "--parity-r-cg-maxit", "1",
    "--speed-target-ratio", "0.25",
    "--speed-stability-margin", "0.03",
    "--speed-stability-boot", "250",
    "--speed-stability-seed", "99"
  ))

  expect_identical(opts$output_dir, "tmp/release")
  expect_equal(opts$parity_scenarios, c("easy", "noisy", "partial"))
  expect_equal(opts$required_quality, c("easy", "partial"))
  expect_equal(opts$pairwise_r_runs, 4L)
  expect_equal(opts$pairwise_py_runs, 5L)
  expect_equal(opts$parity_r_runs, 2L)
  expect_equal(opts$parity_py_runs, 3L)
  expect_equal(opts$parity_r_cg_maxit, 1L)
  expect_equal(opts$speed_target_ratio, 0.25)
  expect_equal(opts$speed_stability_margin, 0.03)
  expect_equal(opts$speed_stability_boot, 250L)
  expect_equal(opts$speed_stability_seed, 99L)
})

test_that("release claim combines speed proof with quality proof", {
  env <- load_release_claim_tools()
  thresholds <- env$benchmark_release_claim_thresholds()

  report <- env$build_release_report(
    pairwise_report = mock_pairwise_report(),
    parity_report = mock_parity_report(),
    thresholds = thresholds,
    root = tempdir(),
    artifacts = list(pairwise_markdown = "pairwise.md", parity_markdown = "parity.md")
  )

  expect_true(isTRUE(report$claim$ok))
  expect_true(isTRUE(report$claim$speed$ok))
  expect_true(isTRUE(report$claim$quality$ok))
  expect_identical(report$claim$speed$source, "parity")
  expect_equal(length(report$claim$quality$checks), 3L)
})

test_that("release claim can require objective parity explicitly", {
  env <- load_release_claim_tools()
  thresholds <- env$benchmark_release_claim_thresholds()
  thresholds$require_objective_parity <- TRUE

  report <- env$build_release_report(
    pairwise_report = mock_pairwise_report(),
    parity_report = mock_parity_report(partial_objective = 4.1),
    thresholds = thresholds,
    root = tempdir()
  )

  expect_false(isTRUE(report$claim$ok))
  partial_chk <- Filter(function(x) identical(x$scenario, "partial"), report$claim$quality$checks)[[1]]
  expect_false(isTRUE(partial_chk$ok))
  expect_true(grepl("objective_abs_gap", partial_chk$reason, fixed = TRUE))
})

test_that("release report writer emits stable latest artifacts", {
  env <- load_release_claim_tools()
  thresholds <- env$benchmark_release_claim_thresholds()
  report <- env$build_release_report(
    pairwise_report = mock_pairwise_report(),
    parity_report = mock_parity_report(),
    thresholds = thresholds,
    root = tempdir(),
    artifacts = list(pairwise_markdown = "pairwise.md", parity_markdown = "parity.md")
  )

  out_dir <- tempfile("release-claim-")
  outputs <- env$write_release_report(report, output_dir = out_dir)

  expect_true(file.exists(outputs$rds))
  expect_true(file.exists(outputs$latest_rds))
  expect_true(file.exists(outputs$markdown))

  md <- readLines(outputs$markdown, warn = FALSE)
  expect_true(any(grepl("Overall claim", md, fixed = TRUE)))
  expect_true(any(grepl("Speed source", md, fixed = TRUE)))
  expect_true(any(grepl("Speed Claim", md, fixed = TRUE)))
  expect_true(any(grepl("Quality Claim", md, fixed = TRUE)))
})
