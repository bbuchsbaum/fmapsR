load_parity_tools <- function() {
  env <- new.env(parent = globalenv())

  root_hints <- unique(Filter(nzchar, c(
    Sys.getenv("FMAPSR_SOURCE_ROOT", unset = ""),
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )))

  candidate_paths <- unlist(lapply(root_hints, function(root) {
    c(
      file.path(root, "tools", "benchmark_parity.R"),
      file.path(root, "fmapsR", "tools", "benchmark_parity.R")
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
    testthat::skip("tools/benchmark_parity.R unavailable outside source checkout")
  }

  source(script, local = env)
  env
}

test_that("parity key-value parser reads script output", {
  env <- load_parity_tools()
  kv <- env$parity_parse_kv(c("status=ok", "runtime_sec=0.123", "note=a=b"))

  expect_identical(kv$status, "ok")
  expect_identical(kv$runtime_sec, "0.123")
  expect_identical(kv$note, "a=b")
})

test_that("parity argument parser supports tuned runtime knobs", {
  env <- load_parity_tools()

  opts_default <- env$parse_args(character())
  expect_equal(opts_default$r_cg_maxit, 3L)
  expect_equal(opts_default$r_refine_nit, 2L)
  expect_equal(opts_default$r_cg_tol, 1e-6)
  expect_identical(opts_default$r_kernel_backend, "auto")
  expect_equal(opts_default$seed, 42L)

  opts <- env$parse_args(c(
    "--r-cg-maxit", "4",
    "--r-cg-tol", "1e-4",
    "--r-refine-nit", "2",
    "--r-kernel-backend", "cpp"
  ))
  expect_equal(opts$r_cg_maxit, 4L)
  expect_equal(opts$r_cg_tol, 1e-4)
  expect_equal(opts$r_refine_nit, 2L)
  expect_identical(opts$r_kernel_backend, "cpp")
})

test_that("parity report writer emits markdown and rds artifacts", {
  env <- load_parity_tools()

  scenario_df <- data.frame(
    scenario = c("easy", "noisy"),
    r_kernel_backend = c("cpp", "cpp"),
    runtime_comparable = c(TRUE, TRUE),
    r_runtime_sec = c(0.1, 0.2),
    py_runtime_sec = c(0.2, 0.3),
    runtime_improvement_ratio = c(0.5, 0.333),
    r_accuracy = c(0.9, 0.8),
    py_accuracy = c(0.85, 0.75),
    accuracy_delta = c(0.05, 0.05),
    r_geodesic_norm = c(0.05, 0.09),
    py_geodesic_norm = c(0.08, 0.12),
    geodesic_norm_delta = c(-0.03, -0.03),
    baseline_available = c(TRUE, TRUE),
    note = c("ok", "ok"),
    stringsAsFactors = FALSE
  )

  report <- list(
    metadata = list(
      timestamp = "2026-02-12 00:00:00 UTC",
      platform = "x86_64",
      scenarios = c("easy", "noisy"),
      r_runs = 1L,
      py_runs = 1L,
      r_cg_maxit = 3L,
      r_cg_tol = 1e-6,
      r_refine_nit = 2L,
      r_kernel_backend = "auto",
      package_load_mode = "pkgload"
    ),
    scenario_results = scenario_df,
    summary = list(
      n_scenarios = 2L,
      n_baseline_available = 2L,
      median_runtime_improvement_ratio = 0.4165,
      mean_accuracy_delta = 0.05,
      mean_geodesic_norm_delta = -0.03
    )
  )

  out_dir <- tempfile("parity-report-")
  outputs <- env$parity_write_report(report, output_dir = out_dir)

  expect_true(file.exists(outputs$rds))
  expect_true(file.exists(outputs$markdown))

  md <- readLines(outputs$markdown, warn = FALSE)
  expect_true(any(grepl("Scenario Results", md, fixed = TRUE)))
  expect_true(any(grepl("Median runtime improvement ratio", md, fixed = TRUE)))
})

test_that("pyFM parity runner is skip-safe when dependencies are unavailable", {
  env <- load_parity_tools()
  cfg <- list(n = 40L, k = 8L, p = 6L, noise = 0.02, basis_noise = 0.0, descriptor_corruption = 0.0, eval_fraction = 1.0)
  res <- env$run_pyfm_parity(cfg = cfg, seed = 1L, n_runs = 1L, maxit = 30L, icp_nit = 1L)

  if (!identical(res$status, "ok")) {
    skip(paste("pyFM parity unavailable:", if (is.null(res$error)) "unknown" else res$error))
  }

  expect_true(is.finite(res$runtime_sec))
  expect_true(is.finite(res$accuracy))
  expect_true(is.finite(res$geodesic_normalized_mean))
})

test_that("compare_scenario computes expected delta fields and unavailable baseline notes", {
  env <- load_parity_tools()
  cfg <- list(n = 20L, k = 6L, p = 4L, noise = 0.01, basis_noise = 0.0, descriptor_corruption = 0.0, eval_fraction = 1.0)

  r_res <- list(
    runtime_sec = 0.2,
    runtime_match_sec = 0.1,
    runtime_refine_sec = 0.1,
    accuracy = 0.9,
    geodesic_normalized_mean = 0.05,
    kernel_backend = "cpp",
    kernel_backend_mixed = FALSE
  )
  py_ok <- list(status = "ok", runtime_sec = 0.4, accuracy = 0.8, geodesic_normalized_mean = 0.07)
  py_missing <- list(status = "error", error = "python_not_found")

  row_ok <- env$compare_scenario("easy", cfg, r_res, py_ok)
  expect_true(isTRUE(row_ok$baseline_available))
  expect_equal(row_ok$runtime_improvement_ratio, 0.5)
  expect_equal(row_ok$accuracy_delta, 0.1)
  expect_equal(row_ok$geodesic_norm_delta, -0.02)
  expect_identical(row_ok$note, "ok")

  row_missing <- env$compare_scenario("easy", cfg, r_res, py_missing)
  expect_false(isTRUE(row_missing$baseline_available))
  expect_true(grepl("pyFM unavailable", row_missing$note, fixed = TRUE))
})

test_that("compare_scenario marks runtime non-comparable on R backend", {
  env <- load_parity_tools()
  cfg <- list(n = 20L, k = 6L, p = 4L, noise = 0.01, basis_noise = 0.0, descriptor_corruption = 0.0, eval_fraction = 1.0)

  r_res <- list(
    runtime_sec = 0.2,
    runtime_match_sec = 0.1,
    runtime_refine_sec = 0.1,
    accuracy = 0.9,
    geodesic_normalized_mean = 0.05,
    kernel_backend = "r",
    kernel_backend_mixed = FALSE
  )
  py_ok <- list(status = "ok", runtime_sec = 0.4, accuracy = 0.8, geodesic_normalized_mean = 0.07)

  row <- env$compare_scenario("easy", cfg, r_res, py_ok)
  expect_true(isTRUE(row$baseline_available))
  expect_false(isTRUE(row$runtime_comparable))
  expect_true(is.na(row$runtime_improvement_ratio))
  expect_true(grepl("runtime not comparable", row$note, fixed = TRUE))
})

test_that("pyFM parity remains mathematically consistent across easy/noisy scenarios", {
  env <- load_parity_tools()
  scenarios <- env$default_scenarios()
  scenario_names <- c("easy", "noisy")

  for (i in seq_along(scenario_names)) {
    nm <- scenario_names[[i]]
    cfg <- scenarios[[nm]]
    seed_i <- 500L + i

    r_res <- env$run_r_parity(
      cfg = cfg,
      seed = seed_i,
      n_runs = 1L,
      cg_maxit = 6L,
      refine_nit = 2L,
      warmup = FALSE
    )
    py_res <- env$run_pyfm_parity(
      cfg = cfg,
      seed = seed_i,
      n_runs = 1L,
      maxit = 80L,
      icp_nit = 2L
    )

    if (!identical(py_res$status, "ok")) {
      skip(paste("pyFM parity unavailable:", if (is.null(py_res$error)) "unknown" else py_res$error))
    }

    row <- env$compare_scenario(nm, cfg, r_res, py_res)
    acc_tol <- if (nm == "easy") 0.10 else 0.15
    geod_tol <- if (nm == "easy") 0.08 else 0.15

    expect_true(isTRUE(row$baseline_available))
    expect_true(row$r_accuracy >= 0 && row$r_accuracy <= 1)
    expect_true(row$py_accuracy >= 0 && row$py_accuracy <= 1)
    expect_gte(row$py_accuracy, 1 / as.numeric(cfg$n))
    expect_gte(row$accuracy_delta, -acc_tol)
    expect_lte(row$geodesic_norm_delta, geod_tol)
  }
})
