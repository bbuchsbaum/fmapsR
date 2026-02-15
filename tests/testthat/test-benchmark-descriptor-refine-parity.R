probe_parse_kv <- function(lines) {
  out <- list()
  for (line in lines) {
    if (!nzchar(line) || !grepl("=", line, fixed = TRUE)) {
      next
    }
    kv <- strsplit(line, "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- paste(kv[-1], collapse = "=")
    out[[key]] <- val
  }
  out
}

probe_parse_num_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    numeric()
  } else {
    as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
  }
}

probe_parse_int_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    integer()
  } else {
    as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
  }
}

probe_parse_matrix <- function(values_csv, dims_csv) {
  dims <- probe_parse_int_csv(dims_csv)
  vals <- probe_parse_num_csv(values_csv)
  matrix(vals, nrow = dims[1], ncol = dims[2], byrow = TRUE)
}

run_descriptor_refine_probe <- function() {
  root_hints <- unique(Filter(nzchar, c(
    Sys.getenv("FMAPSR_SOURCE_ROOT", unset = ""),
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )))

  candidate_scripts <- unlist(lapply(root_hints, function(root) {
    c(
      file.path(root, "tools", "benchmark_pyfm_descriptor_refine.py"),
      file.path(root, "fmapsR", "tools", "benchmark_pyfm_descriptor_refine.py")
    )
  }), use.names = FALSE)

  script <- NULL
  for (p in candidate_scripts) {
    pp <- normalizePath(p, winslash = "/", mustWork = FALSE)
    if (file.exists(pp)) {
      script <- pp
      break
    }
  }
  if (is.null(script)) {
    testthat::skip("tools/benchmark_pyfm_descriptor_refine.py unavailable outside source checkout")
  }

  py_candidates <- unique(Filter(nzchar, c(
    Sys.getenv("PYFM_BENCH_PYTHON", unset = ""),
    file.path(".venv-bench", "bin", "python"),
    file.path("..", ".venv-bench", "bin", "python"),
    file.path("..", "..", ".venv-bench", "bin", "python"),
    Sys.which("python3"),
    Sys.which("python")
  )))

  py <- NULL
  for (cand in py_candidates) {
    if (file.exists(cand)) {
      py <- cand
      break
    }
    wh <- Sys.which(cand)
    if (nzchar(wh)) {
      py <- wh
      break
    }
  }
  if (is.null(py) || !nzchar(py)) {
    testthat::skip("Python interpreter unavailable for pyFM parity probe")
  }

  out <- tryCatch(
    suppressWarnings(system2(py, script, stdout = TRUE, stderr = TRUE)),
    error = function(e) sprintf("status=error\nerror=system2_failure:%s", conditionMessage(e))
  )

  kv <- probe_parse_kv(out)
  if (!identical(kv$status, "ok")) {
    err <- if (is.null(kv$error) || !nzchar(kv$error)) "unknown" else kv$error
    testthat::skip(paste("pyFM descriptor/refine probe unavailable:", err))
  }

  list(
    n_times = as.integer(kv$n_times),
    n_energies = as.integer(kv$n_energies),
    evals = probe_parse_num_csv(kv$evals),
    phi1 = probe_parse_matrix(kv$phi1, kv$phi1_dims),
    phi2 = probe_parse_matrix(kv$phi2, kv$phi2_dims),
    landmarks_1b = probe_parse_int_csv(kv$landmarks_1b),
    hks_time_grid = probe_parse_num_csv(kv$hks_time_grid),
    wks_energy_grid = probe_parse_num_csv(kv$wks_energy_grid),
    wks_sigma = as.numeric(kv$wks_sigma),
    hks = probe_parse_matrix(kv$hks, kv$hks_dims),
    hks_lm = probe_parse_matrix(kv$hks_lm, kv$hks_lm_dims),
    wks = probe_parse_matrix(kv$wks, kv$wks_dims),
    wks_lm = probe_parse_matrix(kv$wks_lm, kv$wks_lm_dims),
    c0 = probe_parse_matrix(kv$c0, kv$c0_dims),
    icp_once = probe_parse_matrix(kv$icp_once, kv$icp_once_dims),
    zoom_once = probe_parse_matrix(kv$zoom_once, kv$zoom_once_dims)
  )
}

test_that("descriptor defaults match pyFM auto signatures on a golden fixture", {
  probe <- run_descriptor_refine_probe()

  domain <- fm_domain_generic(
    n_samples = nrow(probe$phi1),
    basis = list(vectors = probe$phi1, values = probe$evals, k = ncol(probe$phi1))
  )

  expect_equal(
    infer_hks_time_grid(probe$evals, n_times = probe$n_times),
    probe$hks_time_grid,
    tolerance = 1e-12
  )
  expect_equal(
    infer_wks_auto_params(probe$evals, n_energies = probe$n_energies)$energy_grid,
    probe$wks_energy_grid,
    tolerance = 1e-12
  )
  expect_equal(
    infer_wks_auto_params(probe$evals, n_energies = probe$n_energies)$sigma,
    probe$wks_sigma,
    tolerance = 1e-12
  )

  hks_r <- fm_descriptor_hks(domain, n_times = probe$n_times, scaled = TRUE)
  hks_lm_r <- fm_descriptor_hks(domain, n_times = probe$n_times, landmarks = probe$landmarks_1b, scaled = TRUE)
  wks_r <- fm_descriptor_wks(domain, n_energies = probe$n_energies, scaled = TRUE)
  wks_lm_r <- fm_descriptor_wks(domain, n_energies = probe$n_energies, landmarks = probe$landmarks_1b, scaled = TRUE)

  expect_equal(hks_r, probe$hks, tolerance = 1e-8)
  expect_equal(hks_lm_r, probe$hks_lm, tolerance = 1e-8)
  expect_equal(wks_r, probe$wks, tolerance = 1e-8)
  expect_equal(wks_lm_r, probe$wks_lm, tolerance = 1e-8)
})

test_that("single-step refinement operators match pyFM golden outputs", {
  probe <- run_descriptor_refine_probe()

  source <- fm_domain_generic(
    n_samples = nrow(probe$phi1),
    basis = list(vectors = probe$phi1, values = probe$evals, k = ncol(probe$phi1))
  )
  target <- fm_domain_generic(
    n_samples = nrow(probe$phi2),
    basis = list(vectors = probe$phi2, values = probe$evals, k = ncol(probe$phi2))
  )

  sub <- list(source = NULL, target = NULL)
  icp_r <- icp_refine_once(probe$c0, source, target, use_adj = FALSE, subsample = sub)
  zoom_r <- zoomout_refine_once(probe$c0, source, target, step = c(1, 2), subsample = sub)

  expect_equal(icp_r, probe$icp_once, tolerance = 1e-7)
  expect_equal(zoom_r, probe$zoom_once, tolerance = 1e-7)
})
