#!/usr/bin/env Rscript

parse_kv <- function(lines) {
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

or_else <- function(x, y) {
  if (is.null(x) || !nzchar(x)) y else x
}

numeric_or_na <- function(x) {
  if (is.null(x)) return(NA_real_)
  as.numeric(x)
}

parse_args <- function(args) {
  out <- list(
    r_runs = 3L,
    py_runs = 3L,
    target_ratio = 0.30,
    stability_margin = 0.02,
    stability_boot = 500L,
    stability_seed = 2026L,
    output_dir = "benchmarks/pairwise"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NULL

    if (key == "--r-runs" && !is.null(val)) {
      out$r_runs <- as.integer(val)
      i <- i + 2L
    } else if (key == "--py-runs" && !is.null(val)) {
      out$py_runs <- as.integer(val)
      i <- i + 2L
    } else if (key == "--target-ratio" && !is.null(val)) {
      out$target_ratio <- as.numeric(val)
      i <- i + 2L
    } else if (key == "--stability-margin" && !is.null(val)) {
      out$stability_margin <- as.numeric(val)
      i <- i + 2L
    } else if (key == "--stability-boot" && !is.null(val)) {
      out$stability_boot <- as.integer(val)
      i <- i + 2L
    } else if (key == "--stability-seed" && !is.null(val)) {
      out$stability_seed <- as.integer(val)
      i <- i + 2L
    } else if (key == "--output-dir" && !is.null(val)) {
      out$output_dir <- val
      i <- i + 2L
    } else {
      stop(sprintf("Unknown or malformed argument: %s", key), call. = FALSE)
    }
  }

  if (!is.finite(out$target_ratio) || out$target_ratio <= 0 || out$target_ratio >= 1) {
    stop("`--target-ratio` must be in (0, 1)", call. = FALSE)
  }
  if (!is.finite(out$stability_margin) || out$stability_margin < 0 || out$stability_margin >= out$target_ratio) {
    stop("`--stability-margin` must be >= 0 and smaller than target ratio", call. = FALSE)
  }
  if (!is.finite(out$stability_boot) || out$stability_boot < 50L) {
    stop("`--stability-boot` must be >= 50", call. = FALSE)
  }

  out
}

prepare_stability_samples <- function(x, drop_max = TRUE) {
  vals <- as.numeric(x)
  vals <- vals[is.finite(vals) & vals > 0]

  n_raw <- length(vals)
  n_dropped_high <- 0L

  if (isTRUE(drop_max) && length(vals) >= 3L) {
    vals <- vals[-which.max(vals)]
    n_dropped_high <- 1L
  }

  list(
    values = vals,
    n_raw = n_raw,
    n_used = length(vals),
    n_dropped_high = n_dropped_high
  )
}

compute_improvement_samples <- function(r_runtimes, py_runtimes, n_boot = 500L, seed = 2026L) {
  r <- as.numeric(r_runtimes)
  py <- as.numeric(py_runtimes)

  r <- r[is.finite(r) & r > 0]
  py <- py[is.finite(py) & py > 0]

  if (length(r) < 1L || length(py) < 1L) {
    return(numeric())
  }

  n_boot <- max(50L, as.integer(n_boot))
  set.seed(as.integer(seed))

  out <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    r_draw <- sample(r, size = length(r), replace = TRUE)
    py_draw <- sample(py, size = length(py), replace = TRUE)

    py_med <- stats::median(py_draw, na.rm = TRUE)
    r_med <- stats::median(r_draw, na.rm = TRUE)

    out[[i]] <- if (is.finite(py_med) && py_med > 0) {
      (py_med - r_med) / py_med
    } else {
      NA_real_
    }
  }

  out[is.finite(out)]
}

assess_runtime_target <- function(improvement, improvement_samples, target_ratio = 0.30, stability_margin = 0.02) {
  if (!is.finite(improvement)) {
    return(list(
      target_met_raw = NA,
      target_met_stable = NA,
      classification = "unavailable",
      band_low = target_ratio - stability_margin,
      band_high = target_ratio + stability_margin,
      sample_n = 0L,
      sample_median = NA_real_,
      sample_q25 = NA_real_,
      sample_q75 = NA_real_,
      note = "Baseline unavailable; stability assessment skipped"
    ))
  }

  band_low <- target_ratio - stability_margin
  band_high <- target_ratio + stability_margin

  vals <- as.numeric(improvement_samples)
  vals <- vals[is.finite(vals)]

  sample_n <- length(vals)
  sample_median <- if (sample_n > 0L) stats::median(vals, na.rm = TRUE) else NA_real_
  sample_q25 <- if (sample_n > 0L) as.numeric(stats::quantile(vals, probs = 0.25, names = FALSE, na.rm = TRUE)) else NA_real_
  sample_q75 <- if (sample_n > 0L) as.numeric(stats::quantile(vals, probs = 0.75, names = FALSE, na.rm = TRUE)) else NA_real_

  raw <- improvement >= target_ratio

  if (improvement >= band_high) {
    list(
      target_met_raw = raw,
      target_met_stable = TRUE,
      classification = "pass_strong",
      band_low = band_low,
      band_high = band_high,
      sample_n = sample_n,
      sample_median = sample_median,
      sample_q25 = sample_q25,
      sample_q75 = sample_q75,
      note = "Runtime target met outside jitter band"
    )
  } else if (improvement < band_low) {
    list(
      target_met_raw = raw,
      target_met_stable = FALSE,
      classification = "fail_strong",
      band_low = band_low,
      band_high = band_high,
      sample_n = sample_n,
      sample_median = sample_median,
      sample_q25 = sample_q25,
      sample_q75 = sample_q75,
      note = "Runtime target below lower jitter band"
    )
  } else {
    stable <- if (sample_n > 0L) sample_q25 >= band_low else raw
    class <- if (stable) "pass_borderline_stable" else "fail_borderline_unstable"
    note <- if (stable) {
      "Borderline but stable in bootstrap spread"
    } else {
      "Borderline and unstable in bootstrap spread"
    }

    list(
      target_met_raw = raw,
      target_met_stable = stable,
      classification = class,
      band_low = band_low,
      band_high = band_high,
      sample_n = sample_n,
      sample_median = sample_median,
      sample_q25 = sample_q25,
      sample_q75 = sample_q75,
      note = note
    )
  }
}

run_r_pipeline_once <- function(
  seed = 42,
  n = 160,
  k = 32,
  p = 24,
  maxit = 120,
  optimizer = "cg",
  kernel_backend = "auto"
) {
  set.seed(seed)

  make_domain <- function() {
    vals <- seq(0.5, 0.5 + n - 1)
    op <- diag(vals)
    dom <- fm_domain_generic(n_samples = n, operator = op)
    fm_basis(dom, k = k, solver = "rspectra", cache = FALSE, seed = seed)
  }

  source <- make_domain()
  target <- make_domain()

  src_desc <- matrix(rnorm(n * p), nrow = n, ncol = p)
  tgt_desc <- src_desc + matrix(rnorm(n * p, sd = 0.01), nrow = n, ncol = p)

  t_match <- system.time({
    fit <- fm_match(
      source = source,
      target = target,
      descriptors = list(source = src_desc, target = tgt_desc),
      penalties = list(descr = 1e-1, lap = 1e-3, comm = 1),
      init = "identity",
      maxit = maxit,
      optimizer = optimizer,
      kernel_backend = kernel_backend
    )
  })[["elapsed"]]

  t_refine <- system.time({
    fit_refined <- fm_refine(fit, method = "icp", nit = 5, use_adj = FALSE)
  })[["elapsed"]]

  list(
    runtime_sec = t_match + t_refine,
    runtime_match_sec = t_match,
    runtime_refine_sec = t_refine,
    optimizer = fit$diagnostics$optimizer,
    kernel_backend = fit$diagnostics$kernel_backend,
    objective = fit$diagnostics$total_objective,
    refinement_last_delta = tail(fit_refined$diagnostics$refinement$deltas, 1),
    convergence = fit$diagnostics$convergence
  )
}

run_r_benchmark <- function(
  n_runs = 3L,
  optimizer = "cg",
  kernel_backend = "auto",
  seed_offset = 40L,
  measure_memory = TRUE
) {
  invisible(run_r_pipeline_once(seed = 13, optimizer = optimizer, kernel_backend = kernel_backend))

  runs <- lapply(seq_len(n_runs), function(i) {
    run_r_pipeline_once(
      seed = seed_offset + i,
      optimizer = optimizer,
      kernel_backend = kernel_backend
    )
  })

  runtimes <- vapply(runs, `[[`, numeric(1), "runtime_sec")
  objectives <- vapply(runs, `[[`, numeric(1), "objective")

  mem_bytes <- NA_real_
  if (isTRUE(measure_memory) && requireNamespace("bench", quietly = TRUE)) {
    b <- bench::mark(
      run_r_pipeline_once(seed = seed_offset + 37L, optimizer = optimizer, kernel_backend = kernel_backend),
      iterations = 1,
      check = FALSE
    )
    mem_bytes <- as.numeric(b$mem_alloc[[1]])
  }

  list(
    runs = runs,
    median_runtime_sec = median(runtimes),
    mean_runtime_sec = mean(runtimes),
    optimizer = runs[[1]]$optimizer,
    kernel_backend = runs[[1]]$kernel_backend,
    median_objective = median(objectives),
    mem_bytes = mem_bytes
  )
}

run_pyfm_baseline_once <- function(py, n, k, p, maxit, icp_nit, seed) {
  cmd <- c(
    "tools/benchmark_pyfm.py",
    "--n", as.character(n),
    "--k1", as.character(k),
    "--k2", as.character(k),
    "--p", as.character(p),
    "--maxit", as.character(maxit),
    "--icp-nit", as.character(icp_nit),
    "--seed", as.character(seed)
  )

  out <- tryCatch(
    system2(py, cmd, stdout = TRUE, stderr = TRUE),
    error = function(e) paste0("status=error\nerror=system2_failure:", conditionMessage(e))
  )

  kv <- parse_kv(out)
  if (is.null(kv$status)) {
    kv$status <- "error"
    kv$error <- "malformed_output"
    kv$raw <- paste(out, collapse = " | ")
  }
  kv
}

run_pyfm_baseline <- function(
  n = 160L,
  k = 32L,
  p = 24L,
  maxit = 120L,
  icp_nit = 5L,
  seed = 42L,
  n_runs = 3L
) {
  candidates <- c(
    Sys.getenv("PYFM_BENCH_PYTHON", unset = ""),
    file.path(".venv-bench", "bin", "python"),
    Sys.which("python3"),
    Sys.which("python")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  py <- ""
  for (cand in candidates) {
    if (file.exists(cand) || nzchar(Sys.which(cand))) {
      py <- cand
      break
    }
  }

  if (!nzchar(py)) {
    return(list(status = "error", error = "python_not_found"))
  }

  n_runs <- max(1L, as.integer(n_runs))
  runs <- lapply(seq_len(n_runs), function(i) {
    run_pyfm_baseline_once(
      py = py,
      n = n,
      k = k,
      p = p,
      maxit = maxit,
      icp_nit = icp_nit,
      seed = seed + i - 1L
    )
  })

  statuses <- vapply(runs, function(x) or_else(x$status, "error"), character(1))
  if (any(statuses != "ok")) {
    first_err <- runs[[which(statuses != "ok")[1]]]
    return(list(
      status = "error",
      error = or_else(first_err$error, "pyfm_run_failed"),
      py_runs = runs
    ))
  }

  runtimes <- vapply(runs, function(x) numeric_or_na(x$runtime_sec), numeric(1))
  objectives <- vapply(runs, function(x) numeric_or_na(x$objective), numeric(1))
  runtime_opt <- vapply(runs, function(x) numeric_or_na(x$runtime_opt_sec), numeric(1))
  runtime_icp <- vapply(runs, function(x) numeric_or_na(x$runtime_icp_sec), numeric(1))

  list(
    status = "ok",
    mode = or_else(runs[[1]]$mode, "unknown"),
    python = if (is.null(runs[[1]]$python)) NA_character_ else runs[[1]]$python,
    n_runs = n_runs,
    runtime_sec = median(runtimes, na.rm = TRUE),
    runtime_sec_mean = mean(runtimes, na.rm = TRUE),
    runtime_opt_sec = median(runtime_opt, na.rm = TRUE),
    runtime_icp_sec = median(runtime_icp, na.rm = TRUE),
    objective = median(objectives, na.rm = TRUE),
    py_runs = runs
  )
}

write_report <- function(report, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  rds_path <- file.path(output_dir, paste0("pairwise-benchmark-", stamp, ".rds"))
  md_path <- file.path(output_dir, "pairwise-benchmark-latest.md")

  saveRDS(report, rds_path)

  lines <- c(
    "# Pairwise Benchmark",
    "",
    sprintf("- Timestamp: %s", report$metadata$timestamp),
    sprintf("- R version: %s", report$metadata$r_version),
    sprintf("- Platform: %s", report$metadata$platform),
    sprintf("- Target ratio: %s", as.character(report$metadata$target_ratio)),
    sprintf("- Stability margin: %s", as.character(report$metadata$stability_margin)),
    sprintf("- Stability bootstrap samples: %s", as.character(report$metadata$stability_boot)),
    "",
    "## R Pipeline",
    sprintf("- Median runtime (s): %.6f", report$r_pipeline$median_runtime_sec),
    sprintf("- Mean runtime (s): %.6f", report$r_pipeline$mean_runtime_sec),
    sprintf("- Optimizer: %s", report$r_pipeline$optimizer),
    sprintf("- Kernel backend: %s", report$r_pipeline$kernel_backend),
    sprintf("- Median objective: %.6f", report$r_pipeline$median_objective),
    sprintf("- Memory (bytes, optional): %s", as.character(report$r_pipeline$mem_bytes)),
    "",
    "## Backend Comparison",
    sprintf("- Baseline backend: %s", report$backend_comparison$baseline_backend),
    sprintf("- Active backend: %s", report$backend_comparison$active_backend),
    sprintf("- Baseline median runtime (s): %.6f", report$backend_comparison$baseline_runtime_sec),
    sprintf("- Active median runtime (s): %.6f", report$backend_comparison$active_runtime_sec),
    sprintf("- Active vs baseline speedup: %s", as.character(report$backend_comparison$active_vs_baseline_speedup)),
    sprintf("- Objective median gap: %s", as.character(report$backend_comparison$objective_gap)),
    "",
    "## pyFM Baseline",
    sprintf("- Status: %s", report$pyfm_baseline$status),
    sprintf("- Runs: %s", as.character(if (is.null(report$pyfm_baseline$n_runs)) NA_integer_ else report$pyfm_baseline$n_runs)),
    sprintf("- Runtime (s): %s", as.character(report$comparison$pyfm_runtime_sec)),
    sprintf("- Objective: %s", as.character(report$comparison$pyfm_objective)),
    "",
    "## Comparison",
    sprintf("- Runtime improvement ratio: %s", as.character(report$comparison$runtime_improvement_ratio)),
    sprintf("- Runtime improvement sample median: %s", as.character(report$comparison$runtime_improvement_sample_median)),
    sprintf("- Runtime improvement sample q25/q75: %s / %s", as.character(report$comparison$runtime_improvement_sample_q25), as.character(report$comparison$runtime_improvement_sample_q75)),
    sprintf("- Runtime improvement bootstrap samples: %s", as.character(report$comparison$runtime_improvement_boot_n)),
    sprintf("- Stability sample usage (R raw/used, pyFM raw/used): %s/%s, %s/%s",
      as.character(report$comparison$stability_r_runtime_n_raw),
      as.character(report$comparison$stability_r_runtime_n_used),
      as.character(report$comparison$stability_py_runtime_n_raw),
      as.character(report$comparison$stability_py_runtime_n_used)
    ),
    sprintf("- Target ratio: %s", as.character(report$comparison$target_ratio)),
    sprintf("- Stability band: [%s, %s]", as.character(report$comparison$target_band_low), as.character(report$comparison$target_band_high)),
    sprintf("- Target met (raw >= target): %s", as.character(report$comparison$target_met_raw)),
    sprintf("- Stable target met: %s", as.character(report$comparison$target_met_stable)),
    sprintf("- Stability decision: %s", as.character(report$comparison$stability_decision)),
    sprintf("- Baseline available: %s", as.character(report$comparison$baseline_available)),
    sprintf("- Notes: %s", report$comparison$notes),
    ""
  )

  writeLines(lines, md_path)

  list(rds = rds_path, markdown = md_path)
}

main <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    stop("pkgload is required to run benchmark script from source", call. = FALSE)
  }

  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  cpp_available <- exists("fm_match_solve_cg_cpp", mode = "function")
  active_backend <- if (cpp_available) "cpp" else "r"

  r_pipeline <- run_r_benchmark(
    n_runs = opts$r_runs,
    optimizer = "cg",
    kernel_backend = active_backend,
    seed_offset = 40L,
    measure_memory = TRUE
  )
  r_baseline <- run_r_benchmark(
    n_runs = opts$r_runs,
    optimizer = "cg",
    kernel_backend = "r",
    seed_offset = 40L,
    measure_memory = FALSE
  )
  pyfm <- run_pyfm_baseline(
    n = 160L,
    k = 32L,
    p = 24L,
    maxit = 120L,
    icp_nit = 5L,
    seed = 42L,
    n_runs = opts$py_runs
  )

  py_runtime <- numeric_or_na(pyfm$runtime_sec)
  py_objective <- numeric_or_na(pyfm$objective)

  baseline_available <- identical(pyfm$status, "ok") && is.finite(py_runtime)

  improvement <- if (baseline_available) {
    (py_runtime - r_pipeline$median_runtime_sec) / py_runtime
  } else {
    NA_real_
  }

  r_samples_raw <- vapply(r_pipeline$runs, `[[`, numeric(1), "runtime_sec")
  py_samples_raw <- if (baseline_available && !is.null(pyfm$py_runs)) {
    vapply(pyfm$py_runs, function(x) numeric_or_na(x$runtime_sec), numeric(1))
  } else {
    numeric()
  }

  r_samples_meta <- prepare_stability_samples(r_samples_raw, drop_max = TRUE)
  py_samples_meta <- prepare_stability_samples(py_samples_raw, drop_max = TRUE)

  improvement_samples <- compute_improvement_samples(
    r_samples_meta$values,
    py_samples_meta$values,
    n_boot = opts$stability_boot,
    seed = opts$stability_seed
  )

  assess <- assess_runtime_target(
    improvement = improvement,
    improvement_samples = improvement_samples,
    target_ratio = opts$target_ratio,
    stability_margin = opts$stability_margin
  )

  notes <- if (!baseline_available) {
    paste("pyFM baseline unavailable:", or_else(pyfm$error, "unknown"))
  } else {
    assess$note
  }

  report <- list(
    metadata = list(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      r_version = as.character(getRversion()),
      platform = R.version$platform,
      git_commit = tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)[1], error = function(e) NA_character_),
      target_ratio = opts$target_ratio,
      stability_margin = opts$stability_margin,
      stability_boot = opts$stability_boot
    ),
    r_pipeline = r_pipeline,
    r_baseline = r_baseline,
    pyfm_baseline = pyfm,
    backend_comparison = list(
      baseline_backend = r_baseline$kernel_backend,
      active_backend = r_pipeline$kernel_backend,
      baseline_runtime_sec = r_baseline$median_runtime_sec,
      active_runtime_sec = r_pipeline$median_runtime_sec,
      active_vs_baseline_speedup = if (is.finite(r_baseline$median_runtime_sec) && r_baseline$median_runtime_sec > 0) {
        (r_baseline$median_runtime_sec - r_pipeline$median_runtime_sec) / r_baseline$median_runtime_sec
      } else {
        NA_real_
      },
      objective_gap = abs(r_pipeline$median_objective - r_baseline$median_objective)
    ),
    comparison = list(
      baseline_available = baseline_available,
      pyfm_runtime_sec = py_runtime,
      pyfm_objective = py_objective,
      runtime_improvement_ratio = improvement,
      runtime_improvement_sample_median = assess$sample_median,
      runtime_improvement_sample_q25 = assess$sample_q25,
      runtime_improvement_sample_q75 = assess$sample_q75,
      runtime_improvement_boot_n = assess$sample_n,
      stability_r_runtime_n_raw = r_samples_meta$n_raw,
      stability_r_runtime_n_used = r_samples_meta$n_used,
      stability_py_runtime_n_raw = py_samples_meta$n_raw,
      stability_py_runtime_n_used = py_samples_meta$n_used,
      target_ratio = opts$target_ratio,
      target_band_low = assess$band_low,
      target_band_high = assess$band_high,
      target_met_raw = assess$target_met_raw,
      target_met_stable = assess$target_met_stable,
      stability_decision = assess$classification,
      notes = notes
    )
  )

  outputs <- write_report(report, output_dir = opts$output_dir)

  cat("Benchmark complete\n")
  cat("RDS:", outputs$rds, "\n")
  cat("Summary:", outputs$markdown, "\n")
  cat("Target met (raw):", as.character(report$comparison$target_met_raw), "\n")
  cat("Stable target met:", as.character(report$comparison$target_met_stable), "\n")
  cat("Stability decision:", report$comparison$stability_decision, "\n")
  cat("Notes:", report$comparison$notes, "\n")
}

if (sys.nframe() == 0) {
  main()
}
