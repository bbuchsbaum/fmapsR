#!/usr/bin/env Rscript

parity_parse_kv <- function(lines) {
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

parity_numeric_or_na <- function(x) {
  if (is.null(x)) return(NA_real_)
  as.numeric(x)
}

parity_or_else <- function(x, y) {
  if (is.null(x) || !nzchar(x)) y else x
}

parity_tool_path <- function(filename) {
  candidates <- c(
    file.path("tools", filename),
    file.path("..", "..", "tools", filename)
  )
  for (p in candidates) {
    pp <- normalizePath(p, winslash = "/", mustWork = FALSE)
    if (file.exists(pp)) {
      return(pp)
    }
  }
  candidates[[1]]
}

source(parity_tool_path("benchmark_config.R"), local = environment())

parity_python_candidates <- function() {
  raw <- c(
    Sys.getenv("PYFM_BENCH_PYTHON", unset = ""),
    file.path(".venv-bench", "bin", "python"),
    file.path("..", ".venv-bench", "bin", "python"),
    file.path("..", "..", ".venv-bench", "bin", "python"),
    Sys.which("python3"),
    Sys.which("python")
  )

  out <- character()
  for (cand in unique(raw[nzchar(raw)])) {
    if (file.exists(cand)) {
      out <- c(out, cand)
      next
    }

    wh <- Sys.which(cand)
    if (nzchar(wh)) {
      out <- c(out, wh)
    }
  }
  unique(out)
}

parity_find_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- if (length(file_arg) > 0L) {
    normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
  } else {
    ""
  }

  candidates <- unique(c(
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE),
    if (nzchar(script_path)) dirname(dirname(script_path)) else ""
  ))
  candidates <- candidates[nzchar(candidates)]

  for (root in candidates) {
    if (file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "R"))) {
      return(root)
    }
  }

  stop("Could not locate package root (expected DESCRIPTION + R directory)", call. = FALSE)
}

parity_source_package_r_files <- function(root, target_env = globalenv()) {
  r_dir <- file.path(root, "R")
  files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  files <- files[basename(files) != "RcppExports.R"]

  for (path in files) {
    source(path, local = target_env)
  }
}

parity_load_fmaps_package <- function(root = parity_find_repo_root(), target_env = globalenv()) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    load_attempt <- tryCatch(
      {
        pkgload::load_all(root, quiet = TRUE)
        list(ok = TRUE, error = NULL)
      },
      error = function(e) list(ok = FALSE, error = conditionMessage(e))
    )

    if (isTRUE(load_attempt$ok)) {
      return(list(mode = "pkgload", load_error = NULL))
    }
  } else {
    load_attempt <- list(ok = FALSE, error = "pkgload_not_installed")
  }

  parity_source_package_r_files(root = root, target_env = target_env)
  list(mode = "source_fallback", load_error = load_attempt$error)
}

parse_args <- function(args) {
  seed_cfg <- benchmark_seed_registry()

  out <- list(
    scenarios = names(benchmark_parity_scenarios()),
    r_runs = 3L,
    py_runs = 3L,
    r_cg_maxit = 3L,
    r_cg_tol = 1e-6,
    r_refine_nit = 2L,
    r_kernel_backend = "auto",
    py_maxit = 120L,
    py_icp_nit = 3L,
    warmup = TRUE,
    seed = seed_cfg$parity_seed,
    output_dir = "benchmarks/parity"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NULL

    if (key == "--scenarios" && !is.null(val)) {
      out$scenarios <- strsplit(val, ",", fixed = TRUE)[[1]]
      i <- i + 2L
    } else if (key == "--r-runs" && !is.null(val)) {
      out$r_runs <- as.integer(val)
      i <- i + 2L
    } else if (key == "--py-runs" && !is.null(val)) {
      out$py_runs <- as.integer(val)
      i <- i + 2L
    } else if (key == "--r-cg-maxit" && !is.null(val)) {
      out$r_cg_maxit <- as.integer(val)
      i <- i + 2L
    } else if (key == "--r-cg-tol" && !is.null(val)) {
      out$r_cg_tol <- as.numeric(val)
      i <- i + 2L
    } else if (key == "--r-refine-nit" && !is.null(val)) {
      out$r_refine_nit <- as.integer(val)
      i <- i + 2L
    } else if (key == "--r-kernel-backend" && !is.null(val)) {
      out$r_kernel_backend <- val
      i <- i + 2L
    } else if (key == "--py-maxit" && !is.null(val)) {
      out$py_maxit <- as.integer(val)
      i <- i + 2L
    } else if (key == "--py-icp-nit" && !is.null(val)) {
      out$py_icp_nit <- as.integer(val)
      i <- i + 2L
    } else if (key == "--seed" && !is.null(val)) {
      out$seed <- as.integer(val)
      i <- i + 2L
    } else if (key == "--output-dir" && !is.null(val)) {
      out$output_dir <- val
      i <- i + 2L
    } else if (key == "--no-warmup") {
      out$warmup <- FALSE
      i <- i + 1L
    } else {
      stop(sprintf("Unknown or malformed argument: %s", key), call. = FALSE)
    }
  }

  if (!is.finite(out$r_cg_tol) || out$r_cg_tol <= 0) {
    stop("`--r-cg-tol` must be a positive scalar", call. = FALSE)
  }
  out$r_kernel_backend <- match.arg(out$r_kernel_backend, c("auto", "r", "cpp"))
  out
}

default_scenarios <- function() {
  benchmark_parity_scenarios()
}

parity_make_problem <- function(cfg, seed = 42L) {
  set.seed(seed)

  n <- as.integer(cfg$n)
  k <- as.integer(cfg$k)
  p <- as.integer(cfg$p)

  theta <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)]
  coords <- cbind(cos(theta), sin(theta))

  shift <- max(1L, floor(0.15 * n))
  truth <- ((seq_len(n) + shift - 1L) %% n) + 1L # target -> source

  P12 <- matrix(0, nrow = n, ncol = n)
  P12[cbind(seq_len(n), truth)] <- 1

  q <- qr.Q(qr(matrix(rnorm(n * k), nrow = n, ncol = k)))
  phi1 <- q[, seq_len(k), drop = FALSE]
  phi2 <- P12 %*% phi1

  if (cfg$basis_noise > 0) {
    q2 <- qr.Q(qr(phi2 + cfg$basis_noise * matrix(rnorm(n * k), nrow = n, ncol = k)))
    phi2 <- q2[, seq_len(k), drop = FALSE]
  }

  vals <- seq_len(k)
  src_desc <- matrix(rnorm(n * p), nrow = n, ncol = p)
  tgt_desc <- P12 %*% src_desc + matrix(rnorm(n * p, sd = cfg$noise), nrow = n, ncol = p)

  if (cfg$descriptor_corruption > 0) {
    n_bad <- max(1L, as.integer(round(cfg$descriptor_corruption * n)))
    bad_idx <- sample.int(n, size = n_bad, replace = FALSE)
    tgt_desc[bad_idx, ] <- matrix(rnorm(n_bad * p), nrow = n_bad, ncol = p)
  }

  eval_n <- max(1L, as.integer(round(cfg$eval_fraction * n)))
  eval_idx <- sort(sample.int(n, size = eval_n, replace = FALSE))

  source <- fm_domain_generic(
    n_samples = n,
    data = coords,
    basis = list(vectors = phi1, values = vals, k = k)
  )
  target <- fm_domain_generic(
    n_samples = n,
    data = coords[truth, , drop = FALSE],
    basis = list(vectors = phi2, values = vals, k = k)
  )

  list(
    source = source,
    target = target,
    descriptors = list(source = src_desc, target = tgt_desc),
    truth = truth,
    eval_idx = eval_idx,
    source_distance = as.matrix(stats::dist(coords))
  )
}

run_r_parity_once <- function(
  cfg,
  seed,
  cg_maxit = 3L,
  cg_tol = 1e-6,
  refine_nit = 2L,
  kernel_backend = c("auto", "r", "cpp")
) {
  kernel_backend <- match.arg(kernel_backend)
  prob <- parity_make_problem(cfg = cfg, seed = seed)

  t0 <- proc.time()[["elapsed"]]
  fit <- fm_match(
    prob$source,
    prob$target,
    descriptors = prob$descriptors,
    penalties = list(descr = 1e-1, lap = 1e-3, comm = 0.5),
    optimizer = "cg",
    cg_maxit = cg_maxit,
    cg_tol = cg_tol,
    kernel_backend = kernel_backend
  )
  t_match <- proc.time()[["elapsed"]] - t0
  match_objective <- as.numeric(fit$diagnostics$total_objective)

  t1 <- proc.time()[["elapsed"]]
  fit <- fm_refine(fit, method = "icp", nit = refine_nit)
  t_refine <- proc.time()[["elapsed"]] - t1
  t_total <- t_match + t_refine

  C <- as.matrix(fit$C)
  gram <- t(C) %*% C
  I <- diag(1, nrow = nrow(gram), ncol = ncol(gram))

  p2p <- as_p2p(fit)
  metrics <- fm_fit_metrics(
    fit,
    truth = prob$truth,
    source_distance = prob$source_distance,
    target_adjacency = NULL,
    normalize = TRUE
  )

  idx <- prob$eval_idx
  list(
    runtime_sec = t_total,
    runtime_match_sec = t_match,
    runtime_refine_sec = t_refine,
    objective = match_objective,
    map_fro_norm = norm(C, type = "F"),
    map_orth_resid = norm(gram - I, type = "F"),
    accuracy = mean(p2p[idx] == prob$truth[idx]),
    geodesic_mean = mean(metrics$geodesic$per_point[idx], na.rm = TRUE),
    geodesic_normalized_mean = mean(metrics$geodesic$per_point[idx], na.rm = TRUE) / metrics$geodesic$scale,
    coverage_ratio = metrics$coverage$ratio,
    eval_n = length(idx),
    kernel_backend_requested = kernel_backend,
    kernel_backend_used = fit$diagnostics$kernel_backend,
    cg_maxit = cg_maxit,
    cg_tol = cg_tol,
    refine_nit = refine_nit
  )
}

run_r_parity <- function(
  cfg,
  seed = 42L,
  n_runs = 3L,
  cg_maxit = 3L,
  cg_tol = 1e-6,
  refine_nit = 2L,
  kernel_backend = c("auto", "r", "cpp"),
  warmup = TRUE
) {
  kernel_backend <- match.arg(kernel_backend)
  n_runs <- max(1L, as.integer(n_runs))

  warmup_runtime <- NA_real_
  if (isTRUE(warmup)) {
    warm <- run_r_parity_once(
      cfg = cfg,
      seed = seed + 10000L,
      cg_maxit = cg_maxit,
      cg_tol = cg_tol,
      refine_nit = refine_nit,
      kernel_backend = kernel_backend
    )
    warmup_runtime <- warm$runtime_sec
  }

  runs <- lapply(seq_len(n_runs), function(i) {
    run_r_parity_once(
      cfg = cfg,
      seed = seed + i - 1L,
      cg_maxit = cg_maxit,
      cg_tol = cg_tol,
      refine_nit = refine_nit,
      kernel_backend = kernel_backend
    )
  })

  med <- function(key) stats::median(vapply(runs, `[[`, numeric(1), key), na.rm = TRUE)

  list(
    status = "ok",
    n_runs = n_runs,
    warmup_runtime_sec = warmup_runtime,
    runtime_sec = med("runtime_sec"),
    runtime_match_sec = med("runtime_match_sec"),
    runtime_refine_sec = med("runtime_refine_sec"),
    objective = med("objective"),
    map_fro_norm = med("map_fro_norm"),
    map_orth_resid = med("map_orth_resid"),
    accuracy = med("accuracy"),
    geodesic_mean = med("geodesic_mean"),
    geodesic_normalized_mean = med("geodesic_normalized_mean"),
    coverage_ratio = med("coverage_ratio"),
    eval_n = as.integer(round(med("eval_n"))),
    cg_maxit = as.integer(round(med("cg_maxit"))),
    cg_tol = med("cg_tol"),
    refine_nit = as.integer(round(med("refine_nit"))),
    kernel_backend_requested = kernel_backend,
    kernel_backend = runs[[1]]$kernel_backend_used,
    kernel_backend_mixed = length(unique(vapply(runs, `[[`, character(1), "kernel_backend_used"))) > 1L,
    runs = runs
  )
}

run_pyfm_parity_once <- function(py, cfg, seed, maxit = 120L, icp_nit = 3L) {
  cmd <- c(
    parity_tool_path("benchmark_pyfm_quality.py"),
    "--n", as.character(cfg$n),
    "--k", as.character(cfg$k),
    "--p", as.character(cfg$p),
    "--noise", as.character(cfg$noise),
    "--basis-noise", as.character(cfg$basis_noise),
    "--descriptor-corruption", as.character(cfg$descriptor_corruption),
    "--eval-fraction", as.character(cfg$eval_fraction),
    "--seed", as.character(seed),
    "--maxit", as.character(maxit),
    "--icp-nit", as.character(icp_nit)
  )

  out <- tryCatch(
    suppressWarnings(system2(py, cmd, stdout = TRUE, stderr = TRUE)),
    error = function(e) paste0("status=error\nerror=system2_failure:", conditionMessage(e))
  )

  kv <- parity_parse_kv(out)
  if (is.null(kv$status)) {
    kv$status <- "error"
    kv$error <- "malformed_output"
    kv$raw <- paste(out, collapse = " | ")
  }
  kv
}

run_pyfm_parity <- function(cfg, seed = 42L, n_runs = 3L, maxit = 120L, icp_nit = 3L) {
  candidates <- parity_python_candidates()
  py <- if (length(candidates) > 0L) candidates[[1]] else ""

  if (!nzchar(py)) {
    return(list(status = "error", error = "python_not_found"))
  }

  n_runs <- max(1L, as.integer(n_runs))
  runs <- lapply(seq_len(n_runs), function(i) {
    run_pyfm_parity_once(py = py, cfg = cfg, seed = seed + i - 1L, maxit = maxit, icp_nit = icp_nit)
  })

  statuses <- vapply(runs, function(x) if (is.null(x$status)) "error" else x$status, character(1))
  if (any(statuses != "ok")) {
    first_err <- runs[[which(statuses != "ok")[1]]]
    return(list(
      status = "error",
      error = if (is.null(first_err$error)) "pyfm_run_failed" else first_err$error,
      py_runs = runs
    ))
  }

  med <- function(key) stats::median(vapply(runs, function(x) parity_numeric_or_na(x[[key]]), numeric(1)), na.rm = TRUE)

  list(
    status = "ok",
    n_runs = n_runs,
    runtime_sec = med("runtime_sec"),
    objective = med("objective"),
    map_fro_norm = med("map_fro_norm"),
    map_orth_resid = med("map_orth_resid"),
    accuracy = med("accuracy"),
    geodesic_mean = med("geodesic_mean"),
    geodesic_normalized_mean = med("geodesic_normalized_mean"),
    eval_n = as.integer(round(med("eval_n"))),
    py_runs = runs
  )
}

compare_scenario <- function(name, cfg, r_res, py_res) {
  baseline_available <- identical(py_res$status, "ok") && is.finite(parity_numeric_or_na(py_res$runtime_sec))
  r_backend <- parity_or_else(r_res$kernel_backend, "unknown")
  runtime_comparable <- baseline_available && identical(r_backend, "cpp")

  py_runtime <- parity_numeric_or_na(py_res$runtime_sec)
  py_objective <- parity_numeric_or_na(py_res$objective)
  py_map_fro <- parity_numeric_or_na(py_res$map_fro_norm)
  py_map_orth <- parity_numeric_or_na(py_res$map_orth_resid)
  py_acc <- parity_numeric_or_na(py_res$accuracy)
  py_geod_n <- parity_numeric_or_na(py_res$geodesic_normalized_mean)

  runtime_improvement <- if (runtime_comparable) (py_runtime - r_res$runtime_sec) / py_runtime else NA_real_
  accuracy_delta <- if (baseline_available && is.finite(py_acc)) r_res$accuracy - py_acc else NA_real_
  geodesic_norm_delta <- if (baseline_available && is.finite(py_geod_n)) r_res$geodesic_normalized_mean - py_geod_n else NA_real_
  objective_rel_gap <- if (baseline_available && is.finite(py_objective)) {
    abs(r_res$objective - py_objective) / max(abs(py_objective), 1e-12)
  } else {
    NA_real_
  }
  map_fro_norm_delta <- if (baseline_available && is.finite(py_map_fro)) r_res$map_fro_norm - py_map_fro else NA_real_
  map_orth_resid_delta <- if (baseline_available && is.finite(py_map_orth)) r_res$map_orth_resid - py_map_orth else NA_real_

  note <- if (!baseline_available) {
    paste("pyFM unavailable:", parity_or_else(py_res$error, "unknown"))
  } else if (!runtime_comparable) {
    sprintf("runtime not comparable: kernel_backend=%s", r_backend)
  } else {
    "ok"
  }

  data.frame(
    scenario = name,
    n = as.integer(cfg$n),
    k = as.integer(cfg$k),
    p = as.integer(cfg$p),
    noise = as.numeric(cfg$noise),
    basis_noise = as.numeric(cfg$basis_noise),
    descriptor_corruption = as.numeric(cfg$descriptor_corruption),
    eval_fraction = as.numeric(cfg$eval_fraction),
    r_runtime_sec = as.numeric(r_res$runtime_sec),
    r_runtime_match_sec = as.numeric(r_res$runtime_match_sec),
    r_runtime_refine_sec = as.numeric(r_res$runtime_refine_sec),
    r_objective = as.numeric(r_res$objective),
    r_map_fro_norm = as.numeric(r_res$map_fro_norm),
    r_map_orth_resid = as.numeric(r_res$map_orth_resid),
    r_kernel_backend = as.character(r_backend),
    r_kernel_backend_mixed = isTRUE(r_res$kernel_backend_mixed),
    r_accuracy = as.numeric(r_res$accuracy),
    r_geodesic_norm = as.numeric(r_res$geodesic_normalized_mean),
    py_status = as.character(py_res$status),
    py_runtime_sec = py_runtime,
    py_objective = py_objective,
    py_map_fro_norm = py_map_fro,
    py_map_orth_resid = py_map_orth,
    py_accuracy = py_acc,
    py_geodesic_norm = py_geod_n,
    runtime_improvement_ratio = runtime_improvement,
    runtime_comparable = runtime_comparable,
    accuracy_delta = accuracy_delta,
    geodesic_norm_delta = geodesic_norm_delta,
    objective_rel_gap = objective_rel_gap,
    map_fro_norm_delta = map_fro_norm_delta,
    map_orth_resid_delta = map_orth_resid_delta,
    baseline_available = baseline_available,
    note = note,
    stringsAsFactors = FALSE
  )
}

parity_write_report <- function(report, output_dir = "benchmarks/parity") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  rds_path <- file.path(output_dir, paste0("parity-benchmark-", stamp, ".rds"))
  latest_rds_path <- file.path(output_dir, "parity-benchmark-latest.rds")
  md_path <- file.path(output_dir, "parity-benchmark-latest.md")

  saveRDS(report, rds_path)
  saveRDS(report, latest_rds_path)

  rows <- report$scenario_results
  lines <- c(
    "# Parity Benchmark (fmapsR vs pyFM)",
    "",
    sprintf("- Timestamp: %s", report$metadata$timestamp),
    sprintf("- Platform: %s", report$metadata$platform),
    sprintf("- Scenarios: %s", paste(report$metadata$scenarios, collapse = ",")),
    sprintf("- R runs/scenario: %d", report$metadata$r_runs),
    sprintf("- pyFM runs/scenario: %d", report$metadata$py_runs),
    sprintf("- R cg_maxit: %d", report$metadata$r_cg_maxit),
    sprintf("- R cg_tol: %g", report$metadata$r_cg_tol),
    sprintf("- R refine_nit: %d", report$metadata$r_refine_nit),
    sprintf("- R backend request: %s", report$metadata$r_kernel_backend),
    sprintf("- Package load mode: %s", report$metadata$package_load_mode),
    sprintf("- Baseline available scenarios: %d/%d", report$summary$n_baseline_available, report$summary$n_scenarios),
    "",
    "## Scenario Results",
    "| scenario | r_backend | runtime_comparable | r_runtime_sec | py_runtime_sec | runtime_improvement_ratio | r_accuracy | py_accuracy | accuracy_delta | r_geodesic_norm | py_geodesic_norm | geodesic_norm_delta | baseline_available | note |",
    "|---|---|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---|"
  )

  for (i in seq_len(nrow(rows))) {
    r <- rows[i, , drop = FALSE]
    backend <- if ("r_kernel_backend" %in% names(r)) as.character(r$r_kernel_backend) else "unknown"
    runtime_cmp <- if ("runtime_comparable" %in% names(r) && isTRUE(r$runtime_comparable)) "TRUE" else "FALSE"
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %.6f | %s | %s | %.6f | %s | %s | %.6f | %s | %s | %s | %s |",
      r$scenario,
      backend,
      runtime_cmp,
      r$r_runtime_sec,
      as.character(r$py_runtime_sec),
      as.character(r$runtime_improvement_ratio),
      r$r_accuracy,
      as.character(r$py_accuracy),
      as.character(r$accuracy_delta),
      r$r_geodesic_norm,
      as.character(r$py_geodesic_norm),
      as.character(r$geodesic_norm_delta),
      if (isTRUE(r$baseline_available)) "TRUE" else "FALSE",
      r$note
    ))
  }

  lines <- c(
    lines,
    "",
    "## Objective/Map Parity",
    "| scenario | r_objective | py_objective | objective_rel_gap | r_map_fro_norm | py_map_fro_norm | map_fro_norm_delta | r_map_orth_resid | py_map_orth_resid | map_orth_resid_delta |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
  )

  for (i in seq_len(nrow(rows))) {
    r <- rows[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
      r$scenario,
      as.character(r$r_objective),
      as.character(r$py_objective),
      as.character(r$objective_rel_gap),
      as.character(r$r_map_fro_norm),
      as.character(r$py_map_fro_norm),
      as.character(r$map_fro_norm_delta),
      as.character(r$r_map_orth_resid),
      as.character(r$py_map_orth_resid),
      as.character(r$map_orth_resid_delta)
    ))
  }

  lines <- c(
    lines,
    "",
    "## Aggregates",
    sprintf("- Median runtime improvement ratio (available scenarios): %s", as.character(report$summary$median_runtime_improvement_ratio)),
    sprintf("- Mean accuracy delta (R - pyFM): %s", as.character(report$summary$mean_accuracy_delta)),
    sprintf("- Mean geodesic normalized delta (R - pyFM): %s", as.character(report$summary$mean_geodesic_norm_delta)),
    sprintf("- Median objective relative gap: %s", as.character(report$summary$median_objective_rel_gap)),
    sprintf("- Mean map Frobenius norm delta (R - pyFM): %s", as.character(report$summary$mean_map_fro_norm_delta)),
    sprintf("- Mean map orthogonality residual delta (R - pyFM): %s", as.character(report$summary$mean_map_orth_resid_delta)),
    ""
  )

  writeLines(lines, md_path)
  list(rds = rds_path, latest_rds = latest_rds_path, markdown = md_path)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  root <- parity_find_repo_root()
  load_info <- parity_load_fmaps_package(root = root, target_env = globalenv())
  scenarios <- default_scenarios()

  if (!all(opts$scenarios %in% names(scenarios))) {
    stop(sprintf("Unknown scenario(s): %s", paste(setdiff(opts$scenarios, names(scenarios)), collapse = ",")), call. = FALSE)
  }

  rows <- list()
  run_details <- list()

  for (i in seq_along(opts$scenarios)) {
    nm <- opts$scenarios[[i]]
    cfg <- scenarios[[nm]]
    seed_i <- opts$seed + (i - 1L) * benchmark_seed_registry()$parity_scenario_stride

    r_res <- run_r_parity(
      cfg = cfg,
      seed = seed_i,
      n_runs = opts$r_runs,
      cg_maxit = opts$r_cg_maxit,
      cg_tol = opts$r_cg_tol,
      refine_nit = opts$r_refine_nit,
      kernel_backend = opts$r_kernel_backend,
      warmup = opts$warmup
    )
    py_res <- run_pyfm_parity(
      cfg = cfg,
      seed = seed_i,
      n_runs = opts$py_runs,
      maxit = opts$py_maxit,
      icp_nit = opts$py_icp_nit
    )

    rows[[length(rows) + 1L]] <- compare_scenario(nm, cfg, r_res, py_res)
    run_details[[nm]] <- list(r = r_res, py = py_res)
  }

  scenario_df <- do.call(rbind, rows)
  available <- scenario_df$baseline_available %in% TRUE

  report <- list(
    metadata = list(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      platform = R.version$platform,
      scenarios = opts$scenarios,
      r_runs = opts$r_runs,
      py_runs = opts$py_runs,
      r_cg_maxit = opts$r_cg_maxit,
      r_cg_tol = opts$r_cg_tol,
      r_refine_nit = opts$r_refine_nit,
      r_kernel_backend = opts$r_kernel_backend,
      py_maxit = opts$py_maxit,
      py_icp_nit = opts$py_icp_nit,
      warmup = opts$warmup,
      seed = opts$seed,
      package_load_mode = load_info$mode,
      package_load_error = load_info$load_error
    ),
    scenario_results = scenario_df,
    runs = run_details,
    summary = list(
      n_scenarios = nrow(scenario_df),
      n_baseline_available = sum(available),
      median_runtime_improvement_ratio = if (any(available)) stats::median(scenario_df$runtime_improvement_ratio[available], na.rm = TRUE) else NA_real_,
      mean_accuracy_delta = if (any(available)) mean(scenario_df$accuracy_delta[available], na.rm = TRUE) else NA_real_,
      mean_geodesic_norm_delta = if (any(available)) mean(scenario_df$geodesic_norm_delta[available], na.rm = TRUE) else NA_real_,
      median_objective_rel_gap = if (any(available)) stats::median(scenario_df$objective_rel_gap[available], na.rm = TRUE) else NA_real_,
      mean_map_fro_norm_delta = if (any(available)) mean(scenario_df$map_fro_norm_delta[available], na.rm = TRUE) else NA_real_,
      mean_map_orth_resid_delta = if (any(available)) mean(scenario_df$map_orth_resid_delta[available], na.rm = TRUE) else NA_real_
    )
  )

  outputs <- parity_write_report(report, output_dir = opts$output_dir)

  cat("Parity benchmark complete\n")
  cat("RDS:", outputs$rds, "\n")
  cat("Latest RDS:", outputs$latest_rds, "\n")
  cat("Summary:", outputs$markdown, "\n")
  cat(sprintf("Baseline available scenarios: %d/%d\n", report$summary$n_baseline_available, report$summary$n_scenarios))
}

if (sys.nframe() == 0) {
  main()
}
