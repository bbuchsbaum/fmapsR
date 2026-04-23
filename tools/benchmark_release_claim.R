#!/usr/bin/env Rscript

tool_path <- function(filename) {
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

source(tool_path("benchmark_config.R"), local = environment())
source(tool_path("release_claim_helpers.R"), local = environment())

parse_args <- function(args) {
  thresholds <- benchmark_release_claim_thresholds()

  out <- list(
    output_dir = "benchmarks/release-claim",
    parity_scenarios = names(benchmark_parity_scenarios()),
    required_quality = thresholds$required_quality_scenarios,
    pairwise_r_runs = 3L,
    pairwise_py_runs = 3L,
    parity_r_runs = 3L,
    parity_py_runs = 3L,
    parity_r_cg_maxit = thresholds$parity_r_cg_maxit,
    speed_target_ratio = thresholds$speed_target_ratio,
    speed_stability_margin = thresholds$speed_stability_margin,
    speed_stability_boot = 500L,
    speed_stability_seed = benchmark_seed_registry()$network_seed
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NULL

    if (key == "--output-dir" && !is.null(val)) {
      out$output_dir <- val
      i <- i + 2L
    } else if (key == "--parity-scenarios" && !is.null(val)) {
      out$parity_scenarios <- strsplit(val, ",", fixed = TRUE)[[1]]
      i <- i + 2L
    } else if (key == "--required-quality" && !is.null(val)) {
      out$required_quality <- strsplit(val, ",", fixed = TRUE)[[1]]
      i <- i + 2L
    } else if (key == "--pairwise-r-runs" && !is.null(val)) {
      out$pairwise_r_runs <- as.integer(val)
      i <- i + 2L
    } else if (key == "--pairwise-py-runs" && !is.null(val)) {
      out$pairwise_py_runs <- as.integer(val)
      i <- i + 2L
    } else if (key == "--parity-r-runs" && !is.null(val)) {
      out$parity_r_runs <- as.integer(val)
      i <- i + 2L
    } else if (key == "--parity-py-runs" && !is.null(val)) {
      out$parity_py_runs <- as.integer(val)
      i <- i + 2L
    } else if (key == "--parity-r-cg-maxit" && !is.null(val)) {
      out$parity_r_cg_maxit <- as.integer(val)
      i <- i + 2L
    } else if (key == "--speed-target-ratio" && !is.null(val)) {
      out$speed_target_ratio <- as.numeric(val)
      i <- i + 2L
    } else if (key == "--speed-stability-margin" && !is.null(val)) {
      out$speed_stability_margin <- as.numeric(val)
      i <- i + 2L
    } else if (key == "--speed-stability-boot" && !is.null(val)) {
      out$speed_stability_boot <- as.integer(val)
      i <- i + 2L
    } else if (key == "--speed-stability-seed" && !is.null(val)) {
      out$speed_stability_seed <- as.integer(val)
      i <- i + 2L
    } else {
      stop(sprintf("Unknown or malformed argument: %s", key), call. = FALSE)
    }
  }

  if (!all(out$required_quality %in% out$parity_scenarios)) {
    stop("`--required-quality` must be a subset of `--parity-scenarios`", call. = FALSE)
  }

  out
}

release_thresholds_from_opts <- function(opts) {
  thr <- benchmark_release_claim_thresholds()
  thr$required_quality_scenarios <- opts$required_quality
  thr$parity_r_cg_maxit <- opts$parity_r_cg_maxit
  thr$speed_target_ratio <- opts$speed_target_ratio
  thr$speed_stability_margin <- opts$speed_stability_margin
  thr
}

run_benchmark_script <- function(script, args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  out <- tryCatch(
    suppressWarnings(system2(rscript, c(script, args), stdout = TRUE, stderr = TRUE)),
    error = function(e) structure(conditionMessage(e), status = 1L)
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      sprintf("Benchmark script failed: %s\n%s", script, paste(out, collapse = "\n")),
      call. = FALSE
    )
  }
  out
}

build_release_report <- function(pairwise_report, parity_report, thresholds, root = ".", artifacts = list()) {
  vendor <- release_pyfm_vendor_info(root)
  report <- list(
    metadata = list(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      platform = R.version$platform,
      git_commit = tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)[1], error = function(e) NA_character_),
      pyfm_vendor_commit = vendor$commit,
      pyfm_vendor_branch = vendor$branch,
      pyfm_vendor_remote = vendor$remote
    ),
    thresholds = thresholds,
    artifacts = artifacts,
    pairwise = pairwise_report,
    parity = parity_report
  )

  report$claim <- release_claim_checks(report, thresholds = thresholds)
  report
}

write_release_report <- function(report, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  rds_path <- file.path(output_dir, paste0("release-claim-", stamp, ".rds"))
  latest_rds_path <- file.path(output_dir, "release-claim-latest.rds")
  md_path <- file.path(output_dir, "release-claim-latest.md")

  saveRDS(report, rds_path)
  saveRDS(report, latest_rds_path)

  q_checks <- report$claim$quality$checks
  speed <- report$claim$speed
  objective_mode <- if (isTRUE(report$thresholds$require_objective_parity)) "required" else "diagnostic_only"
  pairwise_md <- if (is.null(report$artifacts$pairwise_markdown)) NA_character_ else report$artifacts$pairwise_markdown
  parity_md <- if (is.null(report$artifacts$parity_markdown)) NA_character_ else report$artifacts$parity_markdown

  lines <- c(
    "# Release Claim Benchmark",
    "",
    sprintf("- Timestamp: %s", report$metadata$timestamp),
    sprintf("- Platform: %s", report$metadata$platform),
    sprintf("- Repository commit: %s", as.character(report$metadata$git_commit)),
    sprintf("- Vendored pyFM commit: %s", as.character(report$metadata$pyfm_vendor_commit)),
    sprintf("- Vendored pyFM branch: %s", as.character(report$metadata$pyfm_vendor_branch)),
    sprintf("- Vendored pyFM remote: %s", as.character(report$metadata$pyfm_vendor_remote)),
    sprintf("- Overall claim: %s", if (isTRUE(report$claim$ok)) "PASS" else "FAIL"),
    sprintf("- Objective parity mode: %s", objective_mode),
    "",
    "## Speed Claim",
    sprintf("- Speed source: %s", as.character(speed$source)),
    sprintf("- Pairwise benchmark markdown: %s", as.character(pairwise_md)),
    sprintf("- Baseline available: %s", as.character(speed$baseline_available)),
    sprintf("- Runtime improvement ratio: %s", as.character(speed$runtime_improvement_ratio)),
    sprintf("- Target ratio: %s", as.character(speed$expected_target_ratio)),
    sprintf("- Stable target met: %s", as.character(speed$target_met_stable)),
    sprintf("- Stability decision: %s", as.character(speed$stability_decision)),
    sprintf("- Status: %s", if (isTRUE(speed$ok)) "PASS" else "FAIL"),
    sprintf("- Reason: %s", speed$reason),
    "",
    "## Quality Claim",
    sprintf("- Parity benchmark markdown: %s", as.character(parity_md)),
    "| scenario | status | accuracy_delta | geodesic_norm_delta | map_fro_norm_delta | map_orth_resid_delta | objective_abs_gap | reason |",
    "|---|:---:|---:|---:|---:|---:|---:|---|"
  )

  rows <- report$parity$scenario_results
  for (chk in q_checks) {
    row <- rows[rows$scenario == chk$scenario, , drop = FALSE]
    acc <- if (nrow(row) > 0L) as.character(row$accuracy_delta[[1]]) else NA_character_
    geod <- if (nrow(row) > 0L) as.character(row$geodesic_norm_delta[[1]]) else NA_character_
    map_fro <- if (nrow(row) > 0L) as.character(row$map_fro_norm_delta[[1]]) else NA_character_
    map_orth <- if (nrow(row) > 0L) as.character(row$map_orth_resid_delta[[1]]) else NA_character_
    obj_gap <- if (nrow(row) > 0L) as.character(row$objective_abs_gap[[1]]) else NA_character_
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %s | %s | %s | %s | %s |",
      chk$scenario,
      if (isTRUE(chk$ok)) "PASS" else "FAIL",
      acc,
      geod,
      map_fro,
      map_orth,
      obj_gap,
      chk$reason
    ))
  }

  writeLines(c(lines, ""), md_path)
  list(rds = rds_path, latest_rds = latest_rds_path, markdown = md_path)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  pairwise_dir <- file.path(opts$output_dir, "pairwise")
  parity_dir <- file.path(opts$output_dir, "parity")

  run_benchmark_script(
    tool_path("benchmark_pairwise.R"),
    c(
      "--r-runs", as.character(opts$pairwise_r_runs),
      "--py-runs", as.character(opts$pairwise_py_runs),
      "--target-ratio", as.character(opts$speed_target_ratio),
      "--stability-margin", as.character(opts$speed_stability_margin),
      "--stability-boot", as.character(opts$speed_stability_boot),
      "--stability-seed", as.character(opts$speed_stability_seed),
      "--output-dir", pairwise_dir
    )
  )

  run_benchmark_script(
    tool_path("benchmark_parity.R"),
    c(
      "--scenarios", paste(opts$parity_scenarios, collapse = ","),
      "--r-runs", as.character(opts$parity_r_runs),
      "--py-runs", as.character(opts$parity_py_runs),
      "--r-cg-maxit", as.character(opts$parity_r_cg_maxit),
      "--output-dir", parity_dir
    )
  )

  pairwise_report <- readRDS(file.path(pairwise_dir, "pairwise-benchmark-latest.rds"))
  parity_report <- readRDS(file.path(parity_dir, "parity-benchmark-latest.rds"))

  report <- build_release_report(
    pairwise_report = pairwise_report,
    parity_report = parity_report,
    thresholds = release_thresholds_from_opts(opts),
    root = ".",
    artifacts = list(
      pairwise_rds = normalizePath(file.path(pairwise_dir, "pairwise-benchmark-latest.rds"), winslash = "/", mustWork = FALSE),
      pairwise_markdown = normalizePath(file.path(pairwise_dir, "pairwise-benchmark-latest.md"), winslash = "/", mustWork = FALSE),
      parity_rds = normalizePath(file.path(parity_dir, "parity-benchmark-latest.rds"), winslash = "/", mustWork = FALSE),
      parity_markdown = normalizePath(file.path(parity_dir, "parity-benchmark-latest.md"), winslash = "/", mustWork = FALSE)
    )
  )

  outputs <- write_release_report(report, output_dir = opts$output_dir)

  cat("Release claim benchmark complete\n")
  cat("RDS:", outputs$rds, "\n")
  cat("Latest RDS:", outputs$latest_rds, "\n")
  cat("Summary:", outputs$markdown, "\n")
  cat("Overall claim:", if (isTRUE(report$claim$ok)) "PASS" else "FAIL", "\n")
}

if (sys.nframe() == 0) {
  main()
}
