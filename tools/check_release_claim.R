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
  thr <- benchmark_release_claim_thresholds()
  out <- list(
    report = "benchmarks/release-claim/release-claim-latest.rds",
    required_quality = thr$required_quality_scenarios,
    speed_target_ratio = thr$speed_target_ratio,
    speed_stability_margin = thr$speed_stability_margin,
    require_stable_speed = thr$require_stable_speed,
    require_objective_parity = thr$require_objective_parity
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NULL

    if (key == "--report" && !is.null(val)) {
      out$report <- val
      i <- i + 2L
    } else if (key == "--required-quality" && !is.null(val)) {
      out$required_quality <- strsplit(val, ",", fixed = TRUE)[[1]]
      i <- i + 2L
    } else if (key == "--speed-target-ratio" && !is.null(val)) {
      out$speed_target_ratio <- as.numeric(val)
      i <- i + 2L
    } else if (key == "--speed-stability-margin" && !is.null(val)) {
      out$speed_stability_margin <- as.numeric(val)
      i <- i + 2L
    } else if (key == "--require-objective") {
      out$require_objective_parity <- TRUE
      i <- i + 1L
    } else if (key == "--no-stable-speed") {
      out$require_stable_speed <- FALSE
      i <- i + 1L
    } else {
      stop(sprintf("Unknown or malformed argument: %s", key), call. = FALSE)
    }
  }

  out
}

thresholds_from_args <- function(opts) {
  thr <- benchmark_release_claim_thresholds()
  thr$required_quality_scenarios <- opts$required_quality
  thr$speed_target_ratio <- opts$speed_target_ratio
  thr$speed_stability_margin <- opts$speed_stability_margin
  thr$require_stable_speed <- opts$require_stable_speed
  thr$require_objective_parity <- opts$require_objective_parity
  thr
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!file.exists(opts$report)) {
    stop(sprintf("Release claim report not found: %s", opts$report), call. = FALSE)
  }

  report <- readRDS(opts$report)
  checks <- release_claim_checks(report, thresholds = thresholds_from_args(opts))

  for (chk in checks$quality$checks) {
    cat(sprintf("[release-claim] quality %s: %s\n", chk$scenario, chk$reason))
  }

  cat(sprintf("[release-claim] speed: %s\n", checks$speed$reason))
  cat(sprintf("[release-claim] overall: %s\n", if (isTRUE(checks$ok)) "PASS" else "FAIL"))

  if (!isTRUE(checks$ok)) {
    stop("Release claim gate failed", call. = FALSE)
  }
}

if (sys.nframe() == 0) {
  main()
}
