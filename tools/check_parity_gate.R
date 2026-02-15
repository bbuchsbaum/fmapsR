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

parse_args <- function(args) {
  out <- list(
    report = "benchmarks/parity/parity-benchmark-latest.rds",
    required = benchmark_parity_gate_thresholds()$required_scenarios
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NULL

    if (key == "--report" && !is.null(val)) {
      out$report <- val
      i <- i + 2L
    } else if (key == "--required" && !is.null(val)) {
      out$required <- strsplit(val, ",", fixed = TRUE)[[1]]
      i <- i + 2L
    } else {
      stop(sprintf("Unknown or malformed argument: %s", key), call. = FALSE)
    }
  }

  out
}

scenario_gate_row <- function(row, thresholds) {
  scenario <- as.character(row$scenario[[1]])
  baseline_ok <- isTRUE(row$baseline_available[[1]])

  acc_min <- thresholds$accuracy_delta_min[[scenario]]
  geod_max <- thresholds$geodesic_norm_delta_max[[scenario]]

  if (is.null(acc_min) || is.null(geod_max)) {
    return(list(
      scenario = scenario,
      ok = FALSE,
      reason = "missing_thresholds"
    ))
  }

  acc <- as.numeric(row$accuracy_delta[[1]])
  geod <- as.numeric(row$geodesic_norm_delta[[1]])
  acc_ok <- is.finite(acc) && (acc >= acc_min)
  geod_ok <- is.finite(geod) && (geod <= geod_max)

  reason <- if (!baseline_ok) {
    "baseline_unavailable"
  } else if (!acc_ok) {
    sprintf("accuracy_delta %.6f < %.6f", acc, acc_min)
  } else if (!geod_ok) {
    sprintf("geodesic_norm_delta %.6f > %.6f", geod, geod_max)
  } else {
    "ok"
  }

  list(
    scenario = scenario,
    ok = baseline_ok && acc_ok && geod_ok,
    reason = reason
  )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  if (!file.exists(opts$report)) {
    stop(sprintf("Parity report not found: %s", opts$report), call. = FALSE)
  }

  report <- readRDS(opts$report)
  rows <- report$scenario_results
  if (is.null(rows) || nrow(rows) == 0L) {
    stop("Parity report has no scenario rows", call. = FALSE)
  }

  missing <- setdiff(opts$required, as.character(rows$scenario))
  if (length(missing) > 0L) {
    stop(sprintf("Required scenarios missing from report: %s", paste(missing, collapse = ",")), call. = FALSE)
  }

  thresholds <- benchmark_parity_gate_thresholds()
  checks <- lapply(opts$required, function(scn) {
    row <- rows[rows$scenario == scn, , drop = FALSE]
    scenario_gate_row(row, thresholds = thresholds)
  })

  for (chk in checks) {
    cat(sprintf("[parity-gate] %s: %s\n", chk$scenario, chk$reason))
  }

  ok <- all(vapply(checks, `[[`, logical(1), "ok"))
  cat(sprintf("[parity-gate] overall: %s\n", if (ok) "PASS" else "FAIL"))

  if (!ok) {
    stop("Parity gate failed", call. = FALSE)
  }
}

if (sys.nframe() == 0) {
  main()
}
