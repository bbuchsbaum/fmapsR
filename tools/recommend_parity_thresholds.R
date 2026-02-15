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
    dir = "benchmarks/parity",
    quantile = 0.95,
    safety = 1.25,
    min_n = 3L
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NULL

    if (key == "--dir" && !is.null(val)) {
      out$dir <- val
      i <- i + 2L
    } else if (key == "--quantile" && !is.null(val)) {
      out$quantile <- as.numeric(val)
      i <- i + 2L
    } else if (key == "--safety" && !is.null(val)) {
      out$safety <- as.numeric(val)
      i <- i + 2L
    } else if (key == "--min-n" && !is.null(val)) {
      out$min_n <- as.integer(val)
      i <- i + 2L
    } else {
      stop(sprintf("Unknown or malformed argument: %s", key), call. = FALSE)
    }
  }

  if (!is.finite(out$quantile) || out$quantile <= 0 || out$quantile > 1) {
    stop("`--quantile` must be in (0,1]", call. = FALSE)
  }
  if (!is.finite(out$safety) || out$safety < 1) {
    stop("`--safety` must be >= 1", call. = FALSE)
  }
  if (!is.finite(out$min_n) || out$min_n < 1) {
    stop("`--min-n` must be >= 1", call. = FALSE)
  }

  out
}

collect_parity_rows <- function(dir) {
  files <- list.files(dir, pattern = "^parity-benchmark-.*[.]rds$", full.names = TRUE)
  if (length(files) == 0L) {
    return(data.frame())
  }

  required_cols <- c(
    "scenario",
    "baseline_available",
    "objective_abs_gap",
    "map_fro_norm_delta",
    "map_orth_resid_delta",
    "r_objective",
    "py_objective"
  )

  rows <- lapply(files, function(path) {
    rpt <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(rpt) || is.null(rpt$scenario_results) || !is.data.frame(rpt$scenario_results)) {
      return(NULL)
    }

    d <- rpt$scenario_results
    for (nm in required_cols) {
      if (!nm %in% names(d)) {
        d[[nm]] <- NA
      }
    }
    if (all(!is.finite(as.numeric(d$objective_abs_gap)))) {
      d$objective_abs_gap <- abs(as.numeric(d$r_objective) - as.numeric(d$py_objective))
    }

    d <- d[, required_cols, drop = FALSE]
    d$file <- basename(path)
    d
  })

  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out) || nrow(out) == 0L) {
    return(data.frame())
  }
  out[out$baseline_available %in% TRUE, , drop = FALSE]
}

recommend_threshold <- function(values, current, q, safety, min_n, floor_value) {
  vals <- as.numeric(values)
  vals <- vals[is.finite(vals)]
  if (length(vals) < min_n) {
    return(list(value = current, n = length(vals), reason = "insufficient_data"))
  }

  qv <- as.numeric(stats::quantile(vals, probs = q, names = FALSE, na.rm = TRUE))
  mv <- max(vals, na.rm = TRUE)
  rec <- max(floor_value, safety * qv, 1.05 * mv)
  list(value = rec, n = length(vals), reason = "empirical")
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  rows <- collect_parity_rows(opts$dir)

  if (nrow(rows) == 0L) {
    stop(sprintf("No parity reports found in `%s`", opts$dir), call. = FALSE)
  }

  thr <- benchmark_parity_gate_thresholds()
  scenarios <- names(thr$accuracy_delta_min)

  out <- lapply(scenarios, function(scn) {
    d <- rows[rows$scenario == scn, , drop = FALSE]
    rec_obj <- recommend_threshold(
      values = d$objective_abs_gap,
      current = as.numeric(thr$objective_abs_gap_max[[scn]]),
      q = opts$quantile,
      safety = opts$safety,
      min_n = opts$min_n,
      floor_value = 0.1
    )
    rec_fro <- recommend_threshold(
      values = abs(as.numeric(d$map_fro_norm_delta)),
      current = as.numeric(thr$map_fro_norm_delta_abs_max[[scn]]),
      q = opts$quantile,
      safety = opts$safety,
      min_n = opts$min_n,
      floor_value = 1e-6
    )
    rec_orth <- recommend_threshold(
      values = abs(as.numeric(d$map_orth_resid_delta)),
      current = as.numeric(thr$map_orth_resid_delta_abs_max[[scn]]),
      q = opts$quantile,
      safety = opts$safety,
      min_n = opts$min_n,
      floor_value = 1e-6
    )

    data.frame(
      scenario = scn,
      n = rec_obj$n,
      objective_abs_gap_current = as.numeric(thr$objective_abs_gap_max[[scn]]),
      objective_abs_gap_recommended = rec_obj$value,
      map_fro_abs_current = as.numeric(thr$map_fro_norm_delta_abs_max[[scn]]),
      map_fro_abs_recommended = rec_fro$value,
      map_orth_abs_current = as.numeric(thr$map_orth_resid_delta_abs_max[[scn]]),
      map_orth_abs_recommended = rec_orth$value,
      recommendation_basis = rec_obj$reason,
      stringsAsFactors = FALSE
    )
  })
  rec_df <- do.call(rbind, out)

  cat("# Parity Threshold Recommendations\n\n")
  cat(sprintf("- Input dir: `%s`\n", opts$dir))
  cat(sprintf("- Quantile: %.3f\n", opts$quantile))
  cat(sprintf("- Safety multiplier: %.3f\n", opts$safety))
  cat(sprintf("- Min samples/scenario: %d\n\n", opts$min_n))

  print(rec_df, row.names = FALSE)

  cat("\nSuggested benchmark_parity_gate_thresholds() replacements:\n")
  cat("objective_abs_gap_max = c(\n")
  for (i in seq_len(nrow(rec_df))) {
    row <- rec_df[i, , drop = FALSE]
    cat(sprintf("  %s = %.6f%s\n", row$scenario, row$objective_abs_gap_recommended, if (i < nrow(rec_df)) "," else ""))
  }
  cat("),\n")
  cat("map_fro_norm_delta_abs_max = c(\n")
  for (i in seq_len(nrow(rec_df))) {
    row <- rec_df[i, , drop = FALSE]
    cat(sprintf("  %s = %.6f%s\n", row$scenario, row$map_fro_abs_recommended, if (i < nrow(rec_df)) "," else ""))
  }
  cat("),\n")
  cat("map_orth_resid_delta_abs_max = c(\n")
  for (i in seq_len(nrow(rec_df))) {
    row <- rec_df[i, , drop = FALSE]
    cat(sprintf("  %s = %.6f%s\n", row$scenario, row$map_orth_abs_recommended, if (i < nrow(rec_df)) "," else ""))
  }
  cat(")\n")
}

if (sys.nframe() == 0) {
  main()
}
