release_tool_path <- function(filename) {
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

source(release_tool_path("benchmark_config.R"), local = environment())

release_pyfm_vendor_info <- function(root = ".") {
  pyfm_dir <- file.path(root, "pyFM")
  if (!dir.exists(pyfm_dir)) {
    return(list(
      available = FALSE,
      commit = NA_character_,
      branch = NA_character_,
      remote = NA_character_
    ))
  }

  git_one <- function(args) {
    out <- tryCatch(
      suppressWarnings(system2("git", c("-C", pyfm_dir, args), stdout = TRUE, stderr = TRUE)),
      error = function(e) character()
    )
    out <- trimws(out)
    out[nzchar(out)][1]
  }

  list(
    available = TRUE,
    commit = git_one(c("rev-parse", "--short", "HEAD")),
    branch = git_one(c("branch", "--show-current")),
    remote = git_one(c("remote", "get-url", "origin"))
  )
}

release_quality_check_row <- function(row, thresholds) {
  scenario <- as.character(row$scenario[[1]])
  parity_thresholds <- benchmark_parity_gate_thresholds()
  baseline_ok <- isTRUE(row$baseline_available[[1]])

  acc_min <- parity_thresholds$accuracy_delta_min[[scenario]]
  geod_max <- parity_thresholds$geodesic_norm_delta_max[[scenario]]
  map_fro_abs_max <- parity_thresholds$map_fro_norm_delta_abs_max[[scenario]]
  map_orth_abs_max <- parity_thresholds$map_orth_resid_delta_abs_max[[scenario]]
  objective_required <- isTRUE(thresholds$require_objective_parity)
  obj_abs_max <- parity_thresholds$objective_abs_gap_max[[scenario]]

  acc <- as.numeric(row$accuracy_delta[[1]])
  geod <- as.numeric(row$geodesic_norm_delta[[1]])
  map_fro_abs <- abs(as.numeric(row$map_fro_norm_delta[[1]]))
  map_orth_abs <- abs(as.numeric(row$map_orth_resid_delta[[1]]))
  obj_abs <- as.numeric(row$objective_abs_gap[[1]])

  acc_ok <- is.finite(acc) && (acc >= acc_min)
  geod_ok <- is.finite(geod) && (geod <= geod_max)
  map_fro_ok <- is.finite(map_fro_abs) && (map_fro_abs <= map_fro_abs_max)
  map_orth_ok <- is.finite(map_orth_abs) && (map_orth_abs <= map_orth_abs_max)
  obj_ok <- (!objective_required) || (is.finite(obj_abs) && (obj_abs <= obj_abs_max))

  reason <- if (!baseline_ok) {
    "baseline_unavailable"
  } else if (!acc_ok) {
    sprintf("accuracy_delta %.6f < %.6f", acc, acc_min)
  } else if (!geod_ok) {
    sprintf("geodesic_norm_delta %.6f > %.6f", geod, geod_max)
  } else if (!map_fro_ok) {
    sprintf("|map_fro_norm_delta| %.6f > %.6f", map_fro_abs, map_fro_abs_max)
  } else if (!map_orth_ok) {
    sprintf("|map_orth_resid_delta| %.6f > %.6f", map_orth_abs, map_orth_abs_max)
  } else if (!obj_ok) {
    sprintf("objective_abs_gap %.6f > %.6f", obj_abs, obj_abs_max)
  } else {
    "ok"
  }

  list(
    scenario = scenario,
    ok = baseline_ok && acc_ok && geod_ok && map_fro_ok && map_orth_ok && obj_ok,
    reason = reason
  )
}

release_quality_checks <- function(rows, thresholds = benchmark_release_claim_thresholds()) {
  required <- thresholds$required_quality_scenarios
  missing <- setdiff(required, as.character(rows$scenario))
  if (length(missing) > 0L) {
    return(list(
      ok = FALSE,
      checks = lapply(missing, function(scn) list(scenario = scn, ok = FALSE, reason = "missing_scenario")),
      missing = missing
    ))
  }

  checks <- lapply(required, function(scn) {
    row <- rows[rows$scenario == scn, , drop = FALSE]
    release_quality_check_row(row, thresholds = thresholds)
  })

  list(
    ok = all(vapply(checks, `[[`, logical(1), "ok")),
    checks = checks,
    missing = character()
  )
}

release_speed_check <- function(report, thresholds = benchmark_release_claim_thresholds()) {
  source <- thresholds$speed_source %||% "pairwise"

  if (identical(source, "parity")) {
    rows <- report$parity$scenario_results
    required <- thresholds$required_quality_scenarios
    rows <- rows[rows$scenario %in% required, , drop = FALSE]
    baseline_ok <- nrow(rows) == length(required) &&
      all(rows$baseline_available %in% TRUE) &&
      all(rows$runtime_comparable %in% TRUE)
    improvement <- if (baseline_ok) stats::median(as.numeric(rows$runtime_improvement_ratio), na.rm = TRUE) else NA_real_
    pass <- baseline_ok && is.finite(improvement) && improvement >= as.numeric(thresholds$speed_target_ratio)

    reason <- if (!baseline_ok) {
      "baseline_unavailable_or_runtime_not_comparable"
    } else if (!is.finite(improvement)) {
      "missing_runtime_improvement"
    } else if (!pass) {
      sprintf("parity_speed %.6f below target %.6f", improvement, as.numeric(thresholds$speed_target_ratio))
    } else {
      "ok"
    }

    return(list(
      ok = pass,
      source = "parity",
      reason = reason,
      baseline_available = baseline_ok,
      target_ratio = as.numeric(thresholds$speed_target_ratio),
      expected_target_ratio = as.numeric(thresholds$speed_target_ratio),
      stability_margin = as.numeric(thresholds$speed_stability_margin),
      runtime_improvement_ratio = improvement,
      target_met_raw = pass,
      target_met_stable = NA,
      stability_decision = "not_applicable"
    ))
  }

  cmp <- report$pairwise$comparison
  baseline_ok <- isTRUE(cmp$baseline_available)
  raw_ok <- isTRUE(cmp$target_met_raw)
  stable_ok <- isTRUE(cmp$target_met_stable)
  report_target <- as.numeric(cmp$target_ratio)
  expected_target <- as.numeric(thresholds$speed_target_ratio)
  margin <- as.numeric(thresholds$speed_stability_margin)

  target_match <- is.finite(report_target) && abs(report_target - expected_target) <= 1e-12
  stable_required <- isTRUE(thresholds$require_stable_speed)
  pass <- baseline_ok &&
    target_match &&
    if (stable_required) stable_ok else raw_ok

  reason <- if (!baseline_ok) {
    "baseline_unavailable"
  } else if (!target_match) {
    sprintf("target_ratio %.6f != %.6f", report_target, expected_target)
  } else if (stable_required && !stable_ok) {
    sprintf(
      "stable_speed %.6f below target %.6f (decision: %s)",
      as.numeric(cmp$runtime_improvement_ratio),
      expected_target,
      as.character(cmp$stability_decision)
    )
  } else if (!stable_required && !raw_ok) {
    sprintf(
      "raw_speed %.6f below target %.6f",
      as.numeric(cmp$runtime_improvement_ratio),
      expected_target
    )
  } else {
    "ok"
  }

  list(
    ok = pass,
    source = "pairwise",
    reason = reason,
    baseline_available = baseline_ok,
    target_ratio = report_target,
    expected_target_ratio = expected_target,
    stability_margin = margin,
    runtime_improvement_ratio = as.numeric(cmp$runtime_improvement_ratio),
    target_met_raw = raw_ok,
    target_met_stable = stable_ok,
    stability_decision = as.character(cmp$stability_decision)
  )
}

release_claim_checks <- function(report, thresholds = benchmark_release_claim_thresholds()) {
  quality <- release_quality_checks(report$parity$scenario_results, thresholds = thresholds)
  speed <- release_speed_check(report, thresholds = thresholds)

  list(
    thresholds = thresholds,
    quality = quality,
    speed = speed,
    ok = isTRUE(quality$ok) && isTRUE(speed$ok)
  )
}
