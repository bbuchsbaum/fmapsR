#!/usr/bin/env Rscript

parse_kv <- function(lines) {
  out <- list()
  for (line in lines) {
    if (!grepl("=", line, fixed = TRUE)) next
    kv <- strsplit(line, "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- paste(kv[-1], collapse = "=")
    out[[key]] <- val
  }
  out
}

extract_rss_mb <- function(lines) {
  darwin <- grep("maximum resident set size", lines, value = TRUE)
  if (length(darwin) > 0L) {
    raw <- as.numeric(sub("^\\s*([0-9]+).*$", "\\1", darwin[[1]]))
    return(raw / (1024^2))
  }

  linux <- grep("Maximum resident set size", lines, value = TRUE)
  if (length(linux) > 0L) {
    raw <- as.numeric(sub("^.*: *([0-9]+).*$", "\\1", linux[[1]]))
    return(raw / 1024)
  }

  NA_real_
}

current_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) {
    stop("Unable to resolve script path from commandArgs", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
}

run_child <- function(batch_size = NULL) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required", call. = FALSE)
  }
  pkgload::load_all(".", quiet = TRUE)

  set.seed(123)

  n <- 240
  k <- 96
  p <- 900

  make_domain <- function() {
    op <- diag(seq_len(n))
    fm_basis(
      fm_domain_generic(n_samples = n, operator = op),
      k = k,
      solver = "base",
      cache = FALSE
    )
  }

  source <- make_domain()
  target <- make_domain()

  src_desc <- matrix(rnorm(n * p), nrow = n, ncol = p)
  tgt_desc <- src_desc + matrix(rnorm(n * p, sd = 0.02), nrow = n, ncol = p)

  t <- system.time({
    fit <- fm_match(
      source = source,
      target = target,
      descriptors = list(source = src_desc, target = tgt_desc),
      penalties = list(descr = 1e-1, lap = 1e-3, comm = 1),
      optimizer = "cg",
      kernel_backend = "r",
      cg_maxit = 2,
      descriptor_batch_size = batch_size
    )
  })[["elapsed"]]

  cat(sprintf("status=ok\n"))
  cat(sprintf("runtime_sec=%.6f\n", t))
  cat(sprintf("objective=%.12f\n", fit$diagnostics$total_objective))
  cat(sprintf("commutativity_mode=%s\n", fit$diagnostics$commutativity_mode))
  cat(sprintf("descriptor_batch_size=%s\n", as.character(fit$diagnostics$descriptor_batch_size)))
}

run_case <- function(batch_size = NULL, label = "unbatched") {
  script <- current_script_path()
  rscript <- file.path(R.home("bin"), "Rscript")
  time_bin <- Sys.which("/usr/bin/time")
  if (!nzchar(time_bin)) {
    stop("/usr/bin/time is required for RSS measurement", call. = FALSE)
  }

  time_flag <- if (identical(Sys.info()[["sysname"]], "Darwin")) "-l" else "-v"

  args <- c(time_flag, rscript, script, "--child")
  if (!is.null(batch_size)) {
    args <- c(args, "--batch-size", as.character(batch_size))
  }

  out <- system2(time_bin, args = args, stdout = TRUE, stderr = TRUE)
  kv <- parse_kv(out)

  list(
    label = label,
    batch_size = batch_size,
    runtime_sec = as.numeric(kv$runtime_sec),
    objective = as.numeric(kv$objective),
    commutativity_mode = kv$commutativity_mode,
    rss_mb = extract_rss_mb(out)
  )
}

write_report <- function(report, output_dir = "benchmarks/pairwise") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  rds_path <- file.path(output_dir, paste0("batched-memory-benchmark-", stamp, ".rds"))
  md_path <- file.path(output_dir, "batched-memory-benchmark-latest.md")

  saveRDS(report, rds_path)

  lines <- c(
    "# Batched Memory Benchmark",
    "",
    sprintf("- Timestamp: %s", report$metadata$timestamp),
    sprintf("- Platform: %s", report$metadata$platform),
    sprintf("- R version: %s", report$metadata$r_version),
    "",
    "## Before (Unbatched)",
    sprintf("- Runtime (s): %.6f", report$unbatched$runtime_sec),
    sprintf("- Max RSS (MB): %.2f", report$unbatched$rss_mb),
    sprintf("- Objective: %.6f", report$unbatched$objective),
    "",
    "## After (Descriptor Batching)",
    sprintf("- Runtime (s): %.6f", report$batched$runtime_sec),
    sprintf("- Max RSS (MB): %.2f", report$batched$rss_mb),
    sprintf("- Objective: %.6f", report$batched$objective),
    sprintf("- Batch size: %s", as.character(report$batched$batch_size)),
    "",
    "## Delta",
    sprintf("- RSS reduction ratio: %.6f", report$comparison$rss_reduction_ratio),
    sprintf("- Runtime delta ratio: %.6f", report$comparison$runtime_delta_ratio),
    sprintf("- Objective absolute delta: %.6f", report$comparison$objective_abs_delta),
    sprintf("- Target met (RSS >= 0.30): %s", as.character(report$comparison$target_met)),
    ""
  )

  writeLines(lines, md_path)

  list(rds = rds_path, markdown = md_path)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  is_child <- "--child" %in% args

  if (is_child) {
    batch_size <- NULL
    if ("--batch-size" %in% args) {
      pos <- which(args == "--batch-size")
      if (length(pos) == 1L && pos < length(args)) {
        batch_size <- as.integer(args[[pos + 1L]])
      }
    }
    run_child(batch_size = batch_size)
    return(invisible(NULL))
  }

  unbatched <- run_case(batch_size = NULL, label = "unbatched")
  batched <- run_case(batch_size = 8L, label = "batched")

  rss_reduction <- (unbatched$rss_mb - batched$rss_mb) / unbatched$rss_mb
  runtime_delta <- (batched$runtime_sec - unbatched$runtime_sec) / unbatched$runtime_sec

  report <- list(
    metadata = list(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      platform = R.version$platform,
      r_version = as.character(getRversion())
    ),
    unbatched = unbatched,
    batched = batched,
    comparison = list(
      rss_reduction_ratio = rss_reduction,
      runtime_delta_ratio = runtime_delta,
      objective_abs_delta = abs(unbatched$objective - batched$objective),
      target_met = is.finite(rss_reduction) && rss_reduction >= 0.30
    )
  )

  outputs <- write_report(report)

  cat("Batched memory benchmark complete\n")
  cat("RDS:", outputs$rds, "\n")
  cat("Summary:", outputs$markdown, "\n")
  cat("RSS reduction ratio:", sprintf("%.6f", report$comparison$rss_reduction_ratio), "\n")
  cat("Target met:", as.character(report$comparison$target_met), "\n")
}

if (sys.nframe() == 0) {
  main()
}
