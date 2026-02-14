#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    scales = c(10L, 25L, 50L),
    mode = "robust",
    nit = 8L,
    seed = 2026L,
    output_dir = "benchmarks/network-sync",
    enforce_regression = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NULL

    if (key == "--scales" && !is.null(val)) {
      out$scales <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
      i <- i + 2L
    } else if (key == "--mode" && !is.null(val)) {
      out$mode <- val
      i <- i + 2L
    } else if (key == "--nit" && !is.null(val)) {
      out$nit <- as.integer(val)
      i <- i + 2L
    } else if (key == "--seed" && !is.null(val)) {
      out$seed <- as.integer(val)
      i <- i + 2L
    } else if (key == "--output-dir" && !is.null(val)) {
      out$output_dir <- val
      i <- i + 2L
    } else if (key == "--enforce-regression") {
      out$enforce_regression <- TRUE
      i <- i + 1L
    } else {
      stop(sprintf("Unknown or malformed argument: %s", key), call. = FALSE)
    }
  }

  out
}

mk_domain <- function(k = 4L, n_samples = 10L) {
  op <- diag(seq_len(n_samples))
  fm_basis(
    fm_domain_generic(n_samples = n_samples, operator = op),
    k = k,
    solver = "base",
    cache = FALSE
  )
}

random_orth <- function(k) {
  as.matrix(qr.Q(qr(matrix(rnorm(k * k), nrow = k, ncol = k))))
}

build_synthetic_network <- function(n_nodes, k = 4L, seed = 1L, noise_sd = 0.30, outlier_sd = 1.00, outlier_frac = 0.10) {
  set.seed(seed)

  node_names <- sprintf("n%03d", seq_len(n_nodes))
  domains <- stats::setNames(replicate(n_nodes, mk_domain(k = k), simplify = FALSE), node_names)

  latent <- stats::setNames(lapply(seq_len(n_nodes), function(i) random_orth(k)), node_names)

  net <- fm_network(domains, directed = TRUE)

  for (i in node_names) {
    for (j in node_names) {
      if (identical(i, j)) next
      sd <- if (stats::runif(1) < outlier_frac) outlier_sd else noise_sd
      Cij <- latent[[j]] %*% t(latent[[i]]) + matrix(rnorm(k * k, sd = sd), nrow = k, ncol = k)
      net <- fm_network_add_map(net, i, j, Cij, weight = 1)
    }
  }

  net
}

network_metrics <- function(net) {
  scores <- fm_network_cycle_scores(net)
  stats::median(scores$error, na.rm = TRUE)
}

run_one <- function(n_nodes, mode, nit, seed) {
  net <- build_synthetic_network(
    n_nodes = n_nodes,
    seed = seed,
    noise_sd = 0.30,
    outlier_sd = 1.00,
    outlier_frac = 0.10
  )

  pre <- network_metrics(net)

  t_sync <- system.time({
    fit <- fm_sync(net, mode = mode, nit = nit)
  })[["elapsed"]]

  post <- network_metrics(fit$network)
  reduction <- if (is.finite(pre) && pre > 0 && is.finite(post)) (pre - post) / pre else NA_real_

  list(
    n_nodes = n_nodes,
    n_edges = length(net$maps),
    mode = mode,
    nit = nit,
    runtime_sec = t_sync,
    pre_cycle_median = pre,
    post_cycle_median = post,
    cycle_reduction_ratio = reduction,
    throughput_edges_per_sec = length(net$maps) / t_sync
  )
}

reference_targets <- function() {
  # Conservative quality baselines for regression gating.
  c("10" = 0.60, "25" = 0.55, "50" = 0.50)
}

check_regression <- function(results) {
  targets <- reference_targets()

  rows <- lapply(results, function(x) {
    key <- as.character(x$n_nodes)
    target <- if (key %in% names(targets)) targets[[key]] else NA_real_
    threshold <- 0.90 * target
    ok <- is.finite(x$cycle_reduction_ratio) && is.finite(threshold) && x$cycle_reduction_ratio >= threshold

    list(
      n_nodes = x$n_nodes,
      target = target,
      threshold = threshold,
      observed = x$cycle_reduction_ratio,
      ok = ok
    )
  })

  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

write_report <- function(report, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  rds_path <- file.path(output_dir, paste0("network-sync-benchmark-", stamp, ".rds"))
  md_path <- file.path(output_dir, "network-sync-benchmark-latest.md")

  saveRDS(report, rds_path)

  lines <- c(
    "# Network Sync Benchmark",
    "",
    sprintf("- Timestamp: %s", report$metadata$timestamp),
    sprintf("- Platform: %s", report$metadata$platform),
    sprintf("- Mode: %s", report$metadata$mode),
    sprintf("- Iterations: %d", report$metadata$nit),
    "",
    "## Scale Results",
    "| nodes | edges | runtime_sec | throughput_edges_per_sec | pre_cycle_median | post_cycle_median | cycle_reduction_ratio |",
    "|---:|---:|---:|---:|---:|---:|---:|"
  )

  for (x in report$results) {
    lines <- c(lines, sprintf(
      "| %d | %d | %.6f | %.3f | %.6f | %.6f | %.6f |",
      x$n_nodes,
      x$n_edges,
      x$runtime_sec,
      x$throughput_edges_per_sec,
      x$pre_cycle_median,
      x$post_cycle_median,
      x$cycle_reduction_ratio
    ))
  }

  lines <- c(lines, "", "## Regression Gate")
  lines <- c(lines, "| nodes | target_ratio | threshold_90pct | observed_ratio | pass |")
  lines <- c(lines, "|---:|---:|---:|---:|:---:|")
  for (i in seq_len(nrow(report$regression))) {
    r <- report$regression[i, , drop = FALSE]
    lines <- c(lines, sprintf(
      "| %d | %.6f | %.6f | %.6f | %s |",
      as.integer(r$n_nodes),
      as.numeric(r$target),
      as.numeric(r$threshold),
      as.numeric(r$observed),
      if (isTRUE(r$ok[[1]])) "TRUE" else "FALSE"
    ))
  }

  lines <- c(lines, "", sprintf("- Overall regression pass: %s", as.character(report$overall_pass)), "")

  writeLines(lines, md_path)

  list(rds = rds_path, markdown = md_path)
}

main <- function() {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required to run benchmark script from source", call. = FALSE)
  }
  pkgload::load_all(".", quiet = TRUE)

  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  results <- lapply(seq_along(opts$scales), function(i) {
    n_nodes <- opts$scales[[i]]
    run_one(n_nodes = n_nodes, mode = opts$mode, nit = opts$nit, seed = opts$seed + i - 1L)
  })

  reg <- check_regression(results)
  overall_pass <- all(reg$ok)

  report <- list(
    metadata = list(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      platform = R.version$platform,
      mode = opts$mode,
      nit = opts$nit,
      scales = opts$scales
    ),
    results = results,
    regression = reg,
    overall_pass = overall_pass
  )

  outputs <- write_report(report, output_dir = opts$output_dir)

  cat("Network benchmark complete\n")
  cat("RDS:", outputs$rds, "\n")
  cat("Summary:", outputs$markdown, "\n")
  cat("Regression pass:", as.character(overall_pass), "\n")

  if (isTRUE(opts$enforce_regression) && !isTRUE(overall_pass)) {
    stop("Network sync regression gate failed", call. = FALSE)
  }
}

if (sys.nframe() == 0) {
  main()
}
