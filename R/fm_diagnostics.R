pair_transfer_consistency <- function(network) {
  edge_df <- network_edge_table(network)
  if (nrow(edge_df) == 0L) {
    return(data.frame(i = character(), j = character(), error = numeric(), stringsAsFactors = FALSE))
  }

  rows <- list()
  idx <- 1L

  seen <- character()
  for (r in seq_len(nrow(edge_df))) {
    i <- edge_df$i[[r]]
    j <- edge_df$j[[r]]

    pair_key <- paste(sort(c(i, j)), collapse = "::")
    if (pair_key %in% seen) next
    seen <- c(seen, pair_key)

    Cij <- fm_network_get_map(network, i, j)
    Cji <- fm_network_get_map(network, j, i)
    if (is.null(Cij) || is.null(Cji)) next

    Mij <- as_map_matrix(Cij)
    Mji <- as_map_matrix(Cji)

    if (nrow(Mji) != ncol(Mij)) {
      err <- NA_real_
    } else {
      comp <- Mji %*% Mij
      I <- diag(1, nrow = nrow(comp), ncol = ncol(comp))
      err <- norm(comp - I, type = "F")
    }

    rows[[idx]] <- data.frame(i = i, j = j, error = err, stringsAsFactors = FALSE)
    idx <- idx + 1L
  }

  if (length(rows) == 0L) {
    return(data.frame(i = character(), j = character(), error = numeric(), stringsAsFactors = FALSE))
  }

  do.call(rbind, rows)
}

edge_map_as_fit <- function(network, i, j, map) {
  if (is_fm_fit(map)) {
    return(map)
  }

  C <- as_map_matrix(map)
  src <- network$domains[[i]]
  tgt <- network$domains[[j]]

  structure(
    list(
      C = C,
      k1 = ncol(C),
      k2 = nrow(C),
      source = src,
      target = tgt,
      penalties = list(),
      diagnostics = list(
        convergence = NA_integer_,
        message = "network_matrix",
        counts = c(fn = NA_integer_, gradient = NA_integer_),
        total_objective = NA_real_,
        objective_terms = list(),
        warning = NULL
      )
    ),
    class = "fm_fit"
  )
}

edge_map_quality <- function(network) {
  edge_df <- network_edge_table(network)
  if (nrow(edge_df) == 0L) {
    empty <- data.frame(
      key = character(),
      i = character(),
      j = character(),
      coverage_ratio = numeric(),
      coverage_entropy = numeric(),
      continuity_norm = numeric(),
      stringsAsFactors = FALSE
    )
    return(list(scores = empty, summary = list(n = 0L, n_coverage = 0L, n_continuity = 0L,
      coverage_median = NA_real_, coverage_mean = NA_real_,
      continuity_median = NA_real_, continuity_mean = NA_real_),
      note = "No edges available"))
  }

  rows <- lapply(seq_len(nrow(edge_df)), function(r) {
    key <- edge_df$key[[r]]
    i <- edge_df$i[[r]]
    j <- edge_df$j[[r]]
    map <- network$maps[[key]]

    fit <- edge_map_as_fit(network, i, j, map)

    m <- tryCatch(
      fm_fit_metrics(fit),
      error = function(e) NULL
    )

    cov_ratio <- if (!is.null(m) && !is.null(m$coverage)) as.numeric(m$coverage$ratio) else NA_real_
    cov_entropy <- if (!is.null(m) && !is.null(m$coverage)) as.numeric(m$coverage$entropy) else NA_real_
    cont_norm <- if (!is.null(m) && !is.null(m$continuity)) as.numeric(m$continuity$normalized_mean) else NA_real_

    data.frame(
      key = key,
      i = i,
      j = j,
      coverage_ratio = cov_ratio,
      coverage_entropy = cov_entropy,
      continuity_norm = cont_norm,
      stringsAsFactors = FALSE
    )
  })

  scores <- do.call(rbind, rows)

  cov_vals <- scores$coverage_ratio[is.finite(scores$coverage_ratio)]
  cont_vals <- scores$continuity_norm[is.finite(scores$continuity_norm)]

  note <- if (length(cov_vals) == 0L) {
    "Coverage unavailable for current map/domain combination"
  } else if (length(cont_vals) == 0L) {
    "Coverage available; continuity requires target adjacency"
  } else {
    "Coverage and continuity available"
  }

  list(
    scores = scores,
    summary = list(
      n = nrow(scores),
      n_coverage = length(cov_vals),
      n_continuity = length(cont_vals),
      coverage_median = if (length(cov_vals) > 0L) stats::median(cov_vals, na.rm = TRUE) else NA_real_,
      coverage_mean = if (length(cov_vals) > 0L) mean(cov_vals, na.rm = TRUE) else NA_real_,
      continuity_median = if (length(cont_vals) > 0L) stats::median(cont_vals, na.rm = TRUE) else NA_real_,
      continuity_mean = if (length(cont_vals) > 0L) mean(cont_vals, na.rm = TRUE) else NA_real_
    ),
    note = note
  )
}

#' Compute Network Diagnostic Metrics
#'
#' @param network `fm_network` object.
#'
#' @return List with cycle and transfer consistency diagnostics.
#' @export
fm_network_metrics <- function(network) {
  if (!is_fm_network(network)) {
    stop("`network` must inherit from `fm_network`", call. = FALSE)
  }

  cycles <- fm_network_cycle_scores(network)
  transfer <- pair_transfer_consistency(network)
  quality <- edge_map_quality(network)

  list(
    cycle_scores = cycles,
    cycle_summary = list(
      n = nrow(cycles),
      median = if (nrow(cycles) > 0) stats::median(cycles$error, na.rm = TRUE) else NA_real_,
      mean = if (nrow(cycles) > 0) mean(cycles$error, na.rm = TRUE) else NA_real_
    ),
    transfer_scores = transfer,
    transfer_summary = list(
      n = nrow(transfer),
      median = if (nrow(transfer) > 0) stats::median(transfer$error, na.rm = TRUE) else NA_real_,
      mean = if (nrow(transfer) > 0) mean(transfer$error, na.rm = TRUE) else NA_real_
    ),
    map_quality_scores = quality$scores,
    map_quality_summary = quality$summary,
    coverage = quality$summary$coverage_median,
    coverage_note = quality$note
  )
}

#' Compare Pre/Post Network Diagnostics
#'
#' @param pre Pre-sync `fm_network` object.
#' @param post Post-sync `fm_network` object.
#'
#' @return `fm_network_report` object.
#' @export
fm_network_compare <- function(pre, post) {
  if (!is_fm_network(pre) || !is_fm_network(post)) {
    stop("`pre` and `post` must inherit from `fm_network`", call. = FALSE)
  }

  pre_m <- fm_network_metrics(pre)
  post_m <- fm_network_metrics(post)

  cycle_reduction <- if (
    is.finite(pre_m$cycle_summary$median) && pre_m$cycle_summary$median > 0 &&
      is.finite(post_m$cycle_summary$median)
  ) {
    (pre_m$cycle_summary$median - post_m$cycle_summary$median) / pre_m$cycle_summary$median
  } else {
    NA_real_
  }

  transfer_reduction <- if (
    is.finite(pre_m$transfer_summary$median) && pre_m$transfer_summary$median > 0 &&
      is.finite(post_m$transfer_summary$median)
  ) {
    (pre_m$transfer_summary$median - post_m$transfer_summary$median) / pre_m$transfer_summary$median
  } else {
    NA_real_
  }

  structure(
    list(
      pre = pre_m,
      post = post_m,
      summary = list(
        cycle_reduction_ratio = cycle_reduction,
        transfer_reduction_ratio = transfer_reduction
      )
    ),
    class = "fm_network_report"
  )
}

#' Build Network Report from Sync Output
#'
#' @param x `fm_network_fit`, or `fm_network` if `post` is provided.
#' @param post Optional post-sync `fm_network` when `x` is a pre-sync network.
#'
#' @return `fm_network_report` object.
#' @export
fm_network_report <- function(x, post = NULL) {
  if (inherits(x, "fm_network_fit")) {
    synced <- x$network
    if (!is_fm_network(synced)) {
      stop("`fm_network_fit` does not contain a valid `network`", call. = FALSE)
    }

    pre_median <- x$diagnostics$pre_cycle_median
    report <- fm_network_compare(pre = synced, post = synced)
    report$summary$cycle_reduction_ratio <- x$diagnostics$cycle_reduction_ratio
    report$pre$cycle_summary$median <- pre_median
    return(report)
  }

  if (is_fm_network(x) && is_fm_network(post)) {
    return(fm_network_compare(x, post))
  }

  stop("Provide either `fm_network_fit`, or `pre`/`post` fm_network objects", call. = FALSE)
}

#' Format Network Report as Text Summary
#'
#' @param report `fm_network_report` object.
#'
#' @return Character scalar summary.
#' @export
fm_network_summary <- function(report) {
  if (!inherits(report, "fm_network_report")) {
    stop("`report` must inherit from `fm_network_report`", call. = FALSE)
  }

  paste(
    sprintf("Cycle median (pre/post): %s / %s",
      signif(report$pre$cycle_summary$median, 4),
      signif(report$post$cycle_summary$median, 4)
    ),
    sprintf("Cycle reduction ratio: %s", signif(report$summary$cycle_reduction_ratio, 4)),
    sprintf("Transfer median (pre/post): %s / %s",
      signif(report$pre$transfer_summary$median, 4),
      signif(report$post$transfer_summary$median, 4)
    ),
    sprintf("Transfer reduction ratio: %s", signif(report$summary$transfer_reduction_ratio, 4)),
    sprintf("Coverage median (pre/post): %s / %s",
      signif(report$pre$map_quality_summary$coverage_median, 4),
      signif(report$post$map_quality_summary$coverage_median, 4)
    ),
    sprintf("Continuity median (pre/post): %s / %s",
      signif(report$pre$map_quality_summary$continuity_median, 4),
      signif(report$post$map_quality_summary$continuity_median, 4)
    ),
    sprintf("Coverage note: %s", report$pre$coverage_note),
    sep = "\n"
  )
}

#' @export
print.fm_network_report <- function(x, ...) {
  cat("<fm_network_report>\n")
  cat(fm_network_summary(x), "\n")
  invisible(x)
}
