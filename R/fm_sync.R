network_edge_table <- function(network) {
  keys <- names(network$maps)
  if (is.null(keys) || length(keys) == 0L) {
    return(data.frame(i = character(), j = character(), key = character(), stringsAsFactors = FALSE))
  }

  edges <- lapply(keys, function(k) {
    ij <- parse_edge_key(k)
    data.frame(i = ij[[1]], j = ij[[2]], key = k, stringsAsFactors = FALSE)
  })

  do.call(rbind, edges)
}

truncate_map <- function(M, k) {
  M[seq_len(min(k, nrow(M))), seq_len(min(k, ncol(M))), drop = FALSE]
}

edge_weight_lookup <- function(network, edge_df) {
  if (nrow(edge_df) == 0L) {
    return(numeric())
  }

  vapply(seq_len(nrow(edge_df)), function(ix) {
    ii <- network$index[[edge_df$i[[ix]]]]
    jj <- network$index[[edge_df$j[[ix]]]]
    as.numeric(network$weights[ii, jj])
  }, numeric(1))
}

compute_edge_cycle_stats <- function(network, edge_df) {
  if (nrow(edge_df) == 0L) {
    return(list(
      avg_error = numeric(),
      n_cycles = integer(),
      global_scale = NA_real_
    ))
  }

  cycles <- fm_network_cycles(network)
  if (length(cycles) == 0L) {
    return(list(
      avg_error = rep(1, nrow(edge_df)),
      n_cycles = rep(0L, nrow(edge_df)),
      global_scale = 1
    ))
  }

  cycle_errors <- vapply(cycles, function(cyc) fm_network_cycle_error(network, cyc), numeric(1))
  finite_errors <- cycle_errors[is.finite(cycle_errors)]
  scale <- stats::median(finite_errors)
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }

  edge_err <- rep(0, nrow(edge_df))
  edge_cnt <- rep(0L, nrow(edge_df))

  for (ci in seq_along(cycles)) {
    cyc <- cycles[[ci]]
    err <- cycle_errors[[ci]]
    if (!is.finite(err)) next

    cyc_edges <- list(
      edge_key(cyc[[1]], cyc[[2]]),
      edge_key(cyc[[2]], cyc[[3]]),
      edge_key(cyc[[3]], cyc[[1]])
    )

    for (ek in cyc_edges) {
      idx <- which(edge_df$key == ek)
      if (length(idx) == 1L) {
        edge_err[idx] <- edge_err[idx] + err
        edge_cnt[idx] <- edge_cnt[idx] + 1
      }
    }
  }

  avg_err <- ifelse(edge_cnt > 0, edge_err / edge_cnt, scale)
  list(
    avg_error = avg_err,
    n_cycles = edge_cnt,
    global_scale = scale
  )
}

edge_cycle_penalty <- function(network, edge_df, method = c("exp", "tukey"), tune = 4.685) {
  method <- match.arg(method)
  edge_stats <- compute_edge_cycle_stats(network, edge_df)
  avg_err <- edge_stats$avg_error

  if (length(avg_err) == 0L) {
    return(list(
      weights = numeric(),
      avg_error = numeric(),
      n_cycles = integer(),
      method = method
    ))
  }

  if (method == "exp") {
    scale <- edge_stats$global_scale
    weights <- exp(-avg_err / scale)
  } else {
    # One-sided robust weighting:
    # keep low-error edges near 1 and aggressively downweight high-error outliers.
    center <- stats::median(avg_err)
    scale <- stats::mad(avg_err, center = center, constant = 1.4826)
    if (!is.finite(scale) || scale <= 0) {
      scale <- stats::sd(avg_err)
    }
    if (!is.finite(scale) || scale <= 0) {
      scale <- 1
    }

    u <- pmax(0, (avg_err - center) / scale) / tune
    weights <- ifelse(u < 1, (1 - u^2)^2, 0)
    weights <- pmax(weights, 1e-4)
  }

  list(
    weights = weights,
    avg_error = avg_err,
    n_cycles = edge_stats$n_cycles,
    method = method
  )
}

project_orthogonal <- function(M) {
  s <- svd(M)
  s$u %*% t(s$v)
}

compute_cycle_median <- function(network) {
  scores <- fm_network_cycle_scores(network)
  if (nrow(scores) == 0L) {
    return(NA_real_)
  }
  stats::median(scores$error, na.rm = TRUE)
}

#' Synchronize Maps in a Functional-Map Network
#'
#' @param network `fm_network` object.
#' @param mode Synchronization mode: `"adjacency"`, `"cycle"`, or `"robust"`.
#' @param nit Maximum synchronization iterations.
#' @param tol Convergence tolerance for latent transforms.
#' @param anchor Optional anchor node name (defaults to first domain).
#' @param verbose Emit per-iteration diagnostics.
#'
#' @return `fm_network_fit` object with synchronized maps and diagnostics.
#' @export
fm_sync <- function(
  network,
  mode = c("adjacency", "cycle", "robust"),
  nit = 20,
  tol = 1e-6,
  anchor = NULL,
  verbose = FALSE
) {
  if (!is_fm_network(network)) {
    stop("`network` must inherit from `fm_network`", call. = FALSE)
  }

  mode <- match.arg(mode)

  edge_df <- network_edge_table(network)
  if (nrow(edge_df) == 0L) {
    stop("`network` has no maps to synchronize", call. = FALSE)
  }

  map_dims <- vapply(network$maps, function(m) min(dim(as_map_matrix(m))), integer(1))
  k <- min(map_dims)
  if (k < 1L) {
    stop("Maps must have at least one dimension", call. = FALSE)
  }

  node_names <- names(network$domains)
  anchor <- anchor %||% node_names[[1]]
  if (!(anchor %in% node_names)) {
    stop("`anchor` must be one of the network node names", call. = FALSE)
  }

  W_adj <- edge_weight_lookup(network, edge_df)
  cycle_info <- list(
    weights = rep(1, length(W_adj)),
    avg_error = rep(NA_real_, length(W_adj)),
    n_cycles = rep(0L, length(W_adj)),
    method = "none"
  )

  if (mode == "cycle") {
    cycle_info <- edge_cycle_penalty(network, edge_df, method = "exp")
  } else if (mode == "robust") {
    cycle_info <- edge_cycle_penalty(network, edge_df, method = "tukey")
  }

  if (mode %in% c("cycle", "robust")) {
    W_robust <- cycle_info$weights
  } else {
    W_robust <- rep(1, length(W_adj))
  }
  W <- pmax(W_adj * W_robust, 1e-8)

  X <- stats::setNames(replicate(length(node_names), diag(k), simplify = FALSE), node_names)

  delta_hist <- numeric(nit)

  for (iter in seq_len(nit)) {
    X_prev <- X

    for (node in node_names) {
      if (identical(node, anchor)) {
        X[[node]] <- diag(k)
        next
      }

      proposals <- list()
      prop_w <- numeric()

      for (r in seq_len(nrow(edge_df))) {
        i <- edge_df$i[[r]]
        j <- edge_df$j[[r]]
        Cij <- truncate_map(as_map_matrix(network$maps[[edge_df$key[[r]]]]), k)

        if (j == node) {
          proposals[[length(proposals) + 1L]] <- Cij %*% X_prev[[i]]
          prop_w <- c(prop_w, W[[r]])
        }

        if (i == node) {
          proposals[[length(proposals) + 1L]] <- X_prev[[j]] %*% matrix_pseudoinverse(Cij)
          prop_w <- c(prop_w, W[[r]])
        }
      }

      if (length(proposals) == 0L) {
        next
      }

      avg <- matrix(0, nrow = k, ncol = k)
      for (pi in seq_along(proposals)) {
        avg <- avg + prop_w[[pi]] * proposals[[pi]]
      }
      avg <- avg / sum(prop_w)

      X[[node]] <- project_orthogonal(avg)
    }

    deltas <- vapply(node_names, function(nm) {
      max(abs(X[[nm]] - X_prev[[nm]]))
    }, numeric(1))
    delta_hist[[iter]] <- max(deltas)

    if (isTRUE(verbose)) {
      message(sprintf("[sync:%s] iter=%d delta=%.3e", mode, iter, delta_hist[[iter]]))
    }

    if (delta_hist[[iter]] <= tol) {
      delta_hist <- delta_hist[seq_len(iter)]
      break
    }
  }

  net_sync <- network
  for (r in seq_len(nrow(edge_df))) {
    i <- edge_df$i[[r]]
    j <- edge_df$j[[r]]

    Csync <- X[[j]] %*% matrix_pseudoinverse(X[[i]])
    net_sync$maps[[edge_df$key[[r]]]] <- Csync
  }

  pre_cycle <- compute_cycle_median(network)
  post_cycle <- compute_cycle_median(net_sync)

  diagnostics <- list(
    mode = mode,
    k = k,
    nit_done = length(delta_hist),
    delta = delta_hist,
    edge_weights = data.frame(
      key = edge_df$key,
      adjacency_weight = W_adj,
      cycle_avg_error = cycle_info$avg_error,
      cycle_count = cycle_info$n_cycles,
      robust_weight = W_robust,
      combined_weight = W,
      effective_weight = W,
      robust_method = cycle_info$method,
      stringsAsFactors = FALSE
    ),
    pre_cycle_median = pre_cycle,
    post_cycle_median = post_cycle,
    cycle_reduction_ratio = if (is.finite(pre_cycle) && pre_cycle > 0 && is.finite(post_cycle)) {
      (pre_cycle - post_cycle) / pre_cycle
    } else {
      NA_real_
    }
  )

  structure(
    list(
      network = net_sync,
      transforms = X,
      diagnostics = diagnostics
    ),
    class = "fm_network_fit"
  )
}

#' @export
print.fm_network_fit <- function(x, ...) {
  cat(sprintf("<fm_network_fit> mode=%s nit=%d cycle_reduction=%s\n",
    x$diagnostics$mode,
    x$diagnostics$nit_done,
    as.character(signif(x$diagnostics$cycle_reduction_ratio, 4))
  ))
  invisible(x)
}
