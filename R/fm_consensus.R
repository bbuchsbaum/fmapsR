compute_all_pair_consensus_maps <- function(alignments) {
  nodes <- names(alignments)
  maps <- list()

  for (i in nodes) {
    for (j in nodes) {
      if (identical(i, j)) next
      maps[[edge_key(i, j)]] <- alignments[[j]] %*% matrix_pseudoinverse(alignments[[i]])
    }
  }

  maps
}

edge_confidence_table <- function(net_fit, node_names) {
  keys <- unlist(lapply(node_names, function(i) {
    vapply(node_names[node_names != i], function(j) edge_key(i, j), character(1))
  }), use.names = FALSE)

  edge_weights <- net_fit$diagnostics$edge_weights
  if (is.null(edge_weights) || nrow(edge_weights) == 0L) {
    return(data.frame(
      i = character(),
      j = character(),
      key = character(),
      confidence = numeric(),
      source = character(),
      stringsAsFactors = FALSE
    ))
  }

  conf_col <- if ("effective_weight" %in% names(edge_weights)) "effective_weight" else "combined_weight"
  direct_conf <- stats::setNames(as.numeric(edge_weights[[conf_col]]), edge_weights$key)
  max_conf <- max(direct_conf[is.finite(direct_conf)], na.rm = TRUE)
  if (!is.finite(max_conf) || max_conf <= 0) {
    max_conf <- 1
  }
  direct_conf <- direct_conf / max_conf

  fallback <- stats::median(direct_conf[is.finite(direct_conf)], na.rm = TRUE)
  if (!is.finite(fallback)) {
    fallback <- 1
  }

  rows <- lapply(keys, function(key) {
    ij <- parse_edge_key(key)
    i <- ij[[1]]
    j <- ij[[2]]

    direct_val <- direct_conf[key]
    if (!is.na(direct_val)) {
      return(data.frame(
        i = i,
        j = j,
        key = key,
        confidence = as.numeric(direct_val),
        source = "direct",
        stringsAsFactors = FALSE
      ))
    }

    candidates <- numeric()
    for (k in node_names) {
      if (identical(k, i) || identical(k, j)) next
      c1 <- direct_conf[edge_key(i, k)]
      c2 <- direct_conf[edge_key(k, j)]
      if (is.na(c1) || is.na(c2)) next
      candidates <- c(candidates, min(c1, c2))
    }

    if (length(candidates) > 0L) {
      conf <- max(candidates)
      src <- "two-hop"
    } else {
      conf <- fallback
      src <- "fallback"
    }

    data.frame(
      i = i,
      j = j,
      key = key,
      confidence = as.numeric(conf),
      source = src,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

summarize_confidence <- function(conf_df) {
  if (nrow(conf_df) == 0L) {
    return(list(
      n = 0L,
      median = NA_real_,
      mean = NA_real_,
      ci = c(lower = NA_real_, upper = NA_real_)
    ))
  }

  vals <- conf_df$confidence
  list(
    n = length(vals),
    median = stats::median(vals, na.rm = TRUE),
    mean = mean(vals, na.rm = TRUE),
    ci = stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, names = TRUE)
  )
}

consensus_edge_weights <- function(network, edge_df, edge_weights = NULL) {
  w <- edge_weight_lookup(network, edge_df)

  if (!is.null(edge_weights) && nrow(edge_weights) > 0L && "key" %in% names(edge_weights)) {
    wcol <- if ("effective_weight" %in% names(edge_weights)) {
      "effective_weight"
    } else if ("combined_weight" %in% names(edge_weights)) {
      "combined_weight"
    } else {
      NULL
    }

    if (!is.null(wcol)) {
      lookup <- stats::setNames(as.numeric(edge_weights[[wcol]]), edge_weights$key)
      override <- lookup[edge_df$key]
      keep <- is.finite(override)
      w[keep] <- override[keep]
    }
  }

  w[!is.finite(w) | w <= 0] <- 1
  w
}

consensus_map_objective <- function(network, edge_df, weights, alignments, k) {
  if (nrow(edge_df) == 0L) {
    return(NA_real_)
  }

  val <- 0
  wt <- 0

  for (r in seq_len(nrow(edge_df))) {
    i <- edge_df$i[[r]]
    j <- edge_df$j[[r]]
    Cij <- truncate_map(as_map_matrix(network$maps[[edge_df$key[[r]]]]), k)
    E <- alignments[[j]] %*% matrix_pseudoinverse(alignments[[i]]) - Cij
    w <- weights[[r]]
    val <- val + w * sum(E * E)
    wt <- wt + w
  }

  val / max(wt, .Machine$double.eps)
}

build_consensus_network <- function(network, alignments) {
  maps <- compute_all_pair_consensus_maps(alignments)
  out_net <- fm_network(network$domains, directed = TRUE)

  for (nm in names(maps)) {
    ij <- parse_edge_key(nm)
    out_net <- fm_network_add_map(out_net, ij[[1]], ij[[2]], maps[[nm]], weight = 1)
  }

  list(
    maps = maps,
    network = out_net,
    cycle_median = compute_cycle_median(out_net)
  )
}

consensus_sync_transforms <- function(net_fit, k) {
  alignments <- net_fit$transforms
  if (length(alignments) == 0L) {
    stop("No alignment operators available for consensus", call. = FALSE)
  }

  align_trunc <- lapply(alignments, function(A) A[seq_len(k), seq_len(k), drop = FALSE])
  list(
    alignments = align_trunc,
    diagnostics = list(
      method = "sync_transforms",
      solver = "fm_sync",
      k = k,
      nit_done = as.integer(net_fit$diagnostics$nit_done %||% 0L),
      delta = as.numeric(net_fit$diagnostics$delta %||% numeric()),
      converged = TRUE,
      objective_pre = NA_real_,
      objective_post = NA_real_,
      anchor = names(align_trunc)[[1]],
      eigenvalues = NULL,
      fallback_used = FALSE
    )
  )
}

consensus_latent_clb <- function(network, edge_weights = NULL, k, nit = 10, tol = 1e-6, anchor = NULL) {
  if (!is.numeric(nit) || length(nit) != 1L || nit < 0) {
    stop("`latent_nit` must be a non-negative scalar", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || tol <= 0) {
    stop("`latent_tol` must be a positive scalar", call. = FALSE)
  }

  nit <- as.integer(nit)
  edge_df <- network_edge_table(network)
  if (nrow(edge_df) == 0L) {
    stop("`network` has no maps for latent consensus", call. = FALSE)
  }

  node_names <- names(network$domains)
  n <- length(node_names)
  nk <- n * k

  w <- consensus_edge_weights(network, edge_df, edge_weights = edge_weights)

  M <- matrix(0, nrow = nk, ncol = nk)
  for (r in seq_len(nrow(edge_df))) {
    i <- edge_df$i[[r]]
    j <- edge_df$j[[r]]
    ii <- network$index[[i]]
    jj <- network$index[[j]]

    bi <- ((ii - 1L) * k + 1L):(ii * k)
    bj <- ((jj - 1L) * k + 1L):(jj * k)

    Cij <- truncate_map(as_map_matrix(network$maps[[edge_df$key[[r]]]]), k)
    wr <- w[[r]]

    M[bj, bi] <- M[bj, bi] + wr * Cij
    M[bi, bj] <- M[bi, bj] + wr * t(Cij)
  }

  M <- (M + t(M)) / 2
  eig <- eigen(M, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)[seq_len(k)]
  U <- eig$vectors[, ord, drop = FALSE]

  X <- stats::setNames(vector("list", n), node_names)
  for (idx in seq_along(node_names)) {
    block <- ((idx - 1L) * k + 1L):(idx * k)
    X[[node_names[[idx]]]] <- project_orthogonal(U[block, , drop = FALSE])
  }

  anchor <- anchor %||% node_names[[1]]
  if (!(anchor %in% node_names)) {
    stop("`anchor` must be one of the network node names", call. = FALSE)
  }

  R <- t(X[[anchor]])
  for (nm in node_names) {
    X[[nm]] <- X[[nm]] %*% R
  }
  X[[anchor]] <- diag(k)

  obj_pre <- consensus_map_objective(network, edge_df, w, X, k)

  delta_hist <- numeric(0)
  if (nit > 0L) {
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
            prop_w <- c(prop_w, w[[r]])
          }

          if (i == node) {
            proposals[[length(proposals) + 1L]] <- X_prev[[j]] %*% matrix_pseudoinverse(Cij)
            prop_w <- c(prop_w, w[[r]])
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

      if (delta_hist[[iter]] <= tol) {
        delta_hist <- delta_hist[seq_len(iter)]
        break
      }
    }
  }

  obj_post <- consensus_map_objective(network, edge_df, w, X, k)

  list(
    alignments = X,
    diagnostics = list(
      method = "latent_clb",
      solver = "spectral_procrustes",
      k = k,
      nit_done = length(delta_hist),
      delta = delta_hist,
      converged = if (length(delta_hist) == 0L) TRUE else tail(delta_hist, 1L) <= tol,
      objective_pre = obj_pre,
      objective_post = obj_post,
      anchor = anchor,
      eigenvalues = eig$values[ord],
      fallback_used = FALSE
    )
  )
}

#' Build Consensus Alignment from an FM Network
#'
#' @param x `fm_network` or `fm_network_fit` object.
#' @param mode Sync mode used when `x` is a plain network.
#' @param nit Sync iterations used when `x` is a plain network.
#' @param method Consensus method: `"latent_clb"` (default) or
#'   `"sync_transforms"`.
#' @param latent_nit Iterations for latent Procrustes polishing.
#' @param latent_tol Convergence tolerance for latent polishing.
#' @param anchor Optional anchor node for latent consensus orientation.
#'
#' @return `fm_consensus` object.
#' @export
fm_consensus <- function(
  x,
  mode = "cycle",
  nit = 20,
  method = c("latent_clb", "sync_transforms"),
  latent_nit = 10,
  latent_tol = 1e-6,
  anchor = NULL
) {
  method <- match.arg(method)

  net_fit <- if (inherits(x, "fm_network_fit")) {
    x
  } else if (is_fm_network(x)) {
    fm_sync(x, mode = mode, nit = nit)
  } else {
    stop("`x` must be `fm_network` or `fm_network_fit`", call. = FALSE)
  }

  network <- net_fit$network

  if (length(network$maps) == 0L) {
    stop("`network` has no maps to build consensus", call. = FALSE)
  }

  if (method == "sync_transforms" && length(net_fit$transforms) == 0L) {
    stop("No sync transforms available; use `method = \"latent_clb\"`", call. = FALSE)
  }

  map_dims <- vapply(network$maps, function(m) min(dim(as_map_matrix(m))), integer(1))
  k <- min(map_dims)
  if (k < 1L) {
    stop("Consensus requires maps with at least one spectral dimension", call. = FALSE)
  }

  if (method == "sync_transforms") {
    align_out <- consensus_sync_transforms(net_fit, k = k)
  } else {
    align_out <- consensus_latent_clb(
      network = network,
      edge_weights = net_fit$diagnostics$edge_weights %||% NULL,
      k = k,
      nit = latent_nit,
      tol = latent_tol,
      anchor = anchor
    )

    if (length(net_fit$transforms) > 0L) {
      sync_out <- consensus_sync_transforms(net_fit, k = k)
      latent_cons <- build_consensus_network(network, align_out$alignments)
      sync_cons <- build_consensus_network(network, sync_out$alignments)

      if (is.finite(sync_cons$cycle_median) && is.finite(latent_cons$cycle_median) &&
          latent_cons$cycle_median > sync_cons$cycle_median) {
        align_out <- sync_out
        align_out$diagnostics$method <- "latent_clb"
        align_out$diagnostics$solver <- "sync_transform_fallback"
        align_out$diagnostics$fallback_used <- TRUE
        align_out$diagnostics$objective_pre <- latent_cons$cycle_median
        align_out$diagnostics$objective_post <- sync_cons$cycle_median
      }
    }
  }

  align_trunc <- align_out$alignments
  consensus_obj <- build_consensus_network(network, align_trunc)

  confidence_df <- edge_confidence_table(net_fit, names(align_trunc))
  confidence_summary <- summarize_confidence(confidence_df)

  diagnostics <- list(
    method = method,
    latent = align_out$diagnostics,
    pre_cycle_median = net_fit$diagnostics$pre_cycle_median,
    synced_cycle_median = net_fit$diagnostics$post_cycle_median,
    consensus_cycle_median = consensus_obj$cycle_median,
    confidence_interval = confidence_summary$ci
  )

  structure(
    list(
      latent_basis = diag(k),
      alignments = align_trunc,
      maps = consensus_obj$maps,
      network = consensus_obj$network,
      diagnostics = diagnostics,
      uncertainty = list(
        edge_confidence = confidence_df,
        summary = confidence_summary
      )
    ),
    class = "fm_consensus"
  )
}

#' Save Consensus Object
#'
#' @param object `fm_consensus` object.
#' @param path Output file path.
#'
#' @return Input `path`, invisibly.
#' @export
fm_consensus_save <- function(object, path) {
  if (!inherits(object, "fm_consensus")) {
    stop("`object` must inherit from `fm_consensus`", call. = FALSE)
  }

  saveRDS(object, path)
  invisible(path)
}

#' Load Consensus Object
#'
#' @param path Path created by [fm_consensus_save()].
#'
#' @return `fm_consensus` object.
#' @export
fm_consensus_load <- function(path) {
  obj <- readRDS(path)
  if (!inherits(obj, "fm_consensus")) {
    stop("Serialized object is not an `fm_consensus`", call. = FALSE)
  }
  obj
}

#' Select Consensus Maps by Confidence Threshold
#'
#' @param object `fm_consensus` object.
#' @param min_confidence Minimum confidence threshold in `[0, 1]`.
#'
#' @return Named list of consensus maps passing threshold.
#' @export
fm_consensus_select_maps <- function(object, min_confidence = 0.5) {
  if (!inherits(object, "fm_consensus")) {
    stop("`object` must inherit from `fm_consensus`", call. = FALSE)
  }

  if (!is.numeric(min_confidence) || length(min_confidence) != 1L ||
      min_confidence < 0 || min_confidence > 1) {
    stop("`min_confidence` must be a scalar in [0, 1]", call. = FALSE)
  }

  conf_df <- object$uncertainty$edge_confidence
  if (is.null(conf_df) || nrow(conf_df) == 0L) {
    return(list())
  }

  keep_keys <- conf_df$key[is.finite(conf_df$confidence) & conf_df$confidence >= min_confidence]
  object$maps[keep_keys]
}

#' @export
print.fm_consensus <- function(x, ...) {
  cat(sprintf("<fm_consensus> nodes=%d k=%d\n", length(x$alignments), nrow(x$latent_basis)))
  cat(sprintf("  consensus cycle median: %s\n", as.character(signif(x$diagnostics$consensus_cycle_median, 4))))
  cat(sprintf("  method: %s\n", x$diagnostics$method))
  if (!is.null(x$diagnostics$latent$solver)) {
    cat(sprintf("  latent solver: %s\n", x$diagnostics$latent$solver))
  }
  if (!is.null(x$uncertainty$summary)) {
    cat(sprintf("  confidence median: %s\n", as.character(signif(x$uncertainty$summary$median, 4))))
  }
  invisible(x)
}
