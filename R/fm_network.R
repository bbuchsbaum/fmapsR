edge_key <- function(i, j) {
  paste0(i, "->", j)
}

parse_edge_key <- function(key) {
  parts <- strsplit(key, "->", fixed = TRUE)[[1]]
  if (length(parts) != 2L) {
    stop("Invalid edge key", call. = FALSE)
  }
  c(parts[1], parts[2])
}

as_map_matrix <- function(map) {
  if (is_fm_fit(map)) {
    return(map$C)
  }
  as.matrix(map)
}

ensure_domain_names <- function(domains) {
  if (!is.list(domains) || length(domains) < 2L) {
    stop("`domains` must be a list with at least two domain objects", call. = FALSE)
  }

  if (!all(vapply(domains, is_fm_domain, logical(1)))) {
    stop("All `domains` entries must inherit from `fm_domain`", call. = FALSE)
  }

  nm <- names(domains)
  if (is.null(nm) || any(!nzchar(nm))) {
    nm <- as.character(seq_along(domains))
    names(domains) <- nm
  }

  domains
}

make_empty_weights <- function(n) {
  Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    x = numeric(),
    dims = c(n, n),
    dimnames = list(NULL, NULL)
  )
}

#' Create a Functional-Map Network
#'
#' @param domains Named list of `fm_domain` objects.
#' @param maps Optional named list of maps keyed as `"i->j"`.
#' @param directed Whether to treat edges as directed.
#'
#' @return `fm_network` object.
#' @export
fm_network <- function(domains, maps = NULL, directed = TRUE) {
  domains <- ensure_domain_names(domains)
  n <- length(domains)

  net <- structure(
    list(
      domains = domains,
      maps = list(),
      weights = make_empty_weights(n),
      directed = isTRUE(directed),
      index = stats::setNames(seq_len(n), names(domains))
    ),
    class = "fm_network"
  )

  if (!is.null(maps)) {
    if (is.null(names(maps)) || any(!grepl("->", names(maps), fixed = TRUE))) {
      stop("`maps` must be a named list keyed as 'i->j'", call. = FALSE)
    }

    for (nm in names(maps)) {
      ij <- parse_edge_key(nm)
      net <- fm_network_add_map(net, ij[1], ij[2], maps[[nm]])
    }
  }

  net
}

#' Add a Map to an FM Network
#'
#' @param network `fm_network` object.
#' @param i Source domain name.
#' @param j Target domain name.
#' @param map `fm_fit` object or matrix.
#' @param weight Edge weight.
#' @param add_reverse For undirected networks, add reverse map automatically.
#'
#' @return Updated `fm_network` object.
#' @export
fm_network_add_map <- function(network, i, j, map, weight = 1, add_reverse = !network$directed) {
  if (!inherits(network, "fm_network")) {
    stop("`network` must inherit from `fm_network`", call. = FALSE)
  }

  if (!(i %in% names(network$domains)) || !(j %in% names(network$domains))) {
    stop("`i` and `j` must be valid domain names in network", call. = FALSE)
  }

  key <- edge_key(i, j)
  network$maps[[key]] <- map

  ii <- network$index[[i]]
  jj <- network$index[[j]]
  network$weights[ii, jj] <- as.numeric(weight)

  if (isTRUE(add_reverse) && i != j) {
    rev_key <- edge_key(j, i)
    network$maps[[rev_key]] <- if (is_fm_fit(map)) fm_inverse(map) else t(as.matrix(map))
    network$weights[jj, ii] <- as.numeric(weight)
  }

  network
}

#' Get a Map from an FM Network
#'
#' @param network `fm_network` object.
#' @param i Source domain name.
#' @param j Target domain name.
#'
#' @return `fm_fit` object or matrix map.
#' @export
fm_network_get_map <- function(network, i, j) {
  if (!inherits(network, "fm_network")) {
    stop("`network` must inherit from `fm_network`", call. = FALSE)
  }

  network$maps[[edge_key(i, j)]]
}

normalize_edge_table <- function(domains, edges = NULL, directed = TRUE) {
  nodes <- names(domains)

  if (is.null(edges)) {
    rows <- list()
    idx <- 1L
    for (i in nodes) {
      for (j in nodes) {
        if (identical(i, j)) next
        if (!isTRUE(directed) && which(nodes == i) >= which(nodes == j)) next
        rows[[idx]] <- data.frame(i = i, j = j, stringsAsFactors = FALSE)
        idx <- idx + 1L
      }
    }
    return(do.call(rbind, rows))
  }

  if (is.data.frame(edges) || is.matrix(edges)) {
    edge_df <- as.data.frame(edges, stringsAsFactors = FALSE)
    if (!all(c("i", "j") %in% names(edge_df)) && ncol(edge_df) >= 2L) {
      names(edge_df)[1:2] <- c("i", "j")
    }
  } else if (is.character(edges)) {
    out <- lapply(edges, function(key) {
      ij <- parse_edge_key(key)
      data.frame(i = ij[[1]], j = ij[[2]], stringsAsFactors = FALSE)
    })
    edge_df <- do.call(rbind, out)
  } else {
    stop("`edges` must be NULL, character keys, matrix, or data frame", call. = FALSE)
  }

  if (!all(c("i", "j") %in% names(edge_df))) {
    stop("`edges` must provide `i` and `j` columns", call. = FALSE)
  }

  edge_df$i <- as.character(edge_df$i)
  edge_df$j <- as.character(edge_df$j)

  if (any(!(edge_df$i %in% nodes)) || any(!(edge_df$j %in% nodes))) {
    stop("All edge endpoints must match names in `domains`", call. = FALSE)
  }
  if (any(edge_df$i == edge_df$j)) {
    stop("Self-edges are not supported", call. = FALSE)
  }

  if (!isTRUE(directed)) {
    key <- paste(pmin(edge_df$i, edge_df$j), pmax(edge_df$i, edge_df$j), sep = "::")
    edge_df <- edge_df[!duplicated(key), c("i", "j"), drop = FALSE]
  }

  edge_df
}

normalize_pair_batch_size <- function(pair_batch_size, n_edges) {
  if (is.null(pair_batch_size)) {
    return(n_edges)
  }
  if (!is.numeric(pair_batch_size) || length(pair_batch_size) != 1L || pair_batch_size < 1) {
    stop("`pair_batch_size` must be NULL or a positive scalar", call. = FALSE)
  }
  as.integer(min(pair_batch_size, n_edges))
}

#' Fit Pairwise Maps for a Domain Collection and Build a Network
#'
#' @param domains Named list of `fm_domain` objects.
#' @param descriptors Named list of descriptor matrices keyed by domain name.
#' @param edges Optional edge specification (`NULL`, `"i->j"` keys, or two-column
#'   matrix/data.frame with `i` and `j`).
#' @param directed Whether the resulting network is directed.
#' @param pair_batch_size Optional number of edges processed per batch.
#' @param verbose Emit simple per-batch progress.
#' @param ... Additional arguments passed to [fm_match()] (e.g.
#'   `descriptor_batch_size`, `optimizer`, `cg_maxit`).
#'
#' @return `fm_network` with fitted pairwise maps as edges.
#' @export
fm_network_match <- function(
  domains,
  descriptors,
  edges = NULL,
  directed = TRUE,
  pair_batch_size = NULL,
  verbose = FALSE,
  ...
) {
  domains <- ensure_domain_names(domains)
  node_names <- names(domains)

  if (!is.list(descriptors) || is.null(names(descriptors))) {
    stop("`descriptors` must be a named list keyed by domain name", call. = FALSE)
  }
  if (!all(node_names %in% names(descriptors))) {
    stop("`descriptors` must include entries for every domain", call. = FALSE)
  }

  for (nm in node_names) {
    desc <- as.matrix(descriptors[[nm]])
    if (nrow(desc) != domains[[nm]]$n_samples) {
      stop(sprintf("Descriptor row count mismatch for domain `%s`", nm), call. = FALSE)
    }
    descriptors[[nm]] <- desc
  }

  edge_df <- normalize_edge_table(domains, edges = edges, directed = directed)
  if (nrow(edge_df) == 0L) {
    return(fm_network(domains, directed = directed))
  }

  batch_size <- normalize_pair_batch_size(pair_batch_size, n_edges = nrow(edge_df))
  batches <- split(seq_len(nrow(edge_df)), ceiling(seq_len(nrow(edge_df)) / batch_size))

  net <- fm_network(domains, directed = directed)

  for (bi in seq_along(batches)) {
    idx <- batches[[bi]]
    if (isTRUE(verbose)) {
      message(sprintf("[network-match] batch %d/%d edges=%d", bi, length(batches), length(idx)))
    }

    for (ei in idx) {
      i <- edge_df$i[[ei]]
      j <- edge_df$j[[ei]]
      fit <- fm_match(
        source = domains[[i]],
        target = domains[[j]],
        descriptors = list(
          source = descriptors[[i]],
          target = descriptors[[j]]
        ),
        ...
      )

      net <- fm_network_add_map(net, i, j, fit, weight = 1, add_reverse = !isTRUE(directed))
    }
  }

  net
}

#' Extract 3-Cycles from Network
#'
#' @param network `fm_network` object.
#'
#' @return List of 3-cycles. For directed networks each cycle is ordered `c(i,j,k)`
#'   satisfying edges `i->j`, `j->k`, and `k->i`.
#' @export
fm_network_cycles <- function(network) {
  if (!inherits(network, "fm_network")) {
    stop("`network` must inherit from `fm_network`", call. = FALSE)
  }

  nodes <- names(network$domains)
  cycles <- list()

  if (isTRUE(network$directed)) {
    idx <- 1L
    for (i in seq_along(nodes)) {
      for (j in seq_along(nodes)) {
        for (k in seq_along(nodes)) {
          if (length(unique(c(i, j, k))) < 3L) next

          ni <- nodes[[i]]
          nj <- nodes[[j]]
          nk <- nodes[[k]]

          if (!is.null(fm_network_get_map(network, ni, nj)) &&
              !is.null(fm_network_get_map(network, nj, nk)) &&
              !is.null(fm_network_get_map(network, nk, ni))) {
            if (ni == min(c(ni, nj, nk))) {
              cycles[[idx]] <- c(ni, nj, nk)
              idx <- idx + 1L
            }
          }
        }
      }
    }

    return(cycles)
  }

  idx <- 1L
  for (i in seq_len(length(nodes) - 2L)) {
    for (j in seq.int(i + 1L, length(nodes) - 1L)) {
      for (k in seq.int(j + 1L, length(nodes))) {
        ni <- nodes[[i]]
        nj <- nodes[[j]]
        nk <- nodes[[k]]

        eij <- !is.null(fm_network_get_map(network, ni, nj)) || !is.null(fm_network_get_map(network, nj, ni))
        ejk <- !is.null(fm_network_get_map(network, nj, nk)) || !is.null(fm_network_get_map(network, nk, nj))
        eki <- !is.null(fm_network_get_map(network, nk, ni)) || !is.null(fm_network_get_map(network, ni, nk))

        if (eij && ejk && eki) {
          cycles[[idx]] <- c(ni, nj, nk)
          idx <- idx + 1L
        }
      }
    }
  }

  cycles
}

#' Compute Directed 3-Cycle Error
#'
#' @param network `fm_network` object.
#' @param cycle Character vector of length 3 defining `c(i,j,k)`.
#'
#' @return Frobenius norm cycle error.
#' @export
fm_network_cycle_error <- function(network, cycle) {
  if (!inherits(network, "fm_network")) {
    stop("`network` must inherit from `fm_network`", call. = FALSE)
  }
  if (length(cycle) != 3L) {
    stop("`cycle` must have length 3", call. = FALSE)
  }

  i <- cycle[[1]]
  j <- cycle[[2]]
  k <- cycle[[3]]

  Cij <- fm_network_get_map(network, i, j)
  Cjk <- fm_network_get_map(network, j, k)
  Cki <- fm_network_get_map(network, k, i)

  if (is.null(Cij) || is.null(Cjk) || is.null(Cki)) {
    return(NA_real_)
  }

  Mij <- as_map_matrix(Cij)
  Mjk <- as_map_matrix(Cjk)
  Mki <- as_map_matrix(Cki)

  if (nrow(Mjk) != ncol(Mij) || nrow(Mki) != ncol(Mjk)) {
    return(NA_real_)
  }

  composed <- Mki %*% Mjk %*% Mij
  I <- diag(1, nrow = nrow(composed), ncol = ncol(composed))
  norm(composed - I, type = "F")
}

#' Score All 3-Cycles in a Network
#'
#' @param network `fm_network` object.
#'
#' @return Data frame of cycles and errors.
#' @export
fm_network_cycle_scores <- function(network) {
  cycles <- fm_network_cycles(network)

  if (length(cycles) == 0L) {
    return(data.frame(i = character(), j = character(), k = character(), error = numeric()))
  }

  rows <- lapply(cycles, function(cyc) {
    data.frame(
      i = cyc[[1]],
      j = cyc[[2]],
      k = cyc[[3]],
      error = fm_network_cycle_error(network, cyc),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Check fm_network Class
#'
#' @param x Object to test.
#'
#' @return Logical scalar.
#' @export
is_fm_network <- function(x) {
  inherits(x, "fm_network")
}

#' @export
print.fm_network <- function(x, ...) {
  cat(sprintf("<fm_network> nodes=%d directed=%s edges=%d\n",
    length(x$domains),
    if (isTRUE(x$directed)) "TRUE" else "FALSE",
    length(x$maps)
  ))
  invisible(x)
}
