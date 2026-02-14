operator_eigenpairs <- function(domain) {
  if (!inherits(domain, "fm_domain")) {
    stop("`domain` must inherit from `fm_domain`", call. = FALSE)
  }

  basis <- domain$basis
  if (is.null(basis$vectors) || is.null(basis$values)) {
    stop("Domain basis is required. Run `fm_basis()` first.", call. = FALSE)
  }

  list(
    vectors = as.matrix(basis$vectors),
    values = as.numeric(basis$values)
  )
}

infer_hks_time_grid <- function(evals, n_times, time_range = NULL) {
  if (!is.numeric(n_times) || length(n_times) != 1L || n_times < 1) {
    stop("`n_times` must be a positive scalar", call. = FALSE)
  }

  n_times <- as.integer(n_times)
  lam <- sort(abs(as.numeric(evals)))
  lam <- lam[is.finite(lam) & lam > 0]

  if (length(lam) < 2L) {
    stop("Need at least two positive eigenvalues for HKS time grid", call. = FALSE)
  }

  if (is.null(time_range)) {
    t_min <- 4 * log(10) / max(lam)
    t_max <- 4 * log(10) / min(lam)
  } else {
    if (!is.numeric(time_range) || length(time_range) != 2L || any(time_range <= 0)) {
      stop("`time_range` must be two positive values", call. = FALSE)
    }
    t_min <- min(time_range)
    t_max <- max(time_range)
  }

  exp(seq(log(t_min), log(t_max), length.out = n_times))
}

infer_wks_energy_grid <- function(evals, n_energies, energy_range = NULL) {
  if (!is.numeric(n_energies) || length(n_energies) != 1L || n_energies < 1) {
    stop("`n_energies` must be a positive scalar", call. = FALSE)
  }

  n_energies <- as.integer(n_energies)
  log_ev <- log(abs(as.numeric(evals)))
  log_ev <- log_ev[is.finite(log_ev)]

  if (length(log_ev) < 2L) {
    stop("Need at least two finite log-eigenvalues for WKS energy grid", call. = FALSE)
  }

  if (is.null(energy_range)) {
    e_min <- min(log_ev)
    e_max <- max(log_ev)
  } else {
    if (!is.numeric(energy_range) || length(energy_range) != 2L || any(!is.finite(energy_range))) {
      stop("`energy_range` must contain two finite values", call. = FALSE)
    }
    e_min <- min(energy_range)
    e_max <- max(energy_range)
  }

  seq(e_min, e_max, length.out = n_energies)
}

normalize_descriptor_columns <- function(desc, eps = 1e-12) {
  norms <- sqrt(colSums(desc^2))
  keep <- norms > eps
  if (any(keep)) {
    desc[, keep] <- sweep(desc[, keep, drop = FALSE], 2, norms[keep], "/")
  }
  desc
}

landmark_spectral_response <- function(phi, coeff, landmarks) {
  phi_lm <- phi[landmarks, , drop = FALSE]
  out <- vapply(seq_len(nrow(phi_lm)), function(i) {
    as.numeric(phi %*% (coeff * phi_lm[i, ]))
  }, numeric(nrow(phi)))
  if (is.null(dim(out))) {
    matrix(out, ncol = 1L)
  } else {
    out
  }
}

#' Build Heat Kernel Signature (HKS) Descriptors
#'
#' @param domain `fm_domain` with basis vectors/values.
#' @param n_times Number of diffusion times.
#' @param time_range Optional positive range `c(t_min, t_max)`.
#' @param scaled Whether to normalize each descriptor column.
#' @param landmarks Optional landmark indices for landmark-HKS blocks.
#'
#' @return Numeric matrix with one column per time (or per landmark-time pair).
#' @export
fm_descriptor_hks <- function(domain, n_times = 16, time_range = NULL, scaled = TRUE, landmarks = NULL) {
  eig <- operator_eigenpairs(domain)
  phi <- eig$vectors
  lam <- abs(eig$values)

  if (!is.null(landmarks)) {
    landmarks <- as.integer(landmarks)
    if (any(is.na(landmarks)) || any(landmarks < 1L) || any(landmarks > nrow(phi))) {
      stop("`landmarks` must be valid sample indices", call. = FALSE)
    }
  }

  t_grid <- infer_hks_time_grid(lam, n_times = n_times, time_range = time_range)

  cols <- lapply(t_grid, function(ti) {
    coeff <- exp(-lam * ti)
    if (is.null(landmarks)) {
      matrix(rowSums((phi^2) * rep(coeff, each = nrow(phi))), ncol = 1L)
    } else {
      landmark_spectral_response(phi, coeff = coeff, landmarks = landmarks)
    }
  })

  out <- do.call(cbind, cols)
  if (isTRUE(scaled)) {
    out <- normalize_descriptor_columns(out)
  }
  out
}

#' Build Wave Kernel Signature (WKS) Descriptors
#'
#' @param domain `fm_domain` with basis vectors/values.
#' @param n_energies Number of sampled energies.
#' @param sigma Optional Gaussian width in log-spectrum.
#' @param energy_range Optional energy range on log-spectrum scale.
#' @param scaled Whether to normalize each descriptor column.
#' @param landmarks Optional landmark indices for landmark-WKS blocks.
#'
#' @return Numeric matrix with one column per energy (or per landmark-energy pair).
#' @export
fm_descriptor_wks <- function(domain, n_energies = 16, sigma = NULL, energy_range = NULL, scaled = TRUE, landmarks = NULL) {
  eig <- operator_eigenpairs(domain)
  phi <- eig$vectors
  lam <- abs(eig$values)

  if (!is.null(landmarks)) {
    landmarks <- as.integer(landmarks)
    if (any(is.na(landmarks)) || any(landmarks < 1L) || any(landmarks > nrow(phi))) {
      stop("`landmarks` must be valid sample indices", call. = FALSE)
    }
  }

  e_grid <- infer_wks_energy_grid(lam, n_energies = n_energies, energy_range = energy_range)
  log_lam <- log(pmax(lam, .Machine$double.eps))

  if (is.null(sigma)) {
    sigma <- if (length(e_grid) > 1L) {
      0.5 * abs(e_grid[2] - e_grid[1])
    } else {
      0.5
    }
  }
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("`sigma` must be a positive scalar", call. = FALSE)
  }

  cols <- lapply(e_grid, function(ei) {
    coeff <- exp(-((ei - log_lam)^2) / (2 * sigma^2))
    coeff <- coeff / max(sum(coeff), .Machine$double.eps)

    if (is.null(landmarks)) {
      matrix(rowSums((phi^2) * rep(coeff, each = nrow(phi))), ncol = 1L)
    } else {
      landmark_spectral_response(phi, coeff = coeff, landmarks = landmarks)
    }
  })

  out <- do.call(cbind, cols)
  if (isTRUE(scaled)) {
    out <- normalize_descriptor_columns(out)
  }
  out
}

#' Build Descriptors from a Domain Basis
#'
#' @param domain `fm_domain` with basis/eigenvalues.
#' @param method Descriptor family: `"hks"` or `"wks"`.
#' @param ... Additional method-specific arguments.
#'
#' @return Numeric descriptor matrix.
#' @export
fm_descriptors <- function(domain, method = c("hks", "wks"), ...) {
  method <- match.arg(method)
  switch(
    method,
    hks = fm_descriptor_hks(domain, ...),
    wks = fm_descriptor_wks(domain, ...)
  )
}

mesh_face_adjacency <- function(n, faces) {
  if (is.null(faces)) {
    return(NULL)
  }
  f <- as.matrix(faces)
  if (ncol(f) != 3L) {
    return(NULL)
  }

  A <- matrix(0, nrow = n, ncol = n)
  for (r in seq_len(nrow(f))) {
    tri <- as.integer(f[r, ])
    tri <- tri[tri >= 1L & tri <= n]
    if (length(tri) != 3L) {
      next
    }
    A[tri[1], tri[2]] <- 1
    A[tri[2], tri[1]] <- 1
    A[tri[1], tri[3]] <- 1
    A[tri[3], tri[1]] <- 1
    A[tri[2], tri[3]] <- 1
    A[tri[3], tri[2]] <- 1
  }
  A
}

graph_shortest_paths <- function(adjacency) {
  A <- as.matrix(adjacency)
  if (nrow(A) != ncol(A)) {
    stop("`adjacency` must be square", call. = FALSE)
  }

  n <- nrow(A)
  D <- matrix(Inf, nrow = n, ncol = n)
  diag(D) <- 0

  nz <- which(A > 0, arr.ind = TRUE)
  if (nrow(nz) > 0L) {
    D[nz] <- A[nz]
    D[A > 0 & A < D] <- A[A > 0 & A < D]
  }

  for (k in seq_len(n)) {
    D <- pmin(D, outer(D[, k], D[k, ], "+"))
  }
  D
}

#' Build a Sample Distance Matrix for a Domain
#'
#' @param domain `fm_domain` object.
#' @param method Distance strategy: `"auto"`, `"euclidean"`, or `"graph_shortest"`.
#'
#' @return Dense distance matrix.
#' @export
fm_distance_matrix <- function(domain, method = c("auto", "euclidean", "graph_shortest")) {
  if (!inherits(domain, "fm_domain")) {
    stop("`domain` must inherit from `fm_domain`", call. = FALSE)
  }

  method <- match.arg(method)

  if (method == "auto") {
    method <- switch(
      domain$type,
      mesh = "euclidean",
      pointcloud = "euclidean",
      graph = "graph_shortest",
      if (!is.null(domain$adjacency)) "graph_shortest" else "euclidean"
    )
  }

  if (method == "graph_shortest") {
    A <- domain$adjacency %||% domain$data$adjacency
    if (is.null(A)) {
      stop("No adjacency available for graph shortest-path distance", call. = FALSE)
    }
    return(graph_shortest_paths(A))
  }

  if (domain$type == "mesh") {
    x <- domain$data$vertices
  } else if (domain$type == "pointcloud") {
    x <- domain$data$points
  } else if (!is.null(domain$data) && !is.null(dim(domain$data))) {
    x <- as.matrix(domain$data)
  } else {
    stop("Euclidean distance requires matrix-like sample coordinates/data", call. = FALSE)
  }

  as.matrix(stats::dist(as.matrix(x)))
}

normalize_metric_scale <- function(dist_mat) {
  vals <- as.numeric(dist_mat)
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0L) {
    return(1)
  }
  mean(vals)
}

validate_p2p_pair <- function(p2p, truth, n_source = NULL) {
  p2p <- as.integer(p2p)
  truth <- as.integer(truth)

  if (length(p2p) != length(truth)) {
    stop("`p2p` and `truth` must have the same length", call. = FALSE)
  }
  if (!is.null(n_source)) {
    if (any(p2p < 1L | p2p > n_source) || any(truth < 1L | truth > n_source)) {
      stop("`p2p`/`truth` indices out of source range", call. = FALSE)
    }
  }

  list(p2p = p2p, truth = truth)
}

#' Geodesic-Style Error for Pointwise Maps
#'
#' @param p2p Integer vector mapping target samples to source indices.
#' @param truth Integer vector of ground-truth source indices for each target sample.
#' @param source_distance Source-domain distance matrix (geodesic or surrogate).
#' @param normalize Whether to normalize by mean finite off-diagonal distance.
#'
#' @return List with per-sample errors and summary statistics.
#' @export
fm_eval_geodesic_error <- function(p2p, truth, source_distance, normalize = TRUE) {
  D <- as.matrix(source_distance)
  if (nrow(D) != ncol(D)) {
    stop("`source_distance` must be square", call. = FALSE)
  }

  pair <- validate_p2p_pair(p2p, truth, n_source = nrow(D))
  idx <- cbind(pair$p2p, pair$truth)
  err <- D[idx]

  scale <- if (isTRUE(normalize)) normalize_metric_scale(D) else 1
  out <- list(
    mean = mean(err, na.rm = TRUE),
    median = stats::median(err, na.rm = TRUE),
    max = max(err, na.rm = TRUE),
    normalized_mean = mean(err, na.rm = TRUE) / scale,
    scale = scale,
    per_point = as.numeric(err)
  )
  class(out) <- c("fm_eval_geodesic", class(out))
  out
}

#' Continuity Proxy for Pointwise Maps
#'
#' @param p2p Integer vector mapping target samples to source indices.
#' @param source_distance Source-domain distance matrix.
#' @param target_adjacency Target-domain adjacency matrix.
#' @param normalize Whether to normalize by source-distance scale.
#'
#' @return List with continuity proxy statistics (lower is smoother).
#' @export
fm_eval_continuity <- function(p2p, source_distance, target_adjacency, normalize = TRUE) {
  D <- as.matrix(source_distance)
  A <- as.matrix(target_adjacency)

  if (nrow(D) != ncol(D)) {
    stop("`source_distance` must be square", call. = FALSE)
  }
  if (nrow(A) != ncol(A)) {
    stop("`target_adjacency` must be square", call. = FALSE)
  }

  p2p <- as.integer(p2p)
  if (length(p2p) != nrow(A)) {
    stop("`p2p` length must match target adjacency size", call. = FALSE)
  }
  if (any(p2p < 1L | p2p > nrow(D))) {
    stop("`p2p` indices out of source range", call. = FALSE)
  }

  edges <- which(A > 0, arr.ind = TRUE)
  edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]

  if (nrow(edges) == 0L) {
    out <- list(
      n_edges = 0L,
      mean = NA_real_,
      median = NA_real_,
      normalized_mean = NA_real_,
      scale = if (isTRUE(normalize)) normalize_metric_scale(D) else 1
    )
    class(out) <- c("fm_eval_continuity", class(out))
    return(out)
  }

  dvals <- D[cbind(p2p[edges[, 1]], p2p[edges[, 2]])]
  scale <- if (isTRUE(normalize)) normalize_metric_scale(D) else 1

  out <- list(
    n_edges = nrow(edges),
    mean = mean(dvals, na.rm = TRUE),
    median = stats::median(dvals, na.rm = TRUE),
    normalized_mean = mean(dvals, na.rm = TRUE) / scale,
    scale = scale
  )
  class(out) <- c("fm_eval_continuity", class(out))
  out
}

#' Coverage Metrics for Pointwise Maps
#'
#' @param p2p Integer vector mapping target samples to source indices.
#' @param n_source Source sample count.
#' @param weights Optional non-negative source weights.
#'
#' @return List with support, ratio, weighted ratio, and entropy.
#' @export
fm_eval_coverage <- function(p2p, n_source, weights = NULL) {
  p2p <- as.integer(p2p)
  if (!is.numeric(n_source) || length(n_source) != 1L || n_source < 1) {
    stop("`n_source` must be a positive scalar", call. = FALSE)
  }
  n_source <- as.integer(n_source)

  if (any(p2p < 1L | p2p > n_source)) {
    stop("`p2p` indices out of source range", call. = FALSE)
  }

  tab <- tabulate(p2p, nbins = n_source)
  support <- which(tab > 0)
  ratio <- length(support) / n_source

  if (is.null(weights)) {
    weighted_ratio <- ratio
  } else {
    w <- as.numeric(weights)
    if (length(w) != n_source || any(w < 0) || !any(w > 0)) {
      stop("`weights` must be non-negative with at least one positive value", call. = FALSE)
    }
    weighted_ratio <- sum(w[support]) / sum(w)
  }

  prob <- tab / sum(tab)
  prob <- prob[prob > 0]
  entropy <- -sum(prob * log(prob)) / log(n_source)

  out <- list(
    support = support,
    n_support = length(support),
    ratio = ratio,
    weighted_ratio = weighted_ratio,
    entropy = entropy
  )
  class(out) <- c("fm_eval_coverage", class(out))
  out
}

auto_target_adjacency <- function(domain) {
  if (!inherits(domain, "fm_domain")) {
    return(NULL)
  }

  if (!is.null(domain$adjacency)) {
    return(as.matrix(domain$adjacency))
  }

  if (identical(domain$type, "mesh")) {
    return(mesh_face_adjacency(domain$n_samples, domain$data$faces))
  }

  NULL
}

#' Evaluate a Fitted Functional Map
#'
#' @param fit `fm_fit` object.
#' @param truth Optional ground-truth pointwise map (target -> source).
#' @param source_distance Optional source distance matrix.
#' @param target_adjacency Optional target adjacency matrix.
#' @param normalize Whether to normalize distance-based metrics.
#'
#' @return List with pointwise map and available metric summaries.
#' @export
fm_fit_metrics <- function(
  fit,
  truth = NULL,
  source_distance = NULL,
  target_adjacency = NULL,
  normalize = TRUE
) {
  if (!is_fm_fit(fit)) {
    stop("`fit` must inherit from `fm_fit`", call. = FALSE)
  }

  p2p <- as_p2p(fit)
  source_n <- fit$source$n_samples
  weights <- extract_measure_weights(fit$source)

  out <- list(
    p2p = p2p,
    coverage = fm_eval_coverage(p2p, n_source = source_n, weights = weights)
  )

  need_source_dist <- !is.null(truth) || !is.null(target_adjacency) || !is.null(auto_target_adjacency(fit$target))
  D <- source_distance
  if (is.null(D) && need_source_dist) {
    D <- fm_distance_matrix(fit$source, method = "auto")
  }

  if (!is.null(truth)) {
    out$geodesic <- fm_eval_geodesic_error(
      p2p = p2p,
      truth = truth,
      source_distance = D,
      normalize = normalize
    )
  }

  A_tgt <- target_adjacency
  if (is.null(A_tgt)) {
    A_tgt <- auto_target_adjacency(fit$target)
  }

  if (!is.null(A_tgt) && !is.null(D)) {
    out$continuity <- fm_eval_continuity(
      p2p = p2p,
      source_distance = D,
      target_adjacency = A_tgt,
      normalize = normalize
    )
  }

  out
}
