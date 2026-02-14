#' Check Functional Map Fit Class
#'
#' @param x Object to check.
#'
#' @return Logical scalar.
#' @export
is_fm_fit <- function(x) {
  inherits(x, "fm_fit")
}

nearest_neighbor_index <- function(reference, query) {
  reference <- as.matrix(reference)
  query <- as.matrix(query)

  if (ncol(reference) != ncol(query)) {
    stop("Embedding dimensions must match for nearest-neighbor lookup", call. = FALSE)
  }

  if ((nrow(reference) * nrow(query)) > 5e7) {
    stop(
      "Pointwise conversion is too large for dense nearest-neighbor search. Add an ANN backend.",
      call. = FALSE
    )
  }

  if (exists("fm_nearest_neighbor_index_cpp", mode = "function")) {
    return(as.integer(fm_nearest_neighbor_index_cpp(reference, query)))
  }

  ref_norm <- rowSums(reference * reference)
  query_norm <- rowSums(query * query)
  d2 <- outer(query_norm, ref_norm, "+") - 2 * (query %*% t(reference))

  max.col(-d2)
}

#' Convert Functional Map to Pointwise Map
#'
#' @param fit `fm_fit` object.
#' @param reverse If `TRUE`, return source-to-target mapping.
#' @param use_adj If `TRUE`, use adjoint-style embeddings.
#'
#' @return Integer index vector representing nearest-neighbor pointwise map.
#' @export
as_p2p <- function(fit, reverse = FALSE, use_adj = FALSE) {
  if (!is_fm_fit(fit)) {
    stop("`fit` must inherit from `fm_fit`", call. = FALSE)
  }

  src <- fit$source
  tgt <- fit$target
  C <- fit$C

  if (isTRUE(reverse)) {
    src <- fit$target
    tgt <- fit$source
    C <- t(fit$C)
  }

  if (is.null(src$basis$vectors) || is.null(tgt$basis$vectors)) {
    stop(
      "Pointwise conversion requires basis vectors on both domains. Run fm_basis() first.",
      call. = FALSE
    )
  }

  k2 <- nrow(C)
  k1 <- ncol(C)

  if (k1 > ncol(src$basis$vectors) || k2 > ncol(tgt$basis$vectors)) {
    stop(
      "Functional map dimensions exceed available basis vectors on source/target domains.",
      call. = FALSE
    )
  }

  if (isTRUE(use_adj)) {
    emb_ref <- src$basis$vectors[, seq_len(k1), drop = FALSE]
    emb_query <- tgt$basis$vectors[, seq_len(k2), drop = FALSE] %*% C
  } else {
    emb_ref <- src$basis$vectors[, seq_len(k1), drop = FALSE] %*% t(C)
    emb_query <- tgt$basis$vectors[, seq_len(k2), drop = FALSE]
  }

  nearest_neighbor_index(emb_ref, emb_query)
}

#' Transfer Functions Through a Functional Map
#'
#' @param fit `fm_fit` object.
#' @param x Function values matrix or vector.
#' @param reverse If `TRUE`, transfer target to source.
#'
#' @return Transferred function values on destination domain.
#' @export
fm_transfer <- function(fit, x, reverse = FALSE) {
  if (!is_fm_fit(fit)) {
    stop("`fit` must inherit from `fm_fit`", call. = FALSE)
  }

  src <- if (isTRUE(reverse)) fit$target else fit$source
  tgt <- if (isTRUE(reverse)) fit$source else fit$target
  C <- if (isTRUE(reverse)) t(fit$C) else fit$C

  x_mat <- as.matrix(x)
  if (nrow(x_mat) != src$n_samples) {
    stop("Input function row count must match source domain sample count", call. = FALSE)
  }

  coef_src <- fm_project(src, x_mat)
  coef_tgt <- C %*% coef_src
  y <- fm_unproject(tgt, coef_tgt)

  if (is.vector(x) || is.null(dim(x))) {
    return(as.vector(y))
  }

  y
}

#' Compose Two Functional Fits
#'
#' @param fit_ab `fm_fit` from domain A to B.
#' @param fit_bc `fm_fit` from domain B to C.
#'
#' @return `fm_fit` object mapping A to C.
#' @export
fm_compose <- function(fit_ab, fit_bc) {
  if (!is_fm_fit(fit_ab) || !is_fm_fit(fit_bc)) {
    stop("Both inputs must inherit from `fm_fit`", call. = FALSE)
  }

  if (ncol(fit_bc$C) != nrow(fit_ab$C)) {
    stop("Map dimensions are incompatible for composition", call. = FALSE)
  }

  structure(
    list(
      C = fit_bc$C %*% fit_ab$C,
      k1 = ncol(fit_ab$C),
      k2 = nrow(fit_bc$C),
      source = fit_ab$source,
      target = fit_bc$target,
      penalties = list(),
      diagnostics = list(
        convergence = 0,
        message = "composed",
        counts = c(fn = NA_integer_, gradient = NA_integer_),
        total_objective = NA_real_,
        objective_terms = list(),
        warning = NULL
      )
    ),
    class = "fm_fit"
  )
}

matrix_pseudoinverse <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  keep <- s$d > (tol * max(s$d))

  if (!any(keep)) {
    return(matrix(0, nrow = ncol(A), ncol = nrow(A)))
  }

  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

#' Invert a Functional Fit
#'
#' @param fit `fm_fit` object.
#' @param method One of `"pseudoinverse"` or `"transpose"`.
#'
#' @return `fm_fit` object with swapped source/target direction.
#' @export
fm_inverse <- function(fit, method = c("pseudoinverse", "transpose")) {
  if (!is_fm_fit(fit)) {
    stop("`fit` must inherit from `fm_fit`", call. = FALSE)
  }

  method <- match.arg(method)

  C_inv <- if (method == "transpose") t(fit$C) else matrix_pseudoinverse(fit$C)

  structure(
    list(
      C = C_inv,
      k1 = ncol(C_inv),
      k2 = nrow(C_inv),
      source = fit$target,
      target = fit$source,
      penalties = list(),
      diagnostics = list(
        convergence = 0,
        message = sprintf("inverted:%s", method),
        counts = c(fn = NA_integer_, gradient = NA_integer_),
        total_objective = NA_real_,
        objective_terms = list(),
        warning = NULL
      )
    ),
    class = "fm_fit"
  )
}
