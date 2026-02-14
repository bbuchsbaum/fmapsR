#' Build a Covariance-Derived Sample Operator
#'
#' @param x Numeric matrix with samples in rows and features in columns.
#' @param center Whether to center features.
#' @param scale Whether to scale features to unit variance.
#' @param ridge Diagonal ridge regularization added for numerical stability.
#'
#' @return Symmetric `n_samples x n_samples` operator matrix.
#' @export
fm_operator_covariance <- function(x, center = TRUE, scale = TRUE, ridge = 1e-3) {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`x` must be numeric", call. = FALSE)
  }
  if (nrow(x) < 2L || ncol(x) < 1L) {
    stop("`x` must contain at least two samples and one feature", call. = FALSE)
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || ridge < 0) {
    stop("`ridge` must be a non-negative scalar", call. = FALSE)
  }

  xs <- scale(x, center = center, scale = scale)
  xs[!is.finite(xs)] <- 0

  op <- tcrossprod(xs) / max(1L, ncol(xs))
  op <- (op + t(op)) / 2
  diag(op) <- diag(op) + ridge
  op
}

#' Build a kNN Laplacian Sample Operator
#'
#' @param x Numeric matrix with samples in rows and features in columns.
#' @param k Number of nearest neighbors per sample.
#' @param sigma Optional RBF scale. If `NULL`, estimated from neighbor distances.
#' @param normalized Whether to return normalized Laplacian (`TRUE`) or
#'   combinatorial Laplacian (`FALSE`).
#' @param eps Numerical floor for degree stabilization.
#'
#' @return Symmetric `n_samples x n_samples` Laplacian-like operator matrix.
#' @export
fm_operator_knn_laplacian <- function(x, k = 8, sigma = NULL, normalized = TRUE, eps = 1e-8) {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("`x` must be numeric", call. = FALSE)
  }
  n <- nrow(x)
  if (n < 3L || ncol(x) < 1L) {
    stop("`x` must contain at least three samples and one feature", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || k < 1 || k >= n) {
    stop("`k` must be an integer in [1, n_samples - 1]", call. = FALSE)
  }
  if (!is.null(sigma) && (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0)) {
    stop("`sigma` must be NULL or a positive scalar", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0) {
    stop("`eps` must be a positive scalar", call. = FALSE)
  }

  k <- as.integer(k)
  row_norm <- rowSums(x * x)
  d2 <- outer(row_norm, row_norm, "+") - 2 * (x %*% t(x))
  d2[d2 < 0] <- 0
  diag(d2) <- Inf

  nn_idx <- apply(d2, 1, function(row) order(row)[seq_len(k)])
  if (!is.matrix(nn_idx)) {
    nn_idx <- matrix(nn_idx, nrow = k, ncol = n)
  }

  if (is.null(sigma)) {
    nn_d2 <- vapply(seq_len(n), function(i) d2[i, nn_idx[, i]], numeric(k))
    sigma <- stats::median(sqrt(as.numeric(nn_d2)))
    if (!is.finite(sigma) || sigma <= 0) {
      sigma <- 1
    }
  }

  W <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    idx <- nn_idx[, i]
    wij <- exp(-d2[i, idx] / (2 * sigma^2))
    W[i, idx] <- wij
  }
  W <- pmax(W, t(W))
  diag(W) <- 0

  degree <- rowSums(W)
  D <- diag(degree, nrow = n, ncol = n)

  if (!isTRUE(normalized)) {
    return(D - W)
  }

  inv_sqrt_d <- 1 / sqrt(pmax(degree, eps))
  Dm12 <- diag(inv_sqrt_d, nrow = n, ncol = n)
  I <- diag(1, nrow = n, ncol = n)
  I - Dm12 %*% W %*% Dm12
}

#' Build an Operator for Generic Multidimensional Data
#'
#' @param x Numeric matrix with samples in rows and features in columns.
#' @param method Operator family: `"covariance"` or `"knn_laplacian"`.
#' @param ... Additional method-specific arguments.
#'
#' @return Symmetric sample operator matrix.
#' @export
fm_operator <- function(x, method = c("covariance", "knn_laplacian"), ...) {
  method <- match.arg(method)
  switch(
    method,
    covariance = fm_operator_covariance(x, ...),
    knn_laplacian = fm_operator_knn_laplacian(x, ...)
  )
}
