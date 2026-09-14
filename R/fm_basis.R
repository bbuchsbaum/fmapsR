.fm_basis_cache <- new.env(parent = emptyenv())

operator_signature <- function(operator) {
  raw <- serialize(operator, NULL, version = 2)
  bytes <- as.integer(raw)

  paste(
    length(bytes),
    sum(bytes),
    sum((seq_along(bytes) * bytes) %% 1000000007L),
    sep = "-"
  )
}

basis_cache_key <- function(operator, k, solver, which, control) {
  paste(
    operator_signature(operator),
    as.integer(k),
    solver,
    which,
    paste(names(control), unlist(control), sep = "=", collapse = ";"),
    sep = "|"
  )
}

compute_basis_base <- function(operator, k, which = "SM", sigma = NULL) {
  mat <- as.matrix(operator)
  symmetric <- isTRUE(all.equal(mat, t(mat), tolerance = 1e-10))
  eig <- eigen(mat, symmetric = symmetric)

  magnitude <- if (is.null(sigma)) Mod(eig$values) else 1 / Mod(eig$values - sigma)
  ord <- order(magnitude, decreasing = which == "LM")

  keep <- ord[seq_len(k)]
  list(values = eig$values[keep], vectors = eig$vectors[, keep, drop = FALSE])
}

compute_basis_rspectra <- function(operator, k, which = "SM", control = list()) {
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    return(NULL)
  }

  args <- c(list(A = operator, k = k, which = which), control)

  out <- tryCatch(
    do.call(RSpectra::eigs_sym, args),
    error = function(e) NULL
  )

  if (is.null(out) || !is.matrix(out$vectors) || length(out$values) != k || ncol(out$vectors) != k ||
      any(!is.finite(out$values)) || any(!is.finite(out$vectors))) {
    return(NULL)
  }

  magnitude <- if (is.null(control$sigma)) Mod(out$values) else 1 / Mod(out$values - control$sigma)
  ord <- order(magnitude, decreasing = which == "LM")
  list(values = out$values[ord], vectors = out$vectors[, ord, drop = FALSE])
}

operator_rspectra_hints <- function(operator, k) {
  preferred_which <- attr(operator, "fm_rspectra_preferred_which", exact = TRUE)
  sigma <- attr(operator, "fm_rspectra_sigma", exact = TRUE)
  ncv_min <- attr(operator, "fm_rspectra_ncv_min", exact = TRUE)
  ncv_mult <- attr(operator, "fm_rspectra_ncv_mult", exact = TRUE)

  control <- list()
  if (!is.null(sigma)) {
    control$sigma <- sigma
  }

  if (!is.null(ncv_min) || !is.null(ncv_mult)) {
    min_ncv <- if (is.null(ncv_min)) 0L else as.integer(ncv_min)
    mult_ncv <- if (is.null(ncv_mult)) 0L else as.integer(ncv_mult)
    control$opts <- list(ncv = min(nrow(operator), max(min_ncv, mult_ncv * as.integer(k))))
  }

  list(
    which = preferred_which,
    control = control
  )
}

resolve_basis_solver_controls <- function(operator, k, solver, which, which_missing, control) {
  effective_which <- which
  effective_control <- control

  if (!identical(solver, "rspectra") || !isTRUE(which_missing)) {
    return(list(which = effective_which, control = effective_control))
  }

  hints <- operator_rspectra_hints(operator, k = k)
  if (isTRUE(which_missing) && is.character(hints$which) && length(hints$which) == 1L) {
    effective_which <- hints$which
  }

  if (length(hints$control) > 0L) {
    effective_control <- utils::modifyList(hints$control, effective_control)
  }

  list(which = effective_which, control = effective_control)
}

transform_basis_vectors <- function(operator, basis_vectors) {
  transform_kind <- attr(operator, "fm_basis_transform", exact = TRUE)
  if (is.null(transform_kind)) {
    return(basis_vectors)
  }

  if (identical(transform_kind, "mass_mhalf")) {
    mass <- attr(operator, "fm_basis_mass", exact = TRUE)
    if (is.null(mass) || length(mass) != nrow(basis_vectors)) {
      stop("Operator basis-transform metadata is malformed", call. = FALSE)
    }
    return(sweep(basis_vectors, 1, sqrt(mass), "/"))
  }

  stop("Unsupported operator basis transform", call. = FALSE)
}

#' Compute or Retrieve a Domain Basis
#'
#' @param domain `fm_domain` object with an `operator`.
#' @param k Number of basis vectors.
#' @param solver One of `"rspectra"` or `"base"`.
#' @param which Eigen target used by solver (`"SM"` or `"LM"`).
#' @param cache Whether to use process-level in-memory cache.
#' @param seed Optional RNG seed for reproducibility.
#' @param control Optional arguments passed to `RSpectra::eigs_sym`, for example
#'   `list(sigma = 0, opts = list(ncv = 40))`. When `which` is omitted,
#'   mesh operators can supply shift-invert defaults. An explicit `which`
#'   disables these defaults. Dense fallback preserves shift-invert selection.
#'
#' @return Updated `fm_domain` object with populated `basis`.
#' @export
fm_basis <- function(
  domain,
  k,
  solver = c("rspectra", "base"),
  which = c("SM", "LM"),
  cache = TRUE,
  seed = NULL,
  control = list()
) {
  if (!inherits(domain, "fm_domain")) {
    stop("`domain` must inherit from `fm_domain`", call. = FALSE)
  }

  if (is.null(domain$operator)) {
    stop("`domain$operator` is required before basis computation", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0 || k != floor(k)) {
    stop("`k` must be a positive integer scalar", call. = FALSE)
  }

  k <- as.integer(k)
  if (k > domain$n_samples) {
    stop("`k` cannot exceed `domain$n_samples`", call. = FALSE)
  }

  which_missing <- missing(which)
  solver <- match.arg(solver)
  which <- match.arg(which)
  solver_controls <- resolve_basis_solver_controls(
    operator = domain$operator,
    k = k,
    solver = solver,
    which = which,
    which_missing = which_missing,
    control = control
  )
  which <- match.arg(solver_controls$which, choices = c("SM", "LM"))
  control <- solver_controls$control

  if (!is.null(seed)) {
    set.seed(seed)
  }

  key <- basis_cache_key(domain$operator, k, solver, which, control)
  if (isTRUE(cache) && exists(key, envir = .fm_basis_cache, inherits = FALSE)) {
    basis <- get(key, envir = .fm_basis_cache, inherits = FALSE)
    basis$cached <- TRUE
    domain$basis <- basis
    return(domain)
  }

  use_solver <- solver
  basis_raw <- NULL

  if (solver == "rspectra" && k < domain$n_samples) {
    basis_raw <- compute_basis_rspectra(
      operator = domain$operator,
      k = k,
      which = which,
      control = control
    )

    if (is.null(basis_raw)) {
      use_solver <- "base"
    }
  } else {
    use_solver <- "base"
  }

  if (is.null(basis_raw)) {
    basis_raw <- compute_basis_base(domain$operator, k = k, which = which,
      sigma = if (solver == "rspectra") control$sigma else NULL)
  }

  basis_vectors <- transform_basis_vectors(domain$operator, basis_raw$vectors)
  # Eigenvectors are defined only up to sign. Use a deterministic orientation,
  # including a positive constant mode for Laplacian matching constraints.
  pivots <- max.col(t(abs(basis_vectors)), ties.method = "first")
  signs <- sign(basis_vectors[cbind(pivots, seq_len(k))])
  basis_vectors <- sweep(basis_vectors, 2, signs, "*")

  basis <- list(
    values = basis_raw$values,
    vectors = basis_vectors,
    k = as.integer(k),
    solver = use_solver,
    which = which,
    control = control,
    cached = FALSE
  )

  if (isTRUE(cache)) {
    assign(key, basis, envir = .fm_basis_cache)
  }

  domain$basis <- basis
  domain
}

#' Clear Basis Cache
#'
#' @return `TRUE` invisibly.
#' @export
fm_basis_cache_clear <- function() {
  rm(list = ls(envir = .fm_basis_cache, all.names = TRUE), envir = .fm_basis_cache)
  invisible(TRUE)
}
