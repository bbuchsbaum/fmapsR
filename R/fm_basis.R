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

compute_basis_base <- function(operator, k, which = "SM") {
  mat <- as.matrix(operator)
  symmetric <- isTRUE(all.equal(mat, t(mat), tolerance = 1e-10))
  eig <- eigen(mat, symmetric = symmetric)

  if (which == "LM") {
    ord <- order(Mod(eig$values), decreasing = TRUE)
  } else {
    ord <- order(Mod(eig$values), decreasing = FALSE)
  }

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

  if (is.null(out)) {
    return(NULL)
  }

  list(values = out$values, vectors = out$vectors)
}

#' Compute or Retrieve a Domain Basis
#'
#' @param domain `fm_domain` object with an `operator`.
#' @param k Number of basis vectors.
#' @param solver One of `"rspectra"` or `"base"`.
#' @param which Eigen target used by solver (`"SM"` or `"LM"`).
#' @param cache Whether to use process-level in-memory cache.
#' @param seed Optional RNG seed for reproducibility.
#' @param control Optional solver controls passed to RSpectra.
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

  if (!is.numeric(k) || length(k) != 1L || k <= 0) {
    stop("`k` must be a positive scalar", call. = FALSE)
  }

  k <- as.integer(k)
  if (k > domain$n_samples) {
    stop("`k` cannot exceed `domain$n_samples`", call. = FALSE)
  }

  solver <- match.arg(solver)
  which <- match.arg(which)

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
    basis_raw <- compute_basis_base(domain$operator, k = k, which = which)
  }

  basis <- list(
    values = basis_raw$values,
    vectors = basis_raw$vectors,
    k = as.integer(k),
    solver = use_solver,
    which = which,
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
