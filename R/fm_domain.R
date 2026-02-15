`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

normalize_basis <- function(basis, n_samples) {
  if (is.null(basis)) {
    return(list(values = NULL, vectors = NULL, k = 0L, solver = NULL))
  }

  if (!is.list(basis)) {
    stop("`basis` must be a list", call. = FALSE)
  }

  vectors <- basis$vectors %||% NULL
  values <- basis$values %||% NULL

  if (!is.null(vectors)) {
    if (!is.matrix(vectors)) {
      stop("`basis$vectors` must be a matrix", call. = FALSE)
    }
    if (nrow(vectors) != n_samples) {
      stop("`basis$vectors` row count must match `n_samples`", call. = FALSE)
    }
  }

  if (!is.null(values)) {
    if (!is.numeric(values)) {
      stop("`basis$values` must be numeric", call. = FALSE)
    }
    if (!is.null(vectors) && length(values) != ncol(vectors)) {
      stop("`basis$values` length must match `ncol(basis$vectors)`", call. = FALSE)
    }
  }

  k <- basis$k %||% if (is.null(vectors)) 0L else ncol(vectors)

  list(
    values = values,
    vectors = vectors,
    k = as.integer(k),
    solver = basis$solver %||% NULL
  )
}

normalize_measure <- function(measure, n_samples) {
  if (is.null(measure)) {
    return(rep(1, n_samples))
  }

  if (is.matrix(measure) || inherits(measure, "Matrix")) {
    dims <- dim(measure)
    if (length(dims) != 2 || any(dims != c(n_samples, n_samples))) {
      stop("Matrix `measure` must be `n_samples x n_samples`", call. = FALSE)
    }
    return(measure)
  }

  if (is.numeric(measure) && is.null(dim(measure))) {
    if (length(measure) != n_samples) {
      stop("Numeric `measure` must have length `n_samples`", call. = FALSE)
    }
    return(measure)
  }

  stop("`measure` must be numeric, matrix, sparse Matrix, or NULL", call. = FALSE)
}

validate_square_operator <- function(operator, n_samples, name = "operator") {
  if (is.null(operator)) {
    return(NULL)
  }

  if (!(is.matrix(operator) || inherits(operator, "Matrix"))) {
    stop(sprintf("`%s` must be a matrix or sparse Matrix", name), call. = FALSE)
  }

  dims <- dim(operator)
  if (length(dims) != 2 || any(dims != c(n_samples, n_samples))) {
    stop(
      sprintf("`%s` must be `n_samples x n_samples`", name),
      call. = FALSE
    )
  }

  operator
}

#' Create a New Functional-Map Domain
#'
#' @param type Domain type label.
#' @param n_samples Number of samples in the domain.
#' @param data Backend-specific payload.
#' @param measure Sample measure as weights or matrix.
#' @param operator Domain operator (typically Laplace-like).
#' @param basis Optional basis metadata list.
#' @param adjacency Optional adjacency matrix.
#' @param projector Optional custom projection function.
#' @param unprojector Optional custom reconstruction function.
#' @param metadata Optional named metadata list.
#'
#' @return An `fm_domain` object.
#' @export
fm_new_domain <- function(
  type,
  n_samples,
  data = NULL,
  measure = NULL,
  operator = NULL,
  basis = NULL,
  adjacency = NULL,
  projector = NULL,
  unprojector = NULL,
  metadata = list()
) {
  if (!is.character(type) || length(type) != 1L || !nzchar(type)) {
    stop("`type` must be a non-empty string", call. = FALSE)
  }

  if (!is.numeric(n_samples) || length(n_samples) != 1L || n_samples <= 0) {
    stop("`n_samples` must be a positive scalar", call. = FALSE)
  }

  n_samples <- as.integer(n_samples)
  basis <- normalize_basis(basis, n_samples)
  measure <- normalize_measure(measure, n_samples)
  operator <- validate_square_operator(operator, n_samples, name = "operator")
  adjacency <- validate_square_operator(adjacency, n_samples, name = "adjacency")

  if (!is.list(metadata)) {
    stop("`metadata` must be a list", call. = FALSE)
  }

  default_projector <- function(x, basis_obj, measure_obj) {
    if (is.null(basis_obj$vectors)) {
      stop("No basis vectors available for projection", call. = FALSE)
    }
    x_mat <- as.matrix(x)
    M <- if (is.numeric(measure_obj)) {
      Matrix::Diagonal(x = measure_obj)
    } else {
      measure_obj
    }
    t(basis_obj$vectors) %*% (M %*% x_mat)
  }

  default_unprojector <- function(coef, basis_obj) {
    if (is.null(basis_obj$vectors)) {
      stop("No basis vectors available for unprojection", call. = FALSE)
    }
    basis_obj$vectors %*% as.matrix(coef)
  }

  structure(
    list(
      type = type,
      n_samples = n_samples,
      data = data,
      measure = measure,
      operator = operator,
      basis = basis,
      adjacency = adjacency,
      projector = projector %||% default_projector,
      unprojector = unprojector %||% default_unprojector,
      metadata = metadata
    ),
    class = c(paste0("fm_domain_", type), "fm_domain")
  )
}

#' Unified Domain Constructor
#'
#' @param data Domain payload.
#' @param type One of `"generic"`, `"mesh"`, `"pointcloud"`, or `"graph"`.
#' @param n_samples Number of samples for generic domains.
#' @param ... Adapter-specific arguments.
#'
#' @return An `fm_domain` object.
#' @export
fm_domain <- function(data = NULL, type = c("generic", "mesh", "pointcloud", "graph"), n_samples = NULL, ...) {
  type <- match.arg(type)

  switch(
    type,
    mesh = fm_domain_mesh(vertices = data, ...),
    pointcloud = fm_domain_pointcloud(points = data, ...),
    graph = fm_domain_graph(adjacency = data, ...),
    generic = {
      inferred_n <- n_samples
      if (is.null(inferred_n) && !is.null(data) && !is.null(dim(data))) {
        inferred_n <- nrow(data)
      }
      fm_domain_generic(n_samples = inferred_n, data = data, ...)
    }
  )
}

#' Mesh Domain Adapter
#'
#' @param vertices Vertex matrix.
#' @param faces Optional face index matrix.
#' @param ... Passed to [fm_new_domain()].
#'
#' @return An `fm_domain_mesh` object.
#' @export
fm_domain_mesh <- function(vertices, faces = NULL, ...) {
  if (!is.matrix(vertices)) {
    stop("`vertices` must be a matrix", call. = FALSE)
  }

  if (!is.null(faces) && !is.matrix(faces)) {
    stop("`faces` must be a matrix when provided", call. = FALSE)
  }

  fm_new_domain(
    type = "mesh",
    n_samples = nrow(vertices),
    data = list(vertices = vertices, faces = faces),
    ...
  )
}

#' Point Cloud Domain Adapter
#'
#' @param points Point matrix.
#' @param ... Passed to [fm_new_domain()].
#'
#' @return An `fm_domain_pointcloud` object.
#' @export
fm_domain_pointcloud <- function(points, ...) {
  if (!is.matrix(points)) {
    stop("`points` must be a matrix", call. = FALSE)
  }

  fm_new_domain(
    type = "pointcloud",
    n_samples = nrow(points),
    data = list(points = points),
    ...
  )
}

#' Graph Domain Adapter
#'
#' @param adjacency Adjacency matrix.
#' @param ... Passed to [fm_new_domain()].
#'
#' @return An `fm_domain_graph` object.
#' @export
fm_domain_graph <- function(adjacency, ...) {
  validate_square_operator(adjacency, nrow(adjacency), name = "adjacency")

  fm_new_domain(
    type = "graph",
    n_samples = nrow(adjacency),
    data = list(adjacency = adjacency),
    adjacency = adjacency,
    ...
  )
}

#' Generic Domain Adapter
#'
#' @param n_samples Number of samples in the domain.
#' @param data Optional payload.
#' @param ... Passed to [fm_new_domain()].
#'
#' @return An `fm_domain_generic` object.
#' @export
fm_domain_generic <- function(n_samples, data = NULL, ...) {
  fm_new_domain(
    type = "generic",
    n_samples = n_samples,
    data = data,
    ...
  )
}

#' Number of Samples in a Domain
#'
#' @param domain `fm_domain` object.
#'
#' @return Integer sample count.
#' @export
fm_n_samples <- function(domain) {
  if (!inherits(domain, "fm_domain")) {
    stop("`domain` must inherit from `fm_domain`", call. = FALSE)
  }
  domain$n_samples
}

#' Build a Measure Matrix from Domain Measure
#'
#' @param domain `fm_domain` object.
#'
#' @return Dense or sparse measure matrix.
#' @export
fm_measure_matrix <- function(domain) {
  if (!inherits(domain, "fm_domain")) {
    stop("`domain` must inherit from `fm_domain`", call. = FALSE)
  }

  if (is.numeric(domain$measure) && is.null(dim(domain$measure))) {
    return(Matrix::Diagonal(x = domain$measure))
  }

  domain$measure
}

#' Project Functions on a Domain Basis
#'
#' @param domain `fm_domain` object.
#' @param x Function values on samples.
#' @param basis Optional basis list override.
#'
#' @return Basis coefficients.
#' @export
fm_project <- function(domain, x, basis = NULL) {
  if (!inherits(domain, "fm_domain")) {
    stop("`domain` must inherit from `fm_domain`", call. = FALSE)
  }

  basis_obj <- basis %||% domain$basis
  domain$projector(x, basis_obj, domain$measure)
}

#' Reconstruct Functions from Domain Basis Coefficients
#'
#' @param domain `fm_domain` object.
#' @param coef Basis coefficients.
#' @param basis Optional basis list override.
#'
#' @return Reconstructed function values.
#' @export
fm_unproject <- function(domain, coef, basis = NULL) {
  if (!inherits(domain, "fm_domain")) {
    stop("`domain` must inherit from `fm_domain`", call. = FALSE)
  }

  basis_obj <- basis %||% domain$basis
  domain$unprojector(coef, basis_obj)
}

#' Check Domain Class
#'
#' @param x Object to check.
#'
#' @return Logical scalar.
#' @export
is_fm_domain <- function(x) {
  inherits(x, "fm_domain")
}

#' @export
print.fm_domain <- function(x, ...) {
  cat(sprintf("<fm_domain:%s> n=%d\n", x$type, x$n_samples))
  cat(sprintf("  operator: %s\n", if (is.null(x$operator)) "none" else "set"))
  cat(sprintf("  basis vectors: %s\n", if (is.null(x$basis$vectors)) "none" else ncol(x$basis$vectors)))
  invisible(x)
}
