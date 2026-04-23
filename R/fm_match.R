compute_descriptor_operators <- function(domain, descriptors) {
  phi <- domain$basis$vectors
  if (is.null(phi)) {
    stop("Basis vectors are required to build descriptor operators", call. = FALSE)
  }

  if (is.numeric(domain$measure) && is.null(dim(domain$measure))) {
    pinv <- if (all(domain$measure == 1)) {
      t(phi)
    } else {
      t(phi) * domain$measure
    }
  } else {
    M <- fm_measure_matrix(domain)
    pinv <- t(phi) %*% M
  }

  lapply(seq_len(ncol(descriptors)), function(i) {
    as.matrix(pinv %*% (descriptors[, i] * phi))
  })
}

compute_descriptor_operators_cube <- function(domain, descriptors) {
  phi <- domain$basis$vectors
  if (is.null(phi)) {
    stop("Basis vectors are required to build descriptor operators", call. = FALSE)
  }

  if (exists("fm_compute_descriptor_operators_cpp", mode = "function")) {
    return(fm_compute_descriptor_operators_cpp(
      phi = phi,
      descriptors = as.matrix(descriptors),
      measure = domain$measure
    ))
  }

  M <- fm_measure_matrix(domain)
  pinv <- t(phi) %*% M
  k <- ncol(phi)
  p <- ncol(descriptors)

  out <- array(0, dim = c(k, k, p))
  for (i in seq_len(p)) {
    out[, , i] <- as.matrix(pinv %*% (descriptors[, i] * phi))
  }

  out
}

normalize_descriptor_batch_size <- function(descriptor_batch_size, n_descriptors) {
  if (is.null(descriptor_batch_size)) {
    return(NULL)
  }
  if (!is.numeric(descriptor_batch_size) || length(descriptor_batch_size) != 1L || descriptor_batch_size < 1) {
    stop("`descriptor_batch_size` must be NULL or a positive scalar", call. = FALSE)
  }
  as.integer(min(descriptor_batch_size, n_descriptors))
}

descriptor_batch_indices <- function(n_descriptors, batch_size) {
  if (is.null(batch_size) || batch_size >= n_descriptors) {
    return(list(seq_len(n_descriptors)))
  }
  split(seq_len(n_descriptors), ceiling(seq_len(n_descriptors) / batch_size))
}

make_descriptor_stream <- function(domain, descriptors) {
  phi <- domain$basis$vectors
  if (is.null(phi)) {
    stop("Basis vectors are required to build descriptor operators", call. = FALSE)
  }
  if (is.numeric(domain$measure) && is.null(dim(domain$measure))) {
    pinv <- if (all(domain$measure == 1)) {
      t(phi)
    } else {
      t(phi) * domain$measure
    }
  } else {
    M <- fm_measure_matrix(domain)
    pinv <- t(phi) %*% M
  }
  list(
    phi = phi,
    pinv = pinv,
    descriptors = as.matrix(descriptors)
  )
}

compute_descriptor_operators_stream <- function(stream, indices) {
  phi <- stream$phi
  pinv <- stream$pinv
  desc <- stream$descriptors

  lapply(indices, function(i) {
    as.matrix(pinv %*% (desc[, i] * phi))
  })
}

pack_descriptor_operator_pairs <- function(src_ops, tgt_ops) {
  Map(function(op1, op2) {
    list(
      op1 = op1,
      op2 = op2,
      op1_t = t(op1),
      op2_t = t(op2)
    )
  }, src_ops, tgt_ops)
}

comm_apply <- function(comm_ctx, callback) {
  mode <- comm_ctx$mode %||% "none"

  if (mode == "none") {
    return(invisible(NULL))
  }

  if (mode == "precomputed") {
    for (pair in comm_ctx$list_ops) {
      callback(pair$op1, pair$op2, pair$op1_t, pair$op2_t)
    }
    return(invisible(NULL))
  }

  if (mode == "stream") {
    for (idx in comm_ctx$batches) {
      src_ops <- compute_descriptor_operators_stream(comm_ctx$source_stream, idx)
      tgt_ops <- compute_descriptor_operators_stream(comm_ctx$target_stream, idx)
      for (bi in seq_along(src_ops)) {
        op1 <- src_ops[[bi]]
        op2 <- tgt_ops[[bi]]
        callback(op1, op2, t(op1), t(op2))
      }
    }
    return(invisible(NULL))
  }

  stop("Unknown commutativity context mode", call. = FALSE)
}

compute_comm_diag_terms <- function(comm_ctx, k2, k1) {
  comm_row_diag <- numeric(k2)
  comm_col_diag <- numeric(k1)

  comm_apply(comm_ctx, function(op1, op2, op1_t, op2_t) {
    comm_col_diag <<- comm_col_diag + colSums(op1 * op1)
    comm_row_diag <<- comm_row_diag + colSums(op2 * op2)
  })

  list(row = comm_row_diag, col = comm_col_diag)
}

energy_breakdown <- function(C, A, B, comm_ctx, ev_sqdiff, penalties) {
  e_descr <- 0
  e_lap <- 0
  e_comm <- 0

  if (penalties$descr > 0) {
    residual <- C %*% A - B
    e_descr <- 0.5 * sum(residual * residual)
  }

  if (penalties$lap > 0) {
    e_lap <- 0.5 * sum((C * C) * ev_sqdiff)
  }

  if (penalties$comm > 0 && (comm_ctx$mode %||% "none") != "none") {
    comm_apply(comm_ctx, function(op1, op2, op1_t, op2_t) {
      residual <- C %*% op1 - op2 %*% C
      e_comm <<- e_comm + (0.5 * sum(residual * residual))
    })
  }

  list(
    descr = penalties$descr * e_descr,
    lap = penalties$lap * e_lap,
    comm = penalties$comm * e_comm
  )
}

make_objective_cache <- function(k2, k1, A, B, comm_ctx, ev_sqdiff, penalties) {
  At <- t(A)

  cache <- new.env(parent = emptyenv())
  cache$x <- NULL
  cache$value <- NULL
  cache$grad <- NULL

  eval_fg <- function(x) {
    same_x <- !is.null(cache$x) &&
      length(cache$x) == length(x) &&
      !any(cache$x != x)
    if (same_x) {
      return(list(value = cache$value, grad = cache$grad))
    }

    C <- matrix(x, nrow = k2, ncol = k1)
    g <- matrix(0, nrow = k2, ncol = k1)
    value <- 0

    if (penalties$descr > 0) {
      residual_descr <- C %*% A - B
      value <- value + penalties$descr * (0.5 * sum(residual_descr * residual_descr))
      g <- g + penalties$descr * (residual_descr %*% At)
    }

    if (penalties$lap > 0) {
      value <- value + penalties$lap * (0.5 * sum((C * C) * ev_sqdiff))
      g <- g + penalties$lap * (C * ev_sqdiff)
    }

    if (penalties$comm > 0 && (comm_ctx$mode %||% "none") != "none") {
      comm_apply(comm_ctx, function(op1, op2, op1_t, op2_t) {
        residual <- op2 %*% C - C %*% op1
        value <<- value + penalties$comm * (0.5 * sum(residual * residual))
        g <<- g + penalties$comm * (op2_t %*% residual - residual %*% op1_t)
      })
    }

    cache$x <- x
    cache$value <- value
    cache$grad <- as.vector(g)

    list(value = cache$value, grad = cache$grad)
  }

  list(
    fn = function(x) eval_fg(x)$value,
    gr = function(x) eval_fg(x)$grad
  )
}

make_objective_cache_cpp <- function(
  A,
  B,
  src_ops_cube,
  tgt_ops_cube,
  src_ops_t_cube,
  tgt_ops_t_cube,
  ev_sqdiff,
  penalties
) {
  cache <- new.env(parent = emptyenv())
  cache$x <- NULL
  cache$value <- NULL
  cache$grad <- NULL

  eval_fg <- function(x) {
    same_x <- !is.null(cache$x) &&
      length(cache$x) == length(x) &&
      !any(cache$x != x)
    if (same_x) {
      return(list(value = cache$value, grad = cache$grad))
    }

    C <- matrix(x, nrow = nrow(B), ncol = nrow(A))
    out <- fm_match_value_grad_cpp(
      C = C,
      A = A,
      B = B,
      ev_sqdiff = ev_sqdiff,
      op1_cube = src_ops_cube,
      op2_cube = tgt_ops_cube,
      op1_t_cube = src_ops_t_cube,
      op2_t_cube = tgt_ops_t_cube,
      w_descr = penalties$descr,
      w_lap = penalties$lap,
      w_comm = penalties$comm
    )

    cache$x <- x
    cache$value <- as.numeric(out$value)
    cache$grad <- as.vector(out$grad)

    list(value = cache$value, grad = cache$grad)
  }

  list(
    fn = function(x) eval_fg(x)$value,
    gr = function(x) eval_fg(x)$grad
  )
}

resolve_kernel_backend <- function(kernel_backend) {
  backend <- match.arg(kernel_backend, c("auto", "r", "cpp"))
  has_cpp <- exists("fm_match_solve_cg_cpp", mode = "function") &&
    exists("fm_match_value_grad_cpp", mode = "function") &&
    exists("fm_match_energy_terms_cpp", mode = "function")

  if (backend == "auto") {
    return(if (has_cpp) "cpp" else "r")
  }

  if (backend == "cpp" && !has_cpp) {
    stop(
      "`kernel_backend = \"cpp\"` requested but native kernels are unavailable. Reinstall/build package with compiled code.",
      call. = FALSE
    )
  }

  backend
}

make_operator_apply <- function(A, comm_ctx, ev_sqdiff, penalties) {
  AAt <- A %*% t(A)

  function(C) {
    g <- matrix(0, nrow = nrow(C), ncol = ncol(C))

    if (penalties$descr > 0) {
      g <- g + penalties$descr * (C %*% AAt)
    }

    if (penalties$lap > 0) {
      g <- g + penalties$lap * (C * ev_sqdiff)
    }

    if (penalties$comm > 0 && (comm_ctx$mode %||% "none") != "none") {
      comm_apply(comm_ctx, function(op1, op2, op1_t, op2_t) {
        residual <- op2 %*% C - C %*% op1
        g <<- g + penalties$comm * (op2_t %*% residual - residual %*% op1_t)
      })
    }

    g
  }
}

solve_cg <- function(C0, rhs, apply_operator, maxit, tol, diag_precond = NULL) {
  C <- C0
  R <- rhs - apply_operator(C)
  use_prec <- !is.null(diag_precond)
  if (use_prec) {
    if (!all(dim(diag_precond) == dim(C0))) {
      stop("`diag_precond` dimensions must match `C0`", call. = FALSE)
    }
    Z <- R / diag_precond
  } else {
    Z <- R
  }
  P <- Z
  rz <- sum(R * Z)
  rr <- sum(R * R)
  residual0 <- sqrt(rr)

  if (!is.finite(residual0) || residual0 <= tol) {
    return(list(
      C = C,
      iterations = 0L,
      residual = residual0,
      converged = TRUE,
      breakdown = FALSE
    ))
  }

  converged <- FALSE
  breakdown <- FALSE
  iter_done <- 0L

  for (iter in seq_len(maxit)) {
    HP <- apply_operator(P)
    denom <- sum(P * HP)

    if (!is.finite(denom) || abs(denom) <= .Machine$double.eps) {
      breakdown <- TRUE
      iter_done <- iter - 1L
      break
    }

    alpha <- rz / denom
    C <- C + alpha * P
    R_new <- R - alpha * HP
    rr_new <- sum(R_new * R_new)
    iter_done <- iter

    if (!is.finite(rr_new)) {
      breakdown <- TRUE
      break
    }

    if (sqrt(rr_new) <= tol) {
      R <- R_new
      if (use_prec) {
        Z <- R / diag_precond
        rz <- sum(R * Z)
      } else {
        rz <- rr_new
      }
      converged <- TRUE
      break
    }

    if (use_prec) {
      Z_new <- R_new / diag_precond
      rz_new <- sum(R_new * Z_new)
      beta <- rz_new / rz
      P <- Z_new + beta * P
      Z <- Z_new
      rz <- rz_new
    } else {
      beta <- rr_new / rr
      P <- R_new + beta * P
      rz <- rr_new
    }
    R <- R_new
    rr <- rr_new
  }

  list(
    C = C,
    iterations = as.integer(iter_done),
    residual = sqrt(rr),
    converged = converged,
    breakdown = breakdown
  )
}

as_penalties <- function(penalties) {
  defaults <- list(descr = 1e-1, lap = 1e-3, comm = 1)
  out <- modifyList(defaults, penalties)

  if (any(vapply(out, function(v) !is.numeric(v) || length(v) != 1L || v < 0, logical(1)))) {
    stop("All penalties must be non-negative scalars", call. = FALSE)
  }

  out
}

init_map <- function(k2, k1, init = c("zeros", "identity", "random")) {
  init <- match.arg(init)

  if (init == "identity") {
    return(diag(1, nrow = k2, ncol = k1))
  }

  if (init == "random") {
    return(matrix(stats::runif(k2 * k1, min = -1, max = 1), nrow = k2, ncol = k1))
  }

  matrix(0, nrow = k2, ncol = k1)
}

fixed_first_column <- function(source, target, k2) {
  col <- numeric(k2)
  if (k2 < 1L) {
    return(col)
  }

  source_weights <- extract_measure_weights(source)
  target_weights <- extract_measure_weights(target)
  source_mass <- if (is.null(source_weights)) source$n_samples else sum(source_weights)
  target_mass <- if (is.null(target_weights)) target$n_samples else sum(target_weights)
  area_ratio <- sqrt(target_mass / max(source_mass, 1e-12))
  col[1] <- area_ratio
  col
}

#' Fit a Pairwise Functional Map
#'
#' @param source Source `fm_domain` with basis and descriptors.
#' @param target Target `fm_domain` with basis and descriptors.
#' @param descriptors List with `source` and `target` descriptor matrices.
#' @param penalties Named list with weights for `descr`, `lap`, and `comm`.
#' @param init Initialization strategy (`"zeros"`, `"identity"`, `"random"`).
#' @param maxit Maximum optimization iterations.
#' @param optimizer Optimization backend (`"cg"` or `"lbfgsb"`).
#' @param cg_tol Residual tolerance used by conjugate-gradient solver.
#' @param cg_maxit Optional maximum conjugate-gradient iterations.
#' @param kernel_backend Kernel backend for objective/gradient ops:
#'   `"auto"` (prefer compiled), `"r"`, or `"cpp"`.
#' @param descriptor_batch_size Optional descriptor chunk size for commutativity
#'   operator streaming. Use smaller values to reduce peak memory when many
#'   descriptors are used.
#' @param compute_objective Whether to compute the final objective breakdown
#'   after fitting. Disable to skip post-hoc diagnostics when only the map is
#'   needed.
#' @param trace Whether to print optimizer traces.
#'
#' @return An `fm_fit` object with map, diagnostics, and objective breakdown.
#' @export
fm_match <- function(
  source,
  target,
  descriptors,
  penalties = list(descr = 1e-1, lap = 1e-3, comm = 1),
  init = c("zeros", "identity", "random"),
  maxit = 200,
  optimizer = c("cg", "lbfgsb"),
  cg_tol = 1e-6,
  cg_maxit = NULL,
  kernel_backend = c("auto", "r", "cpp"),
  descriptor_batch_size = NULL,
  compute_objective = TRUE,
  trace = FALSE
) {
  if (!inherits(source, "fm_domain") || !inherits(target, "fm_domain")) {
    stop("`source` and `target` must inherit from `fm_domain`", call. = FALSE)
  }

  if (is.null(source$basis$vectors) || is.null(target$basis$vectors)) {
    stop("Both domains must have basis vectors. Run fm_basis() first.", call. = FALSE)
  }

  if (!is.list(descriptors) || is.null(descriptors$source) || is.null(descriptors$target)) {
    stop("`descriptors` must be a list with `source` and `target` matrices", call. = FALSE)
  }

  src_desc <- as.matrix(descriptors$source)
  tgt_desc <- as.matrix(descriptors$target)

  if (nrow(src_desc) != source$n_samples || nrow(tgt_desc) != target$n_samples) {
    stop("Descriptor row counts must match domain sample counts", call. = FALSE)
  }

  if (ncol(src_desc) != ncol(tgt_desc)) {
    stop("Source and target descriptors must have the same column count", call. = FALSE)
  }

  penalties <- as_penalties(penalties)
  optimizer <- match.arg(optimizer)
  kernel_backend <- resolve_kernel_backend(kernel_backend)
  descriptor_batch_size <- normalize_descriptor_batch_size(descriptor_batch_size, ncol(src_desc))
  compute_objective <- isTRUE(compute_objective)

  k1 <- source$basis$k
  k2 <- target$basis$k

  A <- as.matrix(fm_project(source, src_desc))
  B <- as.matrix(fm_project(target, tgt_desc))

  ev1 <- source$basis$values
  ev2 <- target$basis$values

  if (is.null(ev1) || is.null(ev2) || length(ev1) < k1 || length(ev2) < k2) {
    stop("Both domains must include basis eigenvalues for Laplacian commutativity", call. = FALSE)
  }

  ev_sqdiff <- (matrix(ev2[seq_len(k2)], nrow = k2, ncol = k1) -
    matrix(ev1[seq_len(k1)], nrow = k2, ncol = k1, byrow = TRUE))^2

  if (sum(ev_sqdiff) > 0) {
    ev_sqdiff <- ev_sqdiff / sum(ev_sqdiff)
  }

  comm_ctx <- list(mode = "none", n_batches = 0L)
  src_ops_cpp <- array(numeric(0), dim = c(0L, 0L, 0L))
  tgt_ops_cpp <- array(numeric(0), dim = c(0L, 0L, 0L))
  src_ops_t_cpp <- array(numeric(0), dim = c(0L, 0L, 0L))
  tgt_ops_t_cpp <- array(numeric(0), dim = c(0L, 0L, 0L))
  comm_row_diag <- numeric(k2)
  comm_col_diag <- numeric(k1)
  backend_note <- NULL

  if (penalties$comm > 0) {
    use_stream <- !is.null(descriptor_batch_size) && descriptor_batch_size < ncol(src_desc)
    if (use_stream && kernel_backend == "cpp") {
      kernel_backend <- "r"
      backend_note <- "`descriptor_batch_size` enabled streamed commutativity mode; switched kernel backend from cpp to r for lower memory usage."
    }

    if (kernel_backend == "cpp") {
      src_ops_cpp <- compute_descriptor_operators_cube(source, src_desc)
      tgt_ops_cpp <- compute_descriptor_operators_cube(target, tgt_desc)
      src_ops_t_cpp <- aperm(src_ops_cpp, c(2, 1, 3))
      tgt_ops_t_cpp <- aperm(tgt_ops_cpp, c(2, 1, 3))
      for (si in seq_len(dim(src_ops_cpp)[3])) {
        op_src <- src_ops_cpp[, , si, drop = TRUE]
        op_tgt <- tgt_ops_cpp[, , si, drop = TRUE]
        comm_col_diag <- comm_col_diag + colSums(op_src * op_src)
        comm_row_diag <- comm_row_diag + colSums(op_tgt * op_tgt)
      }
    } else {
      if (use_stream) {
        batches <- descriptor_batch_indices(ncol(src_desc), descriptor_batch_size)
        comm_ctx <- list(
          mode = "stream",
          source_stream = make_descriptor_stream(source, src_desc),
          target_stream = make_descriptor_stream(target, tgt_desc),
          batches = batches,
          n_batches = length(batches)
        )
      } else {
        src_ops <- compute_descriptor_operators(source, src_desc)
        tgt_ops <- compute_descriptor_operators(target, tgt_desc)
        comm_ctx <- list(
          mode = "precomputed",
          list_ops = pack_descriptor_operator_pairs(src_ops, tgt_ops),
          n_batches = 1L
        )
      }

      comm_diag <- compute_comm_diag_terms(comm_ctx, k2 = k2, k1 = k1)
      comm_row_diag <- comm_diag$row
      comm_col_diag <- comm_diag$col
    }
  }

  C0 <- init_map(k2 = k2, k1 = k1, init = init)
  if (k1 >= 1L) {
    C0[, 1] <- fixed_first_column(source, target, k2 = k2)
  }
  if (optimizer == "cg") {
    if (is.null(cg_maxit)) {
      # Keep a practical default budget: enough for moderate convergence
      # without approaching full L-BFGS-B runtime.
      cg_maxit <- max(5L, min(100L, ceiling(maxit / 5)))
    }
    cg_maxit <- as.integer(cg_maxit)
    if (!is.numeric(cg_tol) || length(cg_tol) != 1L || cg_tol <= 0) {
      stop("`cg_tol` must be a positive scalar", call. = FALSE)
    }
    if (!is.numeric(cg_maxit) || length(cg_maxit) != 1L || cg_maxit < 1) {
      stop("`cg_maxit` must be a positive scalar", call. = FALSE)
    }

    rhs <- penalties$descr * (B %*% t(A))
    AAt_diag <- rowSums(A * A)
    diag_precond <- penalties$descr * matrix(AAt_diag, nrow = k2, ncol = k1, byrow = TRUE)
    if (penalties$lap > 0) {
      diag_precond <- diag_precond + penalties$lap * ev_sqdiff
    }
    if (penalties$comm > 0) {
      diag_precond <- diag_precond +
        penalties$comm * (matrix(comm_row_diag, nrow = k2, ncol = k1) +
          matrix(comm_col_diag, nrow = k2, ncol = k1, byrow = TRUE))
    }
    diag_precond <- pmax(diag_precond, 1e-8)

    fixed_col <- if (k1 >= 1L) C0[, 1] else numeric(0)
    fixed_C <- matrix(0, nrow = k2, ncol = k1)
    C0_free <- C0
    AAt <- NULL
    if (kernel_backend == "cpp") {
      AAt <- A %*% t(A)
    }
    if (k1 >= 1L) {
      fixed_C[, 1] <- fixed_col
      C0_free[, 1] <- 0
      rhs <- rhs - {
        if (kernel_backend == "cpp") {
          if (length(fixed_col) == 1L || all(abs(fixed_col[-1]) <= sqrt(.Machine$double.eps))) {
            fm_match_apply_fixed_first_column_cpp(
              fixed_value = fixed_col[[1]],
              AAt = AAt,
              ev_sqdiff = ev_sqdiff,
              op1_cube = src_ops_cpp,
              op2_cube = tgt_ops_cpp,
              w_descr = penalties$descr,
              w_lap = penalties$lap,
              w_comm = penalties$comm
            )
          } else {
            fm_match_apply_operator_cpp(
              C = fixed_C,
              AAt = AAt,
              ev_sqdiff = ev_sqdiff,
              op1_cube = src_ops_cpp,
              op2_cube = tgt_ops_cpp,
              op1_t_cube = src_ops_t_cpp,
              op2_t_cube = tgt_ops_t_cpp,
              w_descr = penalties$descr,
              w_lap = penalties$lap,
              w_comm = penalties$comm
            )
          }
        } else {
          apply_operator_fixed <- make_operator_apply(
            A = A,
            comm_ctx = comm_ctx,
            ev_sqdiff = ev_sqdiff,
            penalties = penalties
          )
          apply_operator_fixed(fixed_C)
        }
      }
      rhs[, 1] <- 0
    }

    if (kernel_backend == "cpp") {
      sol <- fm_match_solve_cg_cpp(
        C0 = C0_free,
        rhs = rhs,
        AAt = AAt,
        ev_sqdiff = ev_sqdiff,
        op1_cube = src_ops_cpp,
        op2_cube = tgt_ops_cpp,
        op1_t_cube = src_ops_t_cpp,
        op2_t_cube = tgt_ops_t_cpp,
        w_descr = penalties$descr,
        w_lap = penalties$lap,
        w_comm = penalties$comm,
        diag_precond = diag_precond,
        maxit = cg_maxit,
        tol = cg_tol
      )
    } else {
      apply_operator <- make_operator_apply(
        A = A,
        comm_ctx = comm_ctx,
        ev_sqdiff = ev_sqdiff,
        penalties = penalties
      )
      if (k1 >= 1L) {
        apply_operator_locked <- function(C) {
          out <- apply_operator(C)
          out[, 1] <- 0
          out
        }
        sol <- solve_cg(
          C0_free,
          rhs,
          apply_operator = apply_operator_locked,
          maxit = cg_maxit,
          tol = cg_tol,
          diag_precond = diag_precond
        )
        sol$C[, 1] <- fixed_col
      } else {
        sol <- solve_cg(
          C0,
          rhs,
          apply_operator = apply_operator,
          maxit = cg_maxit,
          tol = cg_tol,
          diag_precond = diag_precond
        )
      }
    }

    C <- sol$C
    if (k1 >= 1L) {
      C[, 1] <- fixed_col
    }

    if (isTRUE(sol$converged)) {
      conv_code <- 0L
      conv_msg <- sprintf("CG converged in %d iterations", sol$iterations)
    } else if (isTRUE(sol$breakdown)) {
      conv_code <- 2L
      conv_msg <- "CG terminated due to numerical breakdown"
    } else {
      conv_code <- 1L
      conv_msg <- sprintf("CG reached iteration limit (%d)", cg_maxit)
    }

    counts <- stats::setNames(
      c(sol$iterations, sol$iterations),
      c("function", "gradient")
    )
  } else {
    if (kernel_backend == "cpp") {
      objective <- make_objective_cache_cpp(
        A = A,
        B = B,
        src_ops_cube = src_ops_cpp,
        tgt_ops_cube = tgt_ops_cpp,
        src_ops_t_cube = src_ops_t_cpp,
        tgt_ops_t_cube = tgt_ops_t_cpp,
        ev_sqdiff = ev_sqdiff,
        penalties = penalties
      )
    } else {
      objective <- make_objective_cache(
        k2 = k2,
        k1 = k1,
        A = A,
        B = B,
        comm_ctx = comm_ctx,
        ev_sqdiff = ev_sqdiff,
        penalties = penalties
      )
    }

    opt <- stats::optim(
      par = as.vector(C0),
      fn = objective$fn,
      gr = objective$gr,
      method = "L-BFGS-B",
      lower = {
        lo <- rep(-Inf, length(C0))
        if (k1 >= 1L) {
          lo[seq_len(k2)] <- C0[, 1]
        }
        lo
      },
      upper = {
        hi <- rep(Inf, length(C0))
        if (k1 >= 1L) {
          hi[seq_len(k2)] <- C0[, 1]
        }
        hi
      },
      control = list(maxit = maxit, trace = if (isTRUE(trace)) 1 else 0)
    )

    C <- matrix(opt$par, nrow = k2, ncol = k1)
    conv_code <- opt$convergence
    conv_msg <- opt$message
    counts <- opt$counts
  }

  if (compute_objective) {
    if (kernel_backend == "cpp") {
      term_vec <- fm_match_energy_terms_cpp(
        C = C,
        A = A,
        B = B,
        ev_sqdiff = ev_sqdiff,
        op1_cube = src_ops_cpp,
        op2_cube = tgt_ops_cpp,
        w_descr = penalties$descr,
        w_lap = penalties$lap,
        w_comm = penalties$comm
      )
      terms <- list(
        descr = as.numeric(term_vec$descr),
        lap = as.numeric(term_vec$lap),
        comm = as.numeric(term_vec$comm)
      )
    } else {
      terms <- energy_breakdown(C, A, B, comm_ctx, ev_sqdiff, penalties)
    }
  } else {
    terms <- list(descr = NA_real_, lap = NA_real_, comm = NA_real_)
  }

  diagnostics <- list(
    convergence = conv_code,
    message = conv_msg,
    counts = counts,
    total_objective = if (compute_objective) sum(unlist(terms)) else NA_real_,
    objective_terms = terms,
    objective_computed = compute_objective,
    optimizer = optimizer,
    kernel_backend = kernel_backend,
    descriptor_batch_size = descriptor_batch_size,
    commutativity_mode = comm_ctx$mode %||% "none",
    commutativity_batches = as.integer(comm_ctx$n_batches %||% 0L),
    backend_note = backend_note,
    warning = if (conv_code == 0) NULL else {
      sprintf("Optimizer returned convergence code %d", conv_code)
    }
  )

  structure(
    list(
      C = C,
      k1 = k1,
      k2 = k2,
      source = source,
      target = target,
      penalties = penalties,
      diagnostics = diagnostics
    ),
    class = "fm_fit"
  )
}

#' Extract Functional Map Matrix
#'
#' @param fit `fm_fit` object.
#'
#' @return Functional map matrix.
#' @export
fm_fit_matrix <- function(fit) {
  if (!inherits(fit, "fm_fit")) {
    stop("`fit` must inherit from `fm_fit`", call. = FALSE)
  }
  fit$C
}

#' @export
print.fm_fit <- function(x, ...) {
  cat(sprintf("<fm_fit> C: %d x %d\n", nrow(x$C), ncol(x$C)))
  cat(sprintf("  convergence: %d\n", x$diagnostics$convergence))
  cat(sprintf("  total objective: %.6f\n", x$diagnostics$total_objective))
  invisible(x)
}
