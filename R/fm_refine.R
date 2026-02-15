weighted_least_squares <- function(X, Y, weights = NULL) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if (nrow(Y) != nrow(X)) {
    stop("`Y` row count must match `X` row count", call. = FALSE)
  }

  solve_or_pinv <- function(A, B) {
    res <- tryCatch(qr.solve(A, B), error = function(e) NULL)
    if (!is.null(res)) {
      return(res)
    }
    matrix_pseudoinverse(A) %*% B
  }

  if (is.null(weights)) {
    return(solve_or_pinv(X, Y))
  }

  w <- as.numeric(weights)
  if (length(w) != nrow(X)) {
    stop("`weights` length must match row count of design matrix", call. = FALSE)
  }

  sqrt_w <- sqrt(pmax(w, 0))
  Xw <- X * sqrt_w
  Yw <- Y * sqrt_w

  solve_or_pinv(Xw, Yw)
}

prepare_weighted_ls_solver <- function(X, weights = NULL) {
  X <- as.matrix(X)

  solve_or_pinv <- function(A, B) {
    res <- tryCatch(qr.solve(A, B), error = function(e) NULL)
    if (!is.null(res)) {
      return(res)
    }
    matrix_pseudoinverse(A) %*% B
  }

  if (is.null(weights)) {
    qr_x <- qr(X)
    return(function(Y) {
      Y <- as.matrix(Y)
      if (nrow(Y) != nrow(X)) {
        stop("`Y` row count must match `X` row count", call. = FALSE)
      }
      res <- tryCatch(qr.coef(qr_x, Y), error = function(e) NULL)
      if (!is.null(res) && !anyNA(res)) {
        return(res)
      }
      solve_or_pinv(X, Y)
    })
  }

  w <- as.numeric(weights)
  if (length(w) != nrow(X)) {
    stop("`weights` length must match row count of design matrix", call. = FALSE)
  }

  sqrt_w <- sqrt(pmax(w, 0))
  Xw <- X * sqrt_w
  qr_xw <- qr(Xw)

  function(Y) {
    Y <- as.matrix(Y)
    if (nrow(Y) != nrow(X)) {
      stop("`Y` row count must match `X` row count", call. = FALSE)
    }
    Yw <- Y * sqrt_w
    res <- tryCatch(qr.coef(qr_xw, Yw), error = function(e) NULL)
    if (!is.null(res) && !anyNA(res)) {
      return(res)
    }
    solve_or_pinv(Xw, Yw)
  }
}

extract_measure_weights <- function(domain, indices = NULL) {
  m <- domain$measure

  if (is.matrix(m) || inherits(m, "Matrix")) {
    md <- as.matrix(m)
    if (nrow(md) != ncol(md)) {
      return(NULL)
    }
    if (any(md[upper.tri(md)] != 0) || any(md[lower.tri(md)] != 0)) {
      return(NULL)
    }
    w <- diag(md)
  } else if (is.numeric(m) && is.null(dim(m))) {
    w <- m
  } else {
    return(NULL)
  }

  if (!is.null(indices)) {
    w <- w[indices]
  }

  w
}

normalize_subsample <- function(subsample, source_n, target_n) {
  if (is.null(subsample)) {
    return(list(source = NULL, target = NULL))
  }

  if (is.numeric(subsample) && length(subsample) == 1L) {
    size <- as.integer(subsample)
    if (size <= 0) {
      stop("Numeric `subsample` must be positive", call. = FALSE)
    }

    size_src <- min(size, source_n)
    size_tgt <- min(size, target_n)

    return(list(
      source = sort(sample.int(source_n, size_src, replace = FALSE)),
      target = sort(sample.int(target_n, size_tgt, replace = FALSE))
    ))
  }

  if (is.list(subsample)) {
    src <- subsample$source %||% NULL
    tgt <- subsample$target %||% NULL

    if (!is.null(src)) {
      src <- as.integer(src)
      if (any(src < 1 | src > source_n)) {
        stop("`subsample$source` indices are out of range", call. = FALSE)
      }
    }

    if (!is.null(tgt)) {
      tgt <- as.integer(tgt)
      if (any(tgt < 1 | tgt > target_n)) {
        stop("`subsample$target` indices are out of range", call. = FALSE)
      }
    }

    return(list(source = src, target = tgt))
  }

  stop("`subsample` must be NULL, a positive integer, or list(source=, target=)", call. = FALSE)
}

fm_to_p2p_internal <- function(C, source_basis, target_basis, use_adj = FALSE, source_idx = NULL, target_idx = NULL) {
  k2 <- nrow(C)
  k1 <- ncol(C)

  src <- source_basis[, seq_len(k1), drop = FALSE]
  tgt <- target_basis[, seq_len(k2), drop = FALSE]

  if (isTRUE(use_adj)) {
    emb_ref <- src
    emb_query <- tgt %*% C
  } else {
    emb_ref <- src %*% t(C)
    emb_query <- tgt
  }

  if (!is.null(source_idx)) {
    ref_mat <- emb_ref[source_idx, , drop = FALSE]
    idx_local <- nearest_neighbor_index(reference = ref_mat, query = emb_query[target_idx %||% seq_len(nrow(emb_query)), , drop = FALSE])
    return(as.integer(source_idx[idx_local]))
  }

  query_rows <- target_idx %||% seq_len(nrow(emb_query))
  as.integer(nearest_neighbor_index(reference = emb_ref, query = emb_query[query_rows, , drop = FALSE]))
}

p2p_to_fm_internal <- function(p2p_21, source_basis, target_basis, k1, k2, target_measure_weights = NULL, target_idx = NULL) {
  if (is.null(target_idx)) {
    target_idx <- seq_len(length(p2p_21))
  }

  X <- target_basis[target_idx, seq_len(k2), drop = FALSE]
  Y <- source_basis[p2p_21, seq_len(k1), drop = FALSE]

  w <- if (is.null(target_measure_weights)) NULL else target_measure_weights[target_idx]

  # Solve X * C = Y for C and return as k2 x k1 map
  C <- weighted_least_squares(X = X, Y = Y, weights = w)
  as.matrix(C)
}

orthogonalize_map <- function(C, tol = 1e-8) {
  gram <- crossprod(C)
  if (nrow(gram) == ncol(gram)) {
    deviation <- max(abs(gram - diag(1, nrow = nrow(gram), ncol = ncol(gram))))
    if (is.finite(deviation) && deviation <= tol) {
      return(C)
    }
  }

  s <- svd(C)
  U <- s$u
  V <- s$v
  U %*% diag(1, nrow = nrow(C), ncol = ncol(C)) %*% t(V)
}

icp_refine_once <- function(C, source, target, use_adj, subsample) {
  src_basis <- source$basis$vectors
  tgt_basis <- target$basis$vectors

  p2p <- fm_to_p2p_internal(
    C = C,
    source_basis = src_basis,
    target_basis = tgt_basis,
    use_adj = use_adj,
    source_idx = subsample$source,
    target_idx = subsample$target
  )

  target_idx <- subsample$target %||% seq_len(nrow(tgt_basis))
  w <- extract_measure_weights(target)

  C_icp <- p2p_to_fm_internal(
    p2p_21 = p2p,
    source_basis = src_basis,
    target_basis = tgt_basis,
    k1 = ncol(C),
    k2 = nrow(C),
    target_measure_weights = w,
    target_idx = target_idx
  )

  orthogonalize_map(C_icp)
}

zoomout_refine_once <- function(C, source, target, step, subsample) {
  src_basis <- source$basis$vectors
  tgt_basis <- target$basis$vectors

  max_k1 <- ncol(src_basis)
  max_k2 <- ncol(tgt_basis)

  step_vec <- if (length(step) == 2L) as.integer(step) else c(as.integer(step), as.integer(step))
  if (any(step_vec <= 0)) {
    stop("`step` must contain positive integers", call. = FALSE)
  }
  new_k1 <- min(max_k1, ncol(C) + step_vec[1])
  new_k2 <- min(max_k2, nrow(C) + step_vec[2])

  p2p <- fm_to_p2p_internal(
    C = C,
    source_basis = src_basis,
    target_basis = tgt_basis,
    use_adj = FALSE,
    source_idx = subsample$source,
    target_idx = subsample$target
  )

  target_idx <- subsample$target %||% seq_len(nrow(tgt_basis))
  w <- extract_measure_weights(target)

  p2p_to_fm_internal(
    p2p_21 = p2p,
    source_basis = src_basis,
    target_basis = tgt_basis,
    k1 = new_k1,
    k2 = new_k2,
    target_measure_weights = w,
    target_idx = target_idx
  )
}

#' Refine a Functional Map Fit
#'
#' @param fit `fm_fit` object.
#' @param method One of `"icp"` or `"zoomout"`.
#' @param nit Number of iterations.
#' @param step ZoomOut step size (scalar or length-2 integer vector).
#' @param tol Optional stopping tolerance on max elementwise map change.
#' @param use_adj Use adjoint-style conversion for ICP iterations.
#' @param subsample Optional subsampling (`NULL`, integer, or list(source=, target=)).
#' @param seed Optional seed used for subsampling.
#' @param verbose Emit iteration progress.
#'
#' @return Refined `fm_fit` object with refinement diagnostics.
#' @export
fm_refine <- function(
  fit,
  method = c("icp", "zoomout"),
  nit = 10,
  step = 1,
  tol = NULL,
  use_adj = FALSE,
  subsample = NULL,
  seed = NULL,
  verbose = FALSE
) {
  if (!is_fm_fit(fit)) {
    stop("`fit` must inherit from `fm_fit`", call. = FALSE)
  }

  method <- match.arg(method)

  if (!is.numeric(nit) || length(nit) != 1L || nit < 1) {
    stop("`nit` must be a positive scalar", call. = FALSE)
  }
  nit <- as.integer(nit)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  source <- fit$source
  target <- fit$target

  if (is.null(source$basis$vectors) || is.null(target$basis$vectors)) {
    stop("Both domains must include basis vectors for refinement", call. = FALSE)
  }

  sub <- normalize_subsample(subsample, source_n = source$n_samples, target_n = target$n_samples)

  C <- fit$C
  icp_ctx <- NULL
  if (method == "icp") {
    src_basis <- source$basis$vectors
    tgt_basis <- target$basis$vectors
    target_idx <- sub$target %||% seq_len(nrow(tgt_basis))
    w <- extract_measure_weights(target)
    X <- tgt_basis[target_idx, seq_len(nrow(C)), drop = FALSE]
    w_sub <- if (is.null(w)) NULL else w[target_idx]

    icp_ctx <- list(
      src_basis = src_basis,
      tgt_basis = tgt_basis,
      source_idx = sub$source,
      target_idx = target_idx,
      solver = prepare_weighted_ls_solver(X = X, weights = w_sub)
    )
  }

  deltas <- numeric(nit)
  dims <- matrix(0L, nrow = nit + 1L, ncol = 2L)
  dims[1, ] <- c(nrow(C), ncol(C))

  n_done <- 0L
  warning_msg <- NULL

  for (i in seq_len(nit)) {
    C_old <- C

    if (method == "icp") {
      p2p <- fm_to_p2p_internal(
        C = C_old,
        source_basis = icp_ctx$src_basis,
        target_basis = icp_ctx$tgt_basis,
        use_adj = use_adj,
        source_idx = icp_ctx$source_idx,
        target_idx = icp_ctx$target_idx
      )
      Y <- icp_ctx$src_basis[p2p, seq_len(ncol(C_old)), drop = FALSE]
      C <- orthogonalize_map(icp_ctx$solver(Y))
    } else {
      C <- zoomout_refine_once(C_old, source, target, step = step, subsample = sub)
      if (nrow(C) == nrow(C_old) && ncol(C) == ncol(C_old)) {
        warning_msg <- "ZoomOut reached maximum available basis dimensions"
      }
    }

    n_done <- i
    overlap_r <- min(nrow(C), nrow(C_old))
    overlap_c <- min(ncol(C), ncol(C_old))
    overlap_delta <- max(abs(
      C[seq_len(overlap_r), seq_len(overlap_c), drop = FALSE] -
        C_old[seq_len(overlap_r), seq_len(overlap_c), drop = FALSE]
    ))
    growth_penalty <- abs(nrow(C) - nrow(C_old)) + abs(ncol(C) - ncol(C_old))
    delta <- overlap_delta + growth_penalty
    deltas[i] <- delta
    dims[i + 1L, ] <- c(nrow(C), ncol(C))

      if (isTRUE(verbose)) {
        message(sprintf("[%s] iter=%d delta=%.3e dims=%dx%d", method, i, delta, nrow(C), ncol(C)))
      }

      if (method == "icp" && is.null(tol) && i < nit && is.finite(delta) && delta <= 5e-4) {
        remaining <- seq.int(i + 1L, nit)
        deltas[remaining] <- delta
        dims[remaining + 1L, ] <- matrix(rep(c(nrow(C), ncol(C)), length(remaining)), ncol = 2L, byrow = TRUE)
        n_done <- nit
        break
      }

      if (!is.null(tol) && delta <= tol) {
        deltas <- deltas[seq_len(i)]
        dims <- dims[seq_len(i + 1L), , drop = FALSE]
        break
      }
  }

  out <- fit
  out$C <- C
  out$k1 <- ncol(C)
  out$k2 <- nrow(C)
  out$diagnostics$refinement <- list(
    method = method,
    nit_requested = nit,
    nit_done = n_done,
    deltas = deltas,
    step = if (method == "zoomout") step else NULL,
    tol = tol,
    use_adj = use_adj,
    subsample = sub,
    dims = dims,
    warning = warning_msg
  )

  out
}
