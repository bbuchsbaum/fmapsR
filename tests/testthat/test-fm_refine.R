make_manual_fit <- function(n = 8, k_basis = 6, k_map = 3) {
  op <- diag(seq(1, n, by = 1))
  source <- fm_basis(fm_domain_generic(n_samples = n, operator = op), k = k_basis, solver = "base", cache = FALSE)
  target <- fm_basis(fm_domain_generic(n_samples = n, operator = op), k = k_basis, solver = "base", cache = FALSE)

  structure(
    list(
      C = diag(1, nrow = k_map, ncol = k_map),
      k1 = k_map,
      k2 = k_map,
      source = source,
      target = target,
      penalties = list(),
      diagnostics = list(
        convergence = 0,
        message = "manual",
        counts = c(fn = NA_integer_, gradient = NA_integer_),
        total_objective = NA_real_,
        objective_terms = list(),
        warning = NULL
      )
    ),
    class = "fm_fit"
  )
}

test_that("icp refinement returns fm_fit with iteration diagnostics", {
  fit <- make_manual_fit()

  refined <- fm_refine(fit, method = "icp", nit = 3)

  expect_s3_class(refined, "fm_fit")
  expect_identical(dim(refined$C), c(3L, 3L))
  expect_identical(refined$diagnostics$refinement$method, "icp")
  expect_identical(refined$diagnostics$refinement$nit_done, 3L)
  expect_length(refined$diagnostics$refinement$deltas, 3)
})

test_that("zoomout refinement increases map dimensions", {
  fit <- make_manual_fit()

  refined <- fm_refine(fit, method = "zoomout", nit = 2, step = 1)

  expect_identical(dim(refined$C), c(5L, 5L))
  expect_identical(refined$k1, 5L)
  expect_identical(refined$k2, 5L)
})

test_that("refinement supports optional subsampling", {
  fit <- make_manual_fit()

  refined <- fm_refine(fit, method = "icp", nit = 1, subsample = 4, seed = 7)

  expect_length(refined$diagnostics$refinement$subsample$source, 4)
  expect_length(refined$diagnostics$refinement$subsample$target, 4)

  custom_sub <- list(source = 1:3, target = 2:4)
  refined2 <- fm_refine(fit, method = "icp", nit = 1, subsample = custom_sub)
  expect_equal(refined2$diagnostics$refinement$subsample$source, 1:3)
  expect_equal(refined2$diagnostics$refinement$subsample$target, 2:4)

  expect_error(
    fm_refine(fit, method = "icp", nit = 1, subsample = list(source = c(0, 1))),
    "out of range"
  )
})

test_that("zoomout reports warning when basis dimension cannot increase", {
  fit <- make_manual_fit(k_map = 6)

  refined <- fm_refine(fit, method = "zoomout", nit = 2, step = 2)

  expect_identical(dim(refined$C), c(6L, 6L))
  expect_match(refined$diagnostics$refinement$warning, "maximum available basis")
})

test_that("weighted least squares helpers validate dimensions and weights", {
  X <- diag(3)
  Y <- matrix(1:3, ncol = 1)

  expect_error(weighted_least_squares(X, matrix(1:2, ncol = 1)), "row count must match")
  expect_error(weighted_least_squares(X, Y, weights = c(1, 2)), "weights")

  coef <- weighted_least_squares(X, Y, weights = c(1, 4, 9))
  expect_equal(as.numeric(coef), as.numeric(Y))

  expect_error(prepare_weighted_ls_solver(X, weights = c(1, 2)), "weights")
  solver <- prepare_weighted_ls_solver(X)
  expect_error(solver(matrix(1:2, ncol = 1)), "row count must match")

  solver_w <- prepare_weighted_ls_solver(X, weights = c(1, 1, 1))
  expect_equal(as.numeric(solver_w(Y)), as.numeric(Y))

  Q <- qr.Q(qr(matrix(rnorm(18), nrow = 6, ncol = 3)))
  Y2 <- matrix(rnorm(12), nrow = 6, ncol = 2)
  expect_equal(weighted_least_squares(Q, Y2), crossprod(Q, Y2), tolerance = 1e-10)
  solver_q <- prepare_weighted_ls_solver(Q)
  expect_equal(solver_q(Y2), crossprod(Q, Y2), tolerance = 1e-10)
})

test_that("subsampling and measure extraction helpers handle edge cases", {
  fit <- make_manual_fit()

  expect_error(normalize_subsample(0, 5, 5), "must be positive")
  expect_error(normalize_subsample("bad", 5, 5), "`subsample` must be NULL")
  expect_error(normalize_subsample(list(target = 0), 5, 5), "out of range")

  sub <- normalize_subsample(50, source_n = 8, target_n = 6)
  expect_length(sub$source, 8)
  expect_length(sub$target, 6)

  expect_equal(extract_measure_weights(fit$target), rep(1, fit$target$n_samples))
  fit$target$measure <- diag(seq_len(fit$target$n_samples))
  expect_equal(extract_measure_weights(fit$target, indices = 1:3), 1:3)
  fit$target$measure <- matrix(1, nrow = fit$target$n_samples, ncol = fit$target$n_samples)
  expect_null(extract_measure_weights(fit$target))
})

test_that("internal p2p/fm conversions and orthogonalization cover alternate paths", {
  fit <- make_manual_fit(n = 8, k_basis = 6, k_map = 3)
  C <- fit$C
  src <- fit$source$basis$vectors
  tgt <- fit$target$basis$vectors

  p2p_adj <- fm_to_p2p_internal(C, src, tgt, use_adj = TRUE)
  expect_equal(length(p2p_adj), nrow(tgt))

  p2p_sub <- fm_to_p2p_internal(
    C, src, tgt, use_adj = FALSE,
    source_idx = c(1L, 3L, 5L, 7L),
    target_idx = c(2L, 4L)
  )
  expect_equal(length(p2p_sub), 2L)
  expect_true(all(p2p_sub %in% c(1L, 3L, 5L, 7L)))

  C2 <- p2p_to_fm_internal(
    p2p_21 = p2p_adj,
    source_basis = src,
    target_basis = tgt,
    k1 = 3,
    k2 = 3
  )
  expect_equal(dim(C2), c(3L, 3L))

  C_ortho <- orthogonalize_map(diag(3))
  expect_equal(C_ortho, diag(3))
  C_non <- matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), nrow = 3)
  C_fix <- orthogonalize_map(C_non)
  expect_equal(dim(C_fix), c(3L, 3L))
})

test_that("single-step icp/zoomout helpers and fm_refine validations are exercised", {
  fit <- make_manual_fit(n = 8, k_basis = 6, k_map = 3)
  sub <- list(source = 1:4, target = 2:5)

  C_icp <- icp_refine_once(fit$C, fit$source, fit$target, use_adj = FALSE, subsample = sub)
  expect_equal(dim(C_icp), c(3L, 3L))

  C_zoom <- zoomout_refine_once(fit$C, fit$source, fit$target, step = c(1, 2), subsample = sub)
  expect_equal(dim(C_zoom), c(5L, 4L))
  expect_error(
    zoomout_refine_once(fit$C, fit$source, fit$target, step = c(0, 1), subsample = sub),
    "positive integers"
  )

  expect_error(fm_refine(list(), method = "icp"), "`fit` must inherit")
  expect_error(fm_refine(fit, method = "icp", nit = 0), "`nit` must be a positive scalar")

  fit_no_basis <- fit
  fit_no_basis$source$basis$vectors <- NULL
  expect_error(fm_refine(fit_no_basis, method = "icp"), "must include basis vectors")

  refined_tol <- fm_refine(fit, method = "icp", nit = 5, tol = 1e9, use_adj = TRUE, verbose = TRUE)
  expect_lte(length(refined_tol$diagnostics$refinement$deltas), 5L)

  expect_error(fm_refine(fit, method = "icp", nit = 1, subsample = "bad"), "`subsample` must be NULL")
})
