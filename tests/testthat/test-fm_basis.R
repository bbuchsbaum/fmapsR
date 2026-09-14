test_that("fm_basis computes basis with explicit k", {
  fm_basis_cache_clear()

  op <- diag(c(0.1, 0.2, 0.5, 1.0, 2.0))
  domain <- fm_domain_generic(n_samples = 5, operator = op)

  out <- fm_basis(domain, k = 3, solver = "base", cache = FALSE)

  expect_equal(length(out$basis$values), 3)
  expect_identical(dim(out$basis$vectors), c(5L, 3L))
  expect_identical(out$basis$k, 3L)
  expect_identical(out$basis$solver, "base")
})

test_that("fm_basis caching is used on repeated calls", {
  fm_basis_cache_clear()

  op <- diag(c(0.1, 0.2, 0.5, 1.0, 2.0))
  domain <- fm_domain_generic(n_samples = 5, operator = op)

  out1 <- fm_basis(domain, k = 3, solver = "base", cache = TRUE)
  out2 <- fm_basis(domain, k = 3, solver = "base", cache = TRUE)

  expect_false(isTRUE(out1$basis$cached))
  expect_true(isTRUE(out2$basis$cached))
  expect_equal(out1$basis$values, out2$basis$values)
})

test_that("seeded fm_basis calls are reproducible", {
  fm_basis_cache_clear()

  op <- diag(seq(1, 10, by = 1))
  domain <- fm_domain_generic(n_samples = 10, operator = op)

  out1 <- fm_basis(domain, k = 4, solver = "rspectra", seed = 123, cache = FALSE)
  out2 <- fm_basis(domain, k = 4, solver = "rspectra", seed = 123, cache = FALSE)

  expect_equal(out1$basis$values, out2$basis$values, tolerance = 1e-8)

  # Eigenvectors can differ by sign; compare absolute per-column dot products.
  dots <- abs(colSums(out1$basis$vectors * out2$basis$vectors))
  expect_true(all(dots > 0.999))
})

test_that("fm_basis honors operator-provided RSpectra hints for cotangent meshes", {
  skip_if_not_installed("RSpectra")

  pair <- fm_example_mesh_pair()
  op <- fm_operator_mesh(
    pair$source$vertices,
    pair$source$faces,
    method = "cotangent"
  )
  domain <- fm_domain_mesh(pair$source$vertices, pair$source$faces, operator = op)

  out_auto <- fm_basis(domain, k = 16, solver = "rspectra", cache = FALSE)
  expect_identical(out_auto$basis$solver, "rspectra")
  expect_identical(out_auto$basis$which, "LM")
  expect_identical(out_auto$basis$control$sigma, 0)
  expect_gte(out_auto$basis$control$opts$ncv, 40)
  expect_true(all(is.finite(out_auto$basis$values)))

  op_hint <- Matrix::Diagonal(x = seq(1, 24))
  attr(op_hint, "fm_rspectra_preferred_which") <- "LM"
  attr(op_hint, "fm_rspectra_sigma") <- 0
  attr(op_hint, "fm_rspectra_ncv_min") <- 40L
  attr(op_hint, "fm_rspectra_ncv_mult") <- 4L
  hinted_domain <- fm_domain_generic(n_samples = 24, operator = op_hint)

  out_override <- fm_basis(
    hinted_domain,
    k = 16,
    solver = "rspectra",
    which = "SM",
    control = list(opts = list(ncv = 18)),
    cache = FALSE
  )
  expect_identical(out_override$basis$which, "SM")
  expect_identical(out_override$basis$control$opts$ncv, 18)
  expect_null(out_override$basis$control$sigma)
  expect_equal(sort(out_override$basis$values), as.numeric(1:16), tolerance = 1e-8)
})

test_that("sparse fallback preserves shift-invert spectral selection", {
  op <- diag(1:8)
  attr(op, "fm_rspectra_preferred_which") <- "LM"
  attr(op, "fm_rspectra_sigma") <- 0
  d <- fm_domain_generic(8, operator = op)
  local_mocked_bindings(compute_basis_rspectra = function(...) NULL)
  out <- fm_basis(d, 3, cache = FALSE)
  expect_identical(out$basis$solver, "base")
  expect_equal(out$basis$values, as.numeric(1:3))
  expect_equal(fm_basis(d, 3, which = "LM", cache = FALSE)$basis$values, as.numeric(8:6))
  expect_equal(fm_basis(d, 8, cache = FALSE)$basis$values, as.numeric(1:8))
})

test_that("partial sparse convergence is not accepted as a complete basis", {
  local_mocked_bindings(eigs_sym = function(...) list(values = 1, vectors = matrix(1, 8, 1)),
    .package = "RSpectra")
  d <- fm_domain_generic(8, operator = diag(1:8))
  out <- fm_basis(d, 3, cache = FALSE)
  expect_identical(out$basis$solver, "base")
  expect_equal(out$basis$values, as.numeric(1:3))
  expect_identical(dim(out$basis$vectors), c(8L, 3L))
})

test_that("dense and sparse bases put the constant Laplacian mode first", {
  n <- 12L
  A <- matrix(0, n, n)
  A[cbind(1:(n - 1), 2:n)] <- 1
  A <- A + t(A)
  L <- diag(rowSums(A)) - A
  d <- fm_domain_graph(A, operator = L)
  base <- fm_basis(d, 4, solver = "base", cache = FALSE)$basis
  sparse <- fm_basis(d, 4, solver = "rspectra", cache = FALSE)$basis
  expect_identical(sparse$solver, "rspectra")
  expect_equal(sparse$values, base$values, tolerance = 1e-9)
  expect_equal(sparse$vectors[, 1], rep(1 / sqrt(n), n), tolerance = 1e-9)
  expect_equal(base$vectors[, 1], rep(1 / sqrt(n), n), tolerance = 1e-9)
  expect_equal(L %*% sparse$vectors, sparse$vectors %*% diag(sparse$values), tolerance = 1e-9)

  op <- fm_operator_mesh(
    rbind(c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(0, 1, 0)),
    rbind(c(1, 2, 3), c(1, 3, 4))
  )
  mesh <- fm_domain_generic(4, operator = op, measure = attr(op, "fm_measure_weights"))
  out <- fm_basis(mesh, 2, cache = FALSE)
  expect_equal(out$basis$values[1], 0, tolerance = 1e-9)
  expect_equal(out$basis$vectors[, 1], rep(1, 4), tolerance = 1e-9)
})
