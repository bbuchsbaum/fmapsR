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
