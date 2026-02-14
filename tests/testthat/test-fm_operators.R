set.seed(202)

make_operator_fixture <- function(n = 24, p = 10) {
  matrix(rnorm(n * p), nrow = n, ncol = p)
}

test_that("covariance operator is symmetric and basis-ready", {
  X <- make_operator_fixture()
  op <- fm_operator_covariance(X, ridge = 1e-3)

  expect_identical(dim(op), c(nrow(X), nrow(X)))
  expect_equal(op, t(op), tolerance = 1e-10)

  d <- fm_domain_generic(n_samples = nrow(X), data = X, operator = op)
  d <- fm_basis(d, k = 6, solver = "base", cache = FALSE)
  expect_true(is.matrix(d$basis$vectors))
  expect_identical(ncol(d$basis$vectors), 6L)
})

test_that("kNN laplacian operator is symmetric and basis-ready", {
  X <- make_operator_fixture(n = 30, p = 8)
  op <- fm_operator_knn_laplacian(X, k = 6, normalized = TRUE)

  expect_identical(dim(op), c(nrow(X), nrow(X)))
  expect_equal(op, t(op), tolerance = 1e-10)
  expect_true(all(diag(op) >= 0))

  d <- fm_domain_generic(n_samples = nrow(X), data = X, operator = op)
  d <- fm_basis(d, k = 5, solver = "base", cache = FALSE)
  expect_identical(ncol(d$basis$vectors), 5L)
})

test_that("fm_operator dispatches to supported families", {
  X <- make_operator_fixture(n = 18, p = 6)

  op_cov <- fm_operator(X, method = "covariance")
  op_knn <- fm_operator(X, method = "knn_laplacian", k = 4)

  expect_identical(dim(op_cov), c(18L, 18L))
  expect_identical(dim(op_knn), c(18L, 18L))
  expect_false(isTRUE(all.equal(op_cov, op_knn)))
})
