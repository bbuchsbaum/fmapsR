make_identity_fit <- function(n = 6, k = 4) {
  op <- diag(seq(1, n, by = 1))
  src <- fm_basis(fm_domain_generic(n_samples = n, operator = op), k = k, solver = "base", cache = FALSE)
  tgt <- fm_basis(fm_domain_generic(n_samples = n, operator = op), k = k, solver = "base", cache = FALSE)

  structure(
    list(
      C = diag(1, nrow = k, ncol = k),
      k1 = k,
      k2 = k,
      source = src,
      target = tgt,
      penalties = list(),
      diagnostics = list(convergence = 0, message = "manual")
    ),
    class = "fm_fit"
  )
}

test_that("fm_transfer preserves signal under identity fit", {
  fit <- make_identity_fit()
  x <- matrix(seq_len(fit$source$n_samples), ncol = 1)

  y <- fm_transfer(fit, x)
  expected <- fm_unproject(fit$source, fm_project(fit$source, x))

  expect_equal(as.numeric(y), as.numeric(expected), tolerance = 1e-8)
})

test_that("as_p2p returns integer map and actionable errors", {
  fit <- make_identity_fit()
  p2p <- as_p2p(fit)

  expect_type(p2p, "integer")
  expect_length(p2p, fit$target$n_samples)

  fit_bad <- fit
  fit_bad$source$basis$vectors <- NULL
  expect_error(
    as_p2p(fit_bad),
    "requires basis vectors"
  )
})

test_that("fm_compose and fm_inverse enforce shape consistency", {
  fit1 <- make_identity_fit()
  fit2 <- make_identity_fit()

  comp <- fm_compose(fit1, fit2)
  inv <- fm_inverse(fit1)

  expect_s3_class(comp, "fm_fit")
  expect_identical(dim(comp$C), c(4L, 4L))
  expect_s3_class(inv, "fm_fit")
  expect_identical(dim(inv$C), c(4L, 4L))

  fit_bad <- fit2
  fit_bad$C <- diag(1, 3, 3)
  expect_error(fm_compose(fit1, fit_bad), "incompatible")
})

test_that("fm_inverse and composition provide round-trip consistency", {
  fit <- make_identity_fit(n = 8, k = 4)
  inv_pinv <- fm_inverse(fit, method = "pseudoinverse")
  inv_t <- fm_inverse(fit, method = "transpose")

  comp_pinv <- fm_compose(fit, inv_pinv)
  comp_t <- fm_compose(fit, inv_t)

  expect_equal(comp_pinv$C, diag(4), tolerance = 1e-8)
  expect_equal(comp_t$C, diag(4), tolerance = 1e-8)
})

test_that("fm_transfer supports reverse direction with matching dimensions", {
  fit <- make_identity_fit(n = 7, k = 4)
  x <- matrix(seq_len(fit$target$n_samples), ncol = 1)

  y_rev <- fm_transfer(fit, x, reverse = TRUE)
  expected <- fm_unproject(fit$source, fm_project(fit$source, x))

  expect_equal(as.numeric(y_rev), as.numeric(expected), tolerance = 1e-8)
})

test_that("as_p2p guards against oversized dense nearest-neighbor problems", {
  fit <- make_identity_fit(n = 8, k = 4)
  fit$source$n_samples <- 10000L
  fit$target$n_samples <- 6000L
  fit$source$basis$vectors <- matrix(0, nrow = 10000, ncol = 4)
  fit$target$basis$vectors <- matrix(0, nrow = 6000, ncol = 4)

  expect_error(as_p2p(fit), "too large for dense nearest-neighbor search")
})

test_that("native nearest-neighbor kernel matches dense distance argmin", {
  skip_if_not(exists("fm_nearest_neighbor_index_cpp", mode = "function"))

  set.seed(7)
  reference <- matrix(rnorm(40), nrow = 10, ncol = 4)
  query <- matrix(rnorm(24), nrow = 6, ncol = 4)

  d2 <- as.matrix(dist(rbind(query, reference)))^2
  d2 <- d2[seq_len(nrow(query)), nrow(query) + seq_len(nrow(reference)), drop = FALSE]
  expected <- max.col(-d2)

  nn_cpp <- fm_nearest_neighbor_index_cpp(reference, query)
  expect_equal(as.integer(nn_cpp), as.integer(expected))
  expect_equal(nearest_neighbor_index(reference, query), as.integer(expected))
})
