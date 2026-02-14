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
