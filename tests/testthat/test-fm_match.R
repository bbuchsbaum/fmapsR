make_test_domain <- function(n) {
  op <- diag(seq(1, n, by = 1))
  domain <- fm_domain_generic(n_samples = n, operator = op)
  fm_basis(domain, k = min(4, n), solver = "base", cache = FALSE)
}

test_that("fm_match returns fm_fit with diagnostics and objective terms", {
  source <- make_test_domain(6)
  target <- make_test_domain(6)

  src_desc <- cbind(seq_len(6), seq_len(6)^2)
  tgt_desc <- src_desc

  fit <- fm_match(
    source = source,
    target = target,
    descriptors = list(source = src_desc, target = tgt_desc),
    penalties = list(descr = 1, lap = 1e-3, comm = 1e-1),
    init = "identity",
    maxit = 30
  )

  expect_s3_class(fit, "fm_fit")
  expect_identical(dim(fm_fit_matrix(fit)), c(4L, 4L))
  expect_true(is.list(fit$diagnostics$objective_terms))
  expect_true(all(c("descr", "lap", "comm") %in% names(fit$diagnostics$objective_terms)))
})

test_that("fm_match supports disabling individual objective terms", {
  source <- make_test_domain(6)
  target <- make_test_domain(6)

  src_desc <- cbind(seq_len(6), seq_len(6)^2)
  tgt_desc <- src_desc

  fit <- fm_match(
    source = source,
    target = target,
    descriptors = list(source = src_desc, target = tgt_desc),
    penalties = list(descr = 1, lap = 0, comm = 0),
    init = "zeros",
    maxit = 20
  )

  expect_equal(fit$diagnostics$objective_terms$lap, 0)
  expect_equal(fit$diagnostics$objective_terms$comm, 0)
})

test_that("fm_match supports compiled kernel backend with parity", {
  skip_if_not(exists("fm_match_solve_cg_cpp", mode = "function"))

  source <- make_test_domain(10)
  target <- make_test_domain(10)

  set.seed(101)
  src_desc <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
  tgt_desc <- src_desc + matrix(rnorm(10 * 3, sd = 0.01), nrow = 10, ncol = 3)

  fit_r <- fm_match(
    source = source,
    target = target,
    descriptors = list(source = src_desc, target = tgt_desc),
    penalties = list(descr = 1e-1, lap = 1e-3, comm = 1),
    init = "identity",
    optimizer = "cg",
    cg_maxit = 10,
    cg_tol = 1e-4,
    kernel_backend = "r"
  )

  fit_cpp <- fm_match(
    source = source,
    target = target,
    descriptors = list(source = src_desc, target = tgt_desc),
    penalties = list(descr = 1e-1, lap = 1e-3, comm = 1),
    init = "identity",
    optimizer = "cg",
    cg_maxit = 10,
    cg_tol = 1e-4,
    kernel_backend = "cpp"
  )

  expect_identical(fit_r$diagnostics$kernel_backend, "r")
  expect_identical(fit_cpp$diagnostics$kernel_backend, "cpp")
  expect_equal(fit_cpp$diagnostics$total_objective, fit_r$diagnostics$total_objective, tolerance = 1e-8)
  expect_equal(fit_cpp$C, fit_r$C, tolerance = 1e-8)
})

test_that("fm_match supports streamed descriptor batching for commutativity terms", {
  source <- make_test_domain(24)
  target <- make_test_domain(24)

  set.seed(222)
  src_desc <- matrix(rnorm(24 * 48), nrow = 24, ncol = 48)
  tgt_desc <- src_desc + matrix(rnorm(24 * 48, sd = 0.03), nrow = 24, ncol = 48)

  fit_pre <- fm_match(
    source = source,
    target = target,
    descriptors = list(source = src_desc, target = tgt_desc),
    penalties = list(descr = 1e-1, lap = 1e-3, comm = 1),
    init = "identity",
    optimizer = "cg",
    cg_maxit = 6,
    kernel_backend = "r"
  )

  fit_batch <- fm_match(
    source = source,
    target = target,
    descriptors = list(source = src_desc, target = tgt_desc),
    penalties = list(descr = 1e-1, lap = 1e-3, comm = 1),
    init = "identity",
    optimizer = "cg",
    cg_maxit = 6,
    kernel_backend = "r",
    descriptor_batch_size = 8
  )

  expect_identical(fit_batch$diagnostics$commutativity_mode, "stream")
  expect_true(fit_batch$diagnostics$commutativity_batches >= 2L)
  expect_identical(fit_batch$diagnostics$descriptor_batch_size, 8L)
  expect_equal(fit_batch$diagnostics$total_objective, fit_pre$diagnostics$total_objective, tolerance = 1e-6)
  expect_equal(fit_batch$C, fit_pre$C, tolerance = 1e-6)
})
