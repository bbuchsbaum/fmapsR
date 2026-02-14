make_small_generic_domain <- function(n = 10, k = 4, seed = 1) {
  set.seed(seed)
  vals <- seq(1, n, by = 1)
  op <- diag(vals)
  fm_basis(
    fm_domain_generic(n_samples = n, operator = op),
    k = k,
    solver = "base",
    cache = FALSE
  )
}

test_that("end-to-end network workflow produces synchronized consensus and serializes", {
  n <- 10
  p <- 4
  k <- 4

  d1 <- make_small_generic_domain(n = n, k = k, seed = 11)
  d2 <- make_small_generic_domain(n = n, k = k, seed = 12)
  d3 <- make_small_generic_domain(n = n, k = k, seed = 13)

  set.seed(99)
  desc1 <- matrix(rnorm(n * p), nrow = n, ncol = p)
  desc2 <- desc1 + matrix(rnorm(n * p, sd = 0.02), nrow = n, ncol = p)
  desc3 <- desc1 + matrix(rnorm(n * p, sd = 0.03), nrow = n, ncol = p)

  fit12 <- fm_refine(
    fm_match(
      source = d1,
      target = d2,
      descriptors = list(source = desc1, target = desc2),
      optimizer = "cg",
      cg_maxit = 2
    ),
    method = "icp",
    nit = 2
  )
  fit23 <- fm_refine(
    fm_match(
      source = d2,
      target = d3,
      descriptors = list(source = desc2, target = desc3),
      optimizer = "cg",
      cg_maxit = 2
    ),
    method = "icp",
    nit = 2
  )
  fit31 <- fm_refine(
    fm_match(
      source = d3,
      target = d1,
      descriptors = list(source = desc3, target = desc1),
      optimizer = "cg",
      cg_maxit = 2
    ),
    method = "icp",
    nit = 2
  )

  net <- fm_network(list(a = d1, b = d2, c = d3), directed = TRUE)
  net <- fm_network_add_map(net, "a", "b", fit12)
  net <- fm_network_add_map(net, "b", "c", fit23)
  net <- fm_network_add_map(net, "c", "a", fit31)

  synced <- fm_sync(net, mode = "cycle", nit = 10)
  cons <- fm_consensus(synced$network)

  expect_s3_class(synced, "fm_network_fit")
  expect_s3_class(cons, "fm_consensus")
  if (isTRUE(all.equal(synced$diagnostics$pre_cycle_median, 0))) {
    expect_true(is.na(synced$diagnostics$cycle_reduction_ratio))
  } else {
    expect_true(is.finite(synced$diagnostics$cycle_reduction_ratio))
  }
  expect_true(length(cons$alignments) == 3)
  expect_true(length(cons$maps) >= 6)

  path <- tempfile(fileext = ".rds")
  fm_consensus_save(cons, path)
  cons_loaded <- fm_consensus_load(path)
  expect_equal(dim(cons_loaded$latent_basis), dim(cons$latent_basis))
  expect_equal(names(cons_loaded$alignments), names(cons$alignments))
})

test_that("convergence and contract diagnostics are actionable", {
  d1 <- make_small_generic_domain(n = 10, k = 4, seed = 21)
  d2 <- make_small_generic_domain(n = 10, k = 4, seed = 22)

  set.seed(123)
  src_desc <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
  tgt_desc <- src_desc + matrix(rnorm(10 * 3, sd = 0.05), nrow = 10, ncol = 3)

  fit <- fm_match(
    source = d1,
    target = d2,
    descriptors = list(source = src_desc, target = tgt_desc),
    optimizer = "cg",
    cg_maxit = 1,
    cg_tol = 1e-30
  )

  expect_identical(fit$diagnostics$convergence, 1L)
  expect_match(fit$diagnostics$message, "iteration limit")
  expect_match(fit$diagnostics$warning, "convergence code 1")

  bad_desc <- matrix(rnorm(8 * 3), nrow = 8, ncol = 3)
  expect_error(
    fm_match(
      source = d1,
      target = d2,
      descriptors = list(source = bad_desc, target = tgt_desc)
    ),
    "Descriptor row counts must match domain sample counts"
  )
})
