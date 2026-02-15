make_real_domain <- function(x, k = 8L) {
  op <- fm_operator_covariance(x, center = TRUE, scale = TRUE, ridge = 1e-2)
  dom <- fm_domain_generic(n_samples = nrow(x), data = x, operator = op)
  fm_basis(dom, k = min(as.integer(k), nrow(x) - 1L), solver = "base", cache = FALSE)
}

test_that("real-data iris workflow remains stable on identity-style alignment", {
  x <- scale(as.matrix(iris[1:60, 1:4]))
  x <- as.matrix(x)

  source <- make_real_domain(x, k = 8L)
  target <- make_real_domain(x, k = 8L)

  fit <- fm_match(
    source = source,
    target = target,
    descriptors = list(source = x, target = x),
    penalties = list(descr = 1e-1, lap = 1e-3, comm = 1e-1),
    init = "identity",
    optimizer = "cg",
    cg_maxit = 10
  )

  p2p <- as_p2p(fit)
  inv <- fm_inverse(fit, method = "pseudoinverse")
  round_trip <- fm_compose(fit, inv)

  expect_true(is.finite(fit$diagnostics$total_objective))
  expect_gte(mean(p2p == seq_len(nrow(x))), 0.90)
  expect_identical(dim(round_trip$C), dim(fit$C))
  expect_true(all(is.finite(round_trip$C)))
})

test_that("real-data USArrests network sync exposes bounded convergence diagnostics", {
  base_x <- scale(as.matrix(USArrests))
  base_x <- as.matrix(base_x)
  n <- nrow(base_x)

  set.seed(1234)
  x1 <- base_x
  x2 <- base_x + matrix(rnorm(n * ncol(base_x), sd = 0.03), nrow = n)
  x3 <- base_x + matrix(rnorm(n * ncol(base_x), sd = 0.05), nrow = n)

  d1 <- make_real_domain(x1, k = 10L)
  d2 <- make_real_domain(x2, k = 10L)
  d3 <- make_real_domain(x3, k = 10L)

  net <- fm_network_match(
    domains = list(a = d1, b = d2, c = d3),
    descriptors = list(a = x1, b = x2, c = x3),
    edges = c("a->b", "b->c", "c->a"),
    directed = TRUE,
    optimizer = "cg",
    cg_maxit = 4
  )

  synced <- fm_sync(net, mode = "robust", nit = 6)

  expect_s3_class(synced, "fm_network_fit")
  expect_lte(synced$diagnostics$nit_done, 6L)
  expect_equal(nrow(synced$diagnostics$edge_weights), 3L)
  expect_true(is.finite(synced$diagnostics$pre_cycle_median))
  expect_true(is.finite(synced$diagnostics$post_cycle_median))

  for (nm in names(synced$transforms)) {
    X <- synced$transforms[[nm]]
    I <- diag(ncol(X))
    expect_lt(norm(t(X) %*% X - I, type = "F"), 1e-6)
  }
})
