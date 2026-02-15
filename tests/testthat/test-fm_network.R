make_domain_for_network <- function(n = 6, k = 3) {
  op <- diag(seq(1, n, by = 1))
  fm_basis(fm_domain_generic(n_samples = n, operator = op), k = k, solver = "base", cache = FALSE)
}

make_identity_fit <- function(source, target, k = 3) {
  structure(
    list(
      C = diag(1, nrow = k, ncol = k),
      k1 = k,
      k2 = k,
      source = source,
      target = target,
      penalties = list(),
      diagnostics = list(convergence = 0, message = "manual")
    ),
    class = "fm_fit"
  )
}

test_that("network supports directed edges with sparse weight storage", {
  d1 <- make_domain_for_network()
  d2 <- make_domain_for_network()
  d3 <- make_domain_for_network()

  net <- fm_network(list(a = d1, b = d2, c = d3), directed = TRUE)

  fit_ab <- make_identity_fit(d1, d2)
  fit_bc <- make_identity_fit(d2, d3)
  fit_ca <- make_identity_fit(d3, d1)

  net <- fm_network_add_map(net, "a", "b", fit_ab, weight = 0.5)
  net <- fm_network_add_map(net, "b", "c", fit_bc, weight = 0.7)
  net <- fm_network_add_map(net, "c", "a", fit_ca, weight = 0.9)

  expect_true(is_fm_network(net))
  expect_true(inherits(net$weights, "dgCMatrix"))
  expect_equal(length(net$maps), 3)
  expect_equal(as.numeric(net$weights[net$index[["a"]], net$index[["b"]]]), 0.5)
})

test_that("network cycle extraction and error scoring work", {
  d1 <- make_domain_for_network()
  d2 <- make_domain_for_network()
  d3 <- make_domain_for_network()

  net <- fm_network(list(a = d1, b = d2, c = d3), directed = TRUE)
  net <- fm_network_add_map(net, "a", "b", make_identity_fit(d1, d2))
  net <- fm_network_add_map(net, "b", "c", make_identity_fit(d2, d3))
  net <- fm_network_add_map(net, "c", "a", make_identity_fit(d3, d1))

  cycles <- fm_network_cycles(net)
  scores <- fm_network_cycle_scores(net)

  expect_true(length(cycles) >= 1)
  expect_equal(nrow(scores), length(cycles))
  expect_true(all(scores$error < 1e-8, na.rm = TRUE))
})

test_that("undirected mode adds reverse edges automatically", {
  d1 <- make_domain_for_network()
  d2 <- make_domain_for_network()

  net <- fm_network(list(a = d1, b = d2), directed = FALSE)
  net <- fm_network_add_map(net, "a", "b", diag(1, 3, 3), weight = 1)

  expect_equal(length(net$maps), 2)
  expect_false(is.null(fm_network_get_map(net, "a", "b")))
  expect_false(is.null(fm_network_get_map(net, "b", "a")))
  expect_equal(as.numeric(net$weights[net$index[["a"]], net$index[["b"]]]), 1)
  expect_equal(as.numeric(net$weights[net$index[["b"]], net$index[["a"]]]), 1)
})

test_that("network ingests fm_match outputs without reshaping", {
  n <- 6
  op <- diag(seq(1, n, by = 1))

  source <- fm_basis(fm_domain_generic(n_samples = n, operator = op), k = 3, solver = "base", cache = FALSE)
  target <- fm_basis(fm_domain_generic(n_samples = n, operator = op), k = 3, solver = "base", cache = FALSE)

  desc <- cbind(seq_len(n), seq_len(n)^2)

  fit <- fm_match(
    source = source,
    target = target,
    descriptors = list(source = desc, target = desc),
    penalties = list(descr = 1, lap = 0, comm = 0),
    init = "identity",
    maxit = 20
  )

  net <- fm_network(list(src = source, tgt = target), directed = TRUE)
  net <- fm_network_add_map(net, "src", "tgt", fit)

  stored <- fm_network_get_map(net, "src", "tgt")
  expect_s3_class(stored, "fm_fit")
  expect_identical(dim(stored$C), dim(fit$C))
})

test_that("fm_network_match supports pair batching and descriptor batching", {
  n <- 12
  p <- 16
  k <- 4
  op <- diag(seq(1, n, by = 1))

  d1 <- fm_basis(fm_domain_generic(n_samples = n, operator = op), k = k, solver = "base", cache = FALSE)
  d2 <- fm_basis(fm_domain_generic(n_samples = n, operator = op), k = k, solver = "base", cache = FALSE)
  d3 <- fm_basis(fm_domain_generic(n_samples = n, operator = op), k = k, solver = "base", cache = FALSE)

  set.seed(333)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  desc <- list(
    a = x,
    b = x + matrix(rnorm(n * p, sd = 0.03), nrow = n, ncol = p),
    c = x + matrix(rnorm(n * p, sd = 0.04), nrow = n, ncol = p)
  )

  net <- fm_network_match(
    domains = list(a = d1, b = d2, c = d3),
    descriptors = desc,
    edges = c("a->b", "b->c", "c->a"),
    directed = TRUE,
    pair_batch_size = 2,
    descriptor_batch_size = 4,
    optimizer = "cg",
    kernel_backend = "r",
    cg_maxit = 3
  )

  expect_s3_class(net, "fm_network")
  expect_equal(length(net$maps), 3)

  ab <- fm_network_get_map(net, "a", "b")
  expect_s3_class(ab, "fm_fit")
  expect_identical(ab$diagnostics$commutativity_mode, "stream")
  expect_identical(ab$diagnostics$descriptor_batch_size, 4L)
})

test_that("fm_network validates constructor and map keys", {
  d1 <- make_domain_for_network()
  d2 <- make_domain_for_network()

  expect_error(fm_network(list(a = d1)), "at least two")
  expect_error(fm_network(list(a = d1, b = list())), "fm_domain")
  expect_error(
    fm_network(list(a = d1, b = d2), maps = list(ab = diag(3))),
    "keyed as 'i->j'"
  )
})

test_that("fm_network normalizes edge specifications and validates failures", {
  d1 <- make_domain_for_network()
  d2 <- make_domain_for_network()
  d3 <- make_domain_for_network()

  set.seed(901)
  x <- matrix(rnorm(6 * 4), nrow = 6, ncol = 4)
  desc <- list(
    a = x,
    b = x + matrix(rnorm(6 * 4, sd = 0.01), nrow = 6, ncol = 4),
    c = x + matrix(rnorm(6 * 4, sd = 0.02), nrow = 6, ncol = 4)
  )

  # Character edge keys
  net_char <- fm_network_match(
    domains = list(a = d1, b = d2, c = d3),
    descriptors = desc,
    edges = c("a->b", "b->c"),
    directed = TRUE,
    optimizer = "cg",
    cg_maxit = 2
  )
  expect_equal(length(net_char$maps), 2)

  # Matrix/data-frame style edge specs
  edge_mat <- matrix(c("a", "b", "b", "c"), ncol = 2, byrow = TRUE)
  net_mat <- fm_network_match(
    domains = list(a = d1, b = d2, c = d3),
    descriptors = desc,
    edges = edge_mat,
    directed = TRUE,
    optimizer = "cg",
    cg_maxit = 2
  )
  expect_equal(length(net_mat$maps), 2)

  # Undirected duplicate removal path
  edge_df <- data.frame(i = c("a", "b"), j = c("b", "a"), stringsAsFactors = FALSE)
  net_undir <- fm_network_match(
    domains = list(a = d1, b = d2, c = d3),
    descriptors = desc,
    edges = edge_df,
    directed = FALSE,
    optimizer = "cg",
    cg_maxit = 2
  )
  expect_equal(length(net_undir$maps), 2)

  expect_error(
    fm_network_match(
      domains = list(a = d1, b = d2, c = d3),
      descriptors = desc,
      edges = c("bad-key"),
      optimizer = "cg",
      cg_maxit = 2
    ),
    "Invalid edge key"
  )
  expect_error(
    fm_network_match(
      domains = list(a = d1, b = d2, c = d3),
      descriptors = desc,
      edges = data.frame(i = "a", j = "a"),
      optimizer = "cg",
      cg_maxit = 2
    ),
    "Self-edges"
  )
  expect_error(
    fm_network_match(
      domains = list(a = d1, b = d2, c = d3),
      descriptors = desc,
      pair_batch_size = 0,
      optimizer = "cg",
      cg_maxit = 2
    ),
    "pair_batch_size"
  )
})

test_that("cycle diagnostics handle invalid cycle definitions and dimension mismatch", {
  d1 <- make_domain_for_network(n = 6, k = 3)
  d2 <- make_domain_for_network(n = 6, k = 2)
  d3 <- make_domain_for_network(n = 6, k = 4)

  net <- fm_network(list(a = d1, b = d2, c = d3), directed = TRUE)
  net <- fm_network_add_map(net, "a", "b", matrix(1, nrow = 2, ncol = 3))
  net <- fm_network_add_map(net, "b", "c", matrix(1, nrow = 4, ncol = 2))
  net <- fm_network_add_map(net, "c", "a", matrix(1, nrow = 3, ncol = 5))

  expect_error(fm_network_cycle_error(net, c("a", "b")), "length 3")
  expect_true(is.na(fm_network_cycle_error(net, c("a", "b", "c"))))
  expect_error(fm_network_get_map(list(), "a", "b"), "fm_network")
})
