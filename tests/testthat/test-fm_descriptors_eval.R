make_descriptor_domain <- function(n = 12, k = 8) {
  op <- diag(seq_len(n))
  fm_basis(
    fm_domain_generic(n_samples = n, operator = op, data = matrix(rnorm(n * 3), nrow = n)),
    k = k,
    solver = "base",
    cache = FALSE
  )
}

test_that("HKS and WKS descriptor builders return finite matrices with expected shapes", {
  set.seed(1)
  d <- make_descriptor_domain(n = 10, k = 6)

  hks <- fm_descriptor_hks(d, n_times = 7)
  expect_equal(dim(hks), c(10, 7))
  expect_true(all(is.finite(hks)))

  hks_lm <- fm_descriptor_hks(d, n_times = 4, landmarks = c(1, 3))
  expect_equal(dim(hks_lm), c(10, 8))

  wks <- fm_descriptor_wks(d, n_energies = 9)
  expect_equal(dim(wks), c(10, 9))
  expect_true(all(is.finite(wks)))

  wks_lm <- fm_descriptor_wks(d, n_energies = 5, landmarks = c(2, 4, 6))
  expect_equal(dim(wks_lm), c(10, 15))
})

test_that("fm_descriptors dispatches descriptor methods", {
  d <- make_descriptor_domain(n = 8, k = 5)

  h <- fm_descriptors(d, method = "hks", n_times = 6)
  w <- fm_descriptors(d, method = "wks", n_energies = 6)

  expect_equal(dim(h), c(8, 6))
  expect_equal(dim(w), c(8, 6))
})

test_that("distance matrix helper supports graph and euclidean domains", {
  points <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
  d_pc <- fm_domain_pointcloud(points)
  D_pc <- fm_distance_matrix(d_pc, method = "auto")
  expect_equal(dim(D_pc), c(3, 3))
  expect_equal(unname(diag(D_pc)), c(0, 0, 0))

  A <- matrix(0, nrow = 4, ncol = 4)
  A[cbind(1:3, 2:4)] <- 1
  A[cbind(2:4, 1:3)] <- 1
  d_gr <- fm_domain_graph(A)
  D_gr <- fm_distance_matrix(d_gr, method = "graph_shortest")
  expect_equal(D_gr[1, 4], 3)
  expect_equal(D_gr[4, 1], 3)
})

test_that("evaluation utilities compute expected geodesic, continuity, and coverage summaries", {
  D <- as.matrix(stats::dist(matrix(seq_len(5), ncol = 1)))
  truth <- c(1L, 2L, 3L, 4L, 5L)

  geod_exact <- fm_eval_geodesic_error(truth, truth, D)
  expect_equal(geod_exact$mean, 0)
  expect_equal(geod_exact$normalized_mean, 0)

  geod_shift <- fm_eval_geodesic_error(c(2L, 3L, 4L, 5L, 5L), truth, D)
  expect_gt(geod_shift$mean, 0)

  A_tgt <- matrix(0, nrow = 5, ncol = 5)
  A_tgt[cbind(1:4, 2:5)] <- 1
  A_tgt[cbind(2:5, 1:4)] <- 1

  cont_smooth <- fm_eval_continuity(truth, D, A_tgt)
  cont_noisy <- fm_eval_continuity(c(1L, 5L, 1L, 5L, 1L), D, A_tgt)
  expect_gt(cont_noisy$mean, cont_smooth$mean)

  cov <- fm_eval_coverage(c(1L, 1L, 2L, 2L), n_source = 4)
  expect_equal(cov$n_support, 2L)
  expect_equal(cov$ratio, 0.5)

  cov_w <- fm_eval_coverage(c(1L, 1L, 2L, 2L), n_source = 4, weights = c(1, 1, 0, 0))
  expect_equal(cov_w$weighted_ratio, 1)
})

test_that("fm_fit_metrics returns coverage plus optional geodesic and continuity diagnostics", {
  set.seed(42)
  n <- 10

  A <- matrix(0, nrow = n, ncol = n)
  A[cbind(1:(n - 1), 2:n)] <- 1
  A[cbind(2:n, 1:(n - 1))] <- 1
  op <- diag(rowSums(A)) - A

  d1 <- fm_basis(fm_domain_graph(A, operator = op), k = 6, solver = "base", cache = FALSE)
  d2 <- fm_basis(fm_domain_graph(A, operator = op), k = 6, solver = "base", cache = FALSE)

  desc <- diag(n)
  fit <- fm_match(
    d1,
    d2,
    descriptors = list(source = desc, target = desc),
    optimizer = "cg",
    cg_maxit = 4,
    penalties = list(descr = 1, lap = 1e-3, comm = 0.2)
  )

  p2p <- as_p2p(fit)
  m <- fm_fit_metrics(fit, truth = p2p)

  expect_true(is.list(m$coverage))
  expect_true(is.list(m$geodesic))
  expect_true(is.list(m$continuity))
  expect_equal(m$geodesic$mean, 0)
  expect_true(m$coverage$ratio > 0)
})

test_that("normalized geodesic/continuity metrics are scale-invariant", {
  D <- as.matrix(stats::dist(matrix(seq_len(6), ncol = 1)))
  D2 <- 7 * D
  truth <- c(1L, 2L, 3L, 4L, 5L, 6L)
  p2p <- c(2L, 2L, 4L, 4L, 6L, 6L)

  g1 <- fm_eval_geodesic_error(p2p, truth, D, normalize = TRUE)
  g2 <- fm_eval_geodesic_error(p2p, truth, D2, normalize = TRUE)
  expect_equal(g2$mean, 7 * g1$mean, tolerance = 1e-10)
  expect_equal(g2$normalized_mean, g1$normalized_mean, tolerance = 1e-10)

  A_tgt <- matrix(0, nrow = 6, ncol = 6)
  A_tgt[cbind(1:5, 2:6)] <- 1
  A_tgt[cbind(2:6, 1:5)] <- 1

  c1 <- fm_eval_continuity(p2p, D, A_tgt, normalize = TRUE)
  c2 <- fm_eval_continuity(p2p, D2, A_tgt, normalize = TRUE)
  expect_equal(c2$mean, 7 * c1$mean, tolerance = 1e-10)
  expect_equal(c2$normalized_mean, c1$normalized_mean, tolerance = 1e-10)
})

test_that("coverage metrics obey bounds and entropy extremes", {
  n_source <- 10L

  p2p_uniform <- rep(seq_len(n_source), each = 3L)
  cov_uniform <- fm_eval_coverage(p2p_uniform, n_source = n_source)
  expect_equal(cov_uniform$ratio, 1)
  expect_equal(cov_uniform$weighted_ratio, 1)
  expect_true(cov_uniform$entropy >= 0 && cov_uniform$entropy <= 1)
  expect_gte(cov_uniform$entropy, 0.99)

  p2p_concentrated <- rep(1L, 30)
  cov_concentrated <- fm_eval_coverage(p2p_concentrated, n_source = n_source)
  expect_equal(cov_concentrated$ratio, 0.1)
  expect_equal(cov_concentrated$entropy, 0)

  w <- seq_len(n_source)
  cov_weighted <- fm_eval_coverage(c(1L, 1L, n_source, n_source), n_source = n_source, weights = w)
  expect_true(cov_weighted$weighted_ratio > 0)
  expect_true(cov_weighted$weighted_ratio < 1)
})
