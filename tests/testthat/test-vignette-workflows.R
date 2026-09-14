make_vignette_path_adj <- function(n) {
  A <- matrix(0, nrow = n, ncol = n)
  A[cbind(1:(n - 1), 2:n)] <- 1
  A[cbind(2:n, 1:(n - 1))] <- 1
  A
}

make_vignette_pairwise_toy <- function(n = 24, shift = 2, k = 10, noise_sd = 0.005, seed = 7) {
  set.seed(seed)
  A <- make_vignette_path_adj(n)
  op <- diag(rowSums(A)) - A
  theta <- seq(0, 1, length.out = n)

  source_data <- cbind(
    theta,
    sin(2 * pi * theta),
    cos(2 * pi * theta),
    sin(4 * pi * theta),
    cos(6 * pi * theta)
  )

  truth <- c((shift + 1):n, seq_len(shift))
  target_data <- source_data[truth, , drop = FALSE] +
    matrix(rnorm(n * ncol(source_data), sd = noise_sd), nrow = n)

  source <- fm_basis(
    fm_domain_generic(
      n_samples = n,
      data = source_data,
      operator = op,
      adjacency = A
    ),
    k = k,
    solver = "base",
    cache = FALSE
  )

  target <- fm_basis(
    fm_domain_generic(
      n_samples = n,
      data = target_data,
      operator = op[truth, truth],
      adjacency = A[truth, truth]
    ),
    k = k,
    solver = "base",
    cache = FALSE
  )

  list(source = source, target = target, truth = truth)
}

make_vignette_nonshape_toy <- function(n = 30, p = 12, k = 8, seed = 123) {
  set.seed(seed)
  t <- seq(0, 2 * pi, length.out = n)
  latent <- cbind(sin(t), cos(t), sin(2 * t), cos(2 * t))

  make_view <- function(latent_x, weight, noise_sd = 0.02) {
    scale(latent_x %*% weight + matrix(rnorm(nrow(latent_x) * ncol(weight), sd = noise_sd), nrow(latent_x)))
  }

  make_multidomain <- function(X) {
    op <- fm_operator(X, method = "knn_laplacian", k = 7)
    fm_basis(
      fm_domain_generic(n_samples = nrow(X), data = X, operator = op),
      k = k,
      solver = "base",
      cache = FALSE
    )
  }

  W1 <- matrix(rnorm(ncol(latent) * p), ncol = p)
  W2 <- W1 + 0.05 * matrix(rnorm(ncol(latent) * p), ncol = p)
  truth <- c(4:n, 1:3)
  X_A <- make_view(latent, W1)
  X_B <- make_view(latent[truth, , drop = FALSE], W2)

  list(
    source = make_multidomain(X_A),
    target = make_multidomain(X_B),
    desc_source = X_A,
    desc_target = X_B,
    truth = truth,
    shared_signal = latent[, 1] + 0.4 * latent[, 3]
  )
}

make_vignette_noisy_cycle_network <- function(k = 4, noise_sd = 0.35, seed = 42) {
  set.seed(seed)

  make_domain <- function() {
    op <- diag(seq_len(8))
    fm_basis(
      fm_domain_generic(n_samples = 8, operator = op),
      k = k,
      solver = "base",
      cache = FALSE
    )
  }

  random_orth <- function(k_dim) qr.Q(qr(matrix(rnorm(k_dim * k_dim), nrow = k_dim, ncol = k_dim)))

  X1 <- diag(k)
  X2 <- random_orth(k)
  X3 <- random_orth(k)

  net <- fm_network(list(a = make_domain(), b = make_domain(), c = make_domain()), directed = TRUE)
  net <- fm_network_add_map(net, "a", "b", X2 %*% t(X1) + matrix(rnorm(k * k, sd = noise_sd), nrow = k))
  net <- fm_network_add_map(net, "b", "c", X3 %*% t(X2) + matrix(rnorm(k * k, sd = noise_sd), nrow = k))
  net <- fm_network_add_map(net, "c", "a", X1 %*% t(X3) + matrix(rnorm(k * k, sd = noise_sd), nrow = k))
  net
}

test_that("shape vignette example recovers the posed-mesh correspondence", {
  pair <- fm_example_mesh_pair()
  op <- fm_operator_mesh(
    pair$source$vertices,
    pair$source$faces,
    method = "face_graph",
    weight_mode = "binary"
  )

  source <- fm_basis(
    fm_domain_mesh(pair$source$vertices, pair$source$faces, operator = op),
    k = 24,
    solver = "rspectra",
    cache = FALSE
  )
  target <- fm_basis(
    fm_domain_mesh(pair$target$vertices, pair$target$faces, operator = op),
    k = 24,
    solver = "rspectra",
    cache = FALSE
  )

  desc_source <- cbind(
    scale(pair$source$descriptors),
    fm_descriptor_hks(source, n_times = 8, landmarks = pair$landmarks)
  )
  desc_target <- cbind(
    scale(pair$target$descriptors),
    fm_descriptor_hks(target, n_times = 8, landmarks = pair$landmarks)
  )

  fit <- fm_match(
    source,
    target,
    descriptors = list(source = desc_source, target = desc_target),
    penalties = list(descr = 1, lap = 1e-4, comm = 0.01),
    init = "identity",
    optimizer = "cg",
    cg_maxit = 8
  )
  fit <- fm_refine(fit, method = "icp", nit = 4)

  landmark_vertex <- pair$source$vertices[pair$landmarks[[3]], ]
  signal_delta <- sweep(pair$source$vertices, 2, landmark_vertex, "-")
  signal <- exp(-rowSums(signal_delta^2) / (2 * 0.18^2))
  signal_hat <- fm_transfer(fit, signal)
  metrics <- fm_fit_metrics(fit, truth = pair$truth)

  expect_identical(mean(as_p2p(fit) == pair$truth), 1)
  expect_identical(metrics$coverage$ratio, 1)
  expect_identical(metrics$geodesic$normalized_mean, 0)
  expect_lte(mean((signal_hat - signal[pair$truth])^2), 0.01)
})

test_that("shape vignette real TOSCA example clears a smooth-transfer quality floor", {
  skip_on_cran()
  skip_if_not_installed("RSpectra")

  pair <- fm_example_mesh_pair("tosca_cat10")
  op_source <- fm_operator_mesh(
    pair$source$vertices,
    pair$source$faces,
    method = "cotangent"
  )
  op_target <- fm_operator_mesh(
    pair$target$vertices,
    pair$target$faces,
    method = "cotangent"
  )

  source <- fm_basis(
    fm_domain_mesh(pair$source$vertices, pair$source$faces, operator = op_source),
    k = 30,
    solver = "rspectra",
    cache = FALSE
  )
  target <- fm_basis(
    fm_domain_mesh(pair$target$vertices, pair$target$faces, operator = op_target),
    k = 30,
    solver = "rspectra",
    cache = FALSE
  )

  desc_source <- fm_descriptor_hks(source, n_times = 6, landmarks = pair$landmarks)
  desc_target <- fm_descriptor_hks(target, n_times = 6, landmarks = pair$landmarks)

  fit <- fm_match(
    source,
    target,
    descriptors = list(source = desc_source, target = desc_target),
    penalties = list(descr = 1, lap = 1e-4, comm = 0.02),
    init = "identity",
    optimizer = "cg",
    cg_maxit = 300,
    cg_tol = 1e-8
  )
  expect_identical(fit$diagnostics$convergence, 0L)
  fit <- fm_refine(fit, method = "icp", nit = 1)

  source_center <- pair$source$vertices[pair$landmarks[[5]], ]
  source_signal <- exp(-rowSums(sweep(pair$source$vertices, 2, source_center, "-")^2) / (2 * 20^2))
  target_signal <- source_signal[pair$truth]
  signal_hat <- fm_transfer(fit, source_signal)

  overlap <- function(frac) {
    top_n <- max(50L, floor(length(target_signal) * frac))
    length(intersect(
      order(signal_hat, decreasing = TRUE)[seq_len(top_n)],
      order(target_signal, decreasing = TRUE)[seq_len(top_n)]
    )) / top_n
  }

  expect_gte(stats::cor(signal_hat, target_signal), 0.95)
  expect_gte(overlap(0.05), 0.95)
  expect_gte(overlap(0.10), 0.90)
  expect_lte(mean((signal_hat - target_signal)^2), 0.005)
})

test_that("pairwise vignette toy reports real ground-truth metrics", {
  toy <- make_vignette_pairwise_toy()
  desc_source <- cbind(scale(toy$source$data), fm_descriptor_hks(toy$source, n_times = 6))
  desc_target <- cbind(scale(toy$target$data), fm_descriptor_hks(toy$target, n_times = 6))

  fit <- fm_match(
    toy$source,
    toy$target,
    descriptors = list(source = desc_source, target = desc_target),
    optimizer = "cg",
    cg_maxit = 6
  )
  fit <- fm_refine(fit, method = "icp", nit = 5)

  p2p <- as_p2p(fit)
  metrics <- fm_fit_metrics(fit, truth = toy$truth)

  expect_gte(mean(p2p == toy$truth), 0.95)
  expect_gte(metrics$coverage$ratio, 0.95)
  expect_lte(metrics$geodesic$normalized_mean, 1e-6)
  expect_true(is.list(metrics$continuity))
  expect_true(is.finite(metrics$continuity$normalized_mean))
})

test_that("network vignette example shows visible sync improvement", {
  net <- make_vignette_noisy_cycle_network()
  synced <- fm_sync(net, mode = "cycle", nit = 10)
  consensus <- fm_consensus(synced, method = "latent_clb", latent_nit = 8)

  expect_gt(synced$diagnostics$pre_cycle_median, 0)
  expect_lt(synced$diagnostics$post_cycle_median, synced$diagnostics$pre_cycle_median)
  expect_true(is.finite(consensus$diagnostics$consensus_cycle_median))
  expect_lte(consensus$diagnostics$consensus_cycle_median, synced$diagnostics$post_cycle_median + 1e-8)
})

test_that("non-shape vignette example beats a trivial transfer baseline", {
  toy <- make_vignette_nonshape_toy()

  fit <- fm_match(
    toy$source,
    toy$target,
    descriptors = list(source = toy$desc_source, target = toy$desc_target),
    penalties = list(descr = 1, lap = 1e-3, comm = 0.2),
    optimizer = "cg",
    cg_maxit = 6
  )
  fit <- fm_refine(fit, method = "icp", nit = 3)

  target_signal <- toy$shared_signal[toy$truth]
  signal_hat <- fm_transfer(fit, toy$shared_signal)
  transfer_mse <- mean((signal_hat - target_signal)^2)
  baseline_mse <- mean((mean(target_signal) - target_signal)^2)
  metrics <- fm_fit_metrics(fit, truth = toy$truth)

  expect_lte(transfer_mse, baseline_mse)
  expect_lte(metrics$geodesic$normalized_mean, 0.01)
  expect_gte(mean(as_p2p(fit) == toy$truth), 0.9)
})
