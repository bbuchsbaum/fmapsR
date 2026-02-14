make_noisy_cycle_network <- function(k = 4, noise_sd = 0.35, seed = 99) {
  set.seed(seed)

  mk_dom <- function() {
    op <- diag(seq(1, 8, by = 1))
    fm_basis(fm_domain_generic(n_samples = 8, operator = op), k = k, solver = "base", cache = FALSE)
  }

  d1 <- mk_dom(); d2 <- mk_dom(); d3 <- mk_dom()

  random_orth <- function(k) {
    Q <- qr.Q(qr(matrix(rnorm(k * k), k, k)))
    as.matrix(Q)
  }

  X1 <- diag(k)
  X2 <- random_orth(k)
  X3 <- random_orth(k)

  C12 <- X2 %*% t(X1) + matrix(rnorm(k * k, sd = noise_sd), k, k)
  C23 <- X3 %*% t(X2) + matrix(rnorm(k * k, sd = noise_sd), k, k)
  C31 <- X1 %*% t(X3) + matrix(rnorm(k * k, sd = noise_sd), k, k)

  net <- fm_network(list(a = d1, b = d2, c = d3), directed = TRUE)
  net <- fm_network_add_map(net, "a", "b", C12, weight = 1)
  net <- fm_network_add_map(net, "b", "c", C23, weight = 1)
  net <- fm_network_add_map(net, "c", "a", C31, weight = 1)

  net
}

make_outlier_network <- function(k = 4, base_noise = 0.08, outlier_noise = 1.30, seed = 308) {
  set.seed(seed)

  mk_dom <- function() {
    op <- diag(seq(1, 10, by = 1))
    fm_basis(fm_domain_generic(n_samples = 10, operator = op), k = k, solver = "base", cache = FALSE)
  }

  random_orth <- function(k) {
    Q <- qr.Q(qr(matrix(rnorm(k * k), k, k)))
    as.matrix(Q)
  }

  X <- list(
    a = diag(k),
    b = random_orth(k),
    c = random_orth(k),
    d = random_orth(k)
  )

  domains <- list(a = mk_dom(), b = mk_dom(), c = mk_dom(), d = mk_dom())
  net <- fm_network(domains, directed = TRUE)
  nodes <- names(domains)
  truth <- list()

  for (i in nodes) {
    for (j in nodes) {
      if (identical(i, j)) next

      sd <- if (
        (identical(i, "a") && identical(j, "d")) ||
        (identical(i, "d") && identical(j, "a"))
      ) {
        outlier_noise
      } else {
        base_noise
      }

      Ctrue <- X[[j]] %*% t(X[[i]])
      Cij <- Ctrue + matrix(rnorm(k * k, sd = sd), k, k)
      truth[[paste0(i, "->", j)]] <- Ctrue
      net <- fm_network_add_map(net, i, j, Cij, weight = 1)
    }
  }

  list(network = net, truth = truth)
}

median_map_error_to_truth <- function(network, truth) {
  keys <- names(truth)
  errs <- vapply(keys, function(key) {
    ij <- strsplit(key, "->", fixed = TRUE)[[1]]
    map <- as.matrix(fm_network_get_map(network, ij[[1]], ij[[2]]))
    norm(map - truth[[key]], type = "F")
  }, numeric(1))
  stats::median(errs)
}

test_that("fm_sync supports adjacency, cycle, and robust modes with diagnostics", {
  net <- make_noisy_cycle_network(noise_sd = 0.2)

  out_adj <- fm_sync(net, mode = "adjacency", nit = 10)
  out_cyc <- fm_sync(net, mode = "cycle", nit = 10)
  out_rob <- fm_sync(net, mode = "robust", nit = 10)

  expect_s3_class(out_adj, "fm_network_fit")
  expect_s3_class(out_cyc, "fm_network_fit")
  expect_s3_class(out_rob, "fm_network_fit")
  expect_identical(out_adj$diagnostics$mode, "adjacency")
  expect_identical(out_cyc$diagnostics$mode, "cycle")
  expect_identical(out_rob$diagnostics$mode, "robust")
  expect_true(nrow(out_cyc$diagnostics$edge_weights) >= 3)
  expect_true(all(c("effective_weight", "robust_method") %in% names(out_rob$diagnostics$edge_weights)))
})

test_that("cycle-weighted sync improves cycle consistency on noisy benchmark", {
  net <- make_noisy_cycle_network(noise_sd = 0.35, seed = 17)

  pre <- fm_network_cycle_scores(net)
  pre_med <- median(pre$error, na.rm = TRUE)

  out <- fm_sync(net, mode = "cycle", nit = 20)
  post <- fm_network_cycle_scores(out$network)
  post_med <- median(post$error, na.rm = TRUE)

  reduction <- (pre_med - post_med) / pre_med

  expect_true(is.finite(reduction))
  expect_true(reduction >= 0.40)
  expect_equal(out$diagnostics$cycle_reduction_ratio, reduction, tolerance = 1e-8)
})

test_that("robust sync downweights outlier edges and improves over cycle mode", {
  obj <- make_outlier_network(k = 4, base_noise = 0.08, outlier_noise = 1.30, seed = 308)
  net <- obj$network
  truth <- obj$truth

  out_cycle <- fm_sync(net, mode = "cycle", nit = 20)
  out_robust <- fm_sync(net, mode = "robust", nit = 20)

  err_cycle <- median_map_error_to_truth(out_cycle$network, truth)
  err_robust <- median_map_error_to_truth(out_robust$network, truth)
  improvement <- (err_cycle - err_robust) / err_cycle

  expect_true(is.finite(improvement))
  expect_gte(improvement, 0.20)

  ew <- out_robust$diagnostics$edge_weights
  outlier_idx <- ew$key %in% c("a->d", "d->a")
  expect_true(all(ew$robust_weight[outlier_idx] < stats::median(ew$robust_weight[!outlier_idx])))
  expect_identical(unique(ew$robust_method), "tukey")
})
