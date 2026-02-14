make_noisy_network_for_consensus <- function(k = 4, noise_sd = 0.35, seed = 123) {
  set.seed(seed)

  mk_dom <- function() {
    op <- diag(seq(1, 8, by = 1))
    fm_basis(fm_domain_generic(n_samples = 8, operator = op), k = k, solver = "base", cache = FALSE)
  }

  d1 <- mk_dom(); d2 <- mk_dom(); d3 <- mk_dom()

  random_orth <- function(k) qr.Q(qr(matrix(rnorm(k * k), k, k)))

  X1 <- diag(k)
  X2 <- random_orth(k)
  X3 <- random_orth(k)

  C12 <- X2 %*% t(X1) + matrix(rnorm(k * k, sd = noise_sd), k, k)
  C23 <- X3 %*% t(X2) + matrix(rnorm(k * k, sd = noise_sd), k, k)
  C31 <- X1 %*% t(X3) + matrix(rnorm(k * k, sd = noise_sd), k, k)

  net <- fm_network(list(a = d1, b = d2, c = d3), directed = TRUE)
  net <- fm_network_add_map(net, "a", "b", C12)
  net <- fm_network_add_map(net, "b", "c", C23)
  net <- fm_network_add_map(net, "c", "a", C31)
  net
}

test_that("fm_consensus returns latent basis and alignments", {
  net <- make_noisy_network_for_consensus(noise_sd = 0.2)

  cons <- fm_consensus(net)

  expect_s3_class(cons, "fm_consensus")
  expect_true(is.matrix(cons$latent_basis))
  expect_true(length(cons$alignments) == 3)
  expect_true(length(cons$maps) >= 6)
  expect_true(is.list(cons$uncertainty))
  expect_true(all(c("edge_confidence", "summary") %in% names(cons$uncertainty)))
  expect_true(is.finite(cons$uncertainty$summary$median))
  expect_identical(cons$diagnostics$method, "latent_clb")
  expect_true(is.list(cons$diagnostics$latent))
  expect_true(cons$diagnostics$latent$method == "latent_clb")
})

test_that("latent consensus path is not worse than sync-transform baseline", {
  net <- make_noisy_network_for_consensus(noise_sd = 0.35, seed = 404)

  cons_sync <- fm_consensus(net, method = "sync_transforms")
  cons_latent <- fm_consensus(net, method = "latent_clb", latent_nit = 8)

  med_sync <- cons_sync$diagnostics$consensus_cycle_median
  med_latent <- cons_latent$diagnostics$consensus_cycle_median

  expect_true(is.finite(med_sync))
  expect_true(is.finite(med_latent))
  expect_lte(med_latent, med_sync + 1e-10)
  expect_true(!is.null(cons_latent$diagnostics$latent$solver))
})

test_that("consensus maps improve consistency versus unsynchronized baseline", {
  net <- make_noisy_network_for_consensus(noise_sd = 0.35, seed = 55)

  pre <- fm_network_cycle_scores(net)
  pre_med <- median(pre$error, na.rm = TRUE)

  cons <- fm_consensus(net)
  post <- fm_network_cycle_scores(cons$network)
  post_med <- median(post$error, na.rm = TRUE)

  expect_true(is.finite(pre_med))
  expect_true(is.finite(post_med))
  expect_lte(post_med, pre_med)
})

test_that("consensus object round-trips through serialization", {
  net <- make_noisy_network_for_consensus(noise_sd = 0.2)
  cons <- fm_consensus(net)

  path <- tempfile(fileext = ".rds")
  fm_consensus_save(cons, path)
  loaded <- fm_consensus_load(path)

  expect_s3_class(loaded, "fm_consensus")
  expect_equal(names(loaded$alignments), names(cons$alignments))
  expect_equal(dim(loaded$latent_basis), dim(cons$latent_basis))
})

test_that("confidence-aware map selection filters by threshold", {
  net <- make_noisy_network_for_consensus(noise_sd = 0.35, seed = 88)
  cons <- fm_consensus(net)

  selected_lo <- fm_consensus_select_maps(cons, min_confidence = 0.2)
  selected_hi <- fm_consensus_select_maps(cons, min_confidence = 0.8)

  expect_true(length(selected_lo) >= length(selected_hi))
  if (length(selected_hi) > 0L) {
    conf <- cons$uncertainty$edge_confidence
    keep <- conf$key %in% names(selected_hi)
    expect_true(all(conf$confidence[keep] >= 0.8))
  }
})

test_that("confidence summaries decrease as synthetic noise increases", {
  net_low <- make_noisy_network_for_consensus(noise_sd = 0.10, seed = 321)
  net_high <- make_noisy_network_for_consensus(noise_sd = 0.60, seed = 321)

  cons_low <- fm_consensus(net_low)
  cons_high <- fm_consensus(net_high)

  expect_gte(cons_low$uncertainty$summary$median, cons_high$uncertainty$summary$median)
})
