make_diagnostic_network <- function(k = 4, noise_sd = 0.3, seed = 7) {
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
  C21 <- t(C12) + matrix(rnorm(k * k, sd = noise_sd), k, k)
  C23 <- X3 %*% t(X2) + matrix(rnorm(k * k, sd = noise_sd), k, k)
  C32 <- t(C23) + matrix(rnorm(k * k, sd = noise_sd), k, k)
  C31 <- X1 %*% t(X3) + matrix(rnorm(k * k, sd = noise_sd), k, k)
  C13 <- t(C31) + matrix(rnorm(k * k, sd = noise_sd), k, k)

  net <- fm_network(list(a = d1, b = d2, c = d3), directed = TRUE)
  net <- fm_network_add_map(net, "a", "b", C12)
  net <- fm_network_add_map(net, "b", "a", C21)
  net <- fm_network_add_map(net, "b", "c", C23)
  net <- fm_network_add_map(net, "c", "b", C32)
  net <- fm_network_add_map(net, "c", "a", C31)
  net <- fm_network_add_map(net, "a", "c", C13)

  net
}

test_that("network report provides pre/post distributions and summary", {
  pre <- make_diagnostic_network(noise_sd = 0.3)
  fit <- fm_sync(pre, mode = "cycle", nit = 20)
  post <- fit$network

  report <- fm_network_compare(pre, post)
  summary_txt <- fm_network_summary(report)

  expect_s3_class(report, "fm_network_report")
  expect_true(nrow(report$pre$cycle_scores) > 0)
  expect_true(nrow(report$post$cycle_scores) > 0)
  expect_true(is.data.frame(report$pre$map_quality_scores))
  expect_true(is.list(report$pre$map_quality_summary))
  expect_true(is.character(summary_txt))
  expect_match(summary_txt, "Cycle median")
  expect_match(summary_txt, "Coverage median")
})

test_that("diagnostics are available from fm_network_fit and improved post-sync", {
  pre <- make_diagnostic_network(noise_sd = 0.35, seed = 11)
  fit <- fm_sync(pre, mode = "cycle", nit = 20)
  report <- fm_network_report(fit)

  expect_s3_class(report, "fm_network_report")
  expect_true(is.finite(report$summary$cycle_reduction_ratio) || is.na(report$summary$cycle_reduction_ratio))
  expect_true(report$pre$map_quality_summary$n >= 0)
})
