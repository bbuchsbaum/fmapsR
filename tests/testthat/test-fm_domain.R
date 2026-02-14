test_that("mesh adapter creates fm_domain with required contract fields", {
  vertices <- matrix(runif(15), ncol = 3)
  faces <- matrix(c(1, 2, 3), ncol = 3)

  domain <- fm_domain_mesh(vertices, faces = faces)

  expect_true(is_fm_domain(domain))
  expect_s3_class(domain, "fm_domain_mesh")
  expect_identical(fm_n_samples(domain), 5L)
  expect_true(!is.null(domain$measure))
  expect_true(!is.null(domain$projector))
  expect_true(!is.null(domain$unprojector))
})

test_that("pointcloud and graph adapters produce consistent sample counts", {
  points <- matrix(runif(20), ncol = 4)
  adjacency <- matrix(0, nrow = 6, ncol = 6)

  pc <- fm_domain_pointcloud(points)
  gr <- fm_domain_graph(adjacency)

  expect_identical(fm_n_samples(pc), 5L)
  expect_identical(fm_n_samples(gr), 6L)
})

test_that("generic adapter accepts operator and basis metadata", {
  op <- diag(4)
  basis <- list(vectors = diag(4), values = 1:4)

  domain <- fm_domain_generic(
    n_samples = 4,
    operator = op,
    basis = basis,
    metadata = list(source = "test")
  )

  expect_identical(dim(domain$operator), c(4L, 4L))
  expect_identical(domain$basis$k, 4L)
  expect_identical(domain$metadata$source, "test")
})

test_that("domain projection helpers work with provided basis", {
  basis <- list(vectors = diag(3), values = c(1, 2, 3))
  domain <- fm_domain_generic(n_samples = 3, basis = basis)

  x <- matrix(c(1, 2, 3), ncol = 1)
  coef <- fm_project(domain, x)
  recon <- fm_unproject(domain, coef)

  expect_equal(as.numeric(coef), c(1, 2, 3))
  expect_equal(as.numeric(recon), c(1, 2, 3))
})

test_that("new adapter types can be built without changing core constructor", {
  tensor_domain <- fm_new_domain(type = "tensor", n_samples = 3)

  expect_true(is_fm_domain(tensor_domain))
  expect_s3_class(tensor_domain, "fm_domain_tensor")
  expect_identical(tensor_domain$type, "tensor")
})
