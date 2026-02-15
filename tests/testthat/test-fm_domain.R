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

test_that("domain constructor validates core inputs and metadata", {
  expect_error(fm_new_domain(type = "", n_samples = 3), "non-empty string")
  expect_error(fm_new_domain(type = "generic", n_samples = 0), "positive scalar")
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, metadata = "bad"),
    "`metadata` must be a list"
  )
})

test_that("basis and measure validation errors are surfaced", {
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, basis = 1),
    "`basis` must be a list"
  )
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, basis = list(vectors = 1:3)),
    "must be a matrix"
  )
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, basis = list(vectors = diag(2))),
    "row count must match"
  )
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, basis = list(vectors = diag(3), values = "x")),
    "must be numeric"
  )
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, basis = list(vectors = diag(3), values = 1:2)),
    "length must match"
  )

  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, measure = 1:2),
    "must have length"
  )
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, measure = matrix(1, 2, 2)),
    "must be `n_samples x n_samples`"
  )
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, measure = list()),
    "must be numeric, matrix, sparse Matrix, or NULL"
  )
})

test_that("operator and adjacency validation checks dimensions and types", {
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, operator = 1:3),
    "`operator` must be a matrix or sparse Matrix"
  )
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, operator = matrix(1, 2, 2)),
    "`operator` must be `n_samples x n_samples`"
  )
  expect_error(
    fm_new_domain(type = "generic", n_samples = 3, adjacency = matrix(1, 2, 2)),
    "`adjacency` must be `n_samples x n_samples`"
  )
})

test_that("domain adapters validate matrix inputs", {
  expect_error(fm_domain_mesh(vertices = 1:9), "`vertices` must be a matrix")
  expect_error(
    fm_domain_mesh(vertices = matrix(runif(9), ncol = 3), faces = 1:3),
    "`faces` must be a matrix"
  )
  expect_error(fm_domain_pointcloud(points = 1:6), "`points` must be a matrix")
  expect_error(fm_domain_graph(adjacency = matrix(1, nrow = 2, ncol = 3)), "must be `n_samples x n_samples`")
})

test_that("generic domain constructor infers sample count from matrix payload", {
  x <- matrix(runif(12), nrow = 4)
  d <- fm_domain(data = x, type = "generic")
  expect_identical(d$n_samples, 4L)
})

test_that("projection helpers require valid domains and basis vectors", {
  d <- fm_domain_generic(n_samples = 3)
  expect_error(fm_n_samples(list()), "`domain` must inherit from `fm_domain`")
  expect_error(fm_measure_matrix(list()), "`domain` must inherit from `fm_domain`")
  expect_error(fm_project(list(), 1:3), "`domain` must inherit from `fm_domain`")
  expect_error(fm_unproject(list(), 1:3), "`domain` must inherit from `fm_domain`")
  expect_error(fm_project(d, 1:3), "No basis vectors available")
  expect_error(fm_unproject(d, 1:3), "No basis vectors available")
})

test_that("measure matrix and print methods handle matrix and numeric measures", {
  d_num <- fm_domain_generic(n_samples = 3, measure = c(1, 2, 3))
  M <- fm_measure_matrix(d_num)
  expect_true(inherits(M, "Matrix"))
  expect_equal(as.numeric(Matrix::diag(M)), c(1, 2, 3))

  d_mat <- fm_domain_generic(n_samples = 3, measure = diag(3))
  expect_equal(fm_measure_matrix(d_mat), diag(3))

  out <- capture.output(print(d_mat))
  expect_true(any(grepl("<fm_domain:generic>", out, fixed = TRUE)))
})
