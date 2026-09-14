test_that("mesh operator supports native backends", {
  pair <- fm_example_mesh_pair()

  op_native <- fm_operator_mesh(
    pair$source$vertices,
    pair$source$faces,
    method = "face_graph",
    weight_mode = "binary"
  )

  expect_identical(dim(op_native), c(nrow(pair$source$vertices), nrow(pair$source$vertices)))
  expect_equal(as.matrix(op_native), t(as.matrix(op_native)), tolerance = 1e-10)

  op_knn <- fm_operator_mesh(
    pair$source$vertices,
    method = "knn_laplacian",
    nnk = 12
  )

  expect_identical(dim(op_knn), dim(op_native))
  expect_equal(as.matrix(op_knn), t(as.matrix(op_knn)), tolerance = 1e-10)
})

test_that("cotangent mesh operator is symmetric and propagates measure weights", {
  pair <- fm_example_mesh_pair()
  op_cot <- fm_operator_mesh(
    pair$source$vertices,
    pair$source$faces,
    method = "cotangent"
  )

  expect_identical(dim(op_cot), c(nrow(pair$source$vertices), nrow(pair$source$vertices)))
  expect_equal(as.matrix(op_cot), t(as.matrix(op_cot)), tolerance = 1e-8)
  expect_true(is.numeric(attr(op_cot, "fm_measure_weights", exact = TRUE)))

  domain <- fm_domain_mesh(pair$source$vertices, pair$source$faces, operator = op_cot)
  expect_equal(domain$measure, attr(op_cot, "fm_measure_weights", exact = TRUE))
})

test_that("example mesh pair includes topology, descriptors, and landmarks", {
  pair <- fm_example_mesh_pair()

  expect_true(is.list(pair$source))
  expect_true(is.matrix(pair$source$vertices))
  expect_true(is.matrix(pair$source$faces))
  expect_true(is.matrix(pair$source$descriptors))
  expect_identical(nrow(pair$source$vertices), nrow(pair$target$vertices))
  expect_identical(nrow(pair$source$descriptors), nrow(pair$source$vertices))
  expect_identical(length(pair$truth), nrow(pair$source$vertices))
  expect_true(all(pair$landmarks %in% pair$truth))
})

test_that("vendored TOSCA mesh pairs load real cat meshes and thumbnails", {
  pair <- fm_example_mesh_pair("tosca_cat2")

  expect_identical(nrow(pair$source$vertices), 27894L)
  expect_identical(nrow(pair$target$vertices), 27894L)
  expect_identical(nrow(pair$source$faces), 55712L)
  expect_identical(length(pair$truth), 27894L)
  expect_true(all(pair$landmarks %in% pair$truth))
  expect_true(file.exists(pair$metadata$source_thumbnail))
  expect_true(file.exists(pair$metadata$target_thumbnail))
  expect_length(pair$metadata$gallery_thumbnails, 5L)
})

test_that("mesh plotting helper returns plotting metadata and optional plotly widget", {
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  pair <- fm_example_mesh_pair()
  values <- pair$source$descriptors[, "dorsal_tag"]

  out_base <- fm_plot_mesh_scalar(pair$source, values, engine = "base", title = "base mesh")
  expect_identical(out_base$engine, "base")
  expect_identical(nrow(out_base$projected), nrow(pair$source$vertices))
  expect_identical(length(out_base$face_values), nrow(pair$source$faces))

  if (requireNamespace("plotly", quietly = TRUE)) {
    out_plotly <- fm_plot_mesh_scalar(pair$source, values, engine = "plotly", title = "plotly mesh")
    expect_true(inherits(out_plotly, "plotly"))
  }
})

test_that("cotangent mass and stiffness obey analytic triangle and scale contracts", {
  v <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0))
  f <- matrix(1:3, 1)
  expected <- rbind(c(1, -.5, -.5), c(-.5, .5, 0), c(-.5, 0, .5))
  cot <- mesh_cotangent_stiffness(v, f)
  expect_equal(as.matrix(cot$stiffness), expected)
  expect_equal(cot$mass, rep(1 / 6, 3))
  small <- mesh_cotangent_stiffness(v * 1e-5, f)
  expect_equal(as.matrix(small$stiffness), expected)
  expect_equal(small$mass / 1e-10, cot$mass)

  for (eps in c(1e-8, 1)) {
    op <- fm_operator_mesh(v, f, eps = eps)
    d <- fm_basis(fm_domain_mesh(v, f, operator = op), 3, solver = "base", cache = FALSE)
    phi <- d$basis$vectors
    expect_equal(crossprod(phi, phi * d$measure), diag(3), tolerance = 1e-10)
    expect_equal(expected %*% phi, (phi * d$measure) %*% diag(d$basis$values), tolerance = 1e-10)
    expect_equal(fm_unproject(d, fm_project(d, v)), v, tolerance = 1e-10)
  }
})

test_that("malformed and degenerate meshes fail before constructing an operator", {
  v <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0))
  f <- matrix(1:3, 1)
  expect_error(fm_operator_mesh(v, f + .2), "integer")
  expect_error(fm_operator_mesh(v, matrix(c(1, NA, 3), 1)), "finite")
  expect_error(fm_operator_mesh(v, matrix(c(1, 1, 3), 1)), "distinct")
  expect_error(fm_operator_mesh(v, matrix(integer(), 0, 3)), "integer")
  expect_error(fm_operator_mesh(v * NA, f), "finite")
  expect_error(fm_operator_mesh(cbind(v, 0), f), "exactly three")
  expect_error(fm_operator_mesh(v * 0, f), "non-degenerate")
  expect_error(fm_operator_mesh(rbind(v, c(1, 1, 0)), f), "every vertex")
})

test_that("mesh distance auto mode prefers graph distances when topology is available", {
  pair <- fm_example_mesh_pair()
  domain <- fm_domain_mesh(pair$source$vertices, pair$source$faces)

  d_auto <- fm_distance_matrix(domain, method = "auto")
  d_graph <- fm_distance_matrix(domain, method = "graph_shortest")

  expect_equal(d_auto, d_graph)
})

test_that("native mesh workflow supports a stable posed-mesh correspondence demo", {
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

  p2p <- as_p2p(fit)
  metrics <- fm_fit_metrics(fit, truth = pair$truth)

  expect_identical(mean(p2p == pair$truth), 1)
  expect_identical(metrics$coverage$ratio, 1)
  expect_identical(metrics$geodesic$normalized_mean, 0)
})
