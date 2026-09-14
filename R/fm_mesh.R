validate_mesh_vertices <- function(vertices) {
  vertices <- as.matrix(vertices)
  if (!is.numeric(vertices)) {
    stop("`vertices` must be numeric", call. = FALSE)
  }
  if (nrow(vertices) < 3L || ncol(vertices) != 3L || any(!is.finite(vertices))) {
    stop("`vertices` must have at least three rows and exactly three finite columns", call. = FALSE)
  }
  vertices
}

validate_mesh_faces <- function(faces, n_vertices) {
  if (is.null(faces)) {
    return(NULL)
  }

  faces <- as.matrix(faces)
  if (!is.numeric(faces) || ncol(faces) != 3L) {
    stop("`faces` must be a numeric triangle matrix with three columns", call. = FALSE)
  }

  if (nrow(faces) == 0L || any(!is.finite(faces)) || any(faces != floor(faces))) {
    stop("`faces` must contain finite integer triangle indices", call. = FALSE)
  }
  bad <- faces < 1L | faces > n_vertices
  if (any(bad)) {
    stop("`faces` indices must be between 1 and `nrow(vertices)`", call. = FALSE)
  }

  faces <- matrix(as.integer(faces), ncol = 3L)
  if (any(faces[, 1] == faces[, 2] | faces[, 1] == faces[, 3] | faces[, 2] == faces[, 3])) {
    stop("`faces` triangles must contain three distinct vertices", call. = FALSE)
  }

  faces
}

estimate_mesh_sigma <- function(vertices, faces = NULL, nnk = 8L) {
  vertices <- validate_mesh_vertices(vertices)
  n <- nrow(vertices)

  if (!is.null(faces)) {
    faces <- validate_mesh_faces(faces, n)
    edge_pairs <- mesh_edge_pairs(faces, n = n)
    edge_lengths <- sqrt(rowSums((vertices[edge_pairs[, 1], , drop = FALSE] - vertices[edge_pairs[, 2], , drop = FALSE])^2))
    edge_lengths <- edge_lengths[is.finite(edge_lengths) & edge_lengths > 0]
    if (length(edge_lengths) > 0L) {
      return(stats::median(edge_lengths))
    }
  }

  nnk <- max(1L, min(as.integer(nnk), n - 1L))
  row_norm <- rowSums(vertices * vertices)
  d2 <- outer(row_norm, row_norm, "+") - 2 * (vertices %*% t(vertices))
  d2[d2 < 0] <- 0
  diag(d2) <- Inf
  nn_idx <- apply(d2, 1, function(row) order(row)[seq_len(nnk)])
  if (!is.matrix(nn_idx)) {
    nn_idx <- matrix(nn_idx, nrow = nnk, ncol = n)
  }
  vals <- vapply(seq_len(n), function(i) d2[i, nn_idx[, i]], numeric(nnk))
  stats::median(sqrt(as.numeric(vals)))
}

mesh_face_weighted_adjacency <- function(vertices, faces, weight_mode = c("heat", "binary"), sigma = NULL) {
  weight_mode <- match.arg(weight_mode)
  vertices <- validate_mesh_vertices(vertices)
  faces <- validate_mesh_faces(faces, nrow(vertices))
  edges <- mesh_edge_pairs(faces, n = nrow(vertices))
  if (nrow(edges) == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(nrow(vertices), nrow(vertices))))
  }

  if (weight_mode == "binary") {
    return(Matrix::sparseMatrix(
      i = c(edges[, 1], edges[, 2]),
      j = c(edges[, 2], edges[, 1]),
      x = 1,
      dims = c(nrow(vertices), nrow(vertices))
    ))
  }

  sigma <- sigma %||% estimate_mesh_sigma(vertices, faces = faces)
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be a positive scalar", call. = FALSE)
  }

  dist2 <- rowSums((vertices[edges[, 1], , drop = FALSE] - vertices[edges[, 2], , drop = FALSE])^2)
  w <- exp(-dist2 / (2 * sigma^2))

  Matrix::sparseMatrix(
    i = c(edges[, 1], edges[, 2]),
    j = c(edges[, 2], edges[, 1]),
    x = c(w, w),
    dims = c(nrow(vertices), nrow(vertices))
  )
}

mesh_cotangent_terms <- function(vertices, faces) {
  vertices <- validate_mesh_vertices(vertices)
  faces <- validate_mesh_faces(faces, nrow(vertices))

  vi <- vertices[faces[, 1], , drop = FALSE]
  vj <- vertices[faces[, 2], , drop = FALSE]
  vk <- vertices[faces[, 3], , drop = FALSE]

  e_ij <- vj - vi
  e_ik <- vk - vi
  e_ji <- vi - vj
  e_jk <- vk - vj
  e_ki <- vi - vk
  e_kj <- vj - vk

  cross_i <- cbind(
    e_ij[, 2] * e_ik[, 3] - e_ij[, 3] * e_ik[, 2],
    e_ij[, 3] * e_ik[, 1] - e_ij[, 1] * e_ik[, 3],
    e_ij[, 1] * e_ik[, 2] - e_ij[, 2] * e_ik[, 1]
  )
  area2 <- sqrt(rowSums(cross_i^2))
  if (any(!is.finite(area2) | area2 <= 0)) {
    stop("Cotangent operators require non-degenerate triangles with finite area", call. = FALSE)
  }

  cot_i <- rowSums(e_ij * e_ik) / area2
  cot_j <- rowSums(e_ji * e_jk) / area2
  cot_k <- rowSums(e_ki * e_kj) / area2
  area <- 0.5 * area2

  list(
    area = area,
    edge_i = faces[, c(2, 3), drop = FALSE],
    edge_j = faces[, c(1, 3), drop = FALSE],
    edge_k = faces[, c(1, 2), drop = FALSE],
    cot_i = cot_i,
    cot_j = cot_j,
    cot_k = cot_k
  )
}

mesh_cotangent_stiffness <- function(vertices, faces) {
  terms <- mesh_cotangent_terms(vertices, faces)
  n <- nrow(vertices)

  edge_pairs <- rbind(terms$edge_i, terms$edge_j, terms$edge_k)
  weights <- 0.5 * c(terms$cot_i, terms$cot_j, terms$cot_k)

  edge_pairs <- cbind(
    pmin(edge_pairs[, 1], edge_pairs[, 2]),
    pmax(edge_pairs[, 1], edge_pairs[, 2])
  )

  key <- paste(edge_pairs[, 1], edge_pairs[, 2], sep = ":")
  w_sum <- rowsum(weights, group = key, reorder = FALSE)
  edge_unique <- edge_pairs[match(rownames(w_sum), key), , drop = FALSE]
  w_unique <- as.numeric(w_sum[, 1])

  W <- Matrix::sparseMatrix(
    i = c(edge_unique[, 1], edge_unique[, 2]),
    j = c(edge_unique[, 2], edge_unique[, 1]),
    x = c(w_unique, w_unique),
    dims = c(n, n)
  )

  D <- Matrix::Diagonal(x = Matrix::rowSums(W))
  L <- D - W

  mass_idx <- c(faces[, 1], faces[, 2], faces[, 3])
  mass_x <- rep(terms$area / 3, 3)
  mass_accum <- rowsum(mass_x, group = mass_idx, reorder = FALSE)
  mass <- numeric(n)
  mass[as.integer(rownames(mass_accum))] <- as.numeric(mass_accum[, 1])
  if (any(!is.finite(mass) | mass <= 0)) {
    stop("Cotangent operators require every vertex to belong to a triangle", call. = FALSE)
  }

  list(stiffness = L, mass = mass)
}

laplacian_from_adjacency <- function(adjacency, normalized = TRUE, eps = 1e-8) {
  A <- Matrix::Matrix(adjacency, sparse = TRUE)
  A <- (A + Matrix::t(A)) / 2
  Matrix::diag(A) <- 0

  degree <- Matrix::rowSums(A)
  D <- Matrix::Diagonal(x = degree)

  if (!isTRUE(normalized)) {
    return(D - A)
  }

  inv_sqrt <- 1 / sqrt(pmax(degree, eps))
  Dm12 <- Matrix::Diagonal(x = inv_sqrt)
  I <- Matrix::Diagonal(n = nrow(A), x = 1)
  I - Dm12 %*% A %*% Dm12
}

#' Build a Mesh Operator
#'
#' @param vertices Finite numeric vertex matrix with at least three rows and
#'   exactly three coordinate columns.
#' @param faces Optional triangle matrix with 1-based indices.
#' @param method Operator backend: `"auto"`, `"face_graph"`, `"cotangent"`,
#'   or `"knn_laplacian"`.
#' @param weight_mode Edge weighting for graph-based operators: `"heat"` or
#'   `"binary"`.
#' @param normalized Whether to normalize graph-based operators. Cotangent
#'   operators always use the mass normalization selected by `mass_mode`.
#' @param sigma Optional scale used by `"heat"` weighting.
#' @param nnk Number of neighbors used by `"knn_laplacian"` when needed.
#' @param eps Positive numerical stabilization floor for graph degrees or
#'   cotangent vertex masses. Stabilized masses are also used for projection.
#' @param mass_mode Mass normalization for `method = "cotangent"`: `"lumped"`
#'   or `"uniform"`.
#'
#' @return Symmetric operator matrix for use with [fm_basis()].
#' @export
fm_operator_mesh <- function(
  vertices,
  faces = NULL,
  method = c("auto", "face_graph", "cotangent", "knn_laplacian"),
  weight_mode = c("heat", "binary"),
  normalized = TRUE,
  sigma = NULL,
  nnk = 12,
  eps = 1e-8,
  mass_mode = c("lumped", "uniform")
) {
  vertices <- validate_mesh_vertices(vertices)
  faces <- validate_mesh_faces(faces, nrow(vertices))
  weight_mode <- match.arg(weight_mode)
  method <- match.arg(method)
  mass_mode <- match.arg(mass_mode)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a positive finite scalar", call. = FALSE)
  }

  if (method == "auto") {
    if (!is.null(faces)) {
      method <- "cotangent"
    } else {
      method <- "knn_laplacian"
    }
  }

  if (method == "cotangent") {
    if (is.null(faces)) {
      stop("`faces` are required for `method = \"cotangent\"`", call. = FALSE)
    }

    cot <- mesh_cotangent_stiffness(vertices, faces)
    mass <- if (identical(mass_mode, "uniform")) rep(1, nrow(vertices)) else cot$mass
    mass <- pmax(mass, eps)
    inv_sqrt_mass <- 1 / sqrt(mass)
    M_half <- Matrix::Diagonal(x = inv_sqrt_mass)
    op <- M_half %*% cot$stiffness %*% M_half
    attr(op, "fm_measure_weights") <- mass
    attr(op, "fm_basis_transform") <- "mass_mhalf"
    attr(op, "fm_basis_mass") <- mass
    attr(op, "fm_rspectra_preferred_which") <- "LM"
    attr(op, "fm_rspectra_sigma") <- 0
    attr(op, "fm_rspectra_ncv_min") <- 40L
    attr(op, "fm_rspectra_ncv_mult") <- 4L
    return(op)
  }

  if (method == "face_graph") {
    if (is.null(faces)) {
      stop("`faces` are required for `method = \"face_graph\"`", call. = FALSE)
    }
    A <- mesh_face_weighted_adjacency(vertices, faces, weight_mode = weight_mode, sigma = sigma)
    return(laplacian_from_adjacency(A, normalized = normalized, eps = eps))
  }

  fm_operator_knn_laplacian(
    x = vertices,
    k = nnk,
    sigma = sigma,
    normalized = normalized,
    eps = eps
  )
}

mesh_grid_faces <- function(n_length, n_ring) {
  idx <- function(i, j) {
    (i - 1L) * n_ring + ((j - 1L) %% n_ring) + 1L
  }

  faces <- matrix(0L, nrow = 2L * (n_length - 1L) * n_ring, ncol = 3L)
  cursor <- 1L
  for (i in seq_len(n_length - 1L)) {
    for (j in seq_len(n_ring)) {
      j2 <- if (j == n_ring) 1L else j + 1L
      a <- idx(i, j)
      b <- idx(i + 1L, j)
      c <- idx(i + 1L, j2)
      d <- idx(i, j2)
      faces[cursor, ] <- c(a, b, d)
      faces[cursor + 1L, ] <- c(b, c, d)
      cursor <- cursor + 2L
    }
  }
  faces
}

generate_pose_tube_mesh <- function(pose = c("source", "target"), n_length = 36L, n_ring = 18L) {
  pose <- match.arg(pose)
  n_length <- as.integer(max(8L, n_length))
  n_ring <- as.integer(max(6L, n_ring))

  s <- seq(0, 1, length.out = n_length)
  theta <- seq(0, 2 * pi, length.out = n_ring + 1L)[- (n_ring + 1L)]
  faces <- mesh_grid_faces(n_length, n_ring)

  bend <- if (pose == "source") 0.12 else 0.34
  sway <- if (pose == "source") 0.05 else -0.10
  twist <- if (pose == "source") 0.10 else 0.55
  muzzle <- if (pose == "source") 0.02 else 0.05

  center_x <- 2.6 * s - 1.3
  center_y <- 0.16 * sin(2 * pi * (s + sway))
  center_z <- bend * sin(pi * s)^2

  body_radius <- 0.05 +
    0.18 * exp(-((s - 0.38) / 0.18)^2) +
    0.08 * exp(-((s - 0.78) / 0.09)^2)
  flatten <- 0.65 + 0.20 * exp(-((s - 0.40) / 0.18)^2)
  head_profile <- exp(-((s - 0.83) / 0.06)^2)

  vertices <- matrix(0, nrow = n_length * n_ring, ncol = 3L)
  descriptors <- matrix(0, nrow = n_length * n_ring, ncol = 5L)
  colnames(descriptors) <- c("arc_length", "ring_cos", "ring_sin", "dorsal_tag", "cheek_tag")
  cursor <- 1L
  for (i in seq_along(s)) {
    for (j in seq_along(theta)) {
      ang <- theta[[j]] + twist * sin(pi * s[[i]])
      dorsal <- exp(-((s[[i]] - 0.52) / 0.08)^2) * exp(-(atan2(sin(ang - pi / 2), cos(ang - pi / 2)) / 0.35)^2)
      cheek <- exp(-((s[[i]] - 0.80) / 0.05)^2) * exp(-(atan2(sin(ang), cos(ang)) / 0.45)^2)
      belly <- exp(-((s[[i]] - 0.28) / 0.10)^2) * exp(-(atan2(sin(ang + pi / 2), cos(ang + pi / 2)) / 0.40)^2)
      radial_scale <- 1 + 0.30 * dorsal + 0.18 * cheek - 0.10 * belly
      r_y <- radial_scale * body_radius[[i]]
      r_z <- radial_scale * flatten[[i]] * body_radius[[i]]

      vertices[cursor, 1] <- center_x[[i]] + muzzle * cos(2 * ang) * head_profile[[i]] + 0.04 * cheek
      vertices[cursor, 2] <- center_y[[i]] + r_y * cos(ang)
      vertices[cursor, 3] <- center_z[[i]] + r_z * sin(ang) + 0.05 * dorsal
      descriptors[cursor, ] <- c(
        s[[i]],
        cos(theta[[j]]),
        sin(theta[[j]]),
        exp(-((s[[i]] - 0.52) / 0.08)^2) * exp(-(atan2(sin(theta[[j]] - pi / 2), cos(theta[[j]] - pi / 2)) / 0.35)^2),
        exp(-((s[[i]] - 0.80) / 0.05)^2) * exp(-(atan2(sin(theta[[j]]), cos(theta[[j]])) / 0.45)^2)
      )
      cursor <- cursor + 1L
    }
  }

  list(vertices = vertices, faces = faces, descriptors = descriptors)
}

package_extdata_file <- function(...) {
  path <- system.file(..., package = "fmapsR")
  if (nzchar(path)) {
    return(path)
  }

  fallback <- file.path("inst", ...)
  if (file.exists(fallback)) {
    return(normalizePath(fallback, mustWork = TRUE))
  }

  stop("Bundled example asset not found", call. = FALSE)
}

read_ascii_mesh_pair_member <- function(stem, dataset_dir = "tosca_cats") {
  vertices <- as.matrix(utils::read.table(package_extdata_file("extdata", dataset_dir, paste0(stem, ".vert"))))
  faces <- as.matrix(utils::read.table(package_extdata_file("extdata", dataset_dir, paste0(stem, ".tri"))))
  vertices <- validate_mesh_vertices(vertices)
  faces <- validate_mesh_faces(faces, nrow(vertices))

  descriptors <- scale(vertices)
  descriptors[!is.finite(descriptors)] <- 0
  colnames(descriptors) <- c("coord_x", "coord_y", "coord_z")

  list(
    vertices = vertices,
    faces = faces,
    descriptors = descriptors,
    thumbnail = package_extdata_file("extdata", dataset_dir, paste0(stem, ".png"))
  )
}

pick_mesh_landmarks <- function(vertices, n_landmarks = 6L) {
  vertices <- validate_mesh_vertices(vertices)
  n_landmarks <- as.integer(max(1L, n_landmarks))

  idx <- unique(c(
    which.max(vertices[, 1]),
    which.min(vertices[, 1]),
    which.max(vertices[, 2]),
    which.min(vertices[, 2]),
    which.max(vertices[, 3]),
    which.min(vertices[, 3])
  ))

  if (length(idx) < n_landmarks) {
    extra <- unique(round(seq(1, nrow(vertices), length.out = n_landmarks + 2L)))[-c(1, n_landmarks + 2L)]
    idx <- unique(c(idx, extra))
  }

  idx[seq_len(min(n_landmarks, length(idx)))]
}

tosca_gallery_paths <- function() {
  stems <- c("cat0", "cat1", "cat2", "cat6", "cat10")
  stats::setNames(
    vapply(stems, function(stem) package_extdata_file("extdata", "tosca_cats", paste0(stem, ".png")), character(1)),
    stems
  )
}

#' Load a Bundled Example Mesh Pair
#'
#' @param name Example mesh-pair name. Available options are `"posed_tube"`,
#'   `"tosca_cat1"`, `"tosca_cat2"`, `"tosca_cat6"`, and `"tosca_cat10"`.
#'
#' @return List containing `source`, `target`, `truth`, and `metadata`.
#' @export
fm_example_mesh_pair <- function(name = c("posed_tube", "tosca_cat1", "tosca_cat2", "tosca_cat6", "tosca_cat10")) {
  name <- match.arg(name)

  if (identical(name, "posed_tube")) {
    source <- generate_pose_tube_mesh("source")
    target <- generate_pose_tube_mesh("target")
    truth <- seq_len(nrow(source$vertices))
    n_ring <- 18L
    idx <- function(i, j) {
      (i - 1L) * n_ring + ((j - 1L) %% n_ring) + 1L
    }
    landmarks <- c(idx(4, 1), idx(12, 5), idx(19, 10), idx(28, 14), idx(33, 3), idx(35, 9))
    return(list(
      source = source,
      target = target,
      truth = truth,
      landmarks = landmarks,
      metadata = list(
        name = name,
        description = "A tiny posed tube mesh pair with shared topology and known correspondence."
      )
    ))
  }

  target_stem <- switch(
    name,
    tosca_cat1 = "cat1",
    tosca_cat2 = "cat2",
    tosca_cat6 = "cat6",
    tosca_cat10 = "cat10",
    NULL
  )

  if (!is.null(target_stem)) {
    source <- read_ascii_mesh_pair_member("cat0")
    target <- read_ascii_mesh_pair_member(target_stem)
    landmarks <- pick_mesh_landmarks(source$vertices, n_landmarks = 6L)
    gallery <- tosca_gallery_paths()

    return(list(
      source = source,
      target = target,
      truth = seq_len(nrow(source$vertices)),
      landmarks = landmarks,
      metadata = list(
        name = name,
        description = sprintf("A real TOSCA cat pair using cat0 as source and %s as target.", target_stem),
        dataset = "TOSCA",
        source_thumbnail = source$thumbnail,
        target_thumbnail = target$thumbnail,
        gallery_thumbnails = unname(gallery),
        archive_url = "https://web.archive.org/web/20240124102714if_/https://tosca.cs.technion.ac.il/data/toscahires-asci.zip",
        correspondence = "Vertex ordering provides the known target-to-source ground truth."
      )
    ))
  }

  stop("Unknown mesh example", call. = FALSE)
}

resolve_mesh_plot_input <- function(x, faces = NULL) {
  if (inherits(x, "fm_domain_mesh")) {
    return(list(vertices = validate_mesh_vertices(x$data$vertices), faces = validate_mesh_faces(x$data$faces, x$n_samples)))
  }

  if (is.list(x) && all(c("vertices", "faces") %in% names(x))) {
    verts <- validate_mesh_vertices(x$vertices)
    face_mat <- validate_mesh_faces(x$faces, nrow(verts))
    return(list(vertices = verts, faces = face_mat))
  }

  verts <- validate_mesh_vertices(x)
  face_mat <- validate_mesh_faces(faces, nrow(verts))
  list(vertices = verts, faces = face_mat)
}

mesh_view_projection <- function(vertices, view = c("isometric", "xy", "xz", "yz")) {
  view <- match.arg(view)
  verts <- validate_mesh_vertices(vertices)

  rot_x <- function(angle) {
    matrix(c(
      1, 0, 0,
      0, cos(angle), -sin(angle),
      0, sin(angle), cos(angle)
    ), nrow = 3L, byrow = TRUE)
  }

  rot_z <- function(angle) {
    matrix(c(
      cos(angle), -sin(angle), 0,
      sin(angle), cos(angle), 0,
      0, 0, 1
    ), nrow = 3L, byrow = TRUE)
  }

  rotated <- switch(
    view,
    isometric = verts %*% t(rot_x(pi / 7) %*% rot_z(pi / 5)),
    xy = verts,
    xz = verts[, c(1, 3, 2), drop = FALSE],
    yz = verts[, c(2, 3, 1), drop = FALSE]
  )

  rotated
}

mesh_palette <- function(n = 64L, palette = NULL) {
  if (is.function(palette)) {
    return(palette(n))
  }
  if (is.character(palette) && length(palette) == 1L) {
    return(grDevices::hcl.colors(n, palette = palette))
  }
  if (is.character(palette) && length(palette) >= 2L) {
    return(grDevices::colorRampPalette(palette)(n))
  }
  grDevices::hcl.colors(n, palette = "Viridis")
}

scale_to_colors <- function(values, palette = NULL, n = 64L, na_color = "#cccccc", value_range = NULL) {
  pal <- mesh_palette(n = n, palette = palette)
  vals <- as.numeric(values)
  rng <- if (is.null(value_range)) {
    range(vals[is.finite(vals)], na.rm = TRUE)
  } else {
    range(as.numeric(value_range), na.rm = TRUE)
  }
  if (!all(is.finite(rng)) || diff(rng) <= 0) {
    idx <- rep(ceiling(n / 2), length(vals))
  } else {
    idx <- 1L + floor((n - 1L) * (vals - rng[[1]]) / diff(rng))
    idx <- pmin(pmax(idx, 1L), n)
  }
  cols <- pal[idx]
  cols[!is.finite(vals)] <- na_color
  cols
}

shade_hex_colors <- function(cols, shading, na_color = "#cccccc") {
  if (length(cols) == 0L) {
    return(cols)
  }

  shading <- pmin(pmax(as.numeric(shading), 0), 1)
  rgb <- grDevices::col2rgb(cols) / 255
  shaded <- rgb * rep(shading, each = 3L)
  shaded[, !is.finite(shading)] <- grDevices::col2rgb(na_color) / 255
  grDevices::rgb(shaded[1, ], shaded[2, ], shaded[3, ])
}

mesh_face_shading <- function(vertices_proj, faces, light_dir = c(-0.35, -0.25, 1)) {
  tri1 <- vertices_proj[faces[, 1], , drop = FALSE]
  tri2 <- vertices_proj[faces[, 2], , drop = FALSE]
  tri3 <- vertices_proj[faces[, 3], , drop = FALSE]

  edge1 <- tri2 - tri1
  edge2 <- tri3 - tri1
  normals <- cbind(
    edge1[, 2] * edge2[, 3] - edge1[, 3] * edge2[, 2],
    edge1[, 3] * edge2[, 1] - edge1[, 1] * edge2[, 3],
    edge1[, 1] * edge2[, 2] - edge1[, 2] * edge2[, 1]
  )

  norm_len <- sqrt(rowSums(normals^2))
  normals <- sweep(normals, 1, pmax(norm_len, .Machine$double.eps), "/")

  light <- as.numeric(light_dir)
  light <- light / sqrt(sum(light^2))
  diffuse <- pmax(0, as.numeric(normals %*% light))

  0.45 + 0.55 * diffuse
}

plot_mesh_scalar_base <- function(
  vertices,
  faces,
  values,
  main = NULL,
  palette = NULL,
  show_edges = FALSE,
  view = c("isometric", "xy", "xz", "yz"),
  value_range = NULL,
  shading = TRUE
) {
  verts3 <- mesh_view_projection(vertices, view = match.arg(view))
  faces <- validate_mesh_faces(faces, nrow(verts3))
  values <- as.numeric(values)

  if (length(values) != nrow(verts3)) {
    stop("`values` must have length `nrow(vertices)`", call. = FALSE)
  }

  face_values <- rowMeans(matrix(values[faces], ncol = 3L))
  face_colors <- scale_to_colors(face_values, palette = palette, value_range = value_range)
  if (isTRUE(shading)) {
    face_colors <- shade_hex_colors(face_colors, mesh_face_shading(verts3, faces))
  }
  z_face <- rowMeans(matrix(verts3[faces, 3], ncol = 3L))
  ord <- order(z_face, decreasing = FALSE)

  xlim <- range(verts3[, 1])
  ylim <- range(verts3[, 2])
  xpad <- diff(xlim) * 0.05
  ypad <- diff(ylim) * 0.05
  xlim <- xlim + c(-xpad, xpad)
  ylim <- ylim + c(-ypad, ypad)

  graphics::plot(verts3[, 1], verts3[, 2],
    type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
    main = main, xlim = xlim, ylim = ylim
  )

  border <- if (isTRUE(show_edges)) "#1f2937" else NA
  for (ii in ord) {
    tri <- faces[ii, ]
    graphics::polygon(verts3[tri, 1], verts3[tri, 2], col = face_colors[[ii]], border = border)
  }

  invisible(list(
    engine = "base",
    projected = verts3,
    face_values = face_values,
    face_order = ord
  ))
}

#' Plot a Scalar Function on a Mesh
#'
#' @param x Mesh domain, mesh list with `vertices`/`faces`, or vertex matrix.
#' @param values Numeric vector with one scalar per vertex.
#' @param faces Optional face matrix when `x` is a plain vertex matrix.
#' @param engine Plot backend: `"auto"`, `"plotly"`, or `"base"`.
#' @param title Optional plot title.
#' @param palette Optional palette name, vector, or function.
#' @param show_edges Whether to draw triangle edges in the base backend.
#' @param view Projection used by the base backend.
#' @param value_range Optional two-value numeric range used to keep color scaling
#'   consistent across multiple panels.
#' @param shading Whether the base backend should apply simple directional
#'   shading to improve surface legibility.
#'
#' @return Plot object. `plotly` backend returns a plotly widget; `base`
#'   backend returns plotting metadata invisibly.
#' @export
fm_plot_mesh_scalar <- function(
  x,
  values,
  faces = NULL,
  engine = c("auto", "plotly", "base"),
  title = NULL,
  palette = "Viridis",
  show_edges = FALSE,
  view = c("isometric", "xy", "xz", "yz"),
  value_range = NULL,
  shading = TRUE
) {
  mesh <- resolve_mesh_plot_input(x, faces = faces)
  engine <- match.arg(engine)
  view <- match.arg(view)
  values <- as.numeric(values)

  if (length(values) != nrow(mesh$vertices)) {
    stop("`values` must have length `nrow(vertices)`", call. = FALSE)
  }

  if (engine == "auto") {
    engine <- if (requireNamespace("plotly", quietly = TRUE)) "plotly" else "base"
  }

  if (engine == "plotly") {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Install `plotly` to use `engine = \"plotly\"`", call. = FALSE)
    }

    if (is.null(mesh$faces)) {
      plt <- plotly::plot_ly(
      x = mesh$vertices[, 1],
      y = mesh$vertices[, 2],
      z = mesh$vertices[, 3],
      type = "scatter3d",
      mode = "markers",
      marker = list(
        color = values,
        colorscale = palette,
        size = 2,
        cmin = if (is.null(value_range)) NULL else min(value_range),
        cmax = if (is.null(value_range)) NULL else max(value_range)
      )
    )
    return(plotly::layout(plt, title = list(text = title %||% "")))
  }

  plt <- plotly::plot_ly(
      x = mesh$vertices[, 1],
      y = mesh$vertices[, 2],
      z = mesh$vertices[, 3],
      i = mesh$faces[, 1] - 1L,
      j = mesh$faces[, 2] - 1L,
      k = mesh$faces[, 3] - 1L,
      intensity = values,
      intensitymode = "vertex",
      type = "mesh3d",
      colorscale = palette,
      cmin = if (is.null(value_range)) NULL else min(value_range),
      cmax = if (is.null(value_range)) NULL else max(value_range),
      showscale = TRUE
    )
    return(plotly::layout(plt, title = list(text = title %||% "")))
  }

  plot_mesh_scalar_base(
    vertices = mesh$vertices,
    faces = mesh$faces,
    values = values,
    main = title,
    palette = palette,
    show_edges = show_edges,
    view = view,
    value_range = value_range,
    shading = shading
  )
}
