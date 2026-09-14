# Build a Mesh Operator

Build a Mesh Operator

## Usage

``` r
fm_operator_mesh(
  vertices,
  faces = NULL,
  method = c("auto", "face_graph", "cotangent", "knn_laplacian"),
  weight_mode = c("heat", "binary"),
  normalized = TRUE,
  sigma = NULL,
  nnk = 12,
  eps = 1e-08,
  mass_mode = c("lumped", "uniform")
)
```

## Arguments

- vertices:

  Finite numeric vertex matrix with at least three rows and exactly
  three coordinate columns.

- faces:

  Optional triangle matrix with 1-based indices.

- method:

  Operator backend: \`"auto"\`, \`"face_graph"\`, \`"cotangent"\`, or
  \`"knn_laplacian"\`.

- weight_mode:

  Edge weighting for graph-based operators: \`"heat"\` or \`"binary"\`.

- normalized:

  Whether to normalize graph-based operators. Cotangent operators always
  use the mass normalization selected by \`mass_mode\`.

- sigma:

  Optional scale used by \`"heat"\` weighting.

- nnk:

  Number of neighbors used by \`"knn_laplacian"\` when needed.

- eps:

  Positive numerical stabilization floor for graph degrees or cotangent
  vertex masses. Stabilized masses are also used for projection.

- mass_mode:

  Mass normalization for \`method = "cotangent"\`: \`"lumped"\` or
  \`"uniform"\`.

## Value

Symmetric operator matrix for use with \[fm_basis()\].
