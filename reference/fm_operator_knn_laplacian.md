# Build a kNN Laplacian Sample Operator

Build a kNN Laplacian Sample Operator

## Usage

``` r
fm_operator_knn_laplacian(
  x,
  k = 8,
  sigma = NULL,
  normalized = TRUE,
  eps = 1e-08
)
```

## Arguments

- x:

  Numeric matrix with samples in rows and features in columns.

- k:

  Number of nearest neighbors per sample.

- sigma:

  Optional RBF scale. If \`NULL\`, estimated from neighbor distances.

- normalized:

  Whether to return normalized Laplacian (\`TRUE\`) or combinatorial
  Laplacian (\`FALSE\`).

- eps:

  Numerical floor for degree stabilization.

## Value

Symmetric \`n_samples x n_samples\` Laplacian-like operator matrix.
