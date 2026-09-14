# Build an Operator for Generic Multidimensional Data

Build an Operator for Generic Multidimensional Data

## Usage

``` r
fm_operator(x, method = c("covariance", "knn_laplacian"), ...)
```

## Arguments

- x:

  Numeric matrix with samples in rows and features in columns.

- method:

  Operator family: \`"covariance"\` or \`"knn_laplacian"\`.

- ...:

  Additional method-specific arguments.

## Value

Symmetric sample operator matrix.
