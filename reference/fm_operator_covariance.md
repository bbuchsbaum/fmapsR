# Build a Covariance-Derived Sample Operator

Build a Covariance-Derived Sample Operator

## Usage

``` r
fm_operator_covariance(x, center = TRUE, scale = TRUE, ridge = 0.001)
```

## Arguments

- x:

  Numeric matrix with samples in rows and features in columns.

- center:

  Whether to center features.

- scale:

  Whether to scale features to unit variance.

- ridge:

  Diagonal ridge regularization added for numerical stability.

## Value

Symmetric \`n_samples x n_samples\` operator matrix.
