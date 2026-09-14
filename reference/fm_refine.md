# Refine a Functional Map Fit

Refine a Functional Map Fit

## Usage

``` r
fm_refine(
  fit,
  method = c("icp", "zoomout"),
  nit = 10,
  step = 1,
  tol = NULL,
  use_adj = FALSE,
  subsample = NULL,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- fit:

  \`fm_fit\` object.

- method:

  One of \`"icp"\` or \`"zoomout"\`.

- nit:

  Number of iterations.

- step:

  ZoomOut step size (scalar or length-2 integer vector).

- tol:

  Optional stopping tolerance on max elementwise map change.

- use_adj:

  Use adjoint-style conversion for ICP iterations.

- subsample:

  Optional subsampling (\`NULL\`, integer, or list(source=, target=)).

- seed:

  Optional seed used for subsampling.

- verbose:

  Emit iteration progress.

## Value

Refined \`fm_fit\` object with refinement diagnostics.
