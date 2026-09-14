# Create a New Functional-Map Domain

Create a New Functional-Map Domain

## Usage

``` r
fm_new_domain(
  type,
  n_samples,
  data = NULL,
  measure = NULL,
  operator = NULL,
  basis = NULL,
  adjacency = NULL,
  projector = NULL,
  unprojector = NULL,
  metadata = list()
)
```

## Arguments

- type:

  Domain type label.

- n_samples:

  Number of samples in the domain.

- data:

  Backend-specific payload.

- measure:

  Sample measure as weights or matrix.

- operator:

  Domain operator (typically Laplace-like).

- basis:

  Optional basis metadata list.

- adjacency:

  Optional adjacency matrix.

- projector:

  Optional custom projection function.

- unprojector:

  Optional custom reconstruction function.

- metadata:

  Optional named metadata list.

## Value

An \`fm_domain\` object.
