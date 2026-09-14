# Unified Domain Constructor

Unified Domain Constructor

## Usage

``` r
fm_domain(
  data = NULL,
  type = c("generic", "mesh", "pointcloud", "graph"),
  n_samples = NULL,
  ...
)
```

## Arguments

- data:

  Domain payload.

- type:

  One of \`"generic"\`, \`"mesh"\`, \`"pointcloud"\`, or \`"graph"\`.

- n_samples:

  Number of samples for generic domains.

- ...:

  Adapter-specific arguments.

## Value

An \`fm_domain\` object.
