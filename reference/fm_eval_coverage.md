# Coverage Metrics for Pointwise Maps

Coverage Metrics for Pointwise Maps

## Usage

``` r
fm_eval_coverage(p2p, n_source, weights = NULL)
```

## Arguments

- p2p:

  Integer vector mapping target samples to source indices.

- n_source:

  Source sample count.

- weights:

  Optional non-negative source weights.

## Value

List with support, ratio, weighted ratio, and entropy.
