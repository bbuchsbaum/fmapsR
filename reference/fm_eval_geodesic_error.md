# Geodesic-Style Error for Pointwise Maps

Geodesic-Style Error for Pointwise Maps

## Usage

``` r
fm_eval_geodesic_error(p2p, truth, source_distance, normalize = TRUE)
```

## Arguments

- p2p:

  Integer vector mapping target samples to source indices.

- truth:

  Integer vector of ground-truth source indices for each target sample.

- source_distance:

  Source-domain distance matrix (geodesic or surrogate).

- normalize:

  Whether to normalize by mean finite off-diagonal distance.

## Value

List with per-sample errors and summary statistics.
