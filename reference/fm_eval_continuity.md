# Continuity Proxy for Pointwise Maps

Continuity Proxy for Pointwise Maps

## Usage

``` r
fm_eval_continuity(p2p, source_distance, target_adjacency, normalize = TRUE)
```

## Arguments

- p2p:

  Integer vector mapping target samples to source indices.

- source_distance:

  Source-domain distance matrix.

- target_adjacency:

  Target-domain adjacency matrix.

- normalize:

  Whether to normalize by source-distance scale.

## Value

List with continuity proxy statistics (lower is smoother).
