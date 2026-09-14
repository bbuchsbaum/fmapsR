# Build a Sample Distance Matrix for a Domain

Build a Sample Distance Matrix for a Domain

## Usage

``` r
fm_distance_matrix(domain, method = c("auto", "euclidean", "graph_shortest"))
```

## Arguments

- domain:

  \`fm_domain\` object.

- method:

  Distance strategy: \`"auto"\`, \`"euclidean"\`, or
  \`"graph_shortest"\`.

## Value

Dense distance matrix.

## Details

\`"auto"\` uses graph distances for meshes with topology and graphs, and
Euclidean distances for point clouds. Mesh faces supply unit-length
edges, so these distances count edge hops, not physical surface length.
Explicit adjacency entries are treated as edge lengths. The current
dense shortest-path implementation supports at most 2000 samples.
