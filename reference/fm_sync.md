# Synchronize Maps in a Functional-Map Network

Synchronize Maps in a Functional-Map Network

## Usage

``` r
fm_sync(
  network,
  mode = c("adjacency", "cycle", "robust"),
  nit = 20,
  tol = 1e-06,
  anchor = NULL,
  verbose = FALSE
)
```

## Arguments

- network:

  \`fm_network\` object.

- mode:

  Synchronization mode: \`"adjacency"\`, \`"cycle"\`, or \`"robust"\`.

- nit:

  Maximum synchronization iterations.

- tol:

  Convergence tolerance for latent transforms.

- anchor:

  Optional anchor node name (defaults to first domain).

- verbose:

  Emit per-iteration diagnostics.

## Value

\`fm_network_fit\` object with synchronized maps and diagnostics.
