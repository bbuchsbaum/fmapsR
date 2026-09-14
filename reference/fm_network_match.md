# Fit Pairwise Maps for a Domain Collection and Build a Network

Fit Pairwise Maps for a Domain Collection and Build a Network

## Usage

``` r
fm_network_match(
  domains,
  descriptors,
  edges = NULL,
  directed = TRUE,
  pair_batch_size = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- domains:

  Named list of \`fm_domain\` objects.

- descriptors:

  Named list of descriptor matrices keyed by domain name.

- edges:

  Optional edge specification (\`NULL\`, \`"i-\>j"\` keys, or two-column
  matrix/data.frame with \`i\` and \`j\`).

- directed:

  Whether the resulting network is directed.

- pair_batch_size:

  Optional number of edges processed per batch.

- verbose:

  Emit simple per-batch progress.

- ...:

  Additional arguments passed to \[fm_match()\] (e.g.
  \`descriptor_batch_size\`, \`optimizer\`, \`cg_maxit\`).

## Value

\`fm_network\` with fitted pairwise maps as edges.
