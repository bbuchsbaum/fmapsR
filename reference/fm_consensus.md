# Build Consensus Alignment from an FM Network

Build Consensus Alignment from an FM Network

## Usage

``` r
fm_consensus(
  x,
  mode = "cycle",
  nit = 20,
  method = c("latent_clb", "sync_transforms"),
  latent_nit = 10,
  latent_tol = 1e-06,
  anchor = NULL
)
```

## Arguments

- x:

  \`fm_network\` or \`fm_network_fit\` object.

- mode:

  Sync mode used when \`x\` is a plain network.

- nit:

  Sync iterations used when \`x\` is a plain network.

- method:

  Consensus method: \`"latent_clb"\` (default) or \`"sync_transforms"\`.

- latent_nit:

  Iterations for latent Procrustes polishing.

- latent_tol:

  Convergence tolerance for latent polishing.

- anchor:

  Optional anchor node for latent consensus orientation.

## Value

\`fm_consensus\` object.
