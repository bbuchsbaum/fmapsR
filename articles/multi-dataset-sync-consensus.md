# Multi-Dataset Synchronization and Consensus

``` r

library(fmapsR)
set.seed(42)

make_domain <- function(n = 8, k = 4) {
  op <- diag(seq_len(n))
  fm_basis(
    fm_domain_generic(n_samples = n, operator = op),
    k = k,
    solver = "base",
    cache = FALSE
  )
}

random_orth <- function(k) {
  qr.Q(qr(matrix(rnorm(k * k), nrow = k, ncol = k)))
}

make_noisy_cycle_network <- function(k = 4, noise_sd = 0.35) {
  d1 <- make_domain(k = k)
  d2 <- make_domain(k = k)
  d3 <- make_domain(k = k)

  X1 <- diag(k)
  X2 <- random_orth(k)
  X3 <- random_orth(k)

  C12 <- X2 %*% t(X1) + matrix(rnorm(k * k, sd = noise_sd), nrow = k)
  C23 <- X3 %*% t(X2) + matrix(rnorm(k * k, sd = noise_sd), nrow = k)
  C31 <- X1 %*% t(X3) + matrix(rnorm(k * k, sd = noise_sd), nrow = k)

  net <- fm_network(list(a = d1, b = d2, c = d3), directed = TRUE)
  net <- fm_network_add_map(net, "a", "b", C12)
  net <- fm_network_add_map(net, "b", "c", C23)
  net <- fm_network_add_map(net, "c", "a", C31)
  net
}

net <- make_noisy_cycle_network()
```

## End-to-end network sync workflow

This vignette starts after pairwise edge maps already exist. The
synthetic maps below are intentionally inconsistent, so the
synchronization report has something real to show.

``` r

synced <- fm_sync(net, mode = "cycle", nit = 10)
sync_report <- fm_network_report(synced)
consensus <- fm_consensus(synced, method = "latent_clb", latent_nit = 8)

sync_report
#> <fm_network_report>
#> Cycle median (pre/post): 2.461 / 1.447e-15
#> Cycle reduction ratio: 1
#> Transfer median (pre/post): NA / NA
#> Transfer reduction ratio: NA
#> Coverage median (pre/post): 0.375 / 0.375
#> Continuity median (pre/post): NA / NA
#> Coverage note: Coverage available; continuity requires target adjacency

c(
  pre_cycle_median = synced$diagnostics$pre_cycle_median,
  post_cycle_median = synced$diagnostics$post_cycle_median,
  cycle_reduction_ratio = synced$diagnostics$cycle_reduction_ratio,
  consensus_cycle_median = consensus$diagnostics$consensus_cycle_median,
  consensus_confidence_median = consensus$uncertainty$summary$median
)
#>            pre_cycle_median           post_cycle_median 
#>                2.461230e+00                1.446513e-15 
#>       cycle_reduction_ratio      consensus_cycle_median 
#>                1.000000e+00                9.892446e-16 
#> consensus_confidence_median 
#>                1.000000e+00
```

## Notes

- [`fm_sync()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_sync.md)
  is doing visible work here because the pre-sync cycle median starts
  above zero and the report compares it to the synchronized network.
- `fm_consensus(method = "latent_clb")` keeps the network in a
  low-cycle-consistency regime and adds an uncertainty summary for the
  consensus edges.
- If your inputs are raw descriptors rather than pre-fit maps, use
  [`fm_network_match()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_network_match.md)
  to build the initial network first.

## Build the Edge Maps from Raw Descriptors

For bigger runs,
[`fm_network_match()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_network_match.md)
batches the pairwise fits before synchronization:

``` r

d1 <- make_domain()
d2 <- make_domain()
d3 <- make_domain()

desc1 <- matrix(rnorm(d1$n_samples * 4), nrow = d1$n_samples)
desc2 <- desc1 + matrix(rnorm(length(desc1), sd = 0.03), nrow = nrow(desc1))
desc3 <- desc1 + matrix(rnorm(length(desc1), sd = 0.04), nrow = nrow(desc1))

net2 <- fm_network_match(
  domains = list(a = d1, b = d2, c = d3),
  descriptors = list(a = desc1, b = desc2, c = desc3),
  edges = c("a->b", "b->c", "c->a"),
  pair_batch_size = 2,
  descriptor_batch_size = 8,
  optimizer = "cg"
)
```
