# Getting Started with fmapsR

`fmapsR` is organized around a short object flow:

1.  Build a domain with samples, an operator, and optional adjacency.
2.  Add a basis with
    [`fm_basis()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_basis.md).
3.  Fit a map with
    [`fm_match()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_match.md).
4.  Refine, transfer, or score that map with
    [`fm_refine()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_refine.md),
    [`fm_transfer()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_transfer.md),
    [`as_p2p()`](https://bbuchsbaum.github.io/fmapsR/reference/as_p2p.md),
    and
    [`fm_fit_metrics()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_fit_metrics.md).

This vignette shows that flow on a small pairwise example with a known
correspondence, then points to the next articles in the recommended
order.

## First successful fit

``` r

library(fmapsR)
set.seed(7)

make_path_adj <- function(n) {
  A <- matrix(0, nrow = n, ncol = n)
  A[cbind(1:(n - 1), 2:n)] <- 1
  A[cbind(2:n, 1:(n - 1))] <- 1
  A
}

make_toy_pair <- function(n = 24, shift = 2, k = 10, noise_sd = 0.005) {
  A <- make_path_adj(n)
  op <- diag(rowSums(A)) - A

  theta <- seq(0, 1, length.out = n)
  source_data <- cbind(
    theta,
    sin(2 * pi * theta),
    cos(2 * pi * theta),
    sin(4 * pi * theta),
    cos(6 * pi * theta)
  )

  truth <- c((shift + 1):n, seq_len(shift))
  target_data <- source_data[truth, , drop = FALSE] +
    matrix(rnorm(n * ncol(source_data), sd = noise_sd), nrow = n)

  source <- fm_basis(
    fm_domain_generic(
      n_samples = n,
      data = source_data,
      operator = op,
      adjacency = A
    ),
    k = k,
    solver = "base",
    cache = FALSE
  )

  target <- fm_basis(
    fm_domain_generic(
      n_samples = n,
      data = target_data,
      operator = op[truth, truth],
      adjacency = A[truth, truth]
    ),
    k = k,
    solver = "base",
    cache = FALSE
  )

  list(source = source, target = target, truth = truth)
}

toy <- make_toy_pair()
desc_source <- cbind(scale(toy$source$data), fm_descriptor_hks(toy$source, n_times = 8))
desc_target <- cbind(scale(toy$target$data), fm_descriptor_hks(toy$target, n_times = 8))

fit <- fm_match(
  toy$source,
  toy$target,
  descriptors = list(source = desc_source, target = desc_target),
  optimizer = "cg",
  cg_maxit = 6
)
fit <- fm_refine(fit, method = "icp", nit = 5)

p2p <- as_p2p(fit)
metrics <- fm_fit_metrics(fit, truth = toy$truth)

c(
  exact_match_rate = mean(p2p == toy$truth),
  coverage_ratio = metrics$coverage$ratio,
  normalized_geodesic = metrics$geodesic$normalized_mean
)
#>    exact_match_rate      coverage_ratio normalized_geodesic 
#>                   1                   1                   0
```

The fitted object is now the hub for the next operations:

``` r

signal_hat <- fm_transfer(fit, toy$source$data[, 2])

c(
  transfer_mse = mean((signal_hat - toy$target$data[, 2])^2),
  unique_source_matches = length(unique(p2p))
)
#>          transfer_mse unique_source_matches 
#>          3.727986e-04          2.400000e+01
```

## Recommended sequence

1.  `shape-correspondence-posed-meshes`: a visual mesh example with
    transferred heatmaps across two poses.
2.  `pairwise-functional-maps`: the same workflow in more detail,
    including continuity and transfer diagnostics.
3.  `multi-dataset-sync-consensus`: how to combine pairwise maps into a
    network and reduce cycle inconsistency.
4.  `non-shape-multidimensional-data`: how to use the same API on
    generic observation matrices rather than meshes.

If you only remember one object transition, make it this one:
[`fm_basis()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_basis.md)
gives you domains with eigenvectors, and
[`fm_match()`](https://bbuchsbaum.github.io/fmapsR/reference/fm_match.md)
turns those into an `fm_fit` that you can refine, inspect, or apply.
