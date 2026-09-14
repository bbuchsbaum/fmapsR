# Evaluate a Fitted Functional Map

Evaluate a Fitted Functional Map

## Usage

``` r
fm_fit_metrics(
  fit,
  truth = NULL,
  source_distance = NULL,
  target_adjacency = NULL,
  normalize = TRUE
)
```

## Arguments

- fit:

  \`fm_fit\` object.

- truth:

  Optional ground-truth pointwise map (target -\> source).

- source_distance:

  Optional source distance matrix.

- target_adjacency:

  Optional target adjacency matrix.

- normalize:

  Whether to normalize distance-based metrics.

## Value

List with pointwise map and available metric summaries.
