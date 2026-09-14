# Build Heat Kernel Signature (HKS) Descriptors

Build Heat Kernel Signature (HKS) Descriptors

## Usage

``` r
fm_descriptor_hks(
  domain,
  n_times = 16,
  time_range = NULL,
  scaled = TRUE,
  landmarks = NULL
)
```

## Arguments

- domain:

  \`fm_domain\` with basis vectors/values.

- n_times:

  Number of diffusion times.

- time_range:

  Optional positive range \`c(t_min, t_max)\`.

- scaled:

  Whether to normalize each descriptor column.

- landmarks:

  Optional landmark indices for landmark-HKS blocks.

## Value

Numeric matrix with one column per time (or per landmark-time pair).
