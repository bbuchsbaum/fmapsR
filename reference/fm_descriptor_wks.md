# Build Wave Kernel Signature (WKS) Descriptors

Build Wave Kernel Signature (WKS) Descriptors

## Usage

``` r
fm_descriptor_wks(
  domain,
  n_energies = 16,
  sigma = NULL,
  energy_range = NULL,
  scaled = TRUE,
  landmarks = NULL
)
```

## Arguments

- domain:

  \`fm_domain\` with basis vectors/values.

- n_energies:

  Number of sampled energies.

- sigma:

  Optional Gaussian width in log-spectrum.

- energy_range:

  Optional energy range on log-spectrum scale.

- scaled:

  Whether to normalize each descriptor column.

- landmarks:

  Optional landmark indices for landmark-WKS blocks.

## Value

Numeric matrix with one column per energy (or per landmark-energy pair).
