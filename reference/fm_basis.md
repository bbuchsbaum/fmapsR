# Compute or Retrieve a Domain Basis

Compute or Retrieve a Domain Basis

## Usage

``` r
fm_basis(
  domain,
  k,
  solver = c("rspectra", "base"),
  which = c("SM", "LM"),
  cache = TRUE,
  seed = NULL,
  control = list()
)
```

## Arguments

- domain:

  \`fm_domain\` object with an \`operator\`.

- k:

  Number of basis vectors.

- solver:

  One of \`"rspectra"\` or \`"base"\`.

- which:

  Eigen target used by solver (\`"SM"\` or \`"LM"\`).

- cache:

  Whether to use process-level in-memory cache.

- seed:

  Optional RNG seed for reproducibility.

- control:

  Optional arguments passed to \`RSpectra::eigs_sym\`, for example
  \`list(sigma = 0, opts = list(ncv = 40))\`. When \`which\` is omitted,
  mesh operators can supply shift-invert defaults. An explicit \`which\`
  disables these defaults. Dense fallback preserves shift-invert
  selection.

## Value

Updated \`fm_domain\` object with populated \`basis\`.
