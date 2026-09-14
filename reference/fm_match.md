# Fit a Pairwise Functional Map

Fit a Pairwise Functional Map

## Usage

``` r
fm_match(
  source,
  target,
  descriptors,
  penalties = list(descr = 0.1, lap = 0.001, comm = 1),
  init = c("zeros", "identity", "random"),
  maxit = 200,
  optimizer = c("cg", "lbfgsb"),
  cg_tol = 1e-06,
  cg_maxit = NULL,
  kernel_backend = c("auto", "r", "cpp"),
  descriptor_batch_size = NULL,
  compute_objective = TRUE,
  trace = FALSE
)
```

## Arguments

- source:

  Source \`fm_domain\` with basis and descriptors.

- target:

  Target \`fm_domain\` with basis and descriptors.

- descriptors:

  List with \`source\` and \`target\` descriptor matrices.

- penalties:

  Named list with weights for \`descr\`, \`lap\`, and \`comm\`.

- init:

  Initialization strategy (\`"zeros"\`, \`"identity"\`, \`"random"\`).

- maxit:

  Maximum optimization iterations.

- optimizer:

  Optimization backend (\`"cg"\` or \`"lbfgsb"\`).

- cg_tol:

  Residual tolerance used by conjugate-gradient solver.

- cg_maxit:

  Optional maximum conjugate-gradient iterations.

- kernel_backend:

  Kernel backend for objective/gradient ops: \`"auto"\` (prefer
  compiled), \`"r"\`, or \`"cpp"\`.

- descriptor_batch_size:

  Optional descriptor chunk size for commutativity operator streaming.
  Use smaller values to reduce peak memory when many descriptors are
  used.

- compute_objective:

  Whether to compute the final objective breakdown after fitting.
  Disable to skip post-hoc diagnostics when only the map is needed.

- trace:

  Whether to print optimizer traces.

## Value

An \`fm_fit\` object with map, diagnostics, and objective breakdown.
