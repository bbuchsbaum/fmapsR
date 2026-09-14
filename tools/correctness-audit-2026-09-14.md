# Correctness audit, 2026-09-14

This audit reviewed and completed the existing uncommitted mesh implementation,
nearest-neighbor changes, tests, and documentation. It also checked the existing
package suite and local CI smoke gates. It is local validation, not a release or
publication certification.

## Corrections

- Sparse eigenpairs now follow the same spectral ordering as the dense backend.
  Eigenvectors have a consistent sign convention, including a positive constant
  Laplacian mode. An analytic path-graph test verifies that the first mode is
  constant and satisfies the eigen-equation.
- Dense fallback preserves shift-invert selection. Previously, a forced sparse
  failure on a diagonal operator selected eigenvalues 8, 7, 6 instead of 1, 2, 3.
  Partial sparse convergence is rejected. RSpectra subspace-size hints now go
  through `opts$ncv`, with explicit eigen-target overrides respected.
- Cotangent stabilization uses the same vertex masses for the operator,
  eigenvector transformation, and projection. Analytic right-triangle stiffness,
  mass scaling, generalized eigen-equations, and reconstruction are tested.
  Mesh validation rejects fractional, non-finite, and empty triangle indices.
  Cotangent operators also reject degenerate triangles and unused vertices
  instead of silently changing their geometry or mass.
- Graph shortest-path size limits are checked before sparse-to-dense conversion;
  self-distance remains zero even when an adjacency has positive diagonal entries.
  Mesh topology distances are documented as edge hops, not physical geodesics.
- Blocked nearest-neighbor implementations agree with independently computed
  direct distances across block boundaries. Empty queries, invalid inputs, and
  first-index tie handling are covered.
- All five vignette engines are recognized again. API documentation was
  regenerated, including the missing `compute_objective` argument.

## Real-mesh regression

Correcting eigenpair ordering exposed an under-converged six-iteration TOSCA
workflow. Its refined top-5% and top-10% transfer overlaps initially fell to
approximately 0.886 and 0.857, below the existing 0.95 and 0.90 floors.

A convergence probe retained those floors and compared training objectives at
6, 30, 100, and 300 maximum iterations with an absolute CG tolerance of 1e-8.
The fit converged in 124 iterations, with objective 1.139115e-7. After the same
single ICP step, transfer correlation was 0.99495, MSE 0.00066645, and top-5% /
top-10% overlaps were 0.96413 / 0.94120. The vignette and regression now allow
300 iterations and require convergence before refinement. All original transfer
quality thresholds remain unchanged.

This is a supervised, development-consumed example using known landmark
correspondences. Its scores are not evidence of unsupervised or dataset-wide
accuracy. The exact synthetic example shares topology and supplied descriptors.

## Validation

- Full suite: `devtools::test(reporter = "summary", stop_on_failure = TRUE)`
  passed, including the non-CRAN real-mesh regression.
- `covr::package_coverage()`: 86.95%, above the repository's 80% floor.
- `R CMD build ... --no-manual`: passed and built all five vignettes.
- Parity smoke: both `easy` and `noisy` had available Python baselines and passed
  `tools/check_parity_gate.R` with the existing thresholds.
- Network smoke: `--scales 10 --enforce-regression` passed.
- Pairwise smoke: three R and three Python runs completed; the report's stable
  target passed. This concurrent local smoke is not an isolated performance claim.
- Staged whitespace checks passed. Vendored mesh CRLF bytes were preserved via
  path-specific Git attributes.

The initial `R CMD check --as-cran --no-manual` had 0 errors, 3 warnings, and
3 notes. The final check had **0 errors, 0 warnings, and 1 note**, including a
successful rebuild of all vignette outputs. The remaining incoming-feasibility
note concerns new-submission/development-version metadata and the documentation
URL returning HTTP 404. It is not a claim of CRAN submission readiness.

Evidence retained locally:

- Final archive and check log:
  `/private/tmp/fmapsR-verified.tjS9XX/fmapsR_0.0.0.9000.tar.gz` and
  `/private/tmp/fmapsR-verified.tjS9XX/fmapsR.Rcheck/00check.log`.
- Archive SHA-256:
  `e95aeffb4554452e1350cf17cb993578ecedf869f012ac6bad973d29fd67bc2d`.
- Coverage: `/private/tmp/fmapsR-verified.tjS9XX/coverage.rds`.
- Parity: `/private/tmp/fmapsR-final.6Japa1/parity-ordered/`.
- Network and pairwise smoke reports:
  `/private/tmp/fmapsR-verified.tjS9XX/network/` and
  `/private/tmp/fmapsR-verified.tjS9XX/pairwise/`.
- Convergence probe script and snapshots:
  `/private/tmp/fmapsR-verified.tjS9XX/probe-convergence.R` and `probe-*.rds`.
- Initial check log: `/private/tmp/fmapsR-audit.F8CGQq/fmapsR.Rcheck/00check.log`.

## Cleanup and environment

Generated `docs/` files are no longer tracked; the existing GitHub Pages workflow
builds them from source. The local generated site remains on disk. Local theme
backups, Finder metadata, plot output, and development-only directories are
excluded from Git and/or package archives as appropriate.

R checks used R 4.5.1 on arm64 macOS, with `LC_ALL=en_US.UTF-8` and
`R_MAKEVARS_USER=/dev/null`. This selects R's default Apple compiler instead of
the user's Homebrew compiler override, which warned on a pragma in R's own
headers. No global compiler or locale settings were changed. PDF manual
generation, hosted CI, live-site deployment, and release qualification are outside
this validation.
