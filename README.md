# fmapsR

`fmapsR` is an R-first redesign of functional maps for geometric and general multidimensional domains.

## Project Status

Core M1-M3 implementation is in place, including:

- Pairwise functional map estimation and refinement.
- Multi-dataset synchronization and consensus mapping.
- Domain adapters for geometry and generic multidimensional data.
- Basis-derived descriptors (`HKS`, `WKS`) and evaluation metrics (geodesic-style error, continuity, coverage).

## Quick Descriptor + Metrics Example

```r
library(fmapsR)

d <- fm_basis(fm_domain_generic(n_samples = 30, operator = diag(1:30)), k = 12, solver = "base", cache = FALSE)
H <- fm_descriptor_hks(d, n_times = 8)
W <- fm_descriptor_wks(d, n_energies = 8)

# Use descriptors in fm_match(), then score outputs with fm_fit_metrics()
```

## Documentation

Vignettes:

- `getting-started`
- `shape-correspondence-posed-meshes`
- `pairwise-functional-maps`
- `multi-dataset-sync-consensus`
- `non-shape-multidimensional-data`

## Development Setup

Install required R packages:

```bash
Rscript scripts/setup-dev.R
```

Install parity Python deps into `.venv-bench` (optional, required for pyFM parity benchmarks):

```bash
Rscript scripts/setup-dev.R --with-python
```

Then run tests:

```bash
Rscript -e "testthat::test_local('.')"
```

Run a pairwise benchmark report (with configurable stability policy):

```bash
Rscript tools/benchmark_pairwise.R --r-runs 3 --py-runs 3 --target-ratio 0.30 --stability-margin 0.02
```

The pairwise benchmark defaults to the tuned `cg` solver path with `cg_maxit = 1` and performs two unmeasured warmup runs before recording timings, so the report reflects steady-state runtime rather than first-hit overhead.

Pairwise reports include both raw target status and a stability-aware status. Borderline values inside the jitter band (`target ± margin`) are classified using bootstrap spread (`q25`) over independent R/pyFM runtime samples (after dropping one high outlier per side when at least 3 samples are available) to reduce flip-flopping near the threshold.

Run parity benchmark report (runtime + quality vs pyFM, multi-scenario):

```bash
Rscript tools/benchmark_parity.R --scenarios easy,noisy,partial --r-runs 3 --py-runs 3
```

Enforce the parity quality gate (accuracy/geodesic + objective/map deltas):

```bash
Rscript tools/check_parity_gate.R --report benchmarks/parity/parity-benchmark-latest.rds --required easy,noisy
```

Run the combined release claim benchmark (quality parity + pairwise speed proof against vendored `pyFM`):

```bash
Rscript tools/benchmark_release_claim.R --output-dir benchmarks/release-claim
```

Enforce the combined release claim gate:

```bash
Rscript tools/check_release_claim.R --report benchmarks/release-claim/release-claim-latest.rds
```

The release claim uses the multi-scenario parity suite as the blocking speed-and-quality contract, with pairwise benchmarking retained as a diagnostic benchmark in the release report. Objective-gap reporting remains available in parity artifacts, but the default release gate uses accuracy, normalized geodesic error, and map residuals as the required quality metrics.

Recommend updated parity thresholds from historical reports:

```bash
Rscript tools/recommend_parity_thresholds.R --dir benchmarks/parity --quantile 0.95 --safety 1.25 --min-n 3
```

If pyFM dependencies are unavailable, parity tests are skip-safe and parity reports mark unavailable baselines explicitly.

Run the batched-memory benchmark (RSS before/after):

```bash
Rscript tools/benchmark_batched_memory.R
```

Run the multi-scale network sync benchmark and regression gate:

```bash
Rscript tools/benchmark_network_sync.R --scales 10,25,50 --enforce-regression
```

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The defaults are configured via `params$family` and `params$preset` (family = 'teal', preset = 'homage'). The pkgdown site uses `template: { package: albersdown }` together with generated `pkgdown/extra.css` and `pkgdown/extra.js` so the theme is linked and activated on site pages.
<!-- albersdown:theme-note:end -->
