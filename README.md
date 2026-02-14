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

Pairwise reports include both raw target status and a stability-aware status. Borderline values inside the jitter band (`target ± margin`) are classified using bootstrap spread (`q25`) over independent R/pyFM runtime samples (after dropping one high outlier per side when at least 3 samples are available) to reduce flip-flopping near the threshold.

Run parity benchmark report (runtime + quality vs pyFM, multi-scenario):

```bash
Rscript tools/benchmark_parity.R --scenarios easy,noisy,partial --r-runs 3 --py-runs 3
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
