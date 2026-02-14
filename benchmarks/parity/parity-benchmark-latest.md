# Parity Benchmark (fmapsR vs pyFM)

- Timestamp: 2026-02-13 09:36:43 EST
- Platform: aarch64-apple-darwin20
- Scenarios: easy,noisy,partial
- R runs/scenario: 3
- pyFM runs/scenario: 3
- R cg_maxit: 3
- R cg_tol: 1e-06
- R refine_nit: 2
- R backend request: auto
- Package load mode: pkgload
- Baseline available scenarios: 3/3

## Scenario Results
| scenario | r_backend | runtime_comparable | r_runtime_sec | py_runtime_sec | runtime_improvement_ratio | r_accuracy | py_accuracy | accuracy_delta | r_geodesic_norm | py_geodesic_norm | geodesic_norm_delta | baseline_available | note |
|---|---|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---|
| easy | cpp | TRUE | 0.005000 | 0.008522 | 0.413283266838757 | 1.000000 | 1 | 0 | 0.000000 | 0 | 0 | TRUE | ok |
| noisy | cpp | TRUE | 0.005000 | 0.006652084 | 0.248355853594168 | 0.058333 | 0.091666667 | -0.0333333336666667 | 0.899789 | 0.917069477 | -0.0172803802372491 | TRUE | ok |
| partial | cpp | TRUE | 0.005000 | 0.006596792 | 0.242055835624362 | 0.041667 | 0.013888889 | 0.0277777776666667 | 0.991382 | 0.975650728 | 0.0157315151610733 | TRUE | ok |

## Aggregates
- Median runtime improvement ratio (available scenarios): 0.248355853594168
- Mean accuracy delta (R - pyFM): -0.001851852
- Mean geodesic normalized delta (R - pyFM): -0.000516288358725267

