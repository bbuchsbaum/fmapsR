# Pairwise Benchmark

- Timestamp: 2026-02-12 20:06:02 EST
- R version: 4.5.1
- Platform: aarch64-apple-darwin20
- Target ratio: 0.3
- Stability margin: 0.02
- Stability bootstrap samples: 500

## R Pipeline
- Median runtime (s): 0.008000
- Mean runtime (s): 0.034333
- Optimizer: cg
- Kernel backend: cpp
- Median objective: 0.042243
- Memory (bytes, optional): 8494096

## Backend Comparison
- Baseline backend: r
- Active backend: cpp
- Baseline median runtime (s): 0.010000
- Active median runtime (s): 0.008000
- Active vs baseline speedup: 0.19999999999996
- Objective median gap: 0

## pyFM Baseline
- Status: ok
- Runs: 3
- Runtime (s): 0.013264292
- Objective: 117844.188198442

## Comparison
- Runtime improvement ratio: 0.396876968631252
- Runtime improvement sample median: 0.382726956444973
- Runtime improvement sample q25/q75: 0.367897039277139 / 0.396876968631252
- Runtime improvement bootstrap samples: 500
- Stability sample usage (R raw/used, pyFM raw/used): 3/2, 3/2
- Target ratio: 0.3
- Stability band: [0.28, 0.32]
- Target met (raw >= target): TRUE
- Stable target met: TRUE
- Stability decision: pass_strong
- Baseline available: TRUE
- Notes: Runtime target met outside jitter band

