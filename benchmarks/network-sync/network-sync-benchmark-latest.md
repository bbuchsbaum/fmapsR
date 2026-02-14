# Network Sync Benchmark

- Timestamp: 2026-02-12 08:45:38 EST
- Platform: aarch64-apple-darwin20
- Mode: robust
- Iterations: 8

## Scale Results
| nodes | edges | runtime_sec | throughput_edges_per_sec | pre_cycle_median | post_cycle_median | cycle_reduction_ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 90 | 0.177000 | 508.475 | 2.466572 | 0.000000 | 1.000000 |
| 25 | 600 | 2.963000 | 202.497 | 2.693317 | 0.000000 | 1.000000 |
| 50 | 2450 | 35.578000 | 68.863 | 2.573437 | 0.000000 | 1.000000 |

## Regression Gate
| nodes | target_ratio | threshold_90pct | observed_ratio | pass |
|---:|---:|---:|---:|:---:|
| 10 | 0.600000 | 0.540000 | 1.000000 | TRUE |
| 25 | 0.550000 | 0.495000 | 1.000000 | TRUE |
| 50 | 0.500000 | 0.450000 | 1.000000 | TRUE |

- Overall regression pass: TRUE

