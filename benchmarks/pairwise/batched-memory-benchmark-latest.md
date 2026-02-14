# Batched Memory Benchmark

- Timestamp: 2026-02-12 05:12:05 EST
- Platform: aarch64-apple-darwin20
- R version: 4.5.1

## Before (Unbatched)
- Runtime (s): 5.745000
- Max RSS (MB): 1004.78
- Objective: 48.849048

## After (Descriptor Batching)
- Runtime (s): 11.461000
- Max RSS (MB): 603.52
- Objective: 48.849048
- Batch size: 8

## Delta
- RSS reduction ratio: 0.399356
- Runtime delta ratio: 0.994952
- Objective absolute delta: 0.000000
- Target met (RSS >= 0.30): TRUE

