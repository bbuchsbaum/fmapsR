# PRD: fmapsR v1

## 1. Document Control
- Product: `fmapsR`
- Version: `v1.0-draft`
- Date: `2026-02-12`
- Status: `Proposed`
- Audience: Research engineers, package maintainers, applied scientists

## 2. Summary
`fmapsR` will be a high-level, idiomatic R package for functional maps on geometric and non-geometric domains, with first-class support for:
- pairwise map estimation,
- multi-dataset consensus maps,
- map synchronization and cycle consistency optimization,
- fast, scalable execution.

This is not a line-by-line port of `pyFM`. It is a redesign with a cleaner abstraction boundary, stronger API ergonomics, and higher performance targets.

## 3. Problem Statement
Current tooling around functional maps is often:
- mesh-centric to the point of excluding broader multidimensional datasets,
- implementation-heavy for end users,
- weak on multi-dataset synchronization workflows,
- difficult to optimize systematically in R.

Users need one framework that handles shapes, point clouds, graphs, and general multidimensional sampled domains under a common interface.

## 4. Product Goals
1. Deliver an idiomatic R interface that hides low-level plumbing while preserving algorithmic control.
2. Support domain-agnostic functional map workflows (not shape-only).
3. Enable multi-dataset pipelines: consensus bases, synchronized maps, and network optimization.
4. Achieve significant speedups via optimized kernels and sparse linear algebra.
5. Provide diagnostics and uncertainty signals for map quality and synchronization quality.

## 5. Non-Goals (v1)
- Full support for every geometric operator published in all FM papers.
- Perfect one-to-one API compatibility with `pyFM`.
- GPU-first architecture (optional later).
- End-to-end interactive GUI.

## 6. Core Product Principles
- Operator-first design: maps act on functions; geometry adapters provide operators.
- Progressive complexity: simple defaults, expert overrides.
- Composability: pairwise maps, map networks, and consensus objects share common contracts.
- Performance by default: sparse ops, cached eigensystems, vectorized kernels, optional compiled backends.

## 7. Target Users and Primary Jobs
- Method developer: prototype new FM objectives and synchronization terms.
- Applied analyst: compute correspondences and transfer signals across many datasets.
- Pipeline engineer: run reproducible, benchmarkable multi-dataset mapping workflows.

## 8. Scope and Functional Requirements

### FR1. Unified Domain Abstraction
Provide a common `fm_domain` abstraction for:
- triangular meshes,
- point clouds,
- graphs,
- regular multidimensional grids / sampled fields,
- generic `n x d` sampled domains with a user-supplied operator.

Each domain must expose:
- sample count,
- measure/weight matrix or vector,
- basis + eigenvalues,
- projection/unprojection ops,
- optional adjacency/neighborhood.

### FR2. Basis and Operator Pipeline
Support computation and caching of basis/operator pairs:
- Laplace-type operators where defined,
- user-supplied operators for non-geometric data,
- reproducible eigendecompositions with explicit `k` and solver settings.

### FR3. Pairwise Functional Map Estimation
Provide `fm_match()` for pairwise estimation with configurable objective terms:
- descriptor preservation,
- Laplacian commutativity,
- descriptor-operator commutativity,
- optional orientation-like constraints where domain supports them.

### FR4. Refinement
Provide `fm_refine()` with at least:
- ICP-style refinement,
- ZoomOut-style spectral upsampling.

### FR5. Multi-Dataset Map Networks
Provide `fm_network()` to manage many domains and pairwise maps:
- directed/undirected graph of maps,
- weighted edges,
- cycle extraction,
- cycle inconsistency scoring.

### FR6. Consensus and Latent Alignment
Provide consensus mapping utilities:
- canonical/latent basis alignment across datasets,
- consensus map or hub-free shared latent space,
- per-dataset embeddings derived from synchronized maps.

### FR7. Synchronization Optimization
Provide `fm_sync()` with optimization modes:
- adjacency-weighted synchronization,
- cycle-consistency-weighted synchronization,
- robust weighting to downweight noisy edges,
- optional alternating refinement with pairwise map updates.

### FR8. Transfer and Conversion Utilities
Support:
- `as_p2p()` conversion where valid,
- functional transfer (`fm_transfer()`),
- map composition/inversion in functional space,
- accuracy and consistency metrics.

### FR9. Diagnostics and Evaluation
Provide standard reports:
- objective decomposition,
- cycle error distribution,
- coverage/continuity metrics when geometry available,
- synchronization gain vs unsynchronized baseline.

### FR10. Reproducibility
Support seed control, serialized model state, and reproducible benchmark runs.

## 9. API UX Requirements (Idiomatic R)
- Primary API style: function-first with lightweight S3 classes.
- Objects: `fm_domain`, `fm_fit`, `fm_map`, `fm_network_fit`.
- Predictable verbs: `fm_basis()`, `fm_match()`, `fm_refine()`, `fm_sync()`, `fm_consensus()`, `fm_transfer()`, `as_p2p()`.
- Pipe-friendly usage and explicit argument names.
- Clear defaults with structured control lists for advanced settings.

## 10. Performance and Technical Requirements

### NFR1. Runtime
- Pairwise matching + refinement should be at least `30%` faster than a direct `pyFM` reference workflow on matched benchmark tasks (same machine, matched `k`, descriptors, iterations).
- Synchronization on medium graphs (20 datasets, ~80 directed edges) should complete within practical interactive batch times (`< 5 min`) on a modern laptop CPU.

### NFR2. Memory
- For pairwise tasks up to 20k samples/domain and k<=200, memory usage should stay bounded to avoid OOM on 32 GB machines.

### NFR3. Scalability
- Support at least 50 domains in network mode with sparse edge sets.

### NFR4. Reliability
- Numerical routines should return explicit convergence status and diagnostics, not silent failure.

### NFR5. Extensibility
- Adding a new domain backend should not require touching core optimization code.

## 11. Proposed Architecture
- Layer A: Core algebra/optimization kernels (domain-agnostic).
- Layer B: Domain adapters (mesh, cloud, graph, grid, generic operator).
- Layer C: Pairwise FM workflows.
- Layer D: Network synchronization/consensus workflows.
- Layer E: Metrics, diagnostics, and visualization helpers.

Implementation guidance:
- Use `Matrix`, `RSpectra`, and optionally `RcppEigen` for hotspots.
- Keep pure-R fallbacks for portability.
- Cache decompositions and reusable intermediates aggressively.

## 12. Milestones

### M1: Pairwise Core (MVP)
- Domain abstraction + basis pipeline
- Pairwise `fm_match()` and `fm_refine()`
- Basic transfer and p2p conversion
- Baseline benchmark harness

### M2: Network + Synchronization
- `fm_network()`, cycle metrics, `fm_sync()`
- Consensus latent basis and consensus map output
- Synchronization diagnostics

### M3: Performance and Hardening
- Hotspot acceleration
- Robust tests and benchmark thresholds
- Documentation and vignettes

## 13. Acceptance Criteria

### AC1. API Coherence
- End users can complete pairwise workflow in <= 10 lines of R using documented defaults.
- End users can complete multi-dataset sync workflow in <= 20 lines of R.

### AC2. Domain Generality
- Same high-level API works on at least:
  - one triangular mesh dataset,
  - one point cloud dataset,
  - one non-shape multidimensional sampled dataset.

### AC3. Pairwise Quality
- On standard reference shape pairs, accuracy metrics are not worse than `pyFM` by more than `5%` relative (matched settings).

### AC4. Synchronization Quality
- On noisy map-network benchmarks, synchronization reduces median cycle inconsistency by at least `40%` vs pre-sync maps.

### AC5. Consensus Utility
- Consensus latent representation supports downstream transfer/alignment tasks and improves average transfer consistency vs unsynchronized pairwise-only baseline.

### AC6. Performance
- Median runtime improvement >= `30%` vs matched `pyFM` reference on defined benchmark suite.
- Benchmarks and environment metadata are automatically recorded.

### AC7. Reliability and Diagnostics
- All optimization outputs include convergence status, iterations, objective values, and warnings.
- Failures produce actionable errors (invalid domain, unsupported operator, insufficient basis rank).

### AC8. Test Coverage
- Unit + integration tests cover:
  - domain adapters,
  - pairwise optimization,
  - refinement,
  - synchronization,
  - consensus outputs,
  - serialization/reload.
- CI must pass on Linux and macOS.

### AC9. Documentation
- Provide at least three vignettes:
  - pairwise maps,
  - multi-dataset synchronization,
  - non-shape multidimensional workflow.

## 14. Benchmark Plan (Required for Go/No-Go)
- Benchmark suite includes:
  - pairwise mesh benchmark,
  - multi-dataset map-network benchmark,
  - non-shape dataset benchmark.
- Each run logs:
  - runtime,
  - memory,
  - map quality metrics,
  - cycle consistency metrics,
  - solver settings and hardware summary.

## 15. Risks and Mitigations
- Risk: overfitting API to mesh workflows.
  - Mitigation: enforce adapter contract and non-shape acceptance tests from M1.
- Risk: R performance bottlenecks.
  - Mitigation: early profiling; compiled kernels for top hotspots.
- Risk: synchronization objective instability on sparse/noisy graphs.
  - Mitigation: robust weighting, diagnostics, and fallback optimization modes.

## 16. Open Decisions
- S3-only vs selective R6 wrappers for mutable workflows.
- Which compiled backend ships in v1 (`RcppEigen` only vs optional alternatives).
- Default synchronization solver (spectral vs alternating constrained optimization).

## 17. Definition of Done (v1)
`fmapsR v1` is done when AC1-AC9 pass on CI and benchmark targets are met on the defined reference suite.
