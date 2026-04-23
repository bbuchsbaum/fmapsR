benchmark_seed_registry <- function() {
  list(
    pairwise_seed = 42L,
    pairwise_seed_offset = 40L,
    parity_seed = 42L,
    parity_scenario_stride = 100L,
    network_seed = 2026L
  )
}

benchmark_parity_scenarios <- function() {
  list(
    easy = list(
      n = 120L, k = 24L, p = 18L,
      noise = 0.02, basis_noise = 0.00,
      descriptor_corruption = 0.00, eval_fraction = 1.00
    ),
    noisy = list(
      n = 120L, k = 24L, p = 18L,
      noise = 0.12, basis_noise = 0.20,
      descriptor_corruption = 0.30, eval_fraction = 1.00
    ),
    partial = list(
      n = 120L, k = 24L, p = 18L,
      noise = 0.20, basis_noise = 0.40,
      descriptor_corruption = 0.50, eval_fraction = 0.60
    )
  )
}

benchmark_pairwise_problem_spec <- function() {
  list(
    n = 160L,
    k = 32L,
    p = 24L,
    maxit = 120L,
    icp_nit = 5L
  )
}

benchmark_pairwise_solver_defaults <- function() {
  list(
    optimizer = "cg",
    cg_maxit = 1L
  )
}

benchmark_network_defaults <- function() {
  list(
    scales = c(10L, 25L, 50L),
    mode = "robust",
    nit = 8L
  )
}

benchmark_network_reference_targets <- function() {
  c("10" = 0.60, "25" = 0.55, "50" = 0.50)
}

benchmark_parity_gate_thresholds <- function() {
  list(
    required_scenarios = c("easy", "noisy"),
    accuracy_delta_min = c(
      easy = -0.10,
      noisy = -0.15,
      partial = -0.20
    ),
    geodesic_norm_delta_max = c(
      easy = 0.08,
      noisy = 0.15,
      partial = 0.20
    ),
    objective_abs_gap_max = c(
      easy = 1.75,
      noisy = 2.50,
      partial = 3.0
    ),
    map_fro_norm_delta_abs_max = c(
      easy = 0.05,
      noisy = 0.05,
      partial = 0.10
    ),
    map_orth_resid_delta_abs_max = c(
      easy = 0.05,
      noisy = 0.08,
      partial = 0.10
    )
  )
}

benchmark_release_claim_thresholds <- function() {
  list(
    required_quality_scenarios = c("easy", "noisy", "partial"),
    speed_source = "parity",
    speed_target_ratio = 0.30,
    speed_stability_margin = 0.02,
    require_stable_speed = FALSE,
    parity_r_cg_maxit = 1L,
    require_objective_parity = FALSE
  )
}
