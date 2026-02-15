probe2_parse_kv <- function(lines) {
  out <- list()
  for (line in lines) {
    if (!nzchar(line) || !grepl("=", line, fixed = TRUE)) {
      next
    }
    kv <- strsplit(line, "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- paste(kv[-1], collapse = "=")
    out[[key]] <- val
  }
  out
}

probe2_parse_num_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    numeric()
  } else {
    as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
  }
}

probe2_parse_int_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    integer()
  } else {
    as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
  }
}

probe2_parse_matrix <- function(values_csv, dims_csv) {
  dims <- probe2_parse_int_csv(dims_csv)
  vals <- probe2_parse_num_csv(values_csv)
  matrix(vals, nrow = dims[1], ncol = dims[2], byrow = TRUE)
}

run_refine_objective_probe <- function() {
  root_hints <- unique(Filter(nzchar, c(
    Sys.getenv("FMAPSR_SOURCE_ROOT", unset = ""),
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )))

  candidate_scripts <- unlist(lapply(root_hints, function(root) {
    c(
      file.path(root, "tools", "benchmark_pyfm_refine_objective.py"),
      file.path(root, "fmapsR", "tools", "benchmark_pyfm_refine_objective.py")
    )
  }), use.names = FALSE)

  script <- NULL
  for (p in candidate_scripts) {
    pp <- normalizePath(p, winslash = "/", mustWork = FALSE)
    if (file.exists(pp)) {
      script <- pp
      break
    }
  }
  if (is.null(script)) {
    testthat::skip("tools/benchmark_pyfm_refine_objective.py unavailable outside source checkout")
  }

  py_candidates <- unique(Filter(nzchar, c(
    Sys.getenv("PYFM_BENCH_PYTHON", unset = ""),
    file.path(".venv-bench", "bin", "python"),
    file.path("..", ".venv-bench", "bin", "python"),
    file.path("..", "..", ".venv-bench", "bin", "python"),
    Sys.which("python3"),
    Sys.which("python")
  )))

  py <- NULL
  for (cand in py_candidates) {
    if (file.exists(cand)) {
      py <- cand
      break
    }
    wh <- Sys.which(cand)
    if (nzchar(wh)) {
      py <- wh
      break
    }
  }
  if (is.null(py) || !nzchar(py)) {
    testthat::skip("Python interpreter unavailable for pyFM parity probe")
  }

  out <- tryCatch(
    suppressWarnings(system2(py, script, stdout = TRUE, stderr = TRUE)),
    error = function(e) sprintf("status=error\nerror=system2_failure:%s", conditionMessage(e))
  )

  kv <- probe2_parse_kv(out)
  if (!identical(kv$status, "ok")) {
    err <- if (is.null(kv$error) || !nzchar(kv$error)) "unknown" else kv$error
    testthat::skip(paste("pyFM refine/objective probe unavailable:", err))
  }

  n_ops <- as.integer(kv$n_ops)
  op1_list <- vector("list", n_ops)
  op2_list <- vector("list", n_ops)
  for (i in seq_len(n_ops)) {
    op1_list[[i]] <- probe2_parse_matrix(kv[[sprintf("op1_%d", i)]], kv[[sprintf("op1_%d_dims", i)]])
    op2_list[[i]] <- probe2_parse_matrix(kv[[sprintf("op2_%d", i)]], kv[[sprintf("op2_%d_dims", i)]])
  }

  list(
    icp_nit = as.integer(kv$icp_nit),
    zoom_nit = as.integer(kv$zoom_nit),
    step = probe2_parse_int_csv(kv$step),
    phi1 = probe2_parse_matrix(kv$phi1, kv$phi1_dims),
    phi2 = probe2_parse_matrix(kv$phi2, kv$phi2_dims),
    c0 = probe2_parse_matrix(kv$c0, kv$c0_dims),
    sub_source_1b = probe2_parse_int_csv(kv$sub_source_1b),
    sub_target_1b = probe2_parse_int_csv(kv$sub_target_1b),
    icp_sub = probe2_parse_matrix(kv$icp_sub, kv$icp_sub_dims),
    zoom_sub = probe2_parse_matrix(kv$zoom_sub, kv$zoom_sub_dims),
    weights = probe2_parse_num_csv(kv$weights),
    phi1_w = probe2_parse_matrix(kv$phi1_w, kv$phi1_w_dims),
    phi2_w = probe2_parse_matrix(kv$phi2_w, kv$phi2_w_dims),
    c0_w = probe2_parse_matrix(kv$c0_w, kv$c0_w_dims),
    p2p_w_1b = probe2_parse_int_csv(kv$p2p_w_1b),
    fm_weighted = probe2_parse_matrix(kv$fm_weighted, kv$fm_weighted_dims),
    w_descr = as.numeric(kv$w_descr),
    w_lap = as.numeric(kv$w_lap),
    w_comm = as.numeric(kv$w_comm),
    A = probe2_parse_matrix(kv$A, kv$A_dims),
    B = probe2_parse_matrix(kv$B, kv$B_dims),
    C = probe2_parse_matrix(kv$C, kv$C_dims),
    ev_sqdiff = probe2_parse_matrix(kv$ev_sqdiff, kv$ev_sqdiff_dims),
    op1_list = op1_list,
    op2_list = op2_list,
    e_descr = as.numeric(kv$e_descr),
    e_lap = as.numeric(kv$e_lap),
    e_comm = as.numeric(kv$e_comm),
    e_total = as.numeric(kv$e_total),
    grad_raw = probe2_parse_matrix(kv$grad_raw, kv$grad_raw_dims),
    grad_locked = probe2_parse_matrix(kv$grad_locked, kv$grad_locked_dims)
  )
}

make_probe_fit <- function(phi1, phi2, c0) {
  source <- fm_domain_generic(
    n_samples = nrow(phi1),
    basis = list(vectors = phi1, values = seq_len(ncol(phi1)), k = ncol(phi1))
  )
  target <- fm_domain_generic(
    n_samples = nrow(phi2),
    basis = list(vectors = phi2, values = seq_len(ncol(phi2)), k = ncol(phi2))
  )

  structure(
    list(
      C = c0,
      k1 = ncol(c0),
      k2 = nrow(c0),
      source = source,
      target = target,
      penalties = list(),
      diagnostics = list(
        convergence = 0L,
        message = "probe",
        counts = c(fn = NA_integer_, gradient = NA_integer_),
        total_objective = NA_real_,
        objective_terms = list(),
        warning = NULL
      )
    ),
    class = "fm_fit"
  )
}

test_that("multi-iteration ICP and ZoomOut match pyFM under subsampling", {
  probe <- run_refine_objective_probe()
  fit <- make_probe_fit(probe$phi1, probe$phi2, probe$c0)
  sub <- list(source = probe$sub_source_1b, target = probe$sub_target_1b)

  icp_r <- fm_refine(
    fit,
    method = "icp",
    nit = probe$icp_nit,
    tol = -1,
    subsample = sub
  )
  zoom_r <- fm_refine(
    fit,
    method = "zoomout",
    nit = probe$zoom_nit,
    step = probe$step,
    subsample = sub
  )

  expect_equal(icp_r$C, probe$icp_sub, tolerance = 1e-7)
  expect_equal(zoom_r$C, probe$zoom_sub, tolerance = 1e-7)
})

test_that("weighted p2p-to-fm conversion matches pyFM A2 behavior", {
  probe <- run_refine_objective_probe()

  c_r <- p2p_to_fm_internal(
    p2p_21 = probe$p2p_w_1b,
    source_basis = probe$phi1_w,
    target_basis = probe$phi2_w,
    k1 = ncol(probe$c0_w),
    k2 = nrow(probe$c0_w),
    target_measure_weights = probe$weights,
    target_idx = seq_along(probe$p2p_w_1b)
  )

  expect_equal(c_r, probe$fm_weighted, tolerance = 1e-8)
})

test_that("objective term and gradient kernels match pyFM reference", {
  probe <- run_refine_objective_probe()

  comm_ctx <- list(
    mode = "precomputed",
    list_ops = pack_descriptor_operator_pairs(probe$op1_list, probe$op2_list),
    n_batches = 1L
  )
  penalties <- list(descr = probe$w_descr, lap = probe$w_lap, comm = probe$w_comm)

  terms_r <- energy_breakdown(
    C = probe$C,
    A = probe$A,
    B = probe$B,
    comm_ctx = comm_ctx,
    ev_sqdiff = probe$ev_sqdiff,
    penalties = penalties
  )

  expect_equal(terms_r$descr, probe$e_descr, tolerance = 1e-9)
  expect_equal(terms_r$lap, probe$e_lap, tolerance = 1e-9)
  expect_equal(terms_r$comm, probe$e_comm, tolerance = 1e-9)
  expect_equal(sum(unlist(terms_r)), probe$e_total, tolerance = 1e-9)

  obj_r <- make_objective_cache(
    k2 = nrow(probe$C),
    k1 = ncol(probe$C),
    A = probe$A,
    B = probe$B,
    comm_ctx = comm_ctx,
    ev_sqdiff = probe$ev_sqdiff,
    penalties = penalties
  )
  grad_r <- matrix(
    obj_r$gr(as.vector(probe$C)),
    nrow = nrow(probe$C),
    ncol = ncol(probe$C)
  )
  grad_r_locked <- grad_r
  grad_r_locked[, 1] <- 0

  expect_equal(grad_r, probe$grad_raw, tolerance = 1e-9)
  expect_equal(grad_r_locked, probe$grad_locked, tolerance = 1e-9)

  if (exists("fm_match_energy_terms_cpp", mode = "function") &&
      exists("fm_match_value_grad_cpp", mode = "function")) {
    n_ops <- length(probe$op1_list)
    k1 <- ncol(probe$C)
    k2 <- nrow(probe$C)

    op1_cube <- array(0, dim = c(k1, k1, n_ops))
    op2_cube <- array(0, dim = c(k2, k2, n_ops))
    for (i in seq_len(n_ops)) {
      op1_cube[, , i] <- probe$op1_list[[i]]
      op2_cube[, , i] <- probe$op2_list[[i]]
    }
    op1_t_cube <- aperm(op1_cube, c(2, 1, 3))
    op2_t_cube <- aperm(op2_cube, c(2, 1, 3))

    terms_cpp <- fm_match_energy_terms_cpp(
      C = probe$C,
      A = probe$A,
      B = probe$B,
      ev_sqdiff = probe$ev_sqdiff,
      op1_cube = op1_cube,
      op2_cube = op2_cube,
      w_descr = penalties$descr,
      w_lap = penalties$lap,
      w_comm = penalties$comm
    )
    vg_cpp <- fm_match_value_grad_cpp(
      C = probe$C,
      A = probe$A,
      B = probe$B,
      ev_sqdiff = probe$ev_sqdiff,
      op1_cube = op1_cube,
      op2_cube = op2_cube,
      op1_t_cube = op1_t_cube,
      op2_t_cube = op2_t_cube,
      w_descr = penalties$descr,
      w_lap = penalties$lap,
      w_comm = penalties$comm
    )
    grad_cpp <- as.matrix(vg_cpp$grad)

    expect_equal(as.numeric(terms_cpp$descr), probe$e_descr, tolerance = 1e-9)
    expect_equal(as.numeric(terms_cpp$lap), probe$e_lap, tolerance = 1e-9)
    expect_equal(as.numeric(terms_cpp$comm), probe$e_comm, tolerance = 1e-9)
    expect_equal(as.numeric(vg_cpp$value), probe$e_total, tolerance = 1e-9)
    expect_equal(grad_cpp, probe$grad_raw, tolerance = 1e-9)
  }
})
