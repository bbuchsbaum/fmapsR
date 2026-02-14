#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

inline void add_comm_terms(
  const arma::mat& C,
  const arma::cube& op1_cube,
  const arma::cube& op2_cube,
  const arma::cube& op1_t_cube,
  const arma::cube& op2_t_cube,
  const double w_comm,
  double* value,
  arma::mat* grad
) {
  if (w_comm <= 0.0 || op1_cube.n_slices == 0 || op2_cube.n_slices == 0 ||
      op1_t_cube.n_slices == 0 || op2_t_cube.n_slices == 0) {
    return;
  }

  const arma::uword n_ops = std::min(
    std::min(op1_cube.n_slices, op2_cube.n_slices),
    std::min(op1_t_cube.n_slices, op2_t_cube.n_slices)
  );
  for (arma::uword i = 0; i < n_ops; ++i) {
    const arma::mat& op1 = op1_cube.slice(i);
    const arma::mat& op2 = op2_cube.slice(i);
    const arma::mat& op1_t = op1_t_cube.slice(i);
    const arma::mat& op2_t = op2_t_cube.slice(i);
    const arma::mat residual = op2 * C - C * op1;

    if (value != nullptr) {
      *value += w_comm * 0.5 * arma::accu(residual % residual);
    }

    if (grad != nullptr) {
      *grad += w_comm * (op2_t * residual - residual * op1_t);
    }
  }
}

inline arma::mat apply_operator(
  const arma::mat& C,
  const arma::mat& AAt,
  const arma::mat& ev_sqdiff,
  const arma::cube& op1_cube,
  const arma::cube& op2_cube,
  const arma::cube& op1_t_cube,
  const arma::cube& op2_t_cube,
  const double w_descr,
  const double w_lap,
  const double w_comm
) {
  arma::mat out(C.n_rows, C.n_cols, arma::fill::zeros);

  if (w_descr > 0.0) {
    out += w_descr * (C * AAt);
  }

  if (w_lap > 0.0) {
    out += w_lap * (C % ev_sqdiff);
  }

  add_comm_terms(C, op1_cube, op2_cube, op1_t_cube, op2_t_cube, w_comm, nullptr, &out);
  return out;
}

inline arma::uvec nearest_neighbor_index_cpp(
  const arma::mat& reference,
  const arma::mat& query
) {
  if (reference.n_cols != query.n_cols) {
    Rcpp::stop("Embedding dimensions must match for nearest-neighbor lookup");
  }

  arma::uvec out(query.n_rows);
  const arma::vec ref_norm = arma::sum(reference % reference, 1);
  const arma::vec query_norm = arma::sum(query % query, 1);

  arma::mat d2 = -2.0 * (query * reference.t());
  d2.each_col() += query_norm;
  d2.each_row() += ref_norm.t();

  for (arma::uword qi = 0; qi < query.n_rows; ++qi) {
    out[qi] = d2.row(qi).index_min();
  }

  return out;
}

} // namespace

// [[Rcpp::export]]
arma::cube fm_compute_descriptor_operators_cpp(
  const arma::mat& phi,
  const arma::mat& descriptors,
  SEXP measure
) {
  const arma::uword n = phi.n_rows;
  const arma::uword k = phi.n_cols;
  const arma::uword p = descriptors.n_cols;

  if (descriptors.n_rows != n) {
    Rcpp::stop("Descriptor row count must match basis row count");
  }

  arma::mat pinv(k, n, arma::fill::zeros);

  if (Rf_isMatrix(measure)) {
    const arma::mat M = Rcpp::as<arma::mat>(measure);
    if (M.n_rows != n || M.n_cols != n) {
      Rcpp::stop("Measure matrix must be n x n");
    }
    pinv = phi.t() * M;
  } else {
    const arma::vec w = Rcpp::as<arma::vec>(measure);
    if (w.n_elem != n) {
      Rcpp::stop("Measure vector length must equal number of samples");
    }
    pinv = phi.t();
    pinv.each_row() %= w.t();
  }

  arma::cube out(k, k, p, arma::fill::zeros);
  for (arma::uword i = 0; i < p; ++i) {
    const arma::mat Dphi = phi.each_col() % descriptors.col(i);
    out.slice(i) = pinv * Dphi;
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::List fm_match_value_grad_cpp(
  const arma::mat& C,
  const arma::mat& A,
  const arma::mat& B,
  const arma::mat& ev_sqdiff,
  const arma::cube& op1_cube,
  const arma::cube& op2_cube,
  const arma::cube& op1_t_cube,
  const arma::cube& op2_t_cube,
  const double w_descr,
  const double w_lap,
  const double w_comm
) {
  double value = 0.0;
  arma::mat grad(C.n_rows, C.n_cols, arma::fill::zeros);

  if (w_descr > 0.0) {
    const arma::mat residual = C * A - B;
    value += w_descr * 0.5 * arma::accu(residual % residual);
    grad += w_descr * (residual * A.t());
  }

  if (w_lap > 0.0) {
    value += w_lap * 0.5 * arma::accu((C % C) % ev_sqdiff);
    grad += w_lap * (C % ev_sqdiff);
  }

  add_comm_terms(C, op1_cube, op2_cube, op1_t_cube, op2_t_cube, w_comm, &value, &grad);

  return Rcpp::List::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("grad") = grad
  );
}

// [[Rcpp::export]]
Rcpp::List fm_match_energy_terms_cpp(
  const arma::mat& C,
  const arma::mat& A,
  const arma::mat& B,
  const arma::mat& ev_sqdiff,
  const arma::cube& op1_cube,
  const arma::cube& op2_cube,
  const double w_descr,
  const double w_lap,
  const double w_comm
) {
  double e_descr = 0.0;
  double e_lap = 0.0;
  double e_comm = 0.0;

  if (w_descr > 0.0) {
    const arma::mat residual = C * A - B;
    e_descr = w_descr * 0.5 * arma::accu(residual % residual);
  }

  if (w_lap > 0.0) {
    e_lap = w_lap * 0.5 * arma::accu((C % C) % ev_sqdiff);
  }

  if (w_comm > 0.0 && op1_cube.n_slices > 0 && op2_cube.n_slices > 0) {
    const arma::uword n_ops = std::min(op1_cube.n_slices, op2_cube.n_slices);
    for (arma::uword i = 0; i < n_ops; ++i) {
      const arma::mat& op1 = op1_cube.slice(i);
      const arma::mat& op2 = op2_cube.slice(i);
      const arma::mat residual = C * op1 - op2 * C;
      e_comm += w_comm * 0.5 * arma::accu(residual % residual);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("descr") = e_descr,
    Rcpp::Named("lap") = e_lap,
    Rcpp::Named("comm") = e_comm
  );
}

// [[Rcpp::export]]
arma::mat fm_match_apply_operator_cpp(
  const arma::mat& C,
  const arma::mat& AAt,
  const arma::mat& ev_sqdiff,
  const arma::cube& op1_cube,
  const arma::cube& op2_cube,
  const arma::cube& op1_t_cube,
  const arma::cube& op2_t_cube,
  const double w_descr,
  const double w_lap,
  const double w_comm
) {
  return apply_operator(C, AAt, ev_sqdiff, op1_cube, op2_cube, op1_t_cube, op2_t_cube, w_descr, w_lap, w_comm);
}

// [[Rcpp::export]]
Rcpp::List fm_match_solve_cg_cpp(
  const arma::mat& C0,
  const arma::mat& rhs,
  const arma::mat& AAt,
  const arma::mat& ev_sqdiff,
  const arma::cube& op1_cube,
  const arma::cube& op2_cube,
  const arma::cube& op1_t_cube,
  const arma::cube& op2_t_cube,
  const double w_descr,
  const double w_lap,
  const double w_comm,
  const arma::mat& diag_precond,
  const int maxit,
  const double tol
) {
  arma::mat C = C0;
  arma::mat R = rhs - apply_operator(
    C, AAt, ev_sqdiff, op1_cube, op2_cube, op1_t_cube, op2_t_cube, w_descr, w_lap, w_comm
  );
  const bool use_prec = diag_precond.n_elem > 0;
  if (use_prec) {
    if (diag_precond.n_rows != C.n_rows || diag_precond.n_cols != C.n_cols) {
      Rcpp::stop("diag_precond dimensions must match C");
    }
  }

  arma::mat Z = use_prec ? (R / diag_precond) : R;
  arma::mat P = Z;

  double rz = arma::accu(R % Z);
  const double residual0 = std::sqrt(arma::accu(R % R));

  if (!std::isfinite(residual0) || residual0 <= tol) {
    return Rcpp::List::create(
      Rcpp::Named("C") = C,
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("residual") = residual0,
      Rcpp::Named("converged") = true,
      Rcpp::Named("breakdown") = false
    );
  }

  bool converged = false;
  bool breakdown = false;
  int iter_done = 0;

  for (int iter = 1; iter <= maxit; ++iter) {
    const arma::mat HP = apply_operator(
      P, AAt, ev_sqdiff, op1_cube, op2_cube, op1_t_cube, op2_t_cube, w_descr, w_lap, w_comm
    );
    const double denom = arma::accu(P % HP);

    if (!std::isfinite(denom) || std::abs(denom) <= std::numeric_limits<double>::epsilon()) {
      breakdown = true;
      iter_done = iter - 1;
      break;
    }

    const double alpha = rz / denom;
    C += alpha * P;
    const arma::mat R_new = R - alpha * HP;
    const double rr_new = arma::accu(R_new % R_new);
    iter_done = iter;

    if (!std::isfinite(rr_new)) {
      breakdown = true;
      break;
    }

    if (std::sqrt(rr_new) <= tol) {
      R = R_new;
      rz = use_prec ? arma::accu(R_new % (R_new / diag_precond)) : rr_new;
      R = R_new;
      converged = true;
      break;
    }

    const arma::mat Z_new = use_prec ? (R_new / diag_precond) : R_new;
    const double rz_new = arma::accu(R_new % Z_new);
    const double beta = rz_new / rz;
    P = Z_new + beta * P;
    R = R_new;
    Z = Z_new;
    rz = rz_new;
  }

  return Rcpp::List::create(
    Rcpp::Named("C") = C,
    Rcpp::Named("iterations") = iter_done,
    Rcpp::Named("residual") = std::sqrt(arma::accu(R % R)),
    Rcpp::Named("converged") = converged,
    Rcpp::Named("breakdown") = breakdown
  );
}

// [[Rcpp::export]]
Rcpp::IntegerVector fm_nearest_neighbor_index_cpp(
  const arma::mat& reference,
  const arma::mat& query
) {
  const arma::uvec nn = nearest_neighbor_index_cpp(reference, query);
  Rcpp::IntegerVector out(nn.n_elem);
  for (arma::uword i = 0; i < nn.n_elem; ++i) {
    out[static_cast<int>(i)] = static_cast<int>(nn[i] + 1);
  }
  return out;
}
