// Thiago de Paula Oliveira

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "matrixCorr_omp.h"

// Fenwick-tree utilities for prefix sums over ranked y values.
inline void bit_add(std::vector<double>& bit, const int idx, const double val) {
  const int m = static_cast<int>(bit.size()) - 1;
  for (int i = idx; i <= m; i += (i & -i)) bit[static_cast<std::size_t>(i)] += val;
}

inline double bit_sum(const std::vector<double>& bit, const int idx) {
  double out = 0.0;
  for (int i = idx; i > 0; i -= (i & -i)) out += bit[static_cast<std::size_t>(i)];
  return out;
}

inline double dot_ptr(const double* a, const double* b, const int n) {
  double out = 0.0;
  for (int i = 0; i < n; ++i) out += a[i] * b[i];
  return out;
}

inline double dcor_signed_ratio(const double xy, const double x2, const double y2) {
  if (!std::isfinite(x2) || !std::isfinite(y2) || x2 <= 0.0 || y2 <= 0.0) {
    return NA_REAL;
  }
  const double denom = std::sqrt(x2 * y2);
  const double dcor = xy / denom;
  if (!std::isfinite(dcor)) return NA_REAL;
  return dcor;
}

inline double clip_dcor_estimate(const double dcor) {
  if (!std::isfinite(dcor)) return NA_REAL;
  if (dcor < 0.0) return 0.0;
  if (dcor > 1.0) return 1.0;
  return dcor;
}

inline double u_centered_cov_stat(
  const double pair_sum,
  const double row_dot,
  const double total_x,
  const double total_y,
  const int n
) {
  const double dn = static_cast<double>(n);
  const double term1 = pair_sum / (dn * static_cast<double>(n - 3));
  const double term2 = 2.0 * row_dot / (dn * static_cast<double>(n - 2) * static_cast<double>(n - 3));
  const double term3 = (total_x * total_y) /
    (dn * static_cast<double>(n - 1) * static_cast<double>(n - 2) * static_cast<double>(n - 3));
  return term1 - term2 + term3;
}

inline double finalize_dcor(const double xy, const double x2, const double y2) {
  return clip_dcor_estimate(dcor_signed_ratio(xy, x2, y2));
}

inline void dcor_t_test_from_signed(
  const double bc_dcor,
  const int n,
  double& tstat,
  double& df,
  double& p_value
) {
  const double dn = static_cast<double>(n);
  const double M = dn * static_cast<double>(n - 3) / 2.0;
  df = M - 1.0;
  if (!std::isfinite(bc_dcor)) {
    tstat = NA_REAL;
    p_value = NA_REAL;
    return;
  }

  const double bc_clamped = std::max(-1.0, std::min(1.0, bc_dcor));
  const double denom_sq = 1.0 - bc_clamped * bc_clamped;
  if (denom_sq <= 0.0) {
    tstat = (bc_clamped >= 0.0) ? R_PosInf : R_NegInf;
    p_value = (bc_clamped >= 0.0) ? 0.0 : 1.0;
    return;
  }

  tstat = std::sqrt(M - 1.0) * bc_clamped / std::sqrt(denom_sq);
  if (!std::isfinite(tstat)) {
    p_value = NA_REAL;
    return;
  }
  p_value = R::pt(tstat, df, /*lower_tail*/ 0, /*log_p*/ 0);
}

inline void compute_row_sums_and_total_quadratic_ptr(
  const double* x,
  const int n,
  double* row_sums,
  double& total_sum
) {
  std::fill(row_sums, row_sums + n, 0.0);
  total_sum = 0.0;

  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    for (int j = i + 1; j < n; ++j) {
      const double dx = std::abs(xi - x[j]);
      row_sums[i] += dx;
      row_sums[j] += dx;
      total_sum += 2.0 * dx;
    }
  }
}

// O(n log n) row-sum computation in sorted x-order:
// r_i = sum_j |x_i - x_j| and S = sum_i r_i.
inline void compute_row_sums_and_total_fast_ptr(
  const double* x,
  const int n,
  double* row_sums,
  double& total_sum
) {
  static thread_local std::vector<int> order;
  order.resize(static_cast<std::size_t>(n));
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [x](const int a, const int b) {
    if (x[a] < x[b]) return true;
    if (x[a] > x[b]) return false;
    return a < b;
  });

  double grand = 0.0;
  for (int i = 0; i < n; ++i) grand += x[i];

  double prefix = 0.0;
  total_sum = 0.0;
  for (int pos = 0; pos < n; ++pos) {
    const int idx = order[static_cast<std::size_t>(pos)];
    const double xi = x[idx];
    const double ri =
      (2.0 * static_cast<double>(pos) - static_cast<double>(n)) * xi +
      grand - 2.0 * prefix;
    row_sums[idx] = ri;
    total_sum += ri;
    prefix += xi;
  }
}

// Fast O(n log n) evaluation of:
// sum_{i != j} |x_i - x_j| |y_i - y_j|
// using x-order traversal and Fenwick trees over y-ranks.
inline double sum_cross_abs_prod_ptr(const double* x, const double* y, const int n) {
  if (n < 2) return 0.0;

  static thread_local std::vector<int> order;
  static thread_local std::vector<int> y_order;
  static thread_local std::vector<int> rank_y;
  static thread_local std::vector<double> bit_count;
  static thread_local std::vector<double> bit_y;
  static thread_local std::vector<double> bit_x;
  static thread_local std::vector<double> bit_xy;

  order.resize(static_cast<std::size_t>(n));
  y_order.resize(static_cast<std::size_t>(n));
  rank_y.resize(static_cast<std::size_t>(n));
  std::iota(order.begin(), order.end(), 0);
  std::iota(y_order.begin(), y_order.end(), 0);

  std::stable_sort(order.begin(), order.end(), [x, y](const int a, const int b) {
    if (x[a] < x[b]) return true;
    if (x[a] > x[b]) return false;
    if (y[a] < y[b]) return true;
    if (y[a] > y[b]) return false;
    return a < b;
  });
  std::stable_sort(y_order.begin(), y_order.end(), [y](const int a, const int b) {
    if (y[a] < y[b]) return true;
    if (y[a] > y[b]) return false;
    return a < b;
  });

  int m = 0;
  bool first = true;
  double last = 0.0;
  for (int pos = 0; pos < n; ++pos) {
    const int idx = y_order[static_cast<std::size_t>(pos)];
    const double yi = y[idx];
    if (first || yi != last) {
      ++m;
      last = yi;
      first = false;
    }
    rank_y[static_cast<std::size_t>(idx)] = m;
  }

  bit_count.assign(static_cast<std::size_t>(m + 1), 0.0);
  bit_y.assign(static_cast<std::size_t>(m + 1), 0.0);
  bit_x.assign(static_cast<std::size_t>(m + 1), 0.0);
  bit_xy.assign(static_cast<std::size_t>(m + 1), 0.0);

  double total_count = 0.0;
  double total_y = 0.0;
  double total_x = 0.0;
  double total_xy = 0.0;
  double out = 0.0;

  for (int t = 0; t < n; ++t) {
    const int idx = order[static_cast<std::size_t>(t)];
    const int r = rank_y[static_cast<std::size_t>(idx)];
    const double xi = x[idx];
    const double yi = y[idx];

    const double count_le = bit_sum(bit_count, r);
    const double sum_y_le = bit_sum(bit_y, r);
    const double sum_x_le = bit_sum(bit_x, r);
    const double sum_xy_le = bit_sum(bit_xy, r);

    const double count_gt = total_count - count_le;
    const double sum_y_gt = total_y - sum_y_le;
    const double sum_x_gt = total_x - sum_x_le;
    const double sum_xy_gt = total_xy - sum_xy_le;

    const double a_j =
      yi * count_le - sum_y_le + sum_y_gt - yi * count_gt;
    const double b_j =
      yi * (sum_x_le - sum_x_gt) + (sum_xy_gt - sum_xy_le);

    out += xi * a_j - b_j;

    bit_add(bit_count, r, 1.0);
    bit_add(bit_y, r, yi);
    bit_add(bit_x, r, xi);
    bit_add(bit_xy, r, xi * yi);
    total_count += 1.0;
    total_y += yi;
    total_x += xi;
    total_xy += xi * yi;
  }

  return 2.0 * out;
}

inline double ustat_dcor_from_precomputed_quadratic_ptr(
  const double* x,
  const double* y,
  const double* row_sums_x,
  const double* row_sums_y,
  const double total_sum_x,
  const double total_sum_y,
  const int n
) {
  const double inv_nm2 = 1.0 / static_cast<double>(n - 2);
  const double add_x = total_sum_x / static_cast<double>((n - 1) * (n - 2));
  const double add_y = total_sum_y / static_cast<double>((n - 1) * (n - 2));

  double XY = 0.0;
  double X2 = 0.0;
  double Y2 = 0.0;

  for (int i = 0; i < n; ++i) {
    const double rxi = row_sums_x[i];
    const double ryi = row_sums_y[i];
    const double xi = x[i];
    const double yi = y[i];
    for (int j = i + 1; j < n; ++j) {
      const double dx = std::abs(xi - x[j]);
      const double dy = std::abs(yi - y[j]);
      const double ax = dx - (rxi + row_sums_x[j]) * inv_nm2 + add_x;
      const double ay = dy - (ryi + row_sums_y[j]) * inv_nm2 + add_y;
      XY += ax * ay;
      X2 += ax * ax;
      Y2 += ay * ay;
    }
  }

  const double scale = 2.0 / (static_cast<double>(n) * static_cast<double>(n - 3));
  XY *= scale;
  X2 *= scale;
  Y2 *= scale;
  return finalize_dcor(XY, X2, Y2);
}

inline double ustat_dcor_signed_from_precomputed_quadratic_ptr(
  const double* x,
  const double* y,
  const double* row_sums_x,
  const double* row_sums_y,
  const double total_sum_x,
  const double total_sum_y,
  const int n
) {
  const double inv_nm2 = 1.0 / static_cast<double>(n - 2);
  const double add_x = total_sum_x / static_cast<double>((n - 1) * (n - 2));
  const double add_y = total_sum_y / static_cast<double>((n - 1) * (n - 2));

  double XY = 0.0;
  double X2 = 0.0;
  double Y2 = 0.0;

  for (int i = 0; i < n; ++i) {
    const double rxi = row_sums_x[i];
    const double ryi = row_sums_y[i];
    const double xi = x[i];
    const double yi = y[i];
    for (int j = i + 1; j < n; ++j) {
      const double dx = std::abs(xi - x[j]);
      const double dy = std::abs(yi - y[j]);
      const double ax = dx - (rxi + row_sums_x[j]) * inv_nm2 + add_x;
      const double ay = dy - (ryi + row_sums_y[j]) * inv_nm2 + add_y;
      XY += ax * ay;
      X2 += ax * ax;
      Y2 += ay * ay;
    }
  }

  const double scale = 2.0 / (static_cast<double>(n) * static_cast<double>(n - 3));
  XY *= scale;
  X2 *= scale;
  Y2 *= scale;
  return dcor_signed_ratio(XY, X2, Y2);
}

inline double ustat_dcor_quadratic_ptr(const double* x, const double* y, const int n) {
  if (n < 4) Rcpp::stop("Sample size must be at least 4 for unbiased dCor");
  static thread_local std::vector<double> Rx_buf;
  static thread_local std::vector<double> Ry_buf;
  Rx_buf.assign(static_cast<std::size_t>(n), 0.0);
  Ry_buf.assign(static_cast<std::size_t>(n), 0.0);
  double Sx = 0.0;
  double Sy = 0.0;
  compute_row_sums_and_total_quadratic_ptr(x, n, Rx_buf.data(), Sx);
  compute_row_sums_and_total_quadratic_ptr(y, n, Ry_buf.data(), Sy);
  return ustat_dcor_from_precomputed_quadratic_ptr(
    x, y, Rx_buf.data(), Ry_buf.data(), Sx, Sy, n
  );
}

inline double ustat_dcor_quadratic_signed_ptr(const double* x, const double* y, const int n) {
  if (n < 4) Rcpp::stop("Sample size must be at least 4 for unbiased dCor");
  static thread_local std::vector<double> Rx_buf;
  static thread_local std::vector<double> Ry_buf;
  Rx_buf.assign(static_cast<std::size_t>(n), 0.0);
  Ry_buf.assign(static_cast<std::size_t>(n), 0.0);
  double Sx = 0.0;
  double Sy = 0.0;
  compute_row_sums_and_total_quadratic_ptr(x, n, Rx_buf.data(), Sx);
  compute_row_sums_and_total_quadratic_ptr(y, n, Ry_buf.data(), Sy);
  return ustat_dcor_signed_from_precomputed_quadratic_ptr(
    x, y, Rx_buf.data(), Ry_buf.data(), Sx, Sy, n
  );
}

// Fast per-pair path:
// - O(n log n) cross-distance term via Fenwick trees
// - exact U-statistic normalization for unbiased dCor
inline double ustat_dcor_fast_ptr(const double* x, const double* y, const int n) {
  if (n < 4) Rcpp::stop("Sample size must be at least 4 for unbiased dCor");

  static thread_local std::vector<double> Rx_buf;
  static thread_local std::vector<double> Ry_buf;
  Rx_buf.assign(static_cast<std::size_t>(n), 0.0);
  Ry_buf.assign(static_cast<std::size_t>(n), 0.0);

  double Sx = 0.0;
  double Sy = 0.0;
  compute_row_sums_and_total_fast_ptr(x, n, Rx_buf.data(), Sx);
  compute_row_sums_and_total_fast_ptr(y, n, Ry_buf.data(), Sy);

  double sum_x = 0.0, sum_x2 = 0.0;
  double sum_y = 0.0, sum_y2 = 0.0;
  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    const double yi = y[i];
    sum_x += xi;
    sum_x2 += xi * xi;
    sum_y += yi;
    sum_y2 += yi * yi;
  }

  const double Sxy_pair = sum_cross_abs_prod_ptr(x, y, n);
  const double Sxx_pair =
    2.0 * (static_cast<double>(n) * sum_x2 - sum_x * sum_x);
  const double Syy_pair =
    2.0 * (static_cast<double>(n) * sum_y2 - sum_y * sum_y);

  const double row_dot_xy = dot_ptr(Rx_buf.data(), Ry_buf.data(), n);
  const double row_dot_xx = dot_ptr(Rx_buf.data(), Rx_buf.data(), n);
  const double row_dot_yy = dot_ptr(Ry_buf.data(), Ry_buf.data(), n);

  const double XY = u_centered_cov_stat(Sxy_pair, row_dot_xy, Sx, Sy, n);
  const double X2 = u_centered_cov_stat(Sxx_pair, row_dot_xx, Sx, Sx, n);
  const double Y2 = u_centered_cov_stat(Syy_pair, row_dot_yy, Sy, Sy, n);
  return finalize_dcor(XY, X2, Y2);
}

// Dispatch policy:
// - small n uses exact O(n^2) path (lower constant overhead)
// - otherwise use fast O(n log n) path with exact fallback if needed
inline double ustat_dcor_dispatch_ptr(const double* x, const double* y, const int n) {
  if (n < 64) return ustat_dcor_quadratic_ptr(x, y, n);
  const double fast = ustat_dcor_fast_ptr(x, y, n);
  if (std::isfinite(fast)) return fast;
  return ustat_dcor_quadratic_ptr(x, y, n);
}

// Pairwise unbiased distance correlation (U-statistic)
// Székely, Rizzo & Bakirov, 2007
// [[Rcpp::export]]
double ustat_dcor(const arma::vec& x, const arma::vec& y) {
  const int n = static_cast<int>(x.n_elem);
  return ustat_dcor_dispatch_ptr(x.memptr(), y.memptr(), n);
}

// Full matrix of unbiased distance correlations:
// precompute per-column statistics once, then evaluate upper-triangle pairs.
// [[Rcpp::export]]
arma::mat ustat_dcor_matrix_cpp(const arma::mat& X) {
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);
  if (n < 4) Rcpp::stop("Sample size must be at least 4 for unbiased dCor");
  arma::mat R(p, p, arma::fill::eye);

  // Small n: keep simple exact quadratic kernel.
  if (n < 64) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int j = 1; j < p; ++j) {
      for (int i = 0; i < j; ++i) {
        const double* xi = X.colptr(i);
        const double* xj = X.colptr(j);
        const double d = ustat_dcor_quadratic_ptr(xi, xj, n);
        R(i, j) = d;
        R(j, i) = d;
      }
    }
    return R;
  }

  // Memory-guard path for very large n*p: avoid materialising full n x p row_sums.
  if (static_cast<double>(n) * static_cast<double>(p) > 5e7) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int j = 1; j < p; ++j) {
      for (int i = 0; i < j; ++i) {
        const double* xi = X.colptr(i);
        const double* xj = X.colptr(j);
        const double d = ustat_dcor_dispatch_ptr(xi, xj, n);
        R(i, j) = d;
        R(j, i) = d;
      }
    }
    return R;
  }

  arma::mat row_sums(n, p, arma::fill::zeros);
  arma::vec totals(p, arma::fill::zeros);
  arma::vec self_uvar(p, arma::fill::zeros);

  // Precompute per-column row sums, totals, and self U-variance terms.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    const double* xj = X.colptr(j);
    double* rj = row_sums.colptr(j);
    double Sj = 0.0;
    compute_row_sums_and_total_fast_ptr(xj, n, rj, Sj);
    totals[j] = Sj;

    double sum_x = 0.0;
    double sum_x2 = 0.0;
    for (int k = 0; k < n; ++k) {
      const double v = xj[k];
      sum_x += v;
      sum_x2 += v * v;
    }
    const double Sxx_pair =
      2.0 * (static_cast<double>(n) * sum_x2 - sum_x * sum_x);
    const double row_dot = dot_ptr(rj, rj, n);
    self_uvar[j] = u_centered_cov_stat(Sxx_pair, row_dot, Sj, Sj, n);
  }

  // Parallelize over upper-triangular column pairs (uniform work).
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 1; j < p; ++j) {
    const double* xj = X.colptr(j);
    const double* rj = row_sums.colptr(j);
    const double Sj = totals[j];
    const double Y2 = self_uvar[j];
    for (int i = 0; i < j; ++i) {
      const double* xi = X.colptr(i);
      const double* ri = row_sums.colptr(i);
      const double Si = totals[i];
      const double X2 = self_uvar[i];
      double d = NA_REAL;
      if (std::isfinite(X2) && std::isfinite(Y2) && X2 > 0.0 && Y2 > 0.0) {
        const double Sxy_pair = sum_cross_abs_prod_ptr(xi, xj, n);
        const double row_dot_xy = dot_ptr(ri, rj, n);
        const double XY = u_centered_cov_stat(Sxy_pair, row_dot_xy, Si, Sj, n);
        d = finalize_dcor(XY, X2, Y2);
        if (!std::isfinite(d)) d = ustat_dcor_quadratic_ptr(xi, xj, n);
      }
      R(i, j) = d;
      R(j, i) = d;
    }
  }
  return R;
}

// Pairwise matrix path with optional t-test inference.
// [[Rcpp::export]]
Rcpp::List ustat_dcor_matrix_pairwise_cpp(
  const arma::mat& X,
  const bool return_inference = false
) {
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);
  if (n < 4) Rcpp::stop("Sample size must be at least 4 for unbiased dCor");

  arma::mat est(p, p, arma::fill::eye);
  arma::mat n_complete(p, p, arma::fill::value(static_cast<double>(n)));

  if (!return_inference) {
    est = ustat_dcor_matrix_cpp(X);
    return Rcpp::List::create(
      Rcpp::Named("est") = est,
      Rcpp::Named("n_complete") = n_complete
    );
  }

  arma::mat estimate(p, p, arma::fill::eye);
  arma::mat statistic(p, p, arma::fill::value(NA_REAL));
  arma::mat parameter(p, p, arma::fill::value(NA_REAL));
  arma::mat p_value(p, p, arma::fill::value(NA_REAL));

  const double df_value =
    static_cast<double>(n) * static_cast<double>(n - 3) / 2.0 - 1.0;
  for (int j = 0; j < p; ++j) {
    parameter(j, j) = NA_REAL;
  }

  if (n < 64) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int j = 1; j < p; ++j) {
      for (int i = 0; i < j; ++i) {
        const double* xi = X.colptr(i);
        const double* xj = X.colptr(j);
        const double raw = ustat_dcor_quadratic_signed_ptr(xi, xj, n);
        const double clipped = clip_dcor_estimate(raw);
        double tstat = NA_REAL;
        double pval = NA_REAL;
        double df = df_value;
        dcor_t_test_from_signed(raw, n, tstat, df, pval);

        est(i, j) = clipped;
        est(j, i) = clipped;
        estimate(i, j) = raw;
        estimate(j, i) = raw;
        statistic(i, j) = tstat;
        statistic(j, i) = tstat;
        parameter(i, j) = df;
        parameter(j, i) = df;
        p_value(i, j) = pval;
        p_value(j, i) = pval;
      }
    }

    return Rcpp::List::create(
      Rcpp::Named("est") = est,
      Rcpp::Named("n_complete") = n_complete,
      Rcpp::Named("estimate") = estimate,
      Rcpp::Named("statistic") = statistic,
      Rcpp::Named("parameter") = parameter,
      Rcpp::Named("p_value") = p_value
    );
  }

  arma::mat row_sums(n, p, arma::fill::zeros);
  arma::vec totals(p, arma::fill::zeros);
  arma::vec self_uvar(p, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    const double* xj = X.colptr(j);
    double* rj = row_sums.colptr(j);
    double Sj = 0.0;
    compute_row_sums_and_total_fast_ptr(xj, n, rj, Sj);
    totals[j] = Sj;

    double sum_x = 0.0;
    double sum_x2 = 0.0;
    for (int k = 0; k < n; ++k) {
      const double v = xj[k];
      sum_x += v;
      sum_x2 += v * v;
    }
    const double Sxx_pair =
      2.0 * (static_cast<double>(n) * sum_x2 - sum_x * sum_x);
    const double row_dot = dot_ptr(rj, rj, n);
    self_uvar[j] = u_centered_cov_stat(Sxx_pair, row_dot, Sj, Sj, n);
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 1; j < p; ++j) {
    const double* xj = X.colptr(j);
    const double* rj = row_sums.colptr(j);
    const double Sj = totals[j];
    const double Y2 = self_uvar[j];
    for (int i = 0; i < j; ++i) {
      const double* xi = X.colptr(i);
      const double* ri = row_sums.colptr(i);
      const double Si = totals[i];
      const double X2 = self_uvar[i];

      double raw = NA_REAL;
      if (std::isfinite(X2) && std::isfinite(Y2) && X2 > 0.0 && Y2 > 0.0) {
        const double Sxy_pair = sum_cross_abs_prod_ptr(xi, xj, n);
        const double row_dot_xy = dot_ptr(ri, rj, n);
        const double XY = u_centered_cov_stat(Sxy_pair, row_dot_xy, Si, Sj, n);
        raw = dcor_signed_ratio(XY, X2, Y2);
        if (!std::isfinite(raw)) raw = ustat_dcor_quadratic_signed_ptr(xi, xj, n);
      }

      const double clipped = clip_dcor_estimate(raw);
      double tstat = NA_REAL;
      double pval = NA_REAL;
      double df = df_value;
      dcor_t_test_from_signed(raw, n, tstat, df, pval);

      est(i, j) = clipped;
      est(j, i) = clipped;
      estimate(i, j) = raw;
      estimate(j, i) = raw;
      statistic(i, j) = tstat;
      statistic(j, i) = tstat;
      parameter(i, j) = df;
      parameter(j, i) = df;
      p_value(i, j) = pval;
      p_value(j, i) = pval;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("est") = est,
    Rcpp::Named("n_complete") = n_complete,
    Rcpp::Named("estimate") = estimate,
    Rcpp::Named("statistic") = statistic,
    Rcpp::Named("parameter") = parameter,
    Rcpp::Named("p_value") = p_value
  );
}
