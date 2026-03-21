// kendall_corr.cpp
// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "matrixCorr_detail.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using matrixCorr_detail::order_stats::getMs_ll;
using matrixCorr_detail::order_stats::insertion_sort_range;
using matrixCorr_detail::order_stats::inv_count_inplace;
using matrixCorr_detail::order_stats::tau_two_vectors_fast;

namespace {

struct KendallColumnOrder {
  std::vector<int> ord;
  std::vector<std::pair<int, int>> tie_runs;
  long long m1 = 0;
};

inline KendallColumnOrder build_kendall_column_order(const std::vector<long long>& x) {
  const int n = static_cast<int>(x.size());
  KendallColumnOrder state;
  state.ord.resize(n);
  std::iota(state.ord.begin(), state.ord.end(), 0);
  std::sort(state.ord.begin(), state.ord.end(),
            [&](int a, int b) { return x[a] < x[b]; });

  state.tie_runs.reserve(64);
  for (int s = 0; s < n; ) {
    int e = s + 1;
    const long long x0 = x[state.ord[s]];
    while (e < n && x[state.ord[e]] == x0) ++e;
    const int len = e - s;
    if (len > 1) {
      state.tie_runs.emplace_back(s, e);
      state.m1 += 1LL * len * (len - 1) / 2;
    }
    s = e;
  }
  return state;
}

inline double kendall_tau_from_order(const std::vector<long long>& y,
                                     const KendallColumnOrder& state) {
  const int n = static_cast<int>(state.ord.size());
  thread_local std::vector<long long> ybuf, mrg;
  if (static_cast<int>(ybuf.capacity()) < n) {
    ybuf.reserve(n);
    mrg.reserve(n);
  }
  ybuf.resize(n);
  mrg.resize(n);

  for (int k = 0; k < n; ++k) ybuf[k] = y[state.ord[k]];

  long long s_acc = 0;
  for (const auto& run : state.tie_runs) {
    const int s = run.first;
    const int e = run.second;
    const int len = e - s;
    if (len <= 32) insertion_sort_range(ybuf.data(), s, e);
    else           std::sort(ybuf.begin() + s, ybuf.begin() + e);
    s_acc += getMs_ll(ybuf.data() + s, len);
  }

  const long long inv = inv_count_inplace<long long>(ybuf.data(), mrg.data(), n);
  const long long m2  = getMs_ll(ybuf.data(), n);
  const long long n0  = 1LL * n * (n - 1) / 2LL;
  const double s = double(n0) - double(state.m1) - double(m2)
    - 2.0 * double(inv) + double(s_acc);
  const double den1 = double(n0 - state.m1);
  const double den2 = double(n0 - m2);
  if (den1 <= 0.0 || den2 <= 0.0) return NA_REAL;
  return s / std::sqrt(den1 * den2);
}

} // namespace

// ---------------------------- Matrix version ---------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix kendall_matrix_cpp(Rcpp::NumericMatrix mat){
  const int n = mat.nrow();
  const int p = mat.ncol();

  // --- Scalar fast path (p == 2): raw doubles, no discretisation
  if (p == 2) {
    const double* x = &mat(0, 0);
    const double* y = &mat(0, 1);
    const double tau = tau_two_vectors_fast(x, y, n);
    Rcpp::NumericMatrix out(2, 2);
    out(0,0) = out(1,1) = 1.0;
    out(0,1) = out(1,0) = tau;
    return out;
  }

  // --- Discretise once per column for p >= 3 (matrix path)
  std::vector< std::vector<long long> > cols(p, std::vector<long long>(n));
  for (int j = 0; j < p; ++j) {
    const double* cj = &mat(0, j);
    double max_abs = 0.0;
    for (int i = 0; i < n; ++i)
      max_abs = std::max(max_abs, std::abs(cj[i]));

    double scale = 1e8;
    if (max_abs > 0.0) {
      const double max_ll = static_cast<double>(std::numeric_limits<long long>::max()) / 4.0;
      double safe_scale = std::floor(max_ll / max_abs);
      if (!std::isfinite(safe_scale) || safe_scale < 1.0) safe_scale = 1.0;
      scale = std::min(1e8, safe_scale);
    }

    for (int i = 0; i < n; ++i)
      cols[j][i] = static_cast<long long>(std::floor(cj[i] * scale));
  }

  Rcpp::NumericMatrix out(p, p);
  for (int j = 0; j < p; ++j) out(j,j) = 1.0;

  auto fill_row = [&](int i) {
    const KendallColumnOrder state = build_kendall_column_order(cols[i]);
    for (int j = i + 1; j < p; ++j) {
      const double tau = kendall_tau_from_order(cols[j], state);
      out(i, j) = out(j, i) = tau;
    }
  };

  // --- Avoid OpenMP overhead when p is tiny
#if defined(_OPENMP)
  if (p >= 3) {
#pragma omp parallel for schedule(static)
    for (int i = 0; i < p - 1; ++i) {
      fill_row(i);
    }
  } else
#endif
{
  // single-threaded path
  for (int i = 0; i < p - 1; ++i) {
    fill_row(i);
  }
}

  return out;
}

// ---------------------------- Two-vector wrapper ------------------------------
// [[Rcpp::export]]
double kendall_tau2_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  const R_xlen_t n = x.size();
  if (n != y.size() || n < 2) return NA_REAL;
  const double* px = x.begin();
  const double* py = y.begin();
  return matrixCorr_detail::order_stats::tau_two_vectors_fast(px, py, static_cast<int>(n));
}


// [[Rcpp::export]]
double kendall_tau2_from_mat_cpp(Rcpp::NumericMatrix mat) {
  if (mat.ncol() != 2 || mat.nrow() < 2) return NA_REAL;
  const int n = mat.nrow();
  const double* x = &mat(0, 0);
  const double* y = &mat(0, 1);
  return matrixCorr_detail::order_stats::tau_two_vectors_fast(x, y, n);
}
