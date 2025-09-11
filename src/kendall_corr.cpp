// kendall_corr.cpp
// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "matrixCorr_detail.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace matrixCorr_detail;
using matrixCorr_detail::order_stats::discretise;
using matrixCorr_detail::order_stats::compress_ranks_and_T;
using matrixCorr_detail::order_stats::getMs_ll;
using matrixCorr_detail::order_stats::inversion_count;

// ----------------------------- Utilities ------------------------------------

// Non-allocating inversion counter.
// Sorts `a[0..n)` in-place
// Requires external buffer `buf[0..n)`.
// T should be an integer type (e.g., long long).
template <typename T>
static inline long long inv_count_inplace(T* __restrict__ a,
                                          T* __restrict__ buf,
                                          const int n) noexcept {
  if (n < 2) return 0LL;
  long long inv = 0;

  // Tune this threshold for your data/CPU; 64–128 works well.
  constexpr int TH = 96;

  // Insertion sort small blocks to improve locality before merging.
  auto insertion_block = [&](int L, int R) noexcept {
    long long v_inv = 0;
    for (int i = L + 1; i <= R; ++i) {
      T v = a[i]; int j = i - 1;
      while (j >= L && a[j] > v) { a[j + 1] = a[j]; --j; ++v_inv; }
      a[j + 1] = v;
    }
    return v_inv;
  };

  for (int L = 0; L < n; L += TH) {
    const int R = std::min(L + TH - 1, n - 1);
    inv += insertion_block(L, R);
  }

  // Bottom-up merge; count cross inversions.
  for (int width = TH; width < n; width <<= 1) {
    for (int L = 0; L < n; L += 2 * width) {
      const int M = std::min(L + width, n);
      const int R = std::min(L + 2 * width, n);
      int i = L, j = M, k = L;
      while (i < M && j < R) {
        if (a[i] <= a[j]) buf[k++] = a[i++];
        else { buf[k++] = a[j++]; inv += (M - i); }
      }
      while (i < M) buf[k++] = a[i++];
      while (j < R) buf[k++] = a[j++];
      for (int t = L; t < R; ++t) a[t] = buf[t];
    }
  }
  return inv;
}

// Tiny insertion sort for short ranges (used for x-runs)
template <typename T>
static inline void insertion_sort_range(T* a, int s, int e) noexcept {
  for (int i = s + 1; i < e; ++i) {
    T v = a[i]; int j = i - 1;
    while (j >= s && a[j] > v) { a[j + 1] = a[j]; --j; }
    a[j + 1] = v;
  }
}

// ----------------------- Scalar tau-b (arma columns) -------------------------
static inline double kendall_tau_auto_arma(const arma::vec& x,
                                           const arma::vec& y,
                                           const double scale)
{
  const int n = static_cast<int>(x.n_elem);
  if (n != static_cast<int>(y.n_elem) || n < 2) return NA_REAL;

  thread_local std::vector<int>        tl_ord;
  thread_local std::vector<long long>  tl_xv, tl_ys, tl_buf;

  tl_xv.resize(n);
  tl_ys.resize(n);
  tl_ord.resize(n);
  tl_buf.resize(n);

  // Discretise
  for (int i = 0; i < n; ++i) {
    tl_xv[i] = static_cast<long long>(std::floor(x[i] * scale));
    tl_ys[i] = static_cast<long long>(std::floor(y[i] * scale));
    tl_ord[i] = i;
  }

  // Sort indices by x
  std::sort(tl_ord.begin(), tl_ord.end(),
            [&](int a, int b){ return tl_xv[a] < tl_xv[b]; });

  // y reordered by x (into tl_ys via a temporary to avoid aliasing)
            {
              thread_local std::vector<long long> tmp_y;
              tmp_y.resize(n);
              for (int i = 0; i < n; ++i) tmp_y[i] = tl_ys[tl_ord[i]];
              tl_ys.swap(tmp_y);
            }

            // m1 and s_acc from x-runs; sort y within each run to get per-run Ms
            long long m1 = 0, s_acc = 0;
            for (int i = 0; i < n; ) {
              int j = i + 1;
              const long long xi = tl_xv[tl_ord[i]];
              while (j < n && tl_xv[tl_ord[j]] == xi) ++j;
              const int L = j - i;
              if (L > 1) {
                m1 += 1LL * L * (L - 1) / 2;
                if (L <= 32) insertion_sort_range(tl_ys.data(), i, j);
                else         std::sort(tl_ys.begin() + i, tl_ys.begin() + j);
                s_acc += getMs_ll(tl_ys.data() + i, L);
              }
              i = j;
            }

            // Inversion count (in place) + m2 on globally sorted y
            const long long inv = inv_count_inplace(tl_ys.data(), tl_buf.data(), n);
            const long long m2  = getMs_ll(tl_ys.data(), n);

            const long long n0 = 1LL * n * (n - 1) / 2LL;
            long double s = static_cast<long double>(n0);
            s -= static_cast<long double>(m1) + static_cast<long double>(m2);
            s -= 2.0L * static_cast<long double>(inv);
            s += static_cast<long double>(s_acc);

            const long double den1 = static_cast<long double>(n0 - m1);
            const long double den2 = static_cast<long double>(n0 - m2);
            if (den1 <= 0.0L || den2 <= 0.0L) return NA_REAL;

            return static_cast<double>(s / std::sqrt(den1 * den2));
}

// ------------------------------ Exports --------------------------------------

// tau-a: inversion-count formula (ignores tie adjustment in denominator).
// [[Rcpp::export]]
double kendall_tau_a_cpp(Rcpp::NumericVector xr, Rcpp::NumericVector yr, double scale = 1e8){
  const int n = xr.size();
  if (n != yr.size() || n < 2) return NA_REAL;

  arma::vec x(xr.begin(), n, /*copy_aux_mem*/ false, /*strict*/ true);
  arma::vec y(yr.begin(), n, false, true);

  // discretise
  std::vector<long long> xi = discretise(x, scale);
  std::vector<long long> yi = discretise(y, scale);

  // order by x
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b){ return xi[a] < xi[b]; });

  // rank-compress y to small ints (safe for inversion counter)
  auto cy = compress_ranks_and_T(yi);

  // y in x-sorted order (ranks)
  std::vector<int> y_ord(n);
  for (int i = 0; i < n; ++i) y_ord[i] = cy.rank[ord[i]];

  // count inversions (discordant pairs)
  long long discord = inversion_count(y_ord);

  const double n0 = double(n) * double(n - 1) / 2.0;
  return (n0 - 2.0 * double(discord)) / n0;  // tau-a
}

// tau-b (tie adjustment)
// [[Rcpp::export]]
double kendall_tau_b_cpp(Rcpp::NumericVector xr, Rcpp::NumericVector yr, double scale = 1e8){
  arma::vec x(xr.begin(), xr.size(), /*copy_aux_mem*/ false, /*strict*/ true);
  arma::vec y(yr.begin(), yr.size(), false, true);
  return kendall_tau_auto_arma(x, y, scale);
}

// [[Rcpp::export]]
double kendall_tau_auto_cpp(Rcpp::NumericVector xr, Rcpp::NumericVector yr, double scale = 1e8){
  arma::vec x(xr.begin(), xr.size(), /*copy_aux_mem*/ false, /*strict*/ true);
  arma::vec y(yr.begin(), yr.size(), false, true);
  return kendall_tau_auto_arma(x, y, scale);
}

// ---------------------------- Matrix version ---------------------------------

// [[Rcpp::export]]
Rcpp::NumericMatrix kendall_matrix_cpp(Rcpp::NumericMatrix mat){
  const int n = mat.nrow();
  const int p = mat.ncol();
  const double scale = 1e8;

  // Column-wise pre-discretisation
  std::vector< std::vector<long long> > cols(p, std::vector<long long>(n));
  for (int j = 0; j < p; ++j) {
    const double* cj = &mat(0, j);
    for (int i = 0; i < n; ++i)
      cols[j][i] = static_cast<long long>(std::floor(cj[i] * scale));
  }

  Rcpp::NumericMatrix out(p, p);
  for (int j = 0; j < p; ++j) out(j,j) = 1.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < p; ++i) {
    // Precompute ord for column i (indices sorted by x_i)
    std::vector<int> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(),
              [&](int a, int b){ return cols[i][a] < cols[i][b]; });

    // Precompute x-run boundaries for column i
    std::vector<std::pair<int,int>> xruns;
    xruns.reserve(64);
    for (int s = 0; s < n; ) {
      int e = s + 1;
      const long long xi = cols[i][ord[s]];
      while (e < n && cols[i][ord[e]] == xi) ++e;
      xruns.emplace_back(s, e);
      s = e;
    }

    // Helper that uses precomputed ord + xruns and only touches y
    auto tau_from_ord = [&](const std::vector<long long>& y) -> double {
      thread_local std::vector<long long> ybuf, mrg;
      ybuf.resize(n);
      mrg.resize(n);

      // y in x-sorted order
      for (int k = 0; k < n; ++k) ybuf[k] = y[ord[k]];

      // m1 and s_acc from x-runs (sort y within each run, then Ms)
      long long m1 = 0, s_acc = 0;
      for (const auto& run : xruns) {
        const int s = run.first;
        const int e = run.second;
        const int L = e - s;
        if (L > 1) {
          m1 += 1LL * L * (L - 1) / 2;
          if (L <= 32) insertion_sort_range(ybuf.data(), s, e);
          else         std::sort(ybuf.begin() + s, ybuf.begin() + e);
          s_acc += getMs_ll(ybuf.data() + s, L);
        }
      }

      // inversion count (in place) + m2 on globally sorted y
      const long long inv = inv_count_inplace(ybuf.data(), mrg.data(), n);
      const long long m2  = getMs_ll(ybuf.data(), n);

      const long long n0 = 1LL * n * (n - 1) / 2LL;
      long double s = static_cast<long double>(n0);
      s -= static_cast<long double>(m1) + static_cast<long double>(m2);
      s -= 2.0L * static_cast<long double>(inv);
      s += static_cast<long double>(s_acc);

      const long double den1 = static_cast<long double>(n0 - m1);
      const long double den2 = static_cast<long double>(n0 - m2);
      if (den1 <= 0.0L || den2 <= 0.0L) return NA_REAL;
      return static_cast<double>(s / std::sqrt(den1 * den2));
    };

    for (int j = i + 1; j < p; ++j) {
      const double tau = tau_from_ord(cols[j]);
      out(i,j) = out(j,i) = tau;
    }
  }

  return out;
}
