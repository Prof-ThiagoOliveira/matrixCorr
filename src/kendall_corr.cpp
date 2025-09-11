// kendall_corr.cpp
// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include "matrixCorr_detail.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace matrixCorr_detail;
using matrixCorr_detail::order_stats::discretise;
using matrixCorr_detail::order_stats::compress_ranks_and_T;
using matrixCorr_detail::order_stats::inversion_count;

// ---------- core scalar metric (armadillo columns) ----------
static inline double kendall_tau_auto_arma(const arma::vec& x,
                                           const arma::vec& y,
                                           const double scale)
{
  const int n = static_cast<int>(x.n_elem);
  if (n != static_cast<int>(y.n_elem) || n < 2) return NA_REAL;

  thread_local std::vector<int>          tl_ord;
  thread_local std::vector<long long>    tl_xv, tl_ys, tl_buf;

  tl_xv.resize(n);
  tl_ys.resize(n);
  tl_ord.resize(n);

  // Discretise (branch-free)
  for (int i = 0; i < n; ++i) {
    tl_xv[i] = (long long) std::floor(x[i] * scale);
    tl_ys[i] = (long long) std::floor(y[i] * scale);
    tl_ord[i] = i;
  }

  // Sort indices by x
  std::sort(tl_ord.begin(), tl_ord.end(),
            [&](int a, int b){ return tl_xv[a] < tl_xv[b]; });

  // Reorder y in lockstep with x
            {
              thread_local std::vector<long long> tmp_y;
              tmp_y.resize(n);
              for (int i = 0; i < n; ++i) tmp_y[i] = tl_ys[tl_ord[i]];
              tl_ys.swap(tmp_y);
            }

            // Helper: Ms on sorted block
            auto getMs_ll = [](const long long* dat, int len) noexcept -> long long {
              long long Ms = 0, tie = 0;
              for (int i = 1; i < len; ++i) {
                if (dat[i] == dat[i-1]) ++tie;
                else if (tie) { Ms += tie * (tie + 1) / 2; tie = 0; }
              }
              if (tie) Ms += tie * (tie + 1) / 2;
              return Ms;
            };

            // m1 and s_acc; detect blocks via xv[ord[i]] without making xs[]
            long long m1 = 0, s_acc = 0;
            for (int i = 0; i < n; ) {
              int j = i + 1;
              const long long xi = tl_xv[tl_ord[i]];
              while (j < n && tl_xv[tl_ord[j]] == xi) ++j;
              const int L = j - i;
              if (L > 1) {
                m1 += 1LL * L * (L - 1) / 2;
                std::sort(tl_ys.begin() + i, tl_ys.begin() + j);
                s_acc += getMs_ll(tl_ys.data() + i, L);
              }
              i = j;
            }

            // In-place inversion count on tl_ys, with external tl_buf
            auto inversion_count_ll = [&](std::vector<long long>& a) noexcept -> long long {
              const int N = (int)a.size();
              if (N < 2) return 0;
              if ((int)tl_buf.size() < N) tl_buf.resize(N);
              long long inv = 0;

              auto insertion_count = [&](int L, int R) noexcept {
                long long v_inv = 0;
                for (int i = L + 1; i <= R; ++i) {
                  long long v = a[i]; int j = i - 1;
                  while (j >= L && a[j] > v) { a[j + 1] = a[j]; --j; ++v_inv; }
                  a[j + 1] = v;
                }
                return v_inv;
              };

              const int TH = 32;
              for (int L = 0; L < N; L += TH) {
                const int R = std::min(L + TH - 1, N - 1);
                inv += insertion_count(L, R);
              }
              for (int width = TH; width < N; width <<= 1) {
                for (int L = 0; L < N; L += 2 * width) {
                  const int M = std::min(L + width, N);
                  const int R = std::min(L + 2 * width, N);
                  int i = L, j = M, k = L;
                  while (i < M && j < R) {
                    if (a[i] <= a[j]) tl_buf[k++] = a[i++];
                    else { tl_buf[k++] = a[j++]; inv += (M - i); }
                  }
                  while (i < M) tl_buf[k++] = a[i++];
                  while (j < R) tl_buf[k++] = a[j++];
                  for (int t = L; t < R; ++t) a[t] = tl_buf[t];
                }
              }
              return inv;
            };
            // tl_ys becomes globally sorted
            const long long inv = inversion_count_ll(tl_ys);
            // compute ties on sorted tl_ys
            const long long m2  = getMs_ll(tl_ys.data(), n);

            const long long n0 = 1LL * n * (n - 1) / 2LL;
            long double s = (long double)n0;
            s -= (long double)m1 + (long double)m2;
            s -= 2.0L * (long double)inv;
            s += (long double)s_acc;

            const long double den1 = (long double)(n0 - m1);
            const long double den2 = (long double)(n0 - m2);
            if (den1 <= 0.0L || den2 <= 0.0L) return NA_REAL;

            return (double)(s / std::sqrt(den1 * den2));
}


// tau-a: inversion-count formula; computes tau-a even if ties exist
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
  return kendall_tau_auto_arma(x, y, scale); // auto covers ties -> tau-b
}


// [[Rcpp::export]]
double kendall_tau_auto_cpp(Rcpp::NumericVector xr, Rcpp::NumericVector yr, double scale = 1e8){
  arma::vec x(xr.begin(), xr.size(), /*copy_aux_mem*/ false, /*strict*/ true);
  arma::vec y(yr.begin(), yr.size(), false, true);
  return kendall_tau_auto_arma(x, y, scale);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix kendall_matrix_cpp(Rcpp::NumericMatrix mat){
  const int n = mat.nrow();
  const int p = mat.ncol();
  const double scale = 1e8;

  // Column-wise pre-discretisation (stored as long long per column)
  std::vector< std::vector<long long> > cols(p, std::vector<long long>(n));
  for (int j = 0; j < p; ++j) {
    const double* cj = &mat(0, j);
    for (int i = 0; i < n; ++i)
      cols[j][i] = (long long) std::floor(cj[i] * scale);
  }

  // Helper to compute tau-b from two discretised columns (Knight path)
  auto tau_from_discretised = [&](const std::vector<long long>& x,
                                  const std::vector<long long>& y) -> double {
                                    const int n = static_cast<int>(x.size());
                                    if (static_cast<int>(y.size()) != n || n < 2) return NA_REAL;

                                    // Per-thread reusable workspaces to avoid allocations in hot path
                                    thread_local std::vector<int>           ord;
                                    thread_local std::vector<long long>     ybuf;   // y reordered by x
                                    thread_local std::vector<long long>     mrg;    // merge buffer

                                    ord.resize(n);
                                    ybuf.resize(n);
                                    mrg.resize(n);

                                    // order by x (indices only)
                                    std::iota(ord.begin(), ord.end(), 0);
                                    std::sort(ord.begin(), ord.end(),
                                              [&](int a, int b){ return x[a] < x[b]; });

                                    // y in x-sorted order (into workspace buffer)
                                    for (int i = 0; i < n; ++i) ybuf[i] = y[ord[i]];

                                    // inline, non-allocating inversion counter; sorts 'a' into nondecreasing order
                                    auto inversion_count_inplace = [&](long long* a, long long* buf, int N) -> long long {
                                      if (N < 2) return 0LL;
                                      long long inv = 0;

                                      // small-block insertion sort, then bottom-up merges
                                      constexpr int TH = 96; // tune if needed
                                      auto insertion_block = [&](int L, int R) {
                                        long long v_inv = 0;
                                        for (int i = L + 1; i <= R; ++i) {
                                          long long v = a[i]; int j = i - 1;
                                          while (j >= L && a[j] > v) { a[j + 1] = a[j]; --j; ++v_inv; }
                                          a[j + 1] = v;
                                        }
                                        return v_inv;
                                      };

                                      // sort small blocks
                                      for (int L = 0; L < N; L += TH) {
                                        const int R = std::min(L + TH - 1, N - 1);
                                        inv += insertion_block(L, R);
                                      }

                                      // bottom-up merge, counting cross inversions
                                      for (int width = TH; width < N; width <<= 1) {
                                        for (int L = 0; L < N; L += 2 * width) {
                                          const int M = std::min(L + width, N);
                                          const int R = std::min(L + 2 * width, N);
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
                                    };

                                    // m1 and s_acc: scan x-runs without materialising xs[]
                                    long long m1 = 0, s_acc = 0;
                                    for (int i = 0; i < n; ) {
                                      int j = i + 1;
                                      const long long xi = x[ord[i]];
                                      while (j < n && x[ord[j]] == xi) ++j;
                                      const int L = j - i;
                                      if (L > 1) {
                                        // ties in x contribute to m1; for s_acc we need Ms within this block by y
                                        m1 += 1LL * L * (L - 1) / 2;
                                        std::sort(ybuf.begin() + i, ybuf.begin() + j);
                                        s_acc += matrixCorr_detail::order_stats::getMs_ll(ybuf.data() + i, L);
                                      }
                                      i = j;
                                    }

                                    // inversion count on ybuf (in place) + m2 on globally sorted ybuf
                                    const long long inv = inversion_count_inplace(ybuf.data(), mrg.data(), n);
                                    const long long m2  = matrixCorr_detail::order_stats::getMs_ll(ybuf.data(), n);

                                    // assemble tau-b
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

  Rcpp::NumericMatrix out(p, p);
  // Diagonal and symmetry
  for (int j = 0; j < p; ++j) out(j,j) = 1.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      const double tau = tau_from_discretised(cols[i], cols[j]);
      out(i,j) = out(j,i) = tau;
    }
  }
  return out;
}
