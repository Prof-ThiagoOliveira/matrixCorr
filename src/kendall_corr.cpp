// kendall_corr.cpp
// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include "matrixCorr_detail.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace matrixCorr_detail;
using matrixCorr_detail::order_stats::discretise;
using matrixCorr_detail::order_stats::has_ties;
using matrixCorr_detail::order_stats::sort_count;
using matrixCorr_detail::order_stats::Fenwick;
using matrixCorr_detail::order_stats::compress_1based;

// ---------- core scalar metric (armadillo columns) ----------
static inline double kendall_tau_auto_arma(const arma::vec& x, const arma::vec& y, double scale){
  const arma::uword n = x.n_elem;
  if (n != y.n_elem || n < 2) return NA_REAL;

  // discretise to integers (stable equality)
  std::vector<long long> xi = discretise(x, scale);
  std::vector<long long> yi = discretise(y, scale);

  const bool ties_x = has_ties(xi);
  const bool ties_y = has_ties(yi);
  const double n0 = double(n) * double(n - 1) / 2.0;

  // ---------- tau-a (no ties on either margin) ----------
  if (!ties_x && !ties_y){
    // order by x
    std::vector<int> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return xi[a] < xi[b]; });

    // y in x-sorted order
    std::vector<int> y_ord(n), tmp(n);
    for (arma::uword i = 0; i < n; ++i) y_ord[i] = int(yi[ord[i]]);

    long long discord = sort_count(y_ord, tmp, 0, int(n) - 1);
    return (n0 - 2.0 * double(discord)) / n0;  // tau-a
  }

  // ---------- tau-b (ties present) ----------
  // order by x then y
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b){
    if (xi[a] != xi[b]) return xi[a] < xi[b];
    return yi[a] < yi[b];
  });

  // T_x: pairs tied on x
  long long T_x = 0;
  for (arma::uword i = 0; i < n; ){
    arma::uword j = i + 1;
    while (j < n && xi[ord[j]] == xi[ord[i]]) ++j;
    long long g = long(j - i);
    T_x += g * (g - 1) / 2;
    i = j;
  }

  // T_y: pairs tied on y
  long long T_y = 0;
  {
    std::vector<long long> ys = yi;
    std::sort(ys.begin(), ys.end());
    for (size_t i = 0; i < ys.size(); ){
      size_t j = i + 1;
      while (j < ys.size() && ys[j] == ys[i]) ++j;
      long long g = long(j - i);
      T_y += g * (g - 1) / 2;
      i = j;
    }
  }

  // S = C - D via BIT (exclude ties in numerator)
  std::vector<int> ry = compress_1based(yi);
  Fenwick bit(int(*std::max_element(ry.begin(), ry.end())));
  long long processed = 0, S = 0;

  for (arma::uword i = 0; i < n; ){
    arma::uword j = i + 1;
    while (j < n && xi[ord[j]] == xi[ord[i]]) ++j;

    // query
    for (arma::uword k = i; k < j; ++k){
      const int idx = ord[k];
      const int r = ry[idx];
      long long less = bit.sum(r - 1);
      long long leq  = bit.sum(r);
      long long greater = processed - leq;
      S += (less - greater);
    }
    // insert
    for (arma::uword k = i; k < j; ++k){
      bit.add(ry[ord[k]], 1);
    }
    processed += long(j - i);
    i = j;
  }

  const double denom_x = n0 - double(T_x);
  const double denom_y = n0 - double(T_y);
  if (denom_x <= 0.0 || denom_y <= 0.0) return NA_REAL;

  return double(S) / std::sqrt(denom_x * denom_y); // tau-b
}

// tau-a: inversion-count formula; computes tau-a even if ties exist
// [[Rcpp::export]]
double kendall_tau_a_cpp(Rcpp::NumericVector xr, Rcpp::NumericVector yr, double scale = 1e8){
  const int n = xr.size();
  if (n != yr.size() || n < 2) return NA_REAL;

  arma::vec x(xr.begin(), n, /*copy_aux_mem*/ false, /*strict*/ true);
  arma::vec y(yr.begin(), n, false, true);

  // discretise and order by x
  std::vector<long long> xi = discretise(x, scale);
  std::vector<long long> yi = discretise(y, scale);

  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b){ return xi[a] < xi[b]; });

  std::vector<int> y_ord(n), tmp(n);
  for (int i = 0; i < n; ++i) y_ord[i] = int(yi[ord[i]]);
  long long discord = sort_count(y_ord, tmp, 0, n - 1);

  const double n0 = double(n) * double(n - 1) / 2.0;
  return (n0 - 2.0 * double(discord)) / n0;
}

// tau-b: fast BIT path with tie adjustment
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
  const double scale = 1e8;
  arma::mat X(mat.begin(), mat.nrow(), mat.ncol(), /*copy_aux_mem*/ false, /*strict*/ true);

  auto tau_fun = [scale](const arma::vec& a, const arma::vec& b){
    return kendall_tau_auto_arma(a, b, scale);
  };
  arma::mat M = pairwise::apply(X, tau_fun);

  Rcpp::NumericMatrix out(M.n_rows, M.n_cols);
  std::copy(M.begin(), M.end(), out.begin());
  return out;
}
