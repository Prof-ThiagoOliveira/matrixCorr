// Thiago de Paula Oliveira
// src-bicor.cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <limits>     // for std::numeric_limits
#include <cmath>      // for std::floor, std::ceil
#include <algorithm>  // for std::max
#ifdef _OPENMP
#include <omp.h>
#endif

// only what we use from matrixCorr_detail
#include "matrixCorr_detail.h"
using namespace matrixCorr_detail;
using matrixCorr_detail::standardise_bicor::standardise_bicor_column;
using matrixCorr_detail::standardise_bicor::standardise_bicor_column_weighted;

using namespace Rcpp;
using namespace arma;

// Biweight mid-correlation matrix
//
// Computes the biweight mid-correlation for all column pairs of X.
// No missing/Inf values are allowed in X for this fast variant.
// Set `pearson_fallback=1` (individual) to obtain the hybrid correlation when
// MAD==0 in a column.
// X Numeric matrix (rows = observations, cols = variables).
// c_const Tuning constant multiplying raw MAD (default 9.0).
// maxPOutliers Maximum proportion of low *and* high outliers to permit (>0 to cap).
// If < 1, columns are side-rescaled so that these quantiles (if beyond |u>=1) map to |u|=1.
// pearson_fallback 0 = never (returns NA for MAD==0); 1 = individual fallback when needed;
// 2 = force Pearson for all columns (i.e., ordinary Pearson).
// n_threads Number of OpenMP threads (>=1).
// Symmetric p x p matrix of correlations. Diagonals are 1 where defined, NA if
// the column is degenerate.
// [[Rcpp::export]]
arma::mat bicor_matrix_cpp(const arma::mat &X,
                           const double c_const = 9.0,
                           const double maxPOutliers = 1.0,
                           const int pearson_fallback = 1,
                           const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");
  if (!X.is_finite()) stop("X contains NA/NaN/Inf; please handle missingness upstream.");

  // Standardised columns
  arma::mat Z(n, p, fill::zeros);
  std::vector<bool> col_valid(p, false);

  // Standardise each column in parallel
#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    bool ok = false;
    arma::vec zcol(n, fill::zeros);
    standardise_bicor_column(X.col(j), zcol, pearson_fallback, c_const, maxPOutliers, ok);
    Z.col(j) = zcol;
    col_valid[static_cast<std::size_t>(j)] = ok;
  }

  // Correlation matrix R = Z'Z
  arma::mat R(p, p, fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {
      double val = arma::dot(Z.col(j), Z.col(k));
      // Clamp to [-1, 1] if finite
      if (std::isfinite(val)) {
        if (val > 1.0) val = 1.0;
        else if (val < -1.0) val = -1.0;
      }
      R(j, k) = val;
      if (k != static_cast<std::size_t>(j)) R(k, j) = val;
    }
  }

  // Mark invalid columns as NA, others keep unit diagonal
  for (std::size_t j = 0; j < p; ++j) {
    if (!col_valid[j]) {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
    }
    R(j, j) = 1.0;
  }
  return R;
}

// bicor for two vectors (hybrid if one side falls back to Pearson).
// [[Rcpp::export]]
double bicor_vec_cpp(const arma::vec &x, const arma::vec &y,
                     const double c_const = 9.0,
                     const double maxPOutliers = 1.0,
                     const int pearson_fallback = 1) {
  if (x.n_elem != y.n_elem) stop("x and y must have the same length.");
  if (!x.is_finite() || !y.is_finite()) stop("x or y contains NA/NaN/Inf.");

  arma::vec zx(x.n_elem, fill::zeros), zy(y.n_elem, fill::zeros);
  bool okx=false, oky=false;
  standardise_bicor_column(x, zx, pearson_fallback, c_const, maxPOutliers, okx);
  standardise_bicor_column(y, zy, pearson_fallback, c_const, maxPOutliers, oky);

  double val = arma::dot(zx, zy);
  if (std::isfinite(val)) {
    if (val > 1.0) val = 1.0;
    else if (val < -1.0) val = -1.0;
  }
  // If either side is degenerate, return NA
  if (!okx || !oky) return std::numeric_limits<double>::quiet_NaN();
  return val;
}

// [[Rcpp::export]]
arma::mat bicor_matrix_pairwise_cpp(const arma::mat &X,
                                    const double c_const = 9.0,
                                    const double maxPOutliers = 1.0,
                                    const int pearson_fallback = 1,
                                    const int min_n = 5,
                                    const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");

  // initialise with NaN (constructor can't take 'datum::nan')
  arma::mat R(p, p, arma::fill::none);
  R.fill(arma::datum::nan);

  std::vector<arma::uvec> finite_idx(p);
  for (std::size_t j = 0; j < p; ++j) {
    finite_idx[j] = arma::find_finite(X.col(j));
  }

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    const arma::uvec& idx_j = finite_idx[static_cast<std::size_t>(j)];
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {
      if (k == static_cast<std::size_t>(j)) continue;  // handle diagonal later

      const arma::uvec& idx_k = finite_idx[k];
      const std::size_t possible = std::min(idx_j.n_elem, idx_k.n_elem);
      if (possible < static_cast<std::size_t>(min_n)) continue;

      static thread_local std::vector<arma::uword> overlap_idx;
      overlap_idx.clear();
      overlap_idx.reserve(possible);

      arma::uword ia = 0, ib = 0;
      while (ia < idx_j.n_elem && ib < idx_k.n_elem) {
        const arma::uword a_val = idx_j[ia];
        const arma::uword b_val = idx_k[ib];
        if (a_val == b_val) {
          overlap_idx.push_back(a_val);
          ++ia; ++ib;
        } else if (a_val < b_val) {
          ++ia;
        } else {
          ++ib;
        }
      }

      const std::size_t overlap_n = overlap_idx.size();
      if (overlap_n < static_cast<std::size_t>(min_n)) continue;

      static thread_local arma::vec xbuf;
      static thread_local arma::vec ybuf;
      xbuf.set_size(overlap_n);
      ybuf.set_size(overlap_n);

      const double* colj_ptr = X.colptr(static_cast<arma::uword>(j));
      const double* colk_ptr = X.colptr(static_cast<arma::uword>(k));
      for (std::size_t t = 0; t < overlap_n; ++t) {
        const arma::uword row = overlap_idx[t];
        xbuf[t] = colj_ptr[row];
        ybuf[t] = colk_ptr[row];
      }

      static thread_local arma::vec zj_buf;
      static thread_local arma::vec zk_buf;
      bool okj = false, okk = false;
      standardise_bicor_column(xbuf, zj_buf, pearson_fallback, c_const, maxPOutliers, okj);
      standardise_bicor_column(ybuf, zk_buf, pearson_fallback, c_const, maxPOutliers, okk);

      double val = arma::datum::nan;
      if (okj && okk) {
        val = arma::dot(zj_buf, zk_buf);
        if (std::isfinite(val)) {
          if (val > 1.0) val = 1.0;
          else if (val < -1.0) val = -1.0;
        }
      }
      R(j, k) = val;
      R(k, j) = val;
    }
  }

  for (std::size_t j = 0; j < p; ++j) {
    if (finite_idx[j].n_elem < 2) {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
    }
    R(j, j) = 1.0;
  }
  return R;
}


// [[Rcpp::export]]
arma::mat bicor_matrix_weighted_cpp(const arma::mat &X,
                                    const arma::vec &w,
                                    const double c_const = 9.0,
                                    const double maxPOutliers = 1.0,
                                    const int pearson_fallback = 1,
                                    const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");
  if (w.n_elem != n) stop("Length of weights `w` must match nrow(X).");
  if (!X.is_finite()) stop("X contains NA/NaN/Inf; use pairwise kernel or clean data.");
  if (!w.is_finite() || arma::any(w < 0)) stop("Weights must be finite and non-negative.");

  arma::mat Z(n, p, arma::fill::zeros);
  std::vector<bool> col_valid(p, false);

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    bool ok = false;
    arma::vec zcol(n, arma::fill::zeros);
    standardise_bicor_column_weighted(X.col(j), w, zcol, pearson_fallback, c_const, maxPOutliers, ok);
    Z.col(j) = zcol;
    col_valid[static_cast<std::size_t>(j)] = ok;
  }

  arma::mat R(p, p, arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {
      double val = arma::dot(Z.col(j), Z.col(k));
      if (std::isfinite(val)) {
        if (val > 1.0) val = 1.0;
        else if (val < -1.0) val = -1.0;
      }
      R(j, k) = val;
      if (k != static_cast<std::size_t>(j)) R(k, j) = val;
    }
  }
  for (std::size_t j = 0; j < p; ++j) {
    if (col_valid[j]) {
      R(j, j) = 1.0;
    } else {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
      R(j, j) = arma::datum::nan;
    }
  }
  return R;
}

// [[Rcpp::export]]
arma::mat bicor_matrix_weighted_pairwise_cpp(const arma::mat &X,
                                             const arma::vec &w,
                                             const double c_const = 9.0,
                                             const double maxPOutliers = 1.0,
                                             const int pearson_fallback = 1,
                                             const int min_n = 5,
                                             const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");
  if (w.n_elem != n)   stop("Length of weights `w` must match nrow(X).");
  if (!w.is_finite() || arma::any(w < 0)) stop("Weights must be finite and non-negative.");

  arma::mat R(p, p, arma::fill::none);
  R.fill(arma::datum::nan);

  std::vector<arma::uvec> finite_idx(p);
  for (std::size_t j = 0; j < p; ++j) {
    finite_idx[j] = arma::find_finite(X.col(j));
  }

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    const arma::uvec& idx_j = finite_idx[static_cast<std::size_t>(j)];
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {
      if (k == static_cast<std::size_t>(j)) continue;

      const arma::uvec& idx_k = finite_idx[k];
      const std::size_t possible = std::min(idx_j.n_elem, idx_k.n_elem);
      if (possible < static_cast<std::size_t>(min_n)) continue;

      static thread_local std::vector<arma::uword> overlap_idx;
      overlap_idx.clear();
      overlap_idx.reserve(possible);

      arma::uword ia = 0, ib = 0;
      while (ia < idx_j.n_elem && ib < idx_k.n_elem) {
        const arma::uword a_val = idx_j[ia];
        const arma::uword b_val = idx_k[ib];
        if (a_val == b_val) {
          overlap_idx.push_back(a_val);
          ++ia; ++ib;
        } else if (a_val < b_val) {
          ++ia;
        } else {
          ++ib;
        }
      }

      const std::size_t overlap_n = overlap_idx.size();
      if (overlap_n < static_cast<std::size_t>(min_n)) continue;

      static thread_local arma::vec xbuf;
      static thread_local arma::vec ybuf;
      static thread_local arma::vec wbuf;
      xbuf.set_size(overlap_n);
      ybuf.set_size(overlap_n);
      wbuf.set_size(overlap_n);

      const double* colj_ptr = X.colptr(static_cast<arma::uword>(j));
      const double* colk_ptr = X.colptr(static_cast<arma::uword>(k));
      const double* w_ptr    = w.memptr();
      for (std::size_t t = 0; t < overlap_n; ++t) {
        const arma::uword row = overlap_idx[t];
        xbuf[t] = colj_ptr[row];
        ybuf[t] = colk_ptr[row];
        wbuf[t] = w_ptr[row];
      }

      static thread_local arma::vec zj_buf;
      static thread_local arma::vec zk_buf;
      bool okj = false, okk = false;
      standardise_bicor_column_weighted(xbuf, wbuf, zj_buf, pearson_fallback, c_const, maxPOutliers, okj);
      standardise_bicor_column_weighted(ybuf, wbuf, zk_buf, pearson_fallback, c_const, maxPOutliers, okk);

      double val = arma::datum::nan;
      if (okj && okk) {
        val = arma::dot(zj_buf, zk_buf);
        if (std::isfinite(val)) {
          if (val > 1.0) val = 1.0;
          else if (val < -1.0) val = -1.0;
        }
      }
      R(j, k) = val;
      R(k, j) = val;
    }
  }

  for (std::size_t j = 0; j < p; ++j) {
    if (finite_idx[j].n_elem < 2) {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
    }
    R(j, j) = 1.0;
  }
  return R;
}
