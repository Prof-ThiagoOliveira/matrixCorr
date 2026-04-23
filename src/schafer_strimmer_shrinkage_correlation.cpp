// Thiago de Paula Oliveira
// cor_shrink_ss.cpp — Schäfer–Strimmer shrinkage correlation

#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include "matrixCorr_detail.h"
#include "matrixCorr_omp.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using matrixCorr_detail::linalg::centered_crossprod_no_copy;
using matrixCorr_detail::linalg::copy_upper_to_lower_inplace;

// X is n x p (no NA). Returns p x p shrinkage correlation matrix.
// [[Rcpp::export]]
arma::mat sss_cor_cpp(SEXP X_) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");

  const arma::uword n = Rf_nrows(X_);
  const arma::uword p = Rf_ncols(X_);
  if (n < 3 || p < 2) Rcpp::stop("Need >= 3 rows and >= 2 columns.");
  const bool use_omp = (p >= 64u && omp_get_max_threads() > 1);

  // No-copy view
  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  // XtX = (X - mu)'(X - mu), computed without a centered copy of X.
  arma::rowvec mu;
  arma::mat XtX = centered_crossprod_no_copy(X, mu);

  // Sample covariance (unbiased)
  XtX /= static_cast<double>(n - 1);

  // Standard deviations; clamp tiny negatives on diagonal
  arma::vec s = arma::sqrt(
    arma::clamp(XtX.diag(), 0.0, std::numeric_limits<double>::infinity())
  );

  // Build sample correlation on XtX (upper triangle)
  arma::vec inv_s = 1.0 / s;
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static) if (use_omp)
  #endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double sj = inv_s[uj];
    if (!std::isfinite(sj) || sj == 0.0) continue;
    for (arma::uword i = 0; i <= uj; ++i) {
      const double si = inv_s[i];
      if (!std::isfinite(si) || si == 0.0) continue;
      XtX(i, uj) *= (si * sj);
    }
  }
  XtX.diag().ones();

  // Identify zero-variance columns to propagate NA later
  arma::uvec zero = arma::find(s == 0.0);

  // ------------------------------------------------------------
    // Estimate analytic shrinkage intensity lambda
  // Var(r_ij) approx ((1 - r_ij^2)^2) / (n - 1)  for i<j
  // ̂lambda = sum Var(r_ij) / sum r_ij^2, clamped to [0, 1]
  // ------------------------------------------------------------
  long double sum_var = 0.0L;
  long double sum_sq  = 0.0L;

  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:sum_var,sum_sq) schedule(static) if (use_omp)
  #endif
  for (arma::sword j = 1; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    for (arma::uword i = 0; i < uj; ++i) {
      const double rij = XtX(i, uj);
      if (!std::isfinite(rij)) continue;
      const double rij2 = rij * rij;
      sum_sq  += rij2;
      const double one_minus = 1.0 - rij2;
      const double var_ij = (one_minus * one_minus) / static_cast<double>(n - 1);
      sum_var += var_ij;
    }
  }

#ifndef _OPENMP
  for (arma::uword j = 1; j < p; ++j) {
    for (arma::uword i = 0; i < j; ++i) {
      const double rij = XtX(i, j);
      if (!std::isfinite(rij)) continue;
      const double rij2 = rij * rij;
      sum_sq  += rij2;
      const double one_minus = 1.0 - rij2;
      const double var_ij = (one_minus * one_minus) / static_cast<double>(n - 1);
      sum_var += var_ij;
    }
  }
#endif

  // default to full shrinkage if denominator is zero
  double lambda = 1.0;
  if (sum_sq > 0.0L) {
    lambda = static_cast<double>(sum_var / sum_sq);
    lambda = std::max(0.0, std::min(1.0, lambda));
  }
  const double shrink = 1.0 - lambda;

  // ------------------------------------------------------------
  // R_shrunk = (1-lambda)R; diag=1
  // ------------------------------------------------------------
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if (use_omp)
    #endif
    for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
      const arma::uword uj = static_cast<arma::uword>(j);
      XtX(uj, uj) = 1.0;
      for (arma::uword i = 0; i < uj; ++i) {
        const double rij = XtX(i, uj);
        if (std::isfinite(rij)) {
          XtX(i, uj) = shrink * rij;
        }
      }
    }
  copy_upper_to_lower_inplace(XtX);

  // Propagate NA for zero-variance columns
  for (arma::uword k = 0; k < zero.n_elem; ++k) {
    const arma::uword j = zero[k];
    XtX.row(j).fill(arma::datum::nan);
    XtX.col(j).fill(arma::datum::nan);
    XtX(j, j) = arma::datum::nan;
  }

  // shrinkage correlation matrix
  return XtX;
}
