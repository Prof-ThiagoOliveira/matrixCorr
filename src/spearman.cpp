// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include "matrixCorr_detail.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using matrixCorr_detail::ranking::rank_vector;
using matrixCorr_detail::ranking::safe_inv_stddev;

// [[Rcpp::export]]
arma::mat spearman_matrix_cpp(SEXP X_) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");

  const arma::uword n = static_cast<arma::uword>(Rf_nrows(X_));
  const arma::uword p = static_cast<arma::uword>(Rf_ncols(X_));
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  // No-copy view into R memory
  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  // 1) Rank each column (OpenMP-parallel)
  arma::mat R(n, p, arma::fill::none);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    // Create a view onto the j-th column memory to avoid an intermediate vector
    arma::vec out(R.colptr(static_cast<arma::uword>(j)), n, /*copy_aux_mem*/ false, /*strict*/ true);
    rank_vector(X.col(static_cast<arma::uword>(j)), out);  // in-place, fast path when no ties
  }

  // 2) Compute RtR via BLAS SYRK (upper triangle), fallback to GEMM
  arma::mat RtR(p, p);
#if defined(ARMA_USE_BLAS)
{
  RtR.zeros();
  const arma::blas_int N = static_cast<arma::blas_int>(p);
  const arma::blas_int K = static_cast<arma::blas_int>(n);
  const double alpha = 1.0, beta = 0.0;
  const char uplo  = 'U';
  const char trans = 'T'; // C := A' * A
  arma::blas::syrk<double>(&uplo, &trans, &N, &K,
                           &alpha, R.memptr(), &K,
                           &beta,  RtR.memptr(), &N);
}
#else
RtR = R.t() * R;
#endif
RtR = arma::symmatu(RtR);  // Fill lower triangle

// 3) Subtract constant n * mu^2 from all entries; mu = (n+1)/2 for ranks
const double mu = 0.5 * (static_cast<double>(n) + 1.0);
RtR -= static_cast<double>(n) * mu * mu;

// 4) Standard deviations from diagonal, then divide by (n - 1)
arma::vec s2 = RtR.diag() / static_cast<double>(n - 1);
s2 = arma::clamp(s2, 0.0, std::numeric_limits<double>::infinity());
arma::vec s  = arma::sqrt(s2);

// 5) corr = RtR / (n - 1), then D^{-1} * corr * D^{-1}
arma::mat corr = RtR / static_cast<double>(n - 1);
arma::vec inv_s = safe_inv_stddev(s);

// scale columns then rows (equivalent to D^{-1} * corr * D^{-1})
corr.each_row() %= inv_s.t();  // scale columns j by inv_s[j]
corr.each_col() %= inv_s;      // scale rows    i by inv_s[i]

// 6) Set diagonals and mark zero-variance rows/cols as NA
corr.diag().ones();
arma::uvec zero = arma::find(s == 0.0);
for (arma::uword k = 0; k < zero.n_elem; ++k) {
  const arma::uword j = zero[k];
  corr.row(j).fill(arma::datum::nan);
  corr.col(j).fill(arma::datum::nan);
  corr(j, j) = arma::datum::nan;
}

return corr;  // correlation matrix in [-1, 1]
}
