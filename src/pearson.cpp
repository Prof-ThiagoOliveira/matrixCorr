// Thiago de Paula Oliveira
// pearson.cpp -- upper-triangular Pearson kernel with a single mirror pass
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// X is n x p (double, no NA). Returns p x p Pearson correlation matrix.
// [[Rcpp::export]]
arma::mat pearson_matrix_cpp(SEXP X_) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");

  const arma::uword n = Rf_nrows(X_);
  const arma::uword p = Rf_ncols(X_);
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  // No-copy view over the input matrix.
  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  // Column means (1 x p).
  arma::rowvec mu = arma::sum(X, 0) / static_cast<double>(n);

  // XtX := X'X via SYRK (upper triangle); fallback to GEMM.
  arma::mat XtX(p, p, arma::fill::zeros);
#if defined(ARMA_USE_BLAS)
  {
    const arma::blas_int N = static_cast<arma::blas_int>(p);
    const arma::blas_int K = static_cast<arma::blas_int>(n);
    const double alpha = 1.0;
    const double beta = 0.0;
    const char uplo = 'U';
    const char trans = 'T';
    arma::blas::syrk<double>(&uplo, &trans, &N, &K,
                             &alpha, X.memptr(), &K,
                             &beta, XtX.memptr(), &N);
  }
#else
  XtX = X.t() * X;
#endif

  // XtX := XtX - n * mu * mu' on the stored upper triangle only.
  const double neg_n = -static_cast<double>(n);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double fj = neg_n * mu[uj];
    for (arma::uword i = 0; i <= uj; ++i) {
      XtX(i, uj) += fj * mu[i];
    }
  }

  // For correlation, the global covariance factor 1 / (n - 1) cancels.
  // Read the centred cross-product diagonal directly and normalise once.
  arma::vec centered_diag(p);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double d = XtX(uj, uj);
    centered_diag[uj] = d > 0.0 ? d : 0.0;
  }

  arma::vec inv_s = 1.0 / arma::sqrt(centered_diag);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double sj = inv_s[uj];
    if (!std::isfinite(sj) || sj == 0.0) continue;
    for (arma::uword i = 0; i < uj; ++i) {
      const double si = inv_s[i];
      if (!std::isfinite(si) || si == 0.0) continue;
      XtX(i, uj) *= (si * sj);
    }
    XtX(uj, uj) = 1.0;
  }

  XtX = arma::symmatu(XtX);

  // Zero-variance rows/cols are undefined.
  arma::uvec zero = arma::find(centered_diag == 0.0);
  for (arma::uword k = 0; k < zero.n_elem; ++k) {
    const arma::uword j = zero[k];
    XtX.row(j).fill(arma::datum::nan);
    XtX.col(j).fill(arma::datum::nan);
    XtX(j, j) = arma::datum::nan;
  }

  return XtX;
}
