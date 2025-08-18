// Thiago de Paula Oliveira
// pearson.cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Pearson via cross-product
// X is n x p (double, no NA). Returns p x p correlation matrix.
// [[Rcpp::export]]
arma::mat pearson_matrix_cpp(SEXP X_) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");

  const arma::uword n = Rf_nrows(X_);
  const arma::uword p = Rf_ncols(X_);
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  // No-copy view over R memory
  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  // column means: one pass, no big temp
  arma::rowvec mu = arma::sum(X, 0) / static_cast<double>(n);   // 1 x p

  // X'X using SYRK (upper triangle) when BLAS is available; fallback to GEMM
  arma::mat XtX(p, p);
#if defined(ARMA_USE_BLAS)
{
  XtX.zeros();
  const arma::blas_int N = static_cast<arma::blas_int>(p);
  const arma::blas_int K = static_cast<arma::blas_int>(n);
  const double alpha = 1.0, beta = 0.0;
  const char uplo = 'U', trans = 'T';
  // C := alpha*A'*A + beta*C   (upper triangle filled)
  arma::blas::syrk<double>(&uplo, &trans, &N, &K,
                           &alpha, X.memptr(), &K,
                           &beta,  XtX.memptr(), &N);
                           XtX = arma::symmatu(XtX);
}
#else
XtX = X.t() * X;              // GEMM fallback
#endif

// Cov = (X'X - n*mu*mu^T) / (n-1)
arma::mat cov = XtX - static_cast<double>(n) * (mu.t() * mu);
cov /= static_cast<double>(n - 1);

// std devs; clamp to avoid tiny negatives from FP noise
arma::vec s = arma::sqrt( arma::clamp(cov.diag(), 0.0,
                                      std::numeric_limits<double>::infinity()) );

// Normalize to correlation
arma::mat corr = cov;
corr.each_col() /= s;
corr.each_row() /= s.t();

// Handle zero-variance columns like base::cor: all-NA in that row/col
arma::uvec zero = arma::find(s == 0.0);
corr.diag().ones();
for (arma::uword k = 0; k < zero.n_elem; ++k) {
  const arma::uword j = zero[k];
  corr.row(j).fill(arma::datum::nan);
  corr.col(j).fill(arma::datum::nan);
  corr(j, j) = arma::datum::nan;
}
return corr;
}
