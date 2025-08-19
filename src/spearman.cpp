// Thiago de Paula Oliveira
// spearman.cpp — optimised: no column-sum pass; in-place normalisation
#include <RcppArmadillo.h>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

inline arma::vec rank_vector(const arma::vec& x) {
  const arma::uword n = x.n_elem;
  arma::uvec idx = arma::sort_index(x);
  arma::vec ranks(n);

  arma::uword i = 0;
  while (i < n) {
    arma::uword j = i + 1;
    const double xi = x(idx[i]);
    while (j < n && x(idx[j]) == xi) ++j;                 // tie block [i, j)
    const double avg_rank = (static_cast<double>(i + j - 1) * 0.5) + 1.0; // 1-based
    for (arma::uword k = i; k < j; ++k) ranks(idx[k]) = avg_rank;
    i = j;
  }
  return ranks;
}

// [[Rcpp::export]]
arma::mat spearman_matrix_cpp(SEXP X_) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");

  const arma::uword n = static_cast<arma::uword>(Rf_nrows(X_));
  const arma::uword p = static_cast<arma::uword>(Rf_ncols(X_));
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  // No-copy view
  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  // 1) Rank each column (OpenMP-parallel)
  arma::mat R(n, p);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    R.col(static_cast<arma::uword>(j)) = rank_vector(X.col(static_cast<arma::uword>(j)));
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
                           // mirror
                           RtR = arma::symmatu(RtR);
}
#else
RtR = R.t() * R;
#endif

// 3) Subtract constant n * mu^2 from all entries; mu = (n+1)/2 for ranks
const double mu     = 0.5 * (static_cast<double>(n) + 1.0);
const double n_mu2  = static_cast<double>(n) * mu * mu;
RtR -= n_mu2;                         // element-wise: fast constant shift

// 4) Standard deviations from diagonal, then divide by (n - 1)
//    s_j^2 = (RtR_jj) / (n - 1)  after the constant shift above
arma::vec s2 = RtR.diag() / static_cast<double>(n - 1);
s2 = arma::clamp(s2, 0.0, std::numeric_limits<double>::infinity());
arma::vec s  = arma::sqrt(s2);

// 5) corr = RtR / (n - 1), then D^{-1} * corr * D^{-1}
arma::mat corr = RtR;
corr /= static_cast<double>(n - 1);

// Indices with non-zero std dev
arma::uvec nz = arma::find(s > 0.0);
if (!nz.is_empty()) {
  arma::vec s_nz = s.elem(nz);
  // normalise only the nz x nz block to avoid division by zero
  arma::mat sub = corr.submat(nz, nz);
  sub.each_col() /= s_nz;
  sub.each_row() /= s_nz.t();
}

// 6) Set diagonals and mark zero-variance rows/cols as NA
corr.diag().ones();
arma::uvec zero = arma::find(s == 0.0);
for (arma::uword k = 0; k < zero.n_elem; ++k) {
  const arma::uword j = zero[k];
  corr.row(j).fill(arma::datum::nan);
  corr.col(j).fill(arma::datum::nan);
  corr(j, j) = arma::datum::nan;
}

return corr;
}
