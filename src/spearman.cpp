// Thiago de Paula Oliveira
// spearman.cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

inline arma::vec rank_vector(const arma::vec& x) {
  const int n = x.n_elem;
  arma::uvec idx = arma::sort_index(x);     // sorted positions
  arma::vec ranks(n);

  int i = 0;
  while (i < n) {
    int j = i + 1;
    const double xi = x(idx[i]);
    while (j < n && x(idx[j]) == xi) ++j;   // ties
    const double avg_rank = (i + j - 1) * 0.5 + 1.0; // 1-based
    for (int k = i; k < j; ++k) ranks(idx[k]) = avg_rank;
    i = j;
  }
  return ranks;
}

// [[Rcpp::export]]
arma::mat spearman_matrix_cpp(SEXP X_) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");
  const int n  = Rf_nrows(X_);
  const int p  = Rf_ncols(X_);
  arma::mat X(REAL(X_), n, p, false, true);      // no-copy view

  arma::mat ranked_X(n, p);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    ranked_X.col(j) = rank_vector(X.col(j));
  }

  return arma::cor(ranked_X);  // BLAS
}
