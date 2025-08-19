// Partial correlation with sample / ridge / OAS covariance estimators
// Thiago de Paula Oliveira
// partial_correlation.cpp
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

namespace detail {

  // Compute X'X (upper triangle) then mirror
  // Returns p x p matrix XtX = X'X.
  // Uses BLAS SYRK if available.
  inline arma::mat crossprod_no_copy(const arma::mat& X) {
    const arma::uword n = X.n_rows, p = X.n_cols;
    arma::mat XtX(p, p);
    #if defined(ARMA_USE_BLAS)
    {
      XtX.zeros();
      const arma::blas_int N = static_cast<arma::blas_int>(p);
      const arma::blas_int K = static_cast<arma::blas_int>(n);
      const double alpha = 1.0, beta = 0.0;
      const char uplo = 'U', trans = 'T';
      arma::blas::syrk<double>(&uplo, &trans, &N, &K,
                               &alpha, X.memptr(), &K,
                               &beta,  XtX.memptr(), &N);
      XtX = arma::symmatu(XtX);
    }
    #else
    XtX = X.t() * X;
    #endif
    return XtX;
  }

// M := M - n * mu * mu'  (upper) then mirror.
// mu is 1 x p.
inline void subtract_n_outer_mu(arma::mat& M, const arma::rowvec& mu, double n) {
  const arma::uword p = M.n_cols;
  const double scale = -n;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double muj = mu[uj];
    const double fj  = scale * muj;           // = -n * mu_j
    for (arma::uword i = 0; i <= uj; ++i) {
      M(i, uj) += fj * mu[i];                 // -= n * mu_i * mu_j
    }
  }
  M = arma::symmatu(M);
}

// Ensure positive definiteness
inline void make_pd_inplace(arma::mat& S, double& jitter, const double max_jitter = 1e-2) {
  // escalate jitter geometrically until chol succeeds or limit reached
  if (jitter < 0) jitter = 0.0;
  for (;;) {
    bool ok = false;
    arma::mat C;
    // try Cholesky
    ok = arma::chol(C, S, "upper");
    if (ok) { return; }
    if (jitter == 0.0) jitter = 1e-8;
    else jitter *= 10.0;
    if (jitter > max_jitter) {
      Rcpp::stop("Covariance not positive definite; jitter exceeded limit.");
    }
    S.diag() += jitter;
  }
}

// Compute OAS shrinkage of covariance to a scaled identity (Chen–Wiesel–Hero, 2010).
// Input is cov_mle = (1/n) * (X - mu)'(X - mu); n samples, p variables.
// Returns shrunken covariance Sigma = (1 - rho)*S + rho*mu*I, with rho in [0,1].
inline arma::mat oas_shrink(const arma::mat& cov_mle, double n, double& rho_out) {
  const arma::uword p = cov_mle.n_cols;
  const double trS   = arma::trace(cov_mle);
  const double trS2  = arma::accu(cov_mle % cov_mle); // Frobenius^2 = tr(S^2)
  const double mu    = trS / static_cast<double>(p);

  // OAS shrinkage parameter (to identity target mu*I)
  // rho = min(1, max(0, ((1 - 2/p)*tr(S^2) + tr(S)^2) /
  //  ((n + 1 - 2/p)*(tr(S^2) - tr(S)^2/p)) ))
  const double p_d   = static_cast<double>(p);
  const double num   = (1.0 - 2.0 / p_d) * trS2 + trS * trS;
  const double den   = (n + 1.0 - 2.0 / p_d) * (trS2 - (trS * trS) / p_d);

  double rho = (den > 0.0) ? (num / den) : 1.0;
  rho = std::max(0.0, std::min(1.0, rho));
  rho_out = rho;

  arma::mat Sigma = (1.0 - rho) * cov_mle;
  // add rho * mu * I
  Sigma.diag() += rho * mu;
  return Sigma;
}

}


//' Partial correlation matrix with sample / ridge / OAS covariance
//'
//' @param X_ Numeric double matrix (n x p). No NAs.
//' @param method One of "sample", "ridge", "oas". Default "oas" (recommended for p >> n).
//' @param lambda Ridge penalty for "ridge" method (added to diagonal). Ignored otherwise.
//' @param return_cov_precision If TRUE, return covariance and precision matrices.
//' @return A list with elements: \code{pcor}, and optionally \code{cov}, \code{precision},
//'         \code{method}, \code{lambda}, \code{rho} (for OAS).
//' @export
// [[Rcpp::export]]
Rcpp::List partial_correlation_cpp(SEXP X_,
                                   const std::string method = "oas",
                                   const double lambda = 1e-3,
                                   const bool return_cov_precision = true) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");

  const arma::uword n = static_cast<arma::uword>(Rf_nrows(X_));
  const arma::uword p = static_cast<arma::uword>(Rf_ncols(X_));
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  // No-copy
  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  // Column means (1 x p)
  arma::rowvec mu = arma::sum(X, 0) / static_cast<double>(n);

  // X'X, then subtract n * mu * mu' (no centred copy), giving centred
  // cross-product.
  arma::mat XtX = detail::crossprod_no_copy(X);
  detail::subtract_n_outer_mu(XtX, mu, static_cast<double>(n));

  // Two scalings of covariance (MLE and unbiased)
  arma::mat cov_mle      = XtX / static_cast<double>(n);
  arma::mat cov_unbiased = XtX / static_cast<double>(n - 1);

  // Choose estimator for Sigma
  arma::mat Sigma;
  double rho = NA_REAL;         // OAS shrinkage weight (if used)
  if (method == "sample") {
    Sigma = std::move(cov_unbiased);
  } else if (method == "ridge") {
    Sigma = std::move(cov_unbiased);
    if (lambda < 0.0) Rcpp::stop("lambda must be non-negative.");
    if (lambda > 0.0) Sigma.diag() += lambda;
  } else if (method == "oas") {
    Sigma = detail::oas_shrink(cov_mle, static_cast<double>(n), rho);
  } else {
    Rcpp::stop("Unknown method: '%s' (use 'sample', 'ridge', or 'oas').", method.c_str());
  }

  // Ensure positive definiteness (adds minimal jitter if needed)
  double jitter = 0.0;
  detail::make_pd_inplace(Sigma, jitter);

  // Precision matrix: Theta = inv(Sigma)
  arma::mat Theta;
  bool inv_ok = arma::inv_sympd(Theta, Sigma); // SPD path
  if (!inv_ok) {
    // Fallback via solve; add tiny jitter if still unstable
    Sigma.diag() += 1e-12;
    inv_ok = arma::solve(Theta, Sigma, arma::eye<arma::mat>(p, p));
    if (!inv_ok) Rcpp::stop("Failed to invert covariance (even after jitter).");
  }

  // pcor_ij = -Theta_ij / sqrt(Theta_ii * Theta_jj)
  arma::vec d = Theta.diag();
  if (d.min() <= 0.0 || !d.is_finite()) {
    Rcpp::stop("Precision diagonal must be positive and finite.");
  }
  arma::vec inv_sqrt = 1.0 / arma::sqrt(d);

  arma::mat pcor(p, p, arma::fill::ones);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double sj = inv_sqrt[uj];
    for (arma::uword i = 0; i < uj; ++i) {
      const double si = inv_sqrt[i];
      double val = -Theta(i, uj) * (si * sj);
      pcor(i, uj) = val;
      pcor(uj, i) = val;
    }
    pcor(uj, uj) = 1.0;
  }

  if (return_cov_precision) {
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("pcor")      = pcor,
      Rcpp::Named("cov")       = Sigma,
      Rcpp::Named("precision") = Theta,
      Rcpp::Named("method")    = method,
      Rcpp::Named("lambda")    = (method == "ridge" ? lambda : NA_REAL),
      Rcpp::Named("rho")       = (method == "oas"   ? rho    : NA_REAL),
      Rcpp::Named("jitter")    = jitter
    );
    return out;
  } else {
    return Rcpp::List::create(Rcpp::Named("pcor") = pcor);
  }
}
