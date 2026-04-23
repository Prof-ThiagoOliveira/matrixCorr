// Partial correlation with sample / ridge / OAS covariance estimators
// Thiago de Paula Oliveira
// partial_correlation.cpp
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include "matrixCorr_detail.h"
#include "matrixCorr_omp.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// functions we use here:
using namespace matrixCorr_detail;
using matrixCorr_detail::linalg::centered_crossprod_no_copy;
using matrixCorr_detail::linalg::invert_spd_inplace;
using matrixCorr_detail::linalg::precision_to_pcor_inplace;
using matrixCorr_detail::cov_shrinkage::oas_shrink_inplace;
using matrixCorr_detail::sparse_precision::graphical_lasso;

namespace {

arma::mat partial_correlation_p_values(const arma::mat& pcor,
                                       const arma::uword n,
                                       const arma::uword p) {
  if (n <= p) {
    Rcpp::stop("return_p_value requires n > p for the sample partial-correlation test.");
  }
  if (!pcor.is_finite()) {
    Rcpp::stop("Partial correlation matrix must contain only finite values when return_p_value = TRUE.");
  }

  const double df = static_cast<double>(n - p);
  arma::mat p_value(p, p, arma::fill::zeros);
  const bool use_omp = (p >= 64u && omp_get_max_threads() > 1);

 #ifdef _OPENMP
 #pragma omp parallel for schedule(static) if (use_omp)
 #endif
  for (arma::sword j = 1; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    for (arma::uword i = 0; i < uj; ++i) {
      double r = pcor(i, uj);
      r = std::max(-1.0, std::min(1.0, r));
      const double denom = 1.0 - (r * r);
      const double p_ij =
        (denom <= 0.0)
          ? 0.0
          : 2.0 * R::pt(-std::abs(r * std::sqrt(df / denom)), df, /*lower_tail*/ 1, /*log_p*/ 0);

      p_value(i, uj) = p_ij;
      p_value(uj, i) = p_ij;
    }
  }

  return p_value;
}

bool invert_spd_with_jitter_inplace(arma::mat& S,
                                    double& jitter,
                                    const double max_jitter = 1e-2) {
  if (jitter < 0.0) jitter = 0.0;
  for (;;) {
    arma::mat work = S;
    if (invert_spd_inplace(work)) {
      S = std::move(work);
      return true;
    }
    if (jitter == 0.0) jitter = 1e-8; else jitter *= 10.0;
    if (jitter > max_jitter) return false;
    S.diag() += jitter;
  }
}

bool invert_spd_with_jitter_keep_covariance(arma::mat& Sigma,
                                            arma::mat& Theta,
                                            double& jitter,
                                            const double max_jitter = 1e-2) {
  if (jitter < 0.0) jitter = 0.0;
  Theta = Sigma;
  for (;;) {
    arma::mat work = Theta;
    if (invert_spd_inplace(work)) {
      Theta = std::move(work);
      return true;
    }
    if (jitter == 0.0) jitter = 1e-8; else jitter *= 10.0;
    if (jitter > max_jitter) return false;
    Sigma.diag() += jitter;
    Theta.diag() += jitter;
  }
}

} // namespace

// Partial correlation matrix with sample / ridge / OAS / graphical-lasso estimators
// @param X_ Numeric double matrix (n x p). No NAs.
// @param method One of "sample", "ridge", "oas", "glasso". Default "sample".
// @param lambda Penalty parameter for "ridge" (covariance diagonal inflation) or
//        "glasso" (L1 penalty on off-diagonal precision entries). Ignored for "sample" and "oas".
// @param return_cov_precision If TRUE, return covariance and precision matrices.
// @param return_p_value If TRUE, also return two-sided p-values for the sample
//        partial correlations. Supported only for method = "sample" with n > p.
// @return A list with elements: \code{pcor}, and optionally \code{cov}, \code{precision},
//         \code{p_value}, \code{method}, \code{lambda}, \code{rho} (for OAS).
// @export
// [[Rcpp::export]]
 Rcpp::List partial_correlation_cpp(SEXP X_,
                                    const std::string method = "sample",
                                    const double lambda = 1e-3,
                                    const bool return_cov_precision = true,
                                    const bool return_p_value = false) {
   if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
     Rcpp::stop("Numeric double matrix required.");

   const arma::uword n = static_cast<arma::uword>(Rf_nrows(X_));
   const arma::uword p = static_cast<arma::uword>(Rf_ncols(X_));
   if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

   const double* x_ptr = REAL(X_);
   const arma::uword n_elem = n * p;
   for (arma::uword idx = 0; idx < n_elem; ++idx) {
     if (!std::isfinite(x_ptr[idx])) {
       Rcpp::stop("Input matrix must contain only finite values.");
     }
   }

   // No-copy
   arma::mat X(const_cast<double*>(x_ptr), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

   // (X - mu)'(X - mu), computed without forming a centered copy.
   arma::rowvec mu;
   arma::mat Sigma = centered_crossprod_no_copy(X, mu);

   double rho = NA_REAL;         // OAS shrinkage weight (if used)
   bool used_glasso = false;
   if (method == "sample" || method == "ridge") {
     const double inv_unbiased = 1.0 / static_cast<double>(n - 1);
     Sigma *= inv_unbiased;
     if (method == "ridge") {
       if (lambda < 0.0) Rcpp::stop("lambda must be non-negative.");
       if (lambda > 0.0) Sigma.diag() += lambda;
     }
   } else if (method == "oas") {
     Sigma *= 1.0 / static_cast<double>(n);
     oas_shrink_inplace(Sigma, static_cast<double>(n), rho);
   } else if (method == "glasso") {
     if (lambda < 0.0) Rcpp::stop("lambda must be non-negative.");
     Sigma *= 1.0 / static_cast<double>(n);
     used_glasso = true;
   } else {
     Rcpp::stop("Unknown method: '%s' (use 'sample', 'ridge', 'oas', or 'glasso').", method.c_str());
   }

   if (return_p_value && method != "sample") {
     Rcpp::stop("return_p_value is available only for method = 'sample'.");
   }

   double jitter = 0.0;
   arma::mat Theta;

   if (used_glasso) {
     arma::mat Sigma_hat;
     if (lambda == 0.0) {
       if (!invert_spd_with_jitter_keep_covariance(Sigma, Theta, jitter)) {
         Rcpp::stop("Covariance not positive definite; jitter exceeded limit.");
       }
     } else {
       if (!graphical_lasso(Sigma, lambda, Sigma_hat, Theta)) {
         Rcpp::stop("Graphical lasso failed to converge.");
       }
       Sigma = std::move(Sigma_hat);
     }
   }

   if (return_cov_precision) {
     if (!used_glasso) {
       if (!invert_spd_with_jitter_keep_covariance(Sigma, Theta, jitter)) {
         Rcpp::stop("Covariance not positive definite; jitter exceeded limit.");
       }
     }

     arma::mat pcor = Theta;
     if (!precision_to_pcor_inplace(pcor)) {
       Rcpp::stop("Precision diagonal must be positive and finite.");
     }

     arma::mat p_value;
     if (return_p_value) {
       p_value = partial_correlation_p_values(pcor, n, p);
     }

     Rcpp::List out = Rcpp::List::create(
       Rcpp::Named("pcor")      = pcor,
       Rcpp::Named("cov")       = Sigma,
       Rcpp::Named("precision") = Theta,
       Rcpp::Named("p_value")   = return_p_value ? Rcpp::wrap(p_value) : R_NilValue,
       Rcpp::Named("method")    = method,
       Rcpp::Named("lambda")    = (method == "ridge" ? lambda : NA_REAL),
       Rcpp::Named("rho")       = (method == "oas"   ? rho    : NA_REAL),
       Rcpp::Named("jitter")    = jitter
     );
     return out;
   } else {
     if (used_glasso) {
       Sigma = std::move(Theta);
     } else if (!invert_spd_with_jitter_inplace(Sigma, jitter)) {
       Rcpp::stop("Covariance not positive definite; jitter exceeded limit.");
     }
     if (!precision_to_pcor_inplace(Sigma)) {
       Rcpp::stop("Precision diagonal must be positive and finite.");
     }

     arma::mat p_value;
     if (return_p_value) {
       p_value = partial_correlation_p_values(Sigma, n, p);
     }

     return Rcpp::List::create(
       Rcpp::Named("pcor")    = Sigma,
       Rcpp::Named("p_value") = return_p_value ? Rcpp::wrap(p_value) : R_NilValue,
       Rcpp::Named("lambda")  = (method == "ridge" ? lambda : NA_REAL),
       Rcpp::Named("rho")     = (method == "oas"   ? rho    : NA_REAL),
       Rcpp::Named("jitter")  = jitter
     );
   }
 }
