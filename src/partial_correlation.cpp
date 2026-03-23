// Partial correlation with sample / ridge / OAS covariance estimators
// Thiago de Paula Oliveira
// partial_correlation.cpp
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include "matrixCorr_detail.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// functions we use here:
using namespace matrixCorr_detail;
using matrixCorr_detail::linalg::crossprod_no_copy;
using matrixCorr_detail::linalg::subtract_n_outer_mu;
using matrixCorr_detail::linalg::make_pd_inplace;
using matrixCorr_detail::linalg::invert_spd_inplace;
using matrixCorr_detail::linalg::precision_to_pcor_inplace;
using matrixCorr_detail::cov_shrinkage::oas_shrink_inplace;
using matrixCorr_detail::sparse_precision::graphical_lasso;

// Partial correlation matrix with sample / ridge / OAS / graphical-lasso estimators
// @param X_ Numeric double matrix (n x p). No NAs.
// @param method One of "sample", "ridge", "oas", "glasso". Default "sample".
// @param lambda Penalty parameter for "ridge" (covariance diagonal inflation) or
//        "glasso" (L1 penalty on off-diagonal precision entries). Ignored for "sample" and "oas".
// @param return_cov_precision If TRUE, return covariance and precision matrices.
// @return A list with elements: \code{pcor}, and optionally \code{cov}, \code{precision},
//         \code{method}, \code{lambda}, \code{rho} (for OAS).
// @export
// [[Rcpp::export]]
 Rcpp::List partial_correlation_cpp(SEXP X_,
                                    const std::string method = "sample",
                                    const double lambda = 1e-3,
                                    const bool return_cov_precision = true) {
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

   // Column means (1 x p)
   arma::rowvec mu = arma::sum(X, 0) / static_cast<double>(n);

   arma::mat Sigma;
   // X'X, then subtract n * mu * mu' (no centred copy), giving centred cross-product.
   Sigma = crossprod_no_copy(X);
   subtract_n_outer_mu(Sigma, mu, static_cast<double>(n));

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

   double jitter = 0.0;
   arma::mat Theta;

   if (used_glasso) {
     arma::mat Sigma_hat;
     if (lambda == 0.0) {
       make_pd_inplace(Sigma, jitter);
       Theta = Sigma;
       if (!invert_spd_inplace(Theta)) {
         Rcpp::stop("Failed to invert covariance after positive-definite adjustment.");
       }
     } else {
       if (!graphical_lasso(Sigma, lambda, Sigma_hat, Theta)) {
         Rcpp::stop("Graphical lasso failed to converge.");
       }
       Sigma = std::move(Sigma_hat);
     }
   } else {
     // Ensure positive definiteness (adds minimal jitter if needed)
     make_pd_inplace(Sigma, jitter);
   }

   if (return_cov_precision) {
     if (!used_glasso) {
       Theta = Sigma;
       if (!invert_spd_inplace(Theta)) {
         Rcpp::stop("Failed to invert covariance after positive-definite adjustment.");
       }
     }

     arma::mat pcor = Theta;
     if (!precision_to_pcor_inplace(pcor)) {
       Rcpp::stop("Precision diagonal must be positive and finite.");
     }

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
     if (used_glasso) {
       Sigma = std::move(Theta);
     } else if (!invert_spd_inplace(Sigma)) {
       Rcpp::stop("Failed to invert covariance after positive-definite adjustment.");
     }
     if (!precision_to_pcor_inplace(Sigma)) {
       Rcpp::stop("Precision diagonal must be positive and finite.");
     }

     return Rcpp::List::create(
       Rcpp::Named("pcor")   = Sigma,
       Rcpp::Named("lambda") = (method == "ridge" ? lambda : NA_REAL),
       Rcpp::Named("rho")    = (method == "oas"   ? rho    : NA_REAL),
       Rcpp::Named("jitter") = jitter
     );
   }
 }
