// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include <cmath>
#include "matrixCorr_detail.h"
#include "threshold_triplets.h"
using namespace Rcpp;

// only what we use from matrixCorr_detail
using namespace matrixCorr_detail;
using matrixCorr_detail::moments::col_means_vars_pop;
using matrixCorr_detail::moments::cov_xy_pop_manual;
using matrixCorr_detail::ccc_bits::ccc_from_stats_via_r;
using matrixCorr_detail::norm1::qnorm01;
using matrixCorr_detail::ccc_se::se_delta;
using matrixCorr_detail::fisherz::ci_from_z;
using matrixCorr_detail::symm::put;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat ccc_cpp(const arma::mat& X) {
  const int p = X.n_cols;
  arma::mat result(p, p, arma::fill::eye);

  arma::vec means(p), vars(p);
  // Precompute means and population variances (identical scaling)
  col_means_vars_pop(X, means, vars);
  const arma::uword n = static_cast<arma::uword>(X.n_rows);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      const double mean_x = means[i];
      const double mean_y = means[j];
      const double var_x  = vars[i];
      const double var_y  = vars[j];

      // population covariance via manual loop
      const double cov_xy = cov_xy_pop_manual(
        X.colptr(static_cast<arma::uword>(i)),
        X.colptr(static_cast<arma::uword>(j)),
        n,
        mean_x,
        mean_y
      );

      // CCC with identical algebraic order (r -> sxy -> p)
      const double pij = ccc_from_stats_via_r(mean_x, mean_y, var_x, var_y, cov_xy);

      put(result, i, j, pij);
    }
  }

  return result;
}

// [[Rcpp::export]]
List ccc_with_ci_cpp(const arma::mat& X, double conf_level = 0.95) {
  const int n = X.n_rows;
  const int p = X.n_cols;

  arma::mat est(p, p, arma::fill::eye);
  arma::mat lwr(p, p, arma::fill::zeros);
  arma::mat upr(p, p, arma::fill::zeros);

  arma::vec means(p), vars(p);
  // Precompute means and population variances
  col_means_vars_pop(X, means, vars);

  const double alpha  = 1.0 - conf_level;
  const double zcrit  = qnorm01(1.0 - alpha / 2.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      const double mean_x = means[i], mean_y = means[j];
      const double var_x  = vars[i],  var_y  = vars[j];

      const double cov_xy = cov_xy_pop_manual(
        X.colptr(static_cast<arma::uword>(i)),
        X.colptr(static_cast<arma::uword>(j)),
        static_cast<arma::uword>(n),
        mean_x,
        mean_y
      );

      // CCC estimate
      const double p_val = ccc_from_stats_via_r(mean_x, mean_y, var_x, var_y, cov_xy);

      // components for delta-method SE
      const double r   = cov_xy / std::sqrt(var_x * var_y);
      const double u   = (mean_y - mean_x) / std::pow(var_x * var_y, 0.25);
      const double sep = se_delta(r, p_val, u, n);                 // SE on p
      const double set = sep / (1.0 - p_val * p_val);              // SE on Fisher z

      double lci = 0.0, uci = 0.0;
      ci_from_z(p_val, set, zcrit, lci, uci);                      // Fisher-z CI

      put(est, i, j, p_val);
      put(lwr, i, j, lci);
      put(upr, i, j, uci);
    }
  }

  return List::create(
    Named("est")    = est,
    Named("lwr.ci") = lwr,
    Named("upr.ci") = upr
  );
}

// [[Rcpp::export]]
int openmp_threads() {
  int n = 1;
#ifdef _OPENMP
  n = omp_get_max_threads();
#endif
  return n;
}

// Complete-data CCC upper-triplet kernel for thresholded outputs.
// [[Rcpp::export]]
Rcpp::List ccc_threshold_triplets_cpp(const arma::mat& X,
                                      const double threshold = 0.0,
                                      const bool diag = true,
                                      const int block_size = 256) {
  if (!(threshold >= 0.0) || !std::isfinite(threshold)) {
    Rcpp::stop("threshold must be finite and >= 0.");
  }
  if (block_size < 1) {
    Rcpp::stop("block_size must be >= 1.");
  }
  if (X.n_rows < 2 || X.n_cols < 2) {
    Rcpp::stop("Need >= 2 rows and >= 2 columns.");
  }
  if (!X.is_finite()) {
    Rcpp::stop("X contains NA/NaN/Inf; please handle missingness upstream.");
  }

  const std::size_t p = static_cast<std::size_t>(X.n_cols);
  const arma::uword n = static_cast<arma::uword>(X.n_rows);
  const double n_d = static_cast<double>(n);

  arma::vec means(X.n_cols, arma::fill::zeros);
  arma::vec vars(X.n_cols, arma::fill::zeros);
  col_means_vars_pop(X, means, vars);

  const auto trip = matrixCorr_detail::threshold_triplets::collect_upper_triplets(
    p,
    static_cast<std::size_t>(block_size),
    diag,
    threshold,
    [&](std::size_t j0, std::size_t j1, std::size_t k0, std::size_t k1) -> arma::mat {
      arma::mat blk = X.cols(static_cast<arma::uword>(j0), static_cast<arma::uword>(j1 - 1u)).t() *
        X.cols(static_cast<arma::uword>(k0), static_cast<arma::uword>(k1 - 1u));

      const arma::rowvec mu_j = means.subvec(static_cast<arma::uword>(j0), static_cast<arma::uword>(j1 - 1u)).t();
      const arma::rowvec mu_k = means.subvec(static_cast<arma::uword>(k0), static_cast<arma::uword>(k1 - 1u)).t();
      blk -= n_d * (mu_j.t() * mu_k);

      for (std::size_t r = 0u; r < (j1 - j0); ++r) {
        const arma::uword gj = static_cast<arma::uword>(j0 + r);
        for (std::size_t c = 0u; c < (k1 - k0); ++c) {
          const arma::uword gk = static_cast<arma::uword>(k0 + c);
          double val = NA_REAL;
          if (gj == gk) {
            val = 1.0;
          } else {
            const double cov_xy = blk(
              static_cast<arma::uword>(r),
              static_cast<arma::uword>(c)
            ) / n_d;
            val = ccc_from_stats_via_r(
              means[gj],
              means[gk],
              vars[gj],
              vars[gk],
              cov_xy
            );
          }
          blk(static_cast<arma::uword>(r), static_cast<arma::uword>(c)) = val;
        }
      }
      return blk;
    }
  );

  return matrixCorr_detail::threshold_triplets::as_list(trip);
}
