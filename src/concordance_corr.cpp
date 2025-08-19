// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat ccc_cpp(const arma::mat& X) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  arma::mat result(p, p, arma::fill::eye);

  arma::vec means(p);
  arma::vec vars(p);

  // Precompute means and variances
  for (int j = 0; j < p; ++j) {
    arma::vec col = X.col(j);
    double mean_j = arma::mean(col);
    double var_j = arma::var(col) * (n - 1) / n;
    means[j] = mean_j;
    vars[j] = var_j;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      double mean_x = means[i];
      double mean_y = means[j];
      double var_x = vars[i];
      double var_y = vars[j];

      const arma::vec x = X.col(i);
      const arma::vec y = X.col(j);

      double cov_xy = 0.0;
      for (int k = 0; k < n; ++k) {
        cov_xy += (x[k] - mean_x) * (y[k] - mean_y);
      }
      cov_xy /= n;

      double r = cov_xy / std::sqrt(var_x * var_y);
      double sxy = r * std::sqrt(var_x * var_y);
      double p = 2.0 * sxy / (var_x + var_y + std::pow(mean_x - mean_y, 2.0));
      result(i, j) = p;
      result(j, i) = p;
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

  arma::vec means(p);
  arma::vec vars(p);

  // Precompute means and variances
  for (int j = 0; j < p; ++j) {
    arma::vec col = X.col(j);
    means[j] = arma::mean(col);
    vars[j] = arma::var(col) * (n - 1) / n;
  }

  double alpha = 1.0 - conf_level;
  double z = R::qnorm(1.0 - alpha / 2.0, 0.0, 1.0, 1, 0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      const arma::vec x = X.col(i);
      const arma::vec y = X.col(j);
      double mean_x = means[i], mean_y = means[j];
      double var_x = vars[i], var_y = vars[j];

      double cov_xy = arma::as_scalar(arma::cov(x, y)) * (n - 1) / n;
      double r = cov_xy / std::sqrt(var_x * var_y);
      double sxy = r * std::sqrt(var_x * var_y);
      double p_val = 2.0 * sxy / (var_x + var_y + std::pow(mean_x - mean_y, 2.0));

      // Bias correction
      double u = (mean_y - mean_x) / std::pow(var_x * var_y, 0.25);
      double sep = std::sqrt(((1 - r * r) * p_val * p_val * (1 - p_val * p_val) / (r * r) +
                             (2 * std::pow(p_val, 3) * (1 - p_val) * u * u / r) -
                             (0.5 * std::pow(p_val, 4) * std::pow(u, 4) / (r * r))) / (n - 2));

      // z-transform CI
      double t = 0.5 * std::log((1 + p_val) / (1 - p_val));
      double set = sep / (1 - p_val * p_val);
      double llt = t - z * set;
      double ult = t + z * set;
      double lci = (std::exp(2 * llt) - 1) / (std::exp(2 * llt) + 1);
      double uci = (std::exp(2 * ult) - 1) / (std::exp(2 * ult) + 1);

      est(i, j) = est(j, i) = p_val;
      lwr(i, j) = lwr(j, i) = lci;
      upr(i, j) = upr(j, i) = uci;
    }
  }

  return List::create(
    Named("est") = est,
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
