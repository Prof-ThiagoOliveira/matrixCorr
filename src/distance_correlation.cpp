// Thiago de Paula Oliveira

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// Pairwise unbiased distance correlation (U-statistic)
// Székely, Rizzo & Bakirov, 2007
// [[Rcpp::export]]
double ustat_dcor(const arma::vec& x, const arma::vec& y) {
  const int n = x.n_elem;
  if (n < 4) Rcpp::stop("Sample size must be at least 4 for unbiased dCor");

  static thread_local std::vector<double> Rx_buf;
  static thread_local std::vector<double> Ry_buf;
  Rx_buf.assign(static_cast<std::size_t>(n), 0.0);
  Ry_buf.assign(static_cast<std::size_t>(n), 0.0);

  double Sx = 0.0, Sy = 0.0;
  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    const double yi = y[i];
    for (int j = i + 1; j < n; ++j) {
      const double dx = std::abs(xi - x[j]);
      const double dy = std::abs(yi - y[j]);
      Rx_buf[static_cast<std::size_t>(i)] += dx;
      Rx_buf[static_cast<std::size_t>(j)] += dx;
      Ry_buf[static_cast<std::size_t>(i)] += dy;
      Ry_buf[static_cast<std::size_t>(j)] += dy;
      Sx += 2.0 * dx;
      Sy += 2.0 * dy;
    }
  }

  const double inv_nm2 = 1.0 / static_cast<double>(n - 2);
  const double add_x = Sx / static_cast<double>((n - 1) * (n - 2));
  const double add_y = Sy / static_cast<double>((n - 1) * (n - 2));

  double XY = 0.0, X2 = 0.0, Y2 = 0.0;
  for (int i = 0; i < n; ++i) {
    const double rxi = Rx_buf[static_cast<std::size_t>(i)];
    const double ryi = Ry_buf[static_cast<std::size_t>(i)];
    const double xi = x[i];
    const double yi = y[i];
    for (int j = i + 1; j < n; ++j) {
      const double dx = std::abs(xi - x[j]);
      const double dy = std::abs(yi - y[j]);
      const double ax = dx - (rxi + Rx_buf[static_cast<std::size_t>(j)]) * inv_nm2 + add_x;
      const double ay = dy - (ryi + Ry_buf[static_cast<std::size_t>(j)]) * inv_nm2 + add_y;
      XY += ax * ay;
      X2 += ax * ax;
      Y2 += ay * ay;
    }
  }

  const double scale = 2.0 / (static_cast<double>(n) * static_cast<double>(n - 3));
  XY *= scale;
  X2 *= scale;
  Y2 *= scale;

  if (X2 <= 0.0 || Y2 <= 0.0) return NA_REAL;

  const double denom = std::sqrt(X2 * Y2);
  double dcor = XY / denom;
  if (!std::isfinite(dcor)) return NA_REAL;
  if (dcor < 0.0) dcor = 0.0;
  if (dcor > 1.0) dcor = 1.0;
  return dcor;
}

// Full matrix of unbiased distance correlations
// [[Rcpp::export]]
arma::mat ustat_dcor_matrix_cpp(const arma::mat& X) {
  const int p = X.n_cols;
  arma::mat R(p, p, arma::fill::eye);

  // Parallelize over upper-triangular column pairs
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int j = 1; j < p; ++j) {
    for (int i = 0; i < j; ++i) {
      const arma::vec xi = X.col(i);
      const arma::vec xj = X.col(j);
      const double d = ustat_dcor(xi, xj);
      R(i, j) = d;
      R(j, i) = d;
    }
  }
  return R;
}
