// Thiago de Paula Oliveira

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>
#include "matrixCorr_omp.h"

namespace {

constexpr double HSIC_EPS = 1e-12;

inline double clamp01_hsic(const double x) {
  if (!std::isfinite(x)) return NA_REAL;
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

inline double vec_mean(const std::vector<double>& x) {
  double s = 0.0;
  for (double v : x) s += v;
  return s / static_cast<double>(x.size());
}

inline double vec_sd(const std::vector<double>& x) {
  const int n = static_cast<int>(x.size());
  if (n < 2) return 0.0;
  const double m = vec_mean(x);
  double ss = 0.0;
  for (double v : x) {
    const double d = v - m;
    ss += d * d;
  }
  return std::sqrt(ss / static_cast<double>(n - 1));
}

inline double upper_median_pairwise_distance(const std::vector<double>& x) {
  const int n = static_cast<int>(x.size());
  std::vector<double> d;
  d.reserve(static_cast<std::size_t>(n * (n - 1) / 2));
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dij = std::abs(x[i] - x[j]);
      if (std::isfinite(dij)) d.push_back(dij);
    }
  }
  if (d.empty()) return 0.001;
  const std::size_t mid = d.size() / 2U;
  std::nth_element(d.begin(), d.begin() + static_cast<std::ptrdiff_t>(mid), d.end());
  const double sigma = d[mid] / std::sqrt(2.0);
  if (!std::isfinite(sigma) || sigma <= 0.0) return 0.001;
  return sigma;
}

inline double bandwidth_univariate(const std::vector<double>& x, const int bandwidth_code) {
  if (bandwidth_code == 0) {
    return upper_median_pairwise_distance(x);
  }

  const int n = static_cast<int>(x.size());
  const double sd = vec_sd(x);
  if (!std::isfinite(sd) || sd <= HSIC_EPS) return 1.0;
  const double scale = std::pow(static_cast<double>(n), -0.2);
  if (bandwidth_code == 1) return std::max(1.06 * sd * scale, HSIC_EPS);
  if (bandwidth_code == 2) return std::max(sd * scale, HSIC_EPS);
  return std::max(sd, HSIC_EPS);
}

inline arma::mat compute_kernel_univariate(
  const std::vector<double>& x,
  const int kernel_code,
  const int bandwidth_code
) {
  const int n = static_cast<int>(x.size());
  arma::mat K(n, n, arma::fill::zeros);
  const double bw = bandwidth_univariate(x, bandwidth_code);
  const double bw2 = bw * bw;

  for (int i = 0; i < n; ++i) {
    K(i, i) = (kernel_code == 1) ? x[i] * x[i] :
      (kernel_code == 3) ? std::pow(x[i] * x[i] + 1.0, 2.0) : 1.0;
    for (int j = i + 1; j < n; ++j) {
      double kij;
      if (kernel_code == 1) {
        kij = x[i] * x[j];
      } else if (kernel_code == 2) {
        kij = std::exp(-std::abs(x[i] - x[j]) / bw);
      } else if (kernel_code == 3) {
        kij = std::pow(x[i] * x[j] + 1.0, 2.0);
      } else {
        const double d = x[i] - x[j];
        kij = std::exp(-(d * d) / (2.0 * bw2));
      }
      K(i, j) = kij;
      K(j, i) = kij;
    }
  }
  return K;
}

inline void center_gram_inplace(arma::mat& K) {
  const int n = static_cast<int>(K.n_rows);
  arma::vec row_means = arma::mean(K, 1);
  const double grand_mean = arma::mean(row_means);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      K(i, j) = K(i, j) - row_means[i] - row_means[j] + grand_mean;
    }
  }
}

inline void zero_diag_inplace(arma::mat& K) {
  K.diag().zeros();
}

inline double frobenius_inner(const arma::mat& A, const arma::mat& B) {
  return arma::accu(A % B);
}

inline double hsic_biased_from_grams(const arma::mat& Kc, const arma::mat& Lc) {
  const double n = static_cast<double>(Kc.n_rows);
  return frobenius_inner(Kc, Lc) / (n * n);
}

inline double hsic_biased_uncentered_from_stats(
  const double inner,
  const std::vector<double>& rowK,
  const std::vector<double>& rowL,
  const double sumK,
  const double sumL
) {
  const double n = static_cast<double>(rowK.size());
  double rowdot = 0.0;
  for (std::size_t i = 0; i < rowK.size(); ++i) {
    rowdot += rowK[i] * rowL[i];
  }
  const double n2 = n * n;
  const double out =
    inner / n2 +
    (sumK * sumL) / (n2 * n2) -
    2.0 * rowdot / (n2 * n);
  return (std::abs(out) <= HSIC_EPS) ? 0.0 : out;
}

inline double kernel_value_univariate(
  const double x,
  const double y,
  const int kernel_code,
  const double scale
) {
  if (kernel_code == 1) return x * y;
  if (kernel_code == 2) return std::exp(-std::abs(x - y) * scale);
  if (kernel_code == 3) {
    const double z = x * y + 1.0;
    return z * z;
  }
  const double d = x - y;
  return std::exp(-(d * d) * scale);
}

inline double kernel_diag_univariate(
  const double x,
  const int kernel_code
) {
  if (kernel_code == 1) return x * x;
  if (kernel_code == 3) {
    const double z = x * x + 1.0;
    return z * z;
  }
  return 1.0;
}

struct HsicRawVectorPrepared {
  std::vector<double> values;
  double bandwidth = 1.0;
  double scale = 1.0;
};

inline HsicRawVectorPrepared prepare_raw_vector(
  const std::vector<double>& x,
  const int kernel_code,
  const int bandwidth_code
) {
  HsicRawVectorPrepared out;
  out.values = x;
  if (kernel_code == 0 || kernel_code == 2) {
    out.bandwidth = bandwidth_univariate(x, bandwidth_code);
    out.scale = (kernel_code == 0) ?
      1.0 / (2.0 * out.bandwidth * out.bandwidth) :
      1.0 / out.bandwidth;
  }
  return out;
}

inline double gaussian_median_scale_from_distances(std::vector<double>& d) {
  if (d.empty()) return 1.0 / (2.0 * 0.001 * 0.001);
  const std::size_t mid = d.size() / 2U;
  std::nth_element(d.begin(), d.begin() + static_cast<std::ptrdiff_t>(mid), d.end());
  const double med = d[mid];
  if (!std::isfinite(med) || med <= 0.0) return 1.0 / (2.0 * 0.001 * 0.001);
  return 1.0 / (med * med);
}

inline arma::mat hsic_raw_biased_two_column_gaussian_median_impl(const arma::mat& X) {
  const int n = static_cast<int>(X.n_rows);
  if (n < 2) Rcpp::stop("Sample size is too small for the requested HSIC estimator");

  const double* xv = X.colptr(0);
  const double* yv = X.colptr(1);

  std::vector<double> dx;
  std::vector<double> dy;
  dx.reserve(static_cast<std::size_t>(n * (n - 1) / 2));
  dy.reserve(static_cast<std::size_t>(n * (n - 1) / 2));

  for (int i = 0; i < n; ++i) {
    const double xi = xv[i];
    const double yi = yv[i];
    for (int j = i + 1; j < n; ++j) {
      const double xij = std::abs(xi - xv[j]);
      const double yij = std::abs(yi - yv[j]);
      if (std::isfinite(xij)) dx.push_back(xij);
      if (std::isfinite(yij)) dy.push_back(yij);
    }
  }

  const double xscale = gaussian_median_scale_from_distances(dx);
  const double yscale = gaussian_median_scale_from_distances(dy);

  std::vector<double> rowK(static_cast<std::size_t>(n), 1.0);
  std::vector<double> rowL(static_cast<std::size_t>(n), 1.0);
  double innerKL = static_cast<double>(n);
  double innerKK = static_cast<double>(n);
  double innerLL = static_cast<double>(n);
  double sumK = static_cast<double>(n);
  double sumL = static_cast<double>(n);

  double* rowKp = rowK.data();
  double* rowLp = rowL.data();

  for (int i = 0; i < n; ++i) {
    const double xi = xv[i];
    const double yi = yv[i];
    for (int j = i + 1; j < n; ++j) {
      const double xdiff = xi - xv[j];
      const double ydiff = yi - yv[j];
      const double kij = std::exp(-(xdiff * xdiff) * xscale);
      const double lij = std::exp(-(ydiff * ydiff) * yscale);
      innerKL += 2.0 * kij * lij;
      innerKK += 2.0 * kij * kij;
      innerLL += 2.0 * lij * lij;
      rowKp[i] += kij;
      rowKp[j] += kij;
      rowLp[i] += lij;
      rowLp[j] += lij;
      sumK += 2.0 * kij;
      sumL += 2.0 * lij;
    }
  }

  arma::mat out(2, 2);
  out(0, 0) = hsic_biased_uncentered_from_stats(innerKK, rowK, rowK, sumK, sumK);
  out(1, 1) = hsic_biased_uncentered_from_stats(innerLL, rowL, rowL, sumL, sumL);
  const double xy = hsic_biased_uncentered_from_stats(innerKL, rowK, rowL, sumK, sumL);
  out(0, 1) = xy;
  out(1, 0) = xy;
  return out;
}

inline double hsic_biased_raw_pair_direct(
  const HsicRawVectorPrepared& x,
  const HsicRawVectorPrepared& y,
  const int kernel_code
) {
  const int n = static_cast<int>(x.values.size());
  std::vector<double> rowK(static_cast<std::size_t>(n), 0.0);
  std::vector<double> rowL(static_cast<std::size_t>(n), 0.0);
  double inner = 0.0;
  double sumK = 0.0;
  double sumL = 0.0;

  for (int i = 0; i < n; ++i) {
    const double xi = x.values[static_cast<std::size_t>(i)];
    const double yi = y.values[static_cast<std::size_t>(i)];
    const double kii = kernel_diag_univariate(xi, kernel_code);
    const double lii = kernel_diag_univariate(yi, kernel_code);
    inner += kii * lii;
    rowK[static_cast<std::size_t>(i)] += kii;
    rowL[static_cast<std::size_t>(i)] += lii;
    sumK += kii;
    sumL += lii;

    for (int j = i + 1; j < n; ++j) {
      const double kij = kernel_value_univariate(
        xi,
        x.values[static_cast<std::size_t>(j)],
        kernel_code,
        x.scale
      );
      const double lij = kernel_value_univariate(
        yi,
        y.values[static_cast<std::size_t>(j)],
        kernel_code,
        y.scale
      );
      inner += 2.0 * kij * lij;
      rowK[static_cast<std::size_t>(i)] += kij;
      rowK[static_cast<std::size_t>(j)] += kij;
      rowL[static_cast<std::size_t>(i)] += lij;
      rowL[static_cast<std::size_t>(j)] += lij;
      sumK += 2.0 * kij;
      sumL += 2.0 * lij;
    }
  }

  return hsic_biased_uncentered_from_stats(
    inner,
    rowK,
    rowL,
    sumK,
    sumL
  );
}

inline arma::mat hsic_raw_biased_two_column_impl(
  const arma::mat& X,
  const int kernel_code,
  const int bandwidth_code
) {
  const int n = static_cast<int>(X.n_rows);
  if (n < 2) Rcpp::stop("Sample size is too small for the requested HSIC estimator");
  if (kernel_code == 0 && bandwidth_code == 0) {
    return hsic_raw_biased_two_column_gaussian_median_impl(X);
  }

  std::vector<double> x(static_cast<std::size_t>(n));
  std::vector<double> y(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    x[static_cast<std::size_t>(i)] = X(i, 0);
    y[static_cast<std::size_t>(i)] = X(i, 1);
  }

  HsicRawVectorPrepared px;
  HsicRawVectorPrepared py;
  px = prepare_raw_vector(x, kernel_code, bandwidth_code);
  py = prepare_raw_vector(y, kernel_code, bandwidth_code);
  std::vector<double> rowK(static_cast<std::size_t>(n), 0.0);
  std::vector<double> rowL(static_cast<std::size_t>(n), 0.0);
  double innerKL = 0.0;
  double innerKK = 0.0;
  double innerLL = 0.0;
  double sumK = 0.0;
  double sumL = 0.0;

  const double* xv = px.values.data();
  const double* yv = py.values.data();
  double* rowKp = rowK.data();
  double* rowLp = rowL.data();

  if (kernel_code == 0) {
    for (int i = 0; i < n; ++i) {
      innerKL += 1.0;
      innerKK += 1.0;
      innerLL += 1.0;
      rowKp[i] += 1.0;
      rowLp[i] += 1.0;
      sumK += 1.0;
      sumL += 1.0;

      const double xi = xv[i];
      const double yi = yv[i];
      for (int j = i + 1; j < n; ++j) {
        const double dx = xi - xv[j];
        const double dy = yi - yv[j];
        const double kij = std::exp(-(dx * dx) * px.scale);
        const double lij = std::exp(-(dy * dy) * py.scale);
        innerKL += 2.0 * kij * lij;
        innerKK += 2.0 * kij * kij;
        innerLL += 2.0 * lij * lij;
        rowKp[i] += kij;
        rowKp[j] += kij;
        rowLp[i] += lij;
        rowLp[j] += lij;
        sumK += 2.0 * kij;
        sumL += 2.0 * lij;
      }
    }
  } else {
    for (int i = 0; i < n; ++i) {
      const double xi = xv[i];
      const double yi = yv[i];
      const double kii = kernel_diag_univariate(xi, kernel_code);
      const double lii = kernel_diag_univariate(yi, kernel_code);
      innerKL += kii * lii;
      innerKK += kii * kii;
      innerLL += lii * lii;
      rowKp[i] += kii;
      rowLp[i] += lii;
      sumK += kii;
      sumL += lii;

      for (int j = i + 1; j < n; ++j) {
        const double kij = kernel_value_univariate(
          xi,
          xv[j],
          kernel_code,
          px.scale
        );
        const double lij = kernel_value_univariate(
          yi,
          yv[j],
          kernel_code,
          py.scale
        );
        innerKL += 2.0 * kij * lij;
        innerKK += 2.0 * kij * kij;
        innerLL += 2.0 * lij * lij;
        rowKp[i] += kij;
        rowKp[j] += kij;
        rowLp[i] += lij;
        rowLp[j] += lij;
        sumK += 2.0 * kij;
        sumL += 2.0 * lij;
      }
    }
  }

  arma::mat out(2, 2);
  out(0, 0) = hsic_biased_uncentered_from_stats(innerKK, rowK, rowK, sumK, sumK);
  out(1, 1) = hsic_biased_uncentered_from_stats(innerLL, rowL, rowL, sumL, sumL);
  const double xy = hsic_biased_uncentered_from_stats(innerKL, rowK, rowL, sumK, sumL);
  out(0, 1) = xy;
  out(1, 0) = xy;
  return out;
}

inline double hsic_unbiased_from_grams(const arma::mat& K0, const arma::mat& L0) {
  const int n_int = static_cast<int>(K0.n_rows);
  if (n_int < 4) return NA_REAL;
  const double n = static_cast<double>(n_int);
  const double tr = frobenius_inner(K0, L0);
  const arma::vec rK = arma::sum(K0, 1);
  const arma::vec rL = arma::sum(L0, 1);
  const double sumK = arma::accu(rK);
  const double sumL = arma::accu(rL);
  const double rowdot = arma::dot(rK, rL);
  return (
    tr +
      (sumK * sumL) / ((n - 1.0) * (n - 2.0)) -
      2.0 * rowdot / (n - 2.0)
  ) / (n * (n - 3.0));
}

inline double hsic_from_grams(const arma::mat& K, const arma::mat& L, const int estimator_code) {
  if (estimator_code == 1) return hsic_unbiased_from_grams(K, L);
  return hsic_biased_from_grams(K, L);
}

inline double normalise_hsic_value(
  const double raw,
  const double self_x,
  const double self_y,
  const bool normalise
) {
  if (!normalise) return raw;
  if (!(std::isfinite(raw) && std::isfinite(self_x) && std::isfinite(self_y) &&
        self_x > 0.0 && self_y > 0.0)) {
    return NA_REAL;
  }
  return clamp01_hsic(raw / std::sqrt(self_x * self_y));
}

inline arma::mat permute_kernel_rows_cols(
  const arma::mat& K,
  const std::vector<int>& perm
) {
  const int n = static_cast<int>(K.n_rows);
  arma::mat out(n, n);
  for (int i = 0; i < n; ++i) {
    const int pi = perm[static_cast<std::size_t>(i)];
    for (int j = 0; j < n; ++j) {
      out(i, j) = K(pi, perm[static_cast<std::size_t>(j)]);
    }
  }
  return out;
}

struct HsicPrepared {
  arma::mat gram;
  double self = NA_REAL;
  int n_complete = 0;
};

inline HsicPrepared prepare_column(
  const std::vector<double>& x,
  const int kernel_code,
  const int bandwidth_code,
  const int estimator_code
) {
  HsicPrepared out;
  out.n_complete = static_cast<int>(x.size());
  out.gram = compute_kernel_univariate(x, kernel_code, bandwidth_code);
  if (estimator_code == 1) {
    zero_diag_inplace(out.gram);
  } else {
    center_gram_inplace(out.gram);
  }
  out.self = hsic_from_grams(out.gram, out.gram, estimator_code);
  return out;
}

inline Rcpp::List hsic_complete_impl(
  const arma::mat& X,
  const int kernel_code,
  const int bandwidth_code,
  const bool normalise,
  const int estimator_code,
  const bool return_inference,
  const int B,
  const int seed_base
) {
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);
  const int min_n = (estimator_code == 1) ? 4 : 2;
  if (n < min_n) Rcpp::stop("Sample size is too small for the requested HSIC estimator");

  std::vector<HsicPrepared> cols(static_cast<std::size_t>(p));

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    std::vector<double> x(static_cast<std::size_t>(n));
    for (int r = 0; r < n; ++r) x[static_cast<std::size_t>(r)] = X(r, j);
    cols[static_cast<std::size_t>(j)] = prepare_column(
      x, kernel_code, bandwidth_code, estimator_code
    );
  }

  Rcpp::NumericMatrix est(p, p);
  Rcpp::NumericMatrix raw(p, p);
  Rcpp::NumericMatrix n_complete(p, p);
  Rcpp::NumericMatrix p_value(p, p);

  for (int j = 0; j < p; ++j) {
    const double self = cols[static_cast<std::size_t>(j)].self;
    raw(j, j) = self;
    est(j, j) = normalise ? (std::isfinite(self) && self > 0.0 ? 1.0 : NA_REAL) : self;
    n_complete(j, j) = n;
    p_value(j, j) = NA_REAL;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 1; j < p; ++j) {
    for (int i = 0; i < j; ++i) {
      const HsicPrepared& xi = cols[static_cast<std::size_t>(i)];
      const HsicPrepared& xj = cols[static_cast<std::size_t>(j)];
      const double raw_ij = hsic_from_grams(xi.gram, xj.gram, estimator_code);
      const double est_ij = normalise_hsic_value(raw_ij, xi.self, xj.self, normalise);

      raw(i, j) = raw_ij;
      raw(j, i) = raw_ij;
      est(i, j) = est_ij;
      est(j, i) = est_ij;
      n_complete(i, j) = n;
      n_complete(j, i) = n;

      double pval = NA_REAL;
      if (return_inference && std::isfinite(raw_ij) && B > 0) {
        std::mt19937_64 rng(
          static_cast<std::uint64_t>(seed_base) +
            104729ULL * static_cast<std::uint64_t>(i + 1) +
            1000003ULL * static_cast<std::uint64_t>(j + 1)
        );
        std::vector<int> perm(static_cast<std::size_t>(n));
        std::iota(perm.begin(), perm.end(), 0);
        int ge = 0;
        for (int b = 0; b < B; ++b) {
          std::shuffle(perm.begin(), perm.end(), rng);
          const arma::mat permuted = permute_kernel_rows_cols(xj.gram, perm);
          const double stat_b = hsic_from_grams(xi.gram, permuted, estimator_code);
          if (std::isfinite(stat_b) && stat_b >= raw_ij - 1e-15) ++ge;
        }
        pval = (static_cast<double>(ge) + 1.0) / (static_cast<double>(B) + 1.0);
      }
      p_value(i, j) = pval;
      p_value(j, i) = pval;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("est") = est,
    Rcpp::Named("raw") = raw,
    Rcpp::Named("n_complete") = n_complete,
    Rcpp::Named("p_value") = p_value
  );
}

inline Rcpp::List hsic_pairwise_impl(
  const arma::mat& X,
  const int kernel_code,
  const int bandwidth_code,
  const bool normalise,
  const int estimator_code,
  const bool return_inference,
  const int B,
  const int seed_base
) {
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);
  const int min_n = (estimator_code == 1) ? 4 : 2;

  Rcpp::NumericMatrix est(p, p);
  Rcpp::NumericMatrix raw(p, p);
  Rcpp::NumericMatrix n_complete(p, p);
  Rcpp::NumericMatrix p_value(p, p);

  for (int j = 0; j < p; ++j) {
    int nj = 0;
    for (int r = 0; r < n; ++r) {
      if (std::isfinite(X(r, j))) ++nj;
    }
    n_complete(j, j) = nj;
    raw(j, j) = NA_REAL;
    est(j, j) = normalise ? NA_REAL : NA_REAL;
    p_value(j, j) = NA_REAL;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 1; j < p; ++j) {
    for (int i = 0; i < j; ++i) {
      std::vector<double> xi;
      std::vector<double> xj;
      xi.reserve(static_cast<std::size_t>(n));
      xj.reserve(static_cast<std::size_t>(n));
      for (int r = 0; r < n; ++r) {
        const double a = X(r, i);
        const double b = X(r, j);
        if (std::isfinite(a) && std::isfinite(b)) {
          xi.push_back(a);
          xj.push_back(b);
        }
      }

      const int nij = static_cast<int>(xi.size());
      n_complete(i, j) = nij;
      n_complete(j, i) = nij;

      double raw_ij = NA_REAL;
      double est_ij = NA_REAL;
      double pval = NA_REAL;

      if (nij >= min_n) {
        HsicPrepared px = prepare_column(xi, kernel_code, bandwidth_code, estimator_code);
        HsicPrepared py = prepare_column(xj, kernel_code, bandwidth_code, estimator_code);
        raw_ij = hsic_from_grams(px.gram, py.gram, estimator_code);
        est_ij = normalise_hsic_value(raw_ij, px.self, py.self, normalise);

        if (return_inference && std::isfinite(raw_ij) && B > 0) {
          std::mt19937_64 rng(
            static_cast<std::uint64_t>(seed_base) +
              104729ULL * static_cast<std::uint64_t>(i + 1) +
              1000003ULL * static_cast<std::uint64_t>(j + 1)
          );
          std::vector<int> perm(static_cast<std::size_t>(nij));
          std::iota(perm.begin(), perm.end(), 0);
          int ge = 0;
          for (int b = 0; b < B; ++b) {
            std::shuffle(perm.begin(), perm.end(), rng);
            const arma::mat permuted = permute_kernel_rows_cols(py.gram, perm);
            const double stat_b = hsic_from_grams(px.gram, permuted, estimator_code);
            if (std::isfinite(stat_b) && stat_b >= raw_ij - 1e-15) ++ge;
          }
          pval = (static_cast<double>(ge) + 1.0) / (static_cast<double>(B) + 1.0);
        }
      }

      raw(i, j) = raw_ij;
      raw(j, i) = raw_ij;
      est(i, j) = est_ij;
      est(j, i) = est_ij;
      p_value(i, j) = pval;
      p_value(j, i) = pval;
    }
  }

  for (int j = 0; j < p; ++j) {
    if (n_complete(j, j) >= min_n) {
      std::vector<double> xj;
      xj.reserve(static_cast<std::size_t>(n));
      for (int r = 0; r < n; ++r) {
        const double v = X(r, j);
        if (std::isfinite(v)) xj.push_back(v);
      }
      HsicPrepared pj = prepare_column(xj, kernel_code, bandwidth_code, estimator_code);
      raw(j, j) = pj.self;
      est(j, j) = normalise ? (std::isfinite(pj.self) && pj.self > 0.0 ? 1.0 : NA_REAL) : pj.self;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("est") = est,
    Rcpp::Named("raw") = raw,
    Rcpp::Named("n_complete") = n_complete,
    Rcpp::Named("p_value") = p_value
  );
}

inline arma::mat hsic_raw_biased_complete_impl(
  const arma::mat& X,
  const int kernel_code,
  const int bandwidth_code
) {
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);
  if (n < 2) Rcpp::stop("Sample size is too small for the requested HSIC estimator");

  std::vector<HsicRawVectorPrepared> cols(static_cast<std::size_t>(p));

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    std::vector<double> x(static_cast<std::size_t>(n));
    for (int r = 0; r < n; ++r) x[static_cast<std::size_t>(r)] = X(r, j);
    cols[static_cast<std::size_t>(j)] = prepare_raw_vector(
      x, kernel_code, bandwidth_code
    );
  }

  arma::mat out(p, p, arma::fill::zeros);

  for (int j = 0; j < p; ++j) {
    const HsicRawVectorPrepared& xj = cols[static_cast<std::size_t>(j)];
    out(j, j) = hsic_biased_raw_pair_direct(xj, xj, kernel_code);
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 1; j < p; ++j) {
    for (int i = 0; i < j; ++i) {
      const double raw_ij = hsic_biased_raw_pair_direct(
        cols[static_cast<std::size_t>(i)],
        cols[static_cast<std::size_t>(j)],
        kernel_code
      );
      out(i, j) = raw_ij;
      out(j, i) = raw_ij;
    }
  }

  return out;
}

} // namespace

// [[Rcpp::export]]
arma::mat hsic_matrix_cpp(
  const arma::mat& X,
  int kernel_code,
  int bandwidth_code,
  bool normalise,
  int estimator_code
) {
  if (X.n_rows < 2 || X.n_cols < 2) {
    Rcpp::stop("data must have at least two rows and two columns");
  }
  if (!X.is_finite()) {
    Rcpp::stop("data contains missing, NaN, or infinite values");
  }
  if (!normalise && estimator_code == 0) {
    if (X.n_cols == 2) {
      return hsic_raw_biased_two_column_impl(X, kernel_code, bandwidth_code);
    }
    return hsic_raw_biased_complete_impl(X, kernel_code, bandwidth_code);
  }
  Rcpp::List out = hsic_complete_impl(
    X, kernel_code, bandwidth_code, normalise, estimator_code,
    false, 0, 1
  );
  return Rcpp::as<arma::mat>(out["est"]);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix hsic_raw_biased_matrix_object_cpp(
  const arma::mat& X,
  int kernel_code,
  int bandwidth_code,
  std::string kernel,
  std::string bandwidth,
  Rcpp::CharacterVector names
) {
  if (X.n_rows < 2 || X.n_cols < 2) {
    Rcpp::stop("data must have at least two rows and two columns");
  }
  if (!X.is_finite()) {
    Rcpp::stop("data contains missing, NaN, or infinite values");
  }

  arma::mat out;
  if (X.n_cols == 2) {
    out = hsic_raw_biased_two_column_impl(X, kernel_code, bandwidth_code);
  } else {
    out = hsic_raw_biased_complete_impl(X, kernel_code, bandwidth_code);
  }

  Rcpp::NumericMatrix res = Rcpp::wrap(out);
  Rcpp::List dn = Rcpp::List::create(names, Rcpp::clone(names));
  res.attr("dimnames") = dn;

  Rcpp::NumericMatrix raw = Rcpp::clone(res);
  raw.attr("dimnames") = dn;

  res.attr("class") = Rcpp::CharacterVector::create(
    "corr_matrix",
    "hsic",
    "corr_result",
    "matrixCorr_dependence",
    "matrix"
  );
  res.attr("method") = "hsic";
  res.attr("description") = "Pairwise raw Hilbert-Schmidt independence criterion";
  res.attr("package") = "matrixCorr";
  res.attr("output") = "matrix";
  res.attr("threshold") = 0.0;
  res.attr("diag") = true;
  res.attr("corr_result") = true;
  res.attr("corr_output_class") = "corr_matrix";
  res.attr("corr_estimator_class") = "hsic";
  res.attr("corr_dim") = Rcpp::IntegerVector::create(out.n_rows, out.n_cols);
  res.attr("corr_dimnames") = dn;
  res.attr("corr_symmetric") = true;
  res.attr("kernel") = kernel;
  res.attr("bandwidth") = bandwidth;
  res.attr("estimator") = "biased";
  res.attr("normalise") = false;
  res.attr("hsic_raw") = raw;
  return res;
}

// [[Rcpp::export]]
Rcpp::List hsic_matrix_pairwise_cpp(
  const arma::mat& X,
  int kernel_code,
  int bandwidth_code,
  bool normalise,
  int estimator_code,
  bool return_inference,
  int B,
  Rcpp::Nullable<Rcpp::IntegerVector> seed = R_NilValue
) {
  if (X.n_rows < 2 || X.n_cols < 2) {
    Rcpp::stop("data must have at least two rows and two columns");
  }
  if (return_inference && B < 1) {
    Rcpp::stop("B must be positive when permutation inference is requested");
  }

  int seed_base = 1234567;
  if (seed.isNotNull()) {
    Rcpp::IntegerVector seed_vec(seed);
    if (seed_vec.size() > 0 && seed_vec[0] != NA_INTEGER) seed_base = seed_vec[0];
  } else {
    seed_base = static_cast<int>(std::random_device{}());
  }

  if (X.is_finite()) {
    return hsic_complete_impl(
      X, kernel_code, bandwidth_code, normalise, estimator_code,
      return_inference, B, seed_base
    );
  }

  return hsic_pairwise_impl(
    X, kernel_code, bandwidth_code, normalise, estimator_code,
    return_inference, B, seed_base
  );
}
