// Thiago de Paula Oliveira
// src-bicor.cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <limits>     // for std::numeric_limits
#include <cmath>      // for std::floor, std::ceil
#include <algorithm>  // for std::max
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// ---------- helpers: weighted quantile/median on (optionally) sorted inputs ----------

// Weighted quantile for sorted values 's'
// If not sorted, sort both together before calling this function.
inline double weighted_quantile_sorted(const arma::vec &s, const arma::vec &ws, double p) {
  const std::size_t n = s.n_elem;
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  if (p <= 0.0) return s.front();
  if (p >= 1.0) return s.back();
  const double W = arma::accu(ws);
  if (!(W > 0.0)) return std::numeric_limits<double>::quiet_NaN();

  // cumulative weight target
  const double target = p * W;
  double csum = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    csum += ws[i];
    if (csum >= target) {
      // Optional linear interpolation to the previous point
      if (i == 0) return s[0];
      const double csum_prev = csum - ws[i];
      const double wgap = ws[i];
      if (wgap <= 0.0) return s[i];
      const double frac = (target - csum_prev) / wgap; // in [0,1]
      return s[i-1] + frac * (s[i] - s[i-1]);
    }
  }
  return s.back();
}

inline double weighted_median(const arma::vec &x, const arma::vec &w) {
  arma::uvec ord = arma::stable_sort_index(x);       // stable to keep determinism
  arma::vec xs = x.elem(ord), ws = w.elem(ord);
  return weighted_quantile_sorted(xs, ws, 0.5);
}

inline double weighted_mad(const arma::vec &x, const arma::vec &w, double med) {
  arma::vec dev = arma::abs(x - med);
  return weighted_median(dev, w);
}

// Standardise one column with *weights* (no NA, already subset if needed).
static void standardise_bicor_column_weighted(const arma::vec &x,
                                              const arma::vec &wobs,
                                              arma::vec &z,
                                              int pearson_fallback_mode,
                                              double c_const,
                                              double maxPOutliers,
                                              bool &col_is_valid) {
  const std::size_t n = x.n_elem;
  z.zeros();
  col_is_valid = false;

  if (n == 0) { z.fill(arma::datum::nan); return; }

  // If all observation weights are zero → degenerate
  const double Wtot = arma::accu(wobs);
  if (!(Wtot > 0.0)) { z.fill(arma::datum::nan); return; }

  if (pearson_fallback_mode == 2) {
    const double mu = arma::dot(x, wobs) / Wtot;
    arma::vec centered = x - mu;
    const double denom2 = arma::dot(centered % wobs, centered % wobs);
    if (denom2 > 0.0) {
      z = (centered % wobs) / std::sqrt(denom2);
      col_is_valid = true;
    } else {
      z.fill(arma::datum::nan);
    }
    return;
  }

  // Weighted median and MAD
  const double med = weighted_median(x, wobs);
  const double mad = weighted_mad(x, wobs, med);

  if (!(mad > 0.0)) {
    if (pearson_fallback_mode == 0) { z.fill(arma::datum::nan); return; }
    // Pearson fallback with observation weights: mean as weighted mean
    const double mu = arma::dot(x, wobs) / Wtot;
    arma::vec centered = x - mu;
    // Use weight in norm
    double denom2 = arma::dot(centered % wobs, centered % wobs);
    if (!(denom2 > 0.0)) { z.fill(arma::datum::nan); return; }
    z = (centered % wobs) / std::sqrt(denom2);
    col_is_valid = true;
    return;
  }

  // Weighted side-cap quantiles
  double scale_neg = 1.0, scale_pos = 1.0;
  if (maxPOutliers < 1.0) {
    arma::uvec ord = arma::stable_sort_index(x);
    arma::vec xs = x.elem(ord), ws = wobs.elem(ord);
    const double qL = weighted_quantile_sorted(xs, ws, maxPOutliers);
    const double qU = weighted_quantile_sorted(xs, ws, 1.0 - maxPOutliers);
    const double uL = (qL - med) / (c_const * mad);
    const double uU = (qU - med) / (c_const * mad);
    if (std::abs(uL) > 1.0) scale_neg = std::abs(uL);
    if (std::abs(uU) > 1.0) scale_pos = std::abs(uU);
  }

  arma::vec xm = x - med;
  arma::vec u  = xm / (c_const * mad);
  if (maxPOutliers < 1.0 && (scale_neg > 1.0 || scale_pos > 1.0)) {
    for (std::size_t i = 0; i < n; ++i) {
      if (xm[i] < 0.0) u[i] /= scale_neg;
      else if (xm[i] > 0.0) u[i] /= scale_pos;
    }
  }

  arma::vec wt(n, arma::fill::zeros);
  for (std::size_t i = 0; i < n; ++i) {
    const double a = u[i];
    if (std::abs(a) < 1.0) {
      const double t = (1.0 - a*a);
      wt[i] = t * t;
    }
  }
  // Multiply Tukey weights by observation weights
  wt %= wobs;

  arma::vec r = xm % wt;
  const double denom2 = arma::dot(r, r);
  if (!(denom2 > 0.0)) {
    if (pearson_fallback_mode >= 1) {
      // Weighted Pearson fallback
      const double mu = arma::dot(x, wobs) / Wtot;
      arma::vec centered = x - mu;
      double d2 = arma::dot((centered % wobs), (centered % wobs));
      if (d2 > 0.0) {
        z = (centered % wobs) / std::sqrt(d2);
        col_is_valid = true;
        return;
      }
    }
    z.fill(arma::datum::nan);
    return;
  }
  z = r / std::sqrt(denom2);
  col_is_valid = true;
}

// Linear interpolation for quantiles on a sorted vector (p in [0,1]).
inline double quantile_sorted(const arma::vec &s, double p) {
  const std::size_t n = s.n_elem;
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  if (p <= 0.0) return s.front();
  if (p >= 1.0) return s.back();
  double pos = p * (n - 1);
  std::size_t lo = static_cast<std::size_t>(std::floor(pos));
  std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
  double frac = pos - lo;
  return (1.0 - frac) * s[lo] + frac * s[hi];
}

// Implements bicor's weighting
// 0 = robust, 1 = pearson fallback, -1 = degenerate (all NaN)
// Requirements is no NA/Inf in x.
static void standardise_bicor_column(const arma::vec &x,
                                     arma::vec &z,
                                     // 0=none(NA), 1=individual, 2=all
                                     int pearson_fallback_mode,
                                     double c_const,
                                     double maxPOutliers,
                                     bool &col_is_valid,
                                     int n_threads_unused=1) {
  const std::size_t n = x.n_elem;
  col_is_valid = false;
  z.zeros();

  if (pearson_fallback_mode == 2) {
    const double mu = arma::mean(x);
    arma::vec centered = x - mu;
    const double denom2 = arma::dot(centered, centered);
    if (denom2 > 0.0) {
      z = centered / std::sqrt(denom2);
      col_is_valid = true;
    } else {
      z.fill(arma::datum::nan);
    }
    return;
  }

  // Median and MAD
  arma::vec xc = x;
  xc = arma::sort(xc);
  const double med = arma::median(xc);
  arma::vec absdev = arma::abs(x - med);
  const double mad = arma::median(arma::sort(absdev));

  // If MAD == 0, decide Pearson fallback or mark degenerate.
  if (!(mad > 0.0)) {
    if (pearson_fallback_mode == 0) { // none -> NA column
      z.fill(arma::datum::nan);
      return;
    }
    // Pearson standardisation for this column
    const double mu = arma::mean(x);
    arma::vec centered = x - mu;
    double denom2 = arma::dot(centered, centered);
    if (!(denom2 > 0.0)) { // constant column -> NA
      z.fill(arma::datum::nan);
      return;
    }
    z = centered / std::sqrt(denom2);
    col_is_valid = true;
    return;
  }

  // rescale u on each side so that the chosen quantiles
  // just reach |u| = 1 if they would otherwise be >= 1 (Langfelder & Horvath).
  double scale_neg = 1.0, scale_pos = 1.0;
  if (maxPOutliers < 1.0) {
    const double qL = quantile_sorted(xc, maxPOutliers);
    const double qU = quantile_sorted(xc, 1.0 - maxPOutliers);
    const double uL = (qL - med) / (c_const * mad);
    const double uU = (qU - med) / (c_const * mad);
    if (std::abs(uL) > 1.0) scale_neg = std::abs(uL);
    if (std::abs(uU) > 1.0) scale_pos = std::abs(uU);
  }

  // Compute weighted, centred vector r_i = (x_i - med) * w_i.
  arma::vec xm = x - med;
  arma::vec u = xm / (c_const * mad);
  // Side-specific rescaling when requested
  if (maxPOutliers < 1.0 && (scale_neg > 1.0 || scale_pos > 1.0)) {
    for (std::size_t i = 0; i < n; ++i) {
      if (xm[i] < 0.0) u[i] /= scale_neg;
      else if (xm[i] > 0.0) u[i] /= scale_pos;
    }
  }

  arma::vec w(n, fill::zeros);
  for (std::size_t i = 0; i < n; ++i) {
    double a = u[i];
    if (std::abs(a) < 1.0) {
      double t = (1.0 - a*a);
      w[i] = t * t;
    } // else 0
  }
  arma::vec r = xm % w;

  // L2 norm of r
  double denom2 = arma::dot(r, r);
  if (!(denom2 > 0.0)) {
    if (pearson_fallback_mode >= 1) {
      const double mu = arma::mean(x);
      arma::vec centered = x - mu;
      double d2 = arma::dot(centered, centered);
      if (d2 > 0.0) {
        z = centered / std::sqrt(d2);
        col_is_valid = true;
        return;
      }
    }
    z.fill(arma::datum::nan);
    return;
  }

  z = r / std::sqrt(denom2);
  col_is_valid = true;
}

// Biweight mid-correlation matrix
//
// Computes the biweight mid-correlation for all column pairs of X.
// No missing/Inf values are allowed in X for this fast variant.
// Set `pearson_fallback=1` (individual) to obtain the hybrid correlation when
// MAD==0 in a column.
// X Numeric matrix (rows = observations, cols = variables).
// c_const Tuning constant multiplying raw MAD (default 9.0).
// maxPOutliers Maximum proportion of low *and* high outliers to permit (>0 to cap).
// If < 1, columns are side-rescaled so that these quantiles (if beyond |u|>=1) map to |u|=1.
// pearson_fallback 0 = never (returns NA for MAD==0); 1 = individual fallback when needed;
// 2 = force Pearson for all columns (i.e., ordinary Pearson).
// n_threads Number of OpenMP threads (>=1).
// Symmetric p x p matrix of correlations. Diagonals are 1 where defined, NA if
// the column is degenerate.
// [[Rcpp::export]]
arma::mat bicor_matrix_cpp(const arma::mat &X,
                           const double c_const = 9.0,
                           const double maxPOutliers = 1.0,
                           const int pearson_fallback = 1,
                           const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");
  if (!X.is_finite()) stop("X contains NA/NaN/Inf; please handle missingness upstream.");

  // Standardised columns
  arma::mat Z(n, p, fill::zeros);
  std::vector<bool> col_valid(p, false);

  // Standardise each column in parallel
  #ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
  #pragma omp parallel for schedule(static)
  #endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    bool ok = false;
    arma::vec zcol(n, fill::zeros);
    standardise_bicor_column(X.col(j), zcol, pearson_fallback, c_const, maxPOutliers, ok);
    Z.col(j) = zcol;
    col_valid[static_cast<std::size_t>(j)] = ok;
  }

  // Correlation matrix R = Z'Z
  arma::mat R(p, p, fill::zeros);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {
      double val = arma::dot(Z.col(j), Z.col(k));
      // Clamp to [-1, 1] if finite
      if (std::isfinite(val)) {
        if (val > 1.0) val = 1.0;
        else if (val < -1.0) val = -1.0;
      }
      R(j, k) = val;
      if (k != static_cast<std::size_t>(j)) R(k, j) = val;
    }
  }

  // Set diagonals
  R.diag().ones();
  return R;
}

// bicor for two vectors (hybrid if one side falls back to Pearson).
// [[Rcpp::export]]
double bicor_vec_cpp(const arma::vec &x, const arma::vec &y,
                     const double c_const = 9.0,
                     const double maxPOutliers = 1.0,
                     const int pearson_fallback = 1) {
  if (x.n_elem != y.n_elem) stop("x and y must have the same length.");
  if (!x.is_finite() || !y.is_finite()) stop("x or y contains NA/NaN/Inf.");

  arma::vec zx(x.n_elem, fill::zeros), zy(y.n_elem, fill::zeros);
  bool okx=false, oky=false;
  standardise_bicor_column(x, zx, pearson_fallback, c_const, maxPOutliers, okx);
  standardise_bicor_column(y, zy, pearson_fallback, c_const, maxPOutliers, oky);

  double val = arma::dot(zx, zy);
  if (std::isfinite(val)) {
    if (val > 1.0) val = 1.0;
    else if (val < -1.0) val = -1.0;
  }
  // If either side is degenerate, return NA
  if (!okx || !oky) return std::numeric_limits<double>::quiet_NaN();
  return val;
}

// [[Rcpp::export]]
arma::mat bicor_matrix_pairwise_cpp(const arma::mat &X,
                                    const double c_const = 9.0,
                                    const double maxPOutliers = 1.0,
                                    const int pearson_fallback = 1,
                                    const int min_n = 5,
                                    const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");

  // initialise with NaN (constructor can't take 'datum::nan')
  arma::mat R(p, p, arma::fill::none);
  R.fill(arma::datum::nan);

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {

      // rows where both columns are finite
      std::vector<arma::uword> idx; idx.reserve(n);
      for (arma::uword i = 0; i < n; ++i) {
        const double a = X(i, j), b = X(i, k);
        if (std::isfinite(a) && std::isfinite(b)) idx.push_back(i);
      }
      if (idx.size() < static_cast<std::size_t>(min_n)) continue;

      arma::uvec I = arma::conv_to<arma::uvec>::from(idx);

      // materialise columns, then subselect with '.elem()'
      arma::vec colj = X.col(j);
      arma::vec colk = X.col(k);
      arma::vec xj = colj.elem(I);
      arma::vec xk = colk.elem(I);

      arma::vec zj(xj.n_elem, arma::fill::zeros), zk(xk.n_elem, arma::fill::zeros);
      bool okj=false, okk=false;
      standardise_bicor_column(xj, zj, pearson_fallback, c_const, maxPOutliers, okj);
      standardise_bicor_column(xk, zk, pearson_fallback, c_const, maxPOutliers, okk);

      double val = arma::datum::nan;
      if (okj && okk) {
        val = arma::dot(zj, zk);
        if (std::isfinite(val)) {
          if (val > 1.0) val = 1.0;
          else if (val < -1.0) val = -1.0;
        }
      }
      R(j,k) = val;
      if (k != static_cast<std::size_t>(j)) R(k,j) = val;
    }
  }

  // diagonals are defined; set to 1
  for (std::size_t j = 0; j < p; ++j) R(j,j) = 1.0;
  return R;
}


// [[Rcpp::export]]
arma::mat bicor_matrix_weighted_cpp(const arma::mat &X,
                                    const arma::vec &w,
                                    const double c_const = 9.0,
                                    const double maxPOutliers = 1.0,
                                    const int pearson_fallback = 1,
                                    const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");
  if (w.n_elem != n) stop("Length of weights `w` must match nrow(X).");
  if (!X.is_finite()) stop("X contains NA/NaN/Inf; use pairwise kernel or clean data.");
  if (!w.is_finite() || arma::any(w < 0)) stop("Weights must be finite and non-negative.");

  arma::mat Z(n, p, arma::fill::zeros);
  std::vector<bool> col_valid(p, false);

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    bool ok = false;
    arma::vec zcol(n, arma::fill::zeros);
    standardise_bicor_column_weighted(X.col(j), w, zcol, pearson_fallback, c_const, maxPOutliers, ok);
    Z.col(j) = zcol;
    col_valid[static_cast<std::size_t>(j)] = ok;
  }

  arma::mat R(p, p, arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {
      double val = arma::dot(Z.col(j), Z.col(k));
      if (std::isfinite(val)) {
        if (val > 1.0) val = 1.0;
        else if (val < -1.0) val = -1.0;
      }
      R(j, k) = val;
      if (k != static_cast<std::size_t>(j)) R(k, j) = val;
    }
  }
  R.diag().ones();
  return R;
}

// [[Rcpp::export]]
arma::mat bicor_matrix_weighted_pairwise_cpp(const arma::mat &X,
                                             const arma::vec &w,
                                             const double c_const = 9.0,
                                             const double maxPOutliers = 1.0,
                                             const int pearson_fallback = 1,
                                             const int min_n = 5,
                                             const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");
  if (w.n_elem != n)   stop("Length of weights `w` must match nrow(X).");
  if (!w.is_finite() || arma::any(w < 0)) stop("Weights must be finite and non-negative.");

  arma::mat R(p, p, arma::fill::none);
  R.fill(arma::datum::nan);

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {

      std::vector<arma::uword> idx; idx.reserve(n);
      for (arma::uword i = 0; i < n; ++i) {
        const double a = X(i, j), b = X(i, k);
        if (std::isfinite(a) && std::isfinite(b)) idx.push_back(i);
      }
      if (idx.size() < static_cast<std::size_t>(min_n)) continue;

      arma::uvec I = arma::conv_to<arma::uvec>::from(idx);

      arma::vec colj = X.col(j);
      arma::vec colk = X.col(k);
      arma::vec xj = colj.elem(I);
      arma::vec xk = colk.elem(I);
      arma::vec ww = w.elem(I);

      arma::vec zj(xj.n_elem, arma::fill::zeros), zk(xk.n_elem, arma::fill::zeros);
      bool okj=false, okk=false;
      standardise_bicor_column_weighted(xj, ww, zj, pearson_fallback, c_const, maxPOutliers, okj);
      standardise_bicor_column_weighted(xk, ww, zk, pearson_fallback, c_const, maxPOutliers, okk);

      double val = arma::datum::nan;
      if (okj && okk) {
        val = arma::dot(zj, zk);
        if (std::isfinite(val)) {
          if (val > 1.0) val = 1.0;
          else if (val < -1.0) val = -1.0;
        }
      }
      R(j,k) = val;
      if (k != static_cast<std::size_t>(j)) R(k,j) = val;
    }
  }

  for (std::size_t j = 0; j < p; ++j) R(j,j) = 1.0;
  return R;
}

