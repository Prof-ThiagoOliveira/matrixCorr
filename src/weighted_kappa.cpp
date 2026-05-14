#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>
#include "matrixCorr_omp.h"
#include "threshold_triplets.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

namespace {

struct WeightedKappaStats {
  double estimate = NA_REAL;
  int n_complete = 0;
  double observed_agreement = NA_REAL;
  double expected_agreement = NA_REAL;
  double se = NA_REAL;
  double lwr = NA_REAL;
  double upr = NA_REAL;
  double statistic = NA_REAL;
  double p_value = NA_REAL;
};

inline bool valid_code(const int code, const int n_levels) {
  return code >= 1 && code <= n_levels;
}

WeightedKappaStats compute_weighted_kappa_from_counts(
    const std::vector<int>& counts,
    const Rcpp::NumericMatrix& weights,
    const bool return_inference,
    const double conf_level) {
  WeightedKappaStats out;
  const int n_levels = weights.nrow();
  if (n_levels <= 0 || weights.ncol() != n_levels) {
    return out;
  }

  long long n_complete = 0;
  std::vector<double> q(static_cast<std::size_t>(n_levels), 0.0);
  std::vector<double> r(static_cast<std::size_t>(n_levels), 0.0);
  double observed = 0.0;

  for (int a = 0; a < n_levels; ++a) {
    for (int b = 0; b < n_levels; ++b) {
      const std::size_t idx = static_cast<std::size_t>(a * n_levels + b);
      const double cab = static_cast<double>(counts[idx]);
      n_complete += static_cast<long long>(counts[idx]);
      q[static_cast<std::size_t>(a)] += cab;
      r[static_cast<std::size_t>(b)] += cab;
      observed += weights(a, b) * cab;
    }
  }

  out.n_complete = static_cast<int>(n_complete);
  if (n_complete < 2 || n_levels < 1) {
    return out;
  }

  const double n = static_cast<double>(n_complete);
  observed /= n;
  for (int a = 0; a < n_levels; ++a) {
    q[static_cast<std::size_t>(a)] /= n;
    r[static_cast<std::size_t>(a)] /= n;
  }

  double expected = 0.0;
  for (int a = 0; a < n_levels; ++a) {
    for (int b = 0; b < n_levels; ++b) {
      expected += weights(a, b) * q[static_cast<std::size_t>(a)] * r[static_cast<std::size_t>(b)];
    }
  }

  out.observed_agreement = observed;
  out.expected_agreement = expected;

  const double denom = 1.0 - expected;
  if (!std::isfinite(denom) || denom <= 1e-15) {
    return out;
  }

  const double estimate = (observed - expected) / denom;
  if (!std::isfinite(estimate)) {
    return out;
  }
  out.estimate = estimate;

  if (!return_inference) {
    return out;
  }

  std::vector<double> wr(static_cast<std::size_t>(n_levels), 0.0);
  std::vector<double> wtq(static_cast<std::size_t>(n_levels), 0.0);
  for (int a = 0; a < n_levels; ++a) {
    for (int b = 0; b < n_levels; ++b) {
      wr[static_cast<std::size_t>(a)] += weights(a, b) * r[static_cast<std::size_t>(b)];
      wtq[static_cast<std::size_t>(b)] += weights(a, b) * q[static_cast<std::size_t>(a)];
    }
  }

  const double denom2 = denom * denom;
  double mean_g = 0.0;
  double mean_g2 = 0.0;
  for (int a = 0; a < n_levels; ++a) {
    for (int b = 0; b < n_levels; ++b) {
      const std::size_t idx = static_cast<std::size_t>(a * n_levels + b);
      const double pab = static_cast<double>(counts[idx]) / n;
      if (pab <= 0.0) {
        continue;
      }
      const double hab = wr[static_cast<std::size_t>(a)] + wtq[static_cast<std::size_t>(b)];
      const double gab = (denom * weights(a, b) + (observed - 1.0) * hab) / denom2;
      mean_g += pab * gab;
      mean_g2 += pab * gab * gab;
    }
  }

  const double variance = (mean_g2 - mean_g * mean_g) / n;
  if (!std::isfinite(variance)) {
    return out;
  }

  const double se = std::sqrt(std::max(variance, 0.0));
  if (!std::isfinite(se) || se <= 0.0) {
    return out;
  }

  out.se = se;
  out.statistic = estimate / se;
  out.p_value = 2.0 * R::pnorm5(std::fabs(out.statistic), 0.0, 1.0, 0, 0);
  const double z = R::qnorm5(0.5 * (1.0 + conf_level), 0.0, 1.0, 1, 0);
  out.lwr = std::max(-1.0, estimate - z * se);
  out.upr = std::min(1.0, estimate + z * se);
  return out;
}

WeightedKappaStats compute_weighted_kappa_pair(
    const Rcpp::IntegerVector& x,
    const Rcpp::IntegerVector& y,
    const Rcpp::NumericMatrix& weights,
    const bool return_inference,
    const double conf_level) {
  const int n_levels = weights.nrow();
  std::vector<int> counts(static_cast<std::size_t>(std::max(0, n_levels * n_levels)), 0);
  const R_xlen_t n = std::min(x.size(), y.size());
  for (R_xlen_t i = 0; i < n; ++i) {
    const int xi = x[i];
    const int yi = y[i];
    if (xi == NA_INTEGER || yi == NA_INTEGER) {
      continue;
    }
    if (!valid_code(xi, n_levels) || !valid_code(yi, n_levels)) {
      continue;
    }
    counts[static_cast<std::size_t>((xi - 1) * n_levels + (yi - 1))] += 1;
  }
  return compute_weighted_kappa_from_counts(counts, weights, return_inference, conf_level);
}

WeightedKappaStats compute_weighted_kappa_pair_matrix(
    const Rcpp::IntegerMatrix& X,
    const int col_j,
    const int col_k,
    const Rcpp::NumericMatrix& weights,
    const bool return_inference,
    const double conf_level) {
  const int n_levels = weights.nrow();
  std::vector<int> counts(static_cast<std::size_t>(std::max(0, n_levels * n_levels)), 0);
  const int n = X.nrow();
  for (int i = 0; i < n; ++i) {
    const int xij = X(i, col_j);
    const int xik = X(i, col_k);
    if (xij == NA_INTEGER || xik == NA_INTEGER) {
      continue;
    }
    if (!valid_code(xij, n_levels) || !valid_code(xik, n_levels)) {
      continue;
    }
    counts[static_cast<std::size_t>((xij - 1) * n_levels + (xik - 1))] += 1;
  }
  return compute_weighted_kappa_from_counts(counts, weights, return_inference, conf_level);
}

WeightedKappaStats compute_weighted_kappa_diag(
    const Rcpp::IntegerMatrix& X,
    const int col_j,
    const Rcpp::NumericMatrix& weights,
    const bool return_inference) {
  WeightedKappaStats out;
  const int n_levels = weights.nrow();
  if (n_levels <= 0 || weights.ncol() != n_levels) {
    return out;
  }

  std::vector<double> q(static_cast<std::size_t>(n_levels), 0.0);
  const int n = X.nrow();
  for (int i = 0; i < n; ++i) {
    const int code = X(i, col_j);
    if (code == NA_INTEGER || !valid_code(code, n_levels)) {
      continue;
    }
    q[static_cast<std::size_t>(code - 1)] += 1.0;
    out.n_complete += 1;
  }

  if (out.n_complete < 2) {
    return out;
  }

  out.estimate = 1.0;
  out.observed_agreement = 1.0;

  const double n_complete = static_cast<double>(out.n_complete);
  for (int a = 0; a < n_levels; ++a) {
    q[static_cast<std::size_t>(a)] /= n_complete;
  }

  double expected = 0.0;
  for (int a = 0; a < n_levels; ++a) {
    for (int b = 0; b < n_levels; ++b) {
      expected += weights(a, b) * q[static_cast<std::size_t>(a)] * q[static_cast<std::size_t>(b)];
    }
  }
  out.expected_agreement = expected;

  if (return_inference) {
    out.lwr = 1.0;
    out.upr = 1.0;
  }
  return out;
}

inline void store_sym_double(Rcpp::NumericMatrix& M,
                             const int j,
                             const int k,
                             const double value) {
  M(j, k) = value;
  M(k, j) = value;
}

inline void store_sym_int(Rcpp::IntegerMatrix& M,
                          const int j,
                          const int k,
                          const int value) {
  M(j, k) = value;
  M(k, j) = value;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List weighted_kappa_pair_cpp(Rcpp::IntegerVector x,
                                   Rcpp::IntegerVector y,
                                   Rcpp::NumericMatrix weights,
                                   bool return_inference = false,
                                   double conf_level = 0.95) {
  WeightedKappaStats stats = compute_weighted_kappa_pair(
    x,
    y,
    weights,
    return_inference,
    conf_level
  );
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["estimate"] = stats.estimate,
    Rcpp::_["n_complete"] = stats.n_complete,
    Rcpp::_["observed_agreement"] = stats.observed_agreement,
    Rcpp::_["expected_agreement"] = stats.expected_agreement
  );
  if (return_inference) {
    out["se"] = stats.se;
    out["lwr"] = stats.lwr;
    out["upr"] = stats.upr;
    out["statistic"] = stats.statistic;
    out["p_value"] = stats.p_value;
    out["conf_level"] = conf_level;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List weighted_kappa_matrix_cpp(Rcpp::IntegerMatrix X,
                                     Rcpp::NumericMatrix weights,
                                     bool pairwise_complete = false,
                                     bool return_inference = false,
                                     double conf_level = 0.95,
                                     int n_threads = 1) {
  const int p = X.ncol();
  (void) pairwise_complete;

  Rcpp::NumericMatrix est(p, p);
  Rcpp::IntegerMatrix n_complete(p, p);
  Rcpp::NumericMatrix observed_agreement(p, p);
  Rcpp::NumericMatrix expected_agreement(p, p);
  std::fill(est.begin(), est.end(), NA_REAL);
  std::fill(n_complete.begin(), n_complete.end(), NA_INTEGER);
  std::fill(observed_agreement.begin(), observed_agreement.end(), NA_REAL);
  std::fill(expected_agreement.begin(), expected_agreement.end(), NA_REAL);

  Rcpp::NumericMatrix se;
  Rcpp::NumericMatrix lwr;
  Rcpp::NumericMatrix upr;
  Rcpp::NumericMatrix statistic;
  Rcpp::NumericMatrix p_value;
  if (return_inference) {
    se = Rcpp::NumericMatrix(p, p);
    lwr = Rcpp::NumericMatrix(p, p);
    upr = Rcpp::NumericMatrix(p, p);
    statistic = Rcpp::NumericMatrix(p, p);
    p_value = Rcpp::NumericMatrix(p, p);
    std::fill(se.begin(), se.end(), NA_REAL);
    std::fill(lwr.begin(), lwr.end(), NA_REAL);
    std::fill(upr.begin(), upr.end(), NA_REAL);
    std::fill(statistic.begin(), statistic.end(), NA_REAL);
    std::fill(p_value.begin(), p_value.end(), NA_REAL);
  }

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    WeightedKappaStats diag_stats = compute_weighted_kappa_diag(
      X,
      j,
      weights,
      return_inference
    );
    est(j, j) = diag_stats.estimate;
    n_complete(j, j) = diag_stats.n_complete;
    observed_agreement(j, j) = diag_stats.observed_agreement;
    expected_agreement(j, j) = diag_stats.expected_agreement;
    if (return_inference) {
      se(j, j) = diag_stats.se;
      lwr(j, j) = diag_stats.lwr;
      upr(j, j) = diag_stats.upr;
      statistic(j, j) = diag_stats.statistic;
      p_value(j, j) = diag_stats.p_value;
    }

    for (int k = j + 1; k < p; ++k) {
      WeightedKappaStats stats = compute_weighted_kappa_pair_matrix(
        X,
        j,
        k,
        weights,
        return_inference,
        conf_level
      );
      store_sym_double(est, j, k, stats.estimate);
      store_sym_int(n_complete, j, k, stats.n_complete);
      store_sym_double(observed_agreement, j, k, stats.observed_agreement);
      store_sym_double(expected_agreement, j, k, stats.expected_agreement);
      if (return_inference) {
        store_sym_double(se, j, k, stats.se);
        store_sym_double(lwr, j, k, stats.lwr);
        store_sym_double(upr, j, k, stats.upr);
        store_sym_double(statistic, j, k, stats.statistic);
        store_sym_double(p_value, j, k, stats.p_value);
      }
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["est"] = est,
    Rcpp::_["n_complete"] = n_complete,
    Rcpp::_["observed_agreement"] = observed_agreement,
    Rcpp::_["expected_agreement"] = expected_agreement
  );
  if (return_inference) {
    out["se"] = se;
    out["lwr"] = lwr;
    out["upr"] = upr;
    out["statistic"] = statistic;
    out["p_value"] = p_value;
    out["conf_level"] = conf_level;
    out["ci_method"] = "delta";
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List weighted_kappa_threshold_triplets_cpp(Rcpp::IntegerMatrix X,
                                                 Rcpp::NumericMatrix weights,
                                                 double threshold = 0.0,
                                                 bool diag = true,
                                                 int block_size = 256,
                                                 int n_threads = 1) {
  if (!(threshold >= 0.0) || !std::isfinite(threshold)) {
    Rcpp::stop("threshold must be finite and >= 0.");
  }
  if (block_size < 1) {
    Rcpp::stop("block_size must be >= 1.");
  }
  (void) block_size;

  const int p = X.ncol();
  const int n = X.nrow();
  if (p < 1 || n < 1) {
    return matrixCorr_detail::threshold_triplets::as_list(
      matrixCorr_detail::threshold_triplets::TripletBuffer()
    );
  }

  const int max_threads =
#ifdef _OPENMP
    std::max(1, n_threads);
#else
    1;
#endif

  std::vector<matrixCorr_detail::threshold_triplets::TripletBuffer> buffers(
    static_cast<std::size_t>(max_threads)
  );

#ifdef _OPENMP
  omp_set_num_threads(max_threads);
#pragma omp parallel
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    matrixCorr_detail::threshold_triplets::TripletBuffer& local =
      buffers[static_cast<std::size_t>(tid)];

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int j = 0; j < p; ++j) {
      if (diag) {
        WeightedKappaStats d = compute_weighted_kappa_diag(X, j, weights, false);
        if (std::isfinite(d.estimate) && std::fabs(d.estimate) >= threshold) {
          local.i.push_back(j + 1);
          local.j.push_back(j + 1);
          local.x.push_back(d.estimate);
        }
      }
      for (int k = j + 1; k < p; ++k) {
        WeightedKappaStats stats = compute_weighted_kappa_pair_matrix(
          X,
          j,
          k,
          weights,
          false,
          0.95
        );
        if (std::isfinite(stats.estimate) && std::fabs(stats.estimate) >= threshold) {
          local.i.push_back(j + 1);
          local.j.push_back(k + 1);
          local.x.push_back(stats.estimate);
        }
      }
    }
  }

  matrixCorr_detail::threshold_triplets::TripletBuffer out;
  std::size_t total = 0u;
  for (int t = 0; t < max_threads; ++t) {
    total += buffers[static_cast<std::size_t>(t)].i.size();
  }
  out.i.reserve(total);
  out.j.reserve(total);
  out.x.reserve(total);

  for (int t = 0; t < max_threads; ++t) {
    matrixCorr_detail::threshold_triplets::TripletBuffer& buf =
      buffers[static_cast<std::size_t>(t)];
    out.i.insert(out.i.end(), buf.i.begin(), buf.i.end());
    out.j.insert(out.j.end(), buf.j.begin(), buf.j.end());
    out.x.insert(out.x.end(), buf.x.begin(), buf.x.end());
  }

  return matrixCorr_detail::threshold_triplets::as_list(out);
}
