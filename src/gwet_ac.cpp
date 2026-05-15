// gwet.cpp
// Thiago de Paula Oliveira

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "matrixCorr_omp.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

namespace {

struct GwetPairStats {
  double est = NA_REAL;
  int n_complete = 0;
  double po = NA_REAL;
  double pe = NA_REAL;
  double se = NA_REAL;
  double lwr = NA_REAL;
  double upr = NA_REAL;
  double statistic = NA_REAL;
  double p_value = NA_REAL;
};

struct GwetCore {
  double estimate = NA_REAL;
  double observed_agreement = NA_REAL;
  double expected_agreement = NA_REAL;
  int n_items = 0;
  int n_categories = 0;
  int n_ratings_total = 0;
  int n_raters_min = 0;
  int n_raters_max = 0;
  int n_raters = 0;
  bool balanced = true;
  std::vector<int> category_counts;
  std::vector<double> category_proportions;
  std::vector<double> item_agreement;
  std::vector<int> n_raters_by_item;
  std::vector<double> pi_vec;
  arma::mat agree;
};

struct InferenceStats {
  double se = NA_REAL;
  double lwr = NA_REAL;
  double upr = NA_REAL;
  double statistic = NA_REAL;
  double p_value = NA_REAL;
};

inline bool valid_level(const int code, const int n_levels) {
  return code >= 1 && code <= n_levels;
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

inline double weight_sum(const Rcpp::NumericMatrix& weights) {
  return std::accumulate(weights.begin(), weights.end(), 0.0);
}

inline InferenceStats finalize_t_inference(const double estimate,
                                           const double se,
                                           const double conf_level,
                                           const double df) {
  InferenceStats out;
  if (!std::isfinite(estimate) || !std::isfinite(se) || se <= 0.0 ||
      !std::isfinite(df) || df <= 0.0) {
    return out;
  }

  out.se = se;
  out.statistic = estimate / se;
  out.p_value = 2.0 * R::pt(-std::fabs(out.statistic), df, 1, 0);
  const double tcrit = R::qt(1.0 - 0.5 * (1.0 - conf_level), df, 1, 0);
  out.lwr = std::max(-1.0, estimate - tcrit * se);
  out.upr = std::min(1.0, estimate + tcrit * se);
  return out;
}

GwetPairStats compute_gwet_pair_from_counts(const std::vector<int>& counts,
                                            const Rcpp::NumericMatrix& weights,
                                            const bool drop_unused_levels,
                                            const bool return_inference,
                                            const double conf_level) {
  GwetPairStats out;
  const int q = weights.nrow();
  out.n_complete = std::accumulate(counts.begin(), counts.end(), 0);
  if (q < 2 || weights.ncol() != q) {
    return out;
  }

  if (drop_unused_levels) {
    std::vector<double> row_sum(static_cast<std::size_t>(q), 0.0);
    std::vector<double> col_sum(static_cast<std::size_t>(q), 0.0);
    for (int a = 0; a < q; ++a) {
      for (int b = 0; b < q; ++b) {
        const double cab = static_cast<double>(counts[static_cast<std::size_t>(a * q + b)]);
        row_sum[static_cast<std::size_t>(a)] += cab;
        col_sum[static_cast<std::size_t>(b)] += cab;
      }
    }

    std::vector<int> active;
    active.reserve(static_cast<std::size_t>(q));
    for (int k = 0; k < q; ++k) {
      if ((row_sum[static_cast<std::size_t>(k)] + col_sum[static_cast<std::size_t>(k)]) > 0.0) {
        active.push_back(k);
      }
    }

    if (active.size() < 2u) {
      return out;
    }
    if (active.size() < static_cast<std::size_t>(q)) {
      const int q_active = static_cast<int>(active.size());
      std::vector<int> subcounts(static_cast<std::size_t>(q_active * q_active), 0);
      Rcpp::NumericMatrix subweights(q_active, q_active);
      for (int a = 0; a < q_active; ++a) {
        for (int b = 0; b < q_active; ++b) {
          const int ga = active[static_cast<std::size_t>(a)];
          const int gb = active[static_cast<std::size_t>(b)];
          subcounts[static_cast<std::size_t>(a * q_active + b)] =
            counts[static_cast<std::size_t>(ga * q + gb)];
          subweights(a, b) = weights(ga, gb);
        }
      }
      return compute_gwet_pair_from_counts(
        subcounts,
        subweights,
        false,
        return_inference,
        conf_level
      );
    }
  }

  long long n_complete = 0;
  double observed = 0.0;
  std::vector<double> row_prob(static_cast<std::size_t>(q), 0.0);
  std::vector<double> col_prob(static_cast<std::size_t>(q), 0.0);

  for (int a = 0; a < q; ++a) {
    for (int b = 0; b < q; ++b) {
      const std::size_t idx = static_cast<std::size_t>(a * q + b);
      const double cab = static_cast<double>(counts[idx]);
      n_complete += static_cast<long long>(counts[idx]);
      row_prob[static_cast<std::size_t>(a)] += cab;
      col_prob[static_cast<std::size_t>(b)] += cab;
      observed += weights(a, b) * cab;
    }
  }

  out.n_complete = static_cast<int>(n_complete);
  if (n_complete < 2) {
    return out;
  }

  const double n = static_cast<double>(n_complete);
  observed /= n;
  for (int k = 0; k < q; ++k) {
    row_prob[static_cast<std::size_t>(k)] /= n;
    col_prob[static_cast<std::size_t>(k)] /= n;
  }

  std::vector<double> pi_vec(static_cast<std::size_t>(q), 0.0);
  for (int k = 0; k < q; ++k) {
    pi_vec[static_cast<std::size_t>(k)] =
      0.5 * (row_prob[static_cast<std::size_t>(k)] + col_prob[static_cast<std::size_t>(k)]);
  }

  const double tw = weight_sum(weights);
  double expected = 0.0;
  for (int k = 0; k < q; ++k) {
    const double pik = pi_vec[static_cast<std::size_t>(k)];
    expected += pik * (1.0 - pik);
  }
  expected = tw * expected / static_cast<double>(q * (q - 1));

  out.po = observed;
  out.pe = expected;
  const double denom = 1.0 - expected;
  if (!std::isfinite(denom) || denom <= 1e-15) {
    return out;
  }

  const double estimate = (observed - expected) / denom;
  if (!std::isfinite(estimate)) {
    return out;
  }
  out.est = estimate;

  if (!return_inference) {
    return out;
  }

  double sum1 = 0.0;
  for (int a = 0; a < q; ++a) {
    for (int b = 0; b < q; ++b) {
      const double pab = static_cast<double>(counts[static_cast<std::size_t>(a * q + b)]) / n;
      const double term = weights(a, b) -
        2.0 * (1.0 - estimate) * tw *
        (1.0 - (pi_vec[static_cast<std::size_t>(a)] + pi_vec[static_cast<std::size_t>(b)]) / 2.0) /
        static_cast<double>(q * (q - 1));
      sum1 += pab * term * term;
    }
  }

  double variance = (sum1 - std::pow(observed - 2.0 * (1.0 - estimate) * expected, 2.0)) /
    (n * denom * denom);
  variance = std::max(variance, 1e-100);
  const double se = std::sqrt(variance);
  const InferenceStats inf = finalize_t_inference(
    estimate,
    se,
    conf_level,
    static_cast<double>(n_complete - 1)
  );
  out.se = inf.se;
  out.lwr = inf.lwr;
  out.upr = inf.upr;
  out.statistic = inf.statistic;
  out.p_value = inf.p_value;
  return out;
}

GwetPairStats compute_gwet_pair(const Rcpp::IntegerVector& x,
                                const Rcpp::IntegerVector& y,
                                const Rcpp::NumericMatrix& weights,
                                const bool drop_unused_levels,
                                const bool return_inference,
                                const double conf_level) {
  const int q = weights.nrow();
  std::vector<int> counts(static_cast<std::size_t>(std::max(0, q * q)), 0);
  const R_xlen_t n = std::min(x.size(), y.size());
  for (R_xlen_t i = 0; i < n; ++i) {
    const int xi = x[i];
    const int yi = y[i];
    if (xi == NA_INTEGER || yi == NA_INTEGER) {
      continue;
    }
    if (!valid_level(xi, q) || !valid_level(yi, q)) {
      continue;
    }
    counts[static_cast<std::size_t>((xi - 1) * q + (yi - 1))] += 1;
  }
  return compute_gwet_pair_from_counts(
    counts,
    weights,
    drop_unused_levels,
    return_inference,
    conf_level
  );
}

GwetPairStats compute_gwet_pair_matrix(const Rcpp::IntegerMatrix& X,
                                       const int col_j,
                                       const int col_k,
                                       const Rcpp::NumericMatrix& weights,
                                       const bool drop_unused_levels,
                                       const bool return_inference,
                                       const double conf_level) {
  const int q = weights.nrow();
  std::vector<int> counts(static_cast<std::size_t>(std::max(0, q * q)), 0);
  const int n = X.nrow();
  for (int i = 0; i < n; ++i) {
    const int xij = X(i, col_j);
    const int xik = X(i, col_k);
    if (xij == NA_INTEGER || xik == NA_INTEGER) {
      continue;
    }
    if (!valid_level(xij, q) || !valid_level(xik, q)) {
      continue;
    }
    counts[static_cast<std::size_t>((xij - 1) * q + (xik - 1))] += 1;
  }
  return compute_gwet_pair_from_counts(
    counts,
    weights,
    drop_unused_levels,
    return_inference,
    conf_level
  );
}

GwetPairStats compute_gwet_diag(const Rcpp::IntegerMatrix& X,
                                const int col_j,
                                const Rcpp::NumericMatrix& weights,
                                const bool return_inference,
                                const double conf_level) {
  GwetPairStats out;
  const int q = weights.nrow();
  if (q < 2 || weights.ncol() != q) {
    return out;
  }

  std::vector<int> freq(static_cast<std::size_t>(q), 0);
  const int n = X.nrow();
  for (int i = 0; i < n; ++i) {
    const int code = X(i, col_j);
    if (code == NA_INTEGER || !valid_level(code, q)) {
      continue;
    }
    freq[static_cast<std::size_t>(code - 1)] += 1;
    out.n_complete += 1;
  }

  if (out.n_complete < 2) {
    return out;
  }

  out.est = 1.0;
  out.po = 1.0;
  const double n_complete = static_cast<double>(out.n_complete);
  const double tw = weight_sum(weights);
  double expected = 0.0;
  for (int k = 0; k < q; ++k) {
    const double pk = static_cast<double>(freq[static_cast<std::size_t>(k)]) / n_complete;
    expected += pk * (1.0 - pk);
  }
  out.pe = tw * expected / static_cast<double>(q * (q - 1));

  if (return_inference) {
    out.lwr = 1.0;
    out.upr = 1.0;
    out.se = NA_REAL;
    out.statistic = NA_REAL;
    out.p_value = NA_REAL;
    if (!std::isfinite(conf_level)) {
      out.lwr = NA_REAL;
      out.upr = NA_REAL;
    }
  }
  return out;
}

GwetCore compute_gwet_from_counts(const Rcpp::IntegerMatrix& counts,
                                  const Rcpp::NumericMatrix& weights) {
  GwetCore out;
  const int n_items = counts.nrow();
  const int q = counts.ncol();
  out.n_items = n_items;
  out.n_categories = q;
  out.category_counts.assign(static_cast<std::size_t>(std::max(0, q)), 0);
  out.category_proportions.assign(static_cast<std::size_t>(std::max(0, q)), NA_REAL);
  out.item_agreement.assign(static_cast<std::size_t>(std::max(0, n_items)), NA_REAL);
  out.n_raters_by_item.assign(static_cast<std::size_t>(std::max(0, n_items)), 0);
  out.pi_vec.assign(static_cast<std::size_t>(std::max(0, q)), NA_REAL);
  out.agree = arma::mat(std::max(0, n_items), std::max(0, q), arma::fill::zeros);

  if (n_items <= 0 || q < 2 || weights.nrow() != q || weights.ncol() != q) {
    return out;
  }

  const double tw = weight_sum(weights);
  bool first_row = true;
  double observed_sum = 0.0;
  std::vector<double> pi_sum(static_cast<std::size_t>(q), 0.0);
  int n2more = 0;

  arma::mat W = Rcpp::as<arma::mat>(weights);

  for (int i = 0; i < n_items; ++i) {
    double row_sum = 0.0;
    for (int k = 0; k < q; ++k) {
      const int rik = counts(i, k);
      out.agree(i, k) = static_cast<double>(rik);
      out.category_counts[static_cast<std::size_t>(k)] += rik;
      row_sum += rik;
    }

    const int ri = static_cast<int>(row_sum);
    out.n_raters_by_item[static_cast<std::size_t>(i)] = ri;
    out.n_ratings_total += ri;
    if (first_row) {
      out.n_raters_min = ri;
      out.n_raters_max = ri;
      first_row = false;
    } else {
      out.n_raters_min = std::min(out.n_raters_min, ri);
      out.n_raters_max = std::max(out.n_raters_max, ri);
    }

    if (ri >= 2) {
      arma::rowvec agree_i = out.agree.row(i);
      arma::rowvec agree_w = agree_i * W;
      const double sum_q = arma::accu(agree_i % (agree_w - 1.0));
      const double item_agreement = sum_q / (static_cast<double>(ri) * static_cast<double>(ri - 1));
      out.item_agreement[static_cast<std::size_t>(i)] = item_agreement;
      observed_sum += item_agreement;
      n2more += 1;
    }

    if (ri > 0) {
      for (int k = 0; k < q; ++k) {
        pi_sum[static_cast<std::size_t>(k)] += out.agree(i, k) / static_cast<double>(ri);
      }
    }
  }

  out.balanced = (out.n_raters_min == out.n_raters_max);
  if (out.balanced) {
    out.n_raters = out.n_raters_max;
  }
  if (n2more < 1) {
    return out;
  }

  out.observed_agreement = observed_sum / static_cast<double>(n2more);
  for (int k = 0; k < q; ++k) {
    out.category_proportions[static_cast<std::size_t>(k)] =
      static_cast<double>(out.category_counts[static_cast<std::size_t>(k)]) /
      static_cast<double>(std::max(1, out.n_ratings_total));
    out.pi_vec[static_cast<std::size_t>(k)] = pi_sum[static_cast<std::size_t>(k)] / static_cast<double>(n_items);
  }

  double expected = 0.0;
  for (int k = 0; k < q; ++k) {
    const double pik = out.pi_vec[static_cast<std::size_t>(k)];
    expected += pik * (1.0 - pik);
  }
  out.expected_agreement = tw * expected / static_cast<double>(q * (q - 1));

  const double denom = 1.0 - out.expected_agreement;
  if (!std::isfinite(denom) || denom <= 1e-15) {
    return out;
  }
  out.estimate = (out.observed_agreement - out.expected_agreement) / denom;
  if (!std::isfinite(out.estimate)) {
    out.estimate = NA_REAL;
  }
  return out;
}

InferenceStats compute_asymptotic_se_counts(const GwetCore& full,
                                            const Rcpp::NumericMatrix& weights,
                                            const double conf_level) {
  InferenceStats out;
  const int n = full.n_items;
  const int q = full.n_categories;
  if (n < 2 || q < 2 || !std::isfinite(full.estimate)) {
    return out;
  }

  const double tw = weight_sum(weights);
  const double denom = 1.0 - full.expected_agreement;
  arma::vec one_minus_pi(q);
  for (int k = 0; k < q; ++k) {
    one_minus_pi[k] = 1.0 - full.pi_vec[static_cast<std::size_t>(k)];
  }

  arma::vec ri_vec(n);
  arma::vec acx(n);
  for (int i = 0; i < n; ++i) {
    const double ri = static_cast<double>(full.n_raters_by_item[static_cast<std::size_t>(i)]);
    ri_vec[i] = ri;
    const double den_i = ri * (ri - 1.0);
    const double den_safe = (den_i == 0.0) ? -1.0 : den_i;
    arma::rowvec agree_i = full.agree.row(i);
    arma::rowvec agree_w = agree_i * Rcpp::as<arma::mat>(weights);
    const double sum_q = arma::accu(agree_i % (agree_w - 1.0));
    const double pa_i = sum_q / den_safe;
    const double ac_i = (pa_i - full.expected_agreement) / denom;
    const double pe_i = (tw / static_cast<double>(q * (q - 1))) *
      arma::dot(agree_i.t(), one_minus_pi) / ri;
    acx[i] = ac_i - 2.0 * (1.0 - full.estimate) * (pe_i - full.expected_agreement) / denom;
  }

  const double variance = arma::accu(arma::square(acx - full.estimate)) /
    (static_cast<double>(n) * static_cast<double>(n - 1));
  if (!std::isfinite(variance)) {
    return out;
  }

  return finalize_t_inference(
    full.estimate,
    std::sqrt(std::max(variance, 1e-100)),
    conf_level,
    static_cast<double>(n - 1)
  );
}

InferenceStats compute_jackknife_se_counts(const Rcpp::IntegerMatrix& counts,
                                           const Rcpp::NumericMatrix& weights,
                                           const double conf_level,
                                           const int n_threads) {
  InferenceStats out;
  const int n_items = counts.nrow();
  if (n_items < 3) {
    return out;
  }

  const GwetCore full = compute_gwet_from_counts(counts, weights);
  if (!std::isfinite(full.estimate)) {
    return out;
  }

  std::vector<double> theta_minus(static_cast<std::size_t>(n_items), NA_REAL);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) if (n_items >= 32)
#endif
  for (int i = 0; i < n_items; ++i) {
    Rcpp::IntegerMatrix reduced(n_items - 1, counts.ncol());
    int rr = 0;
    for (int s = 0; s < n_items; ++s) {
      if (s == i) {
        continue;
      }
      for (int k = 0; k < counts.ncol(); ++k) {
        reduced(rr, k) = counts(s, k);
      }
      rr += 1;
    }
    const GwetCore minus = compute_gwet_from_counts(reduced, weights);
    theta_minus[static_cast<std::size_t>(i)] = minus.estimate;
  }

  for (int i = 0; i < n_items; ++i) {
    if (!std::isfinite(theta_minus[static_cast<std::size_t>(i)])) {
      return out;
    }
  }

  const double theta_bar = std::accumulate(theta_minus.begin(), theta_minus.end(), 0.0) /
    static_cast<double>(n_items);
  double sum_sq = 0.0;
  for (int i = 0; i < n_items; ++i) {
    const double diff = theta_minus[static_cast<std::size_t>(i)] - theta_bar;
    sum_sq += diff * diff;
  }
  const double variance =
    (static_cast<double>(n_items - 1) / static_cast<double>(n_items)) * sum_sq;
  if (!std::isfinite(variance)) {
    return out;
  }

  return finalize_t_inference(
    full.estimate,
    std::sqrt(std::max(variance, 0.0)),
    conf_level,
    static_cast<double>(n_items - 1)
  );
}

Rcpp::IntegerMatrix collapse_counts_for_category(const Rcpp::IntegerMatrix& counts,
                                                 const int category_index) {
  const int n_items = counts.nrow();
  const int q = counts.ncol();
  Rcpp::IntegerMatrix out(n_items, 2);
  for (int i = 0; i < n_items; ++i) {
    int total = 0;
    for (int k = 0; k < q; ++k) {
      total += counts(i, k);
    }
    out(i, 0) = counts(i, category_index);
    out(i, 1) = total - counts(i, category_index);
  }
  return out;
}

Rcpp::IntegerMatrix ratings_to_counts(const Rcpp::IntegerMatrix& ratings,
                                      const int n_levels,
                                      const int na_code,
                                      const int min_raters) {
  const int n_items = ratings.nrow();
  const int n_raters = ratings.ncol();
  std::vector<int> keep_rows;
  keep_rows.reserve(static_cast<std::size_t>(std::max(0, n_items)));

  for (int i = 0; i < n_items; ++i) {
    int observed = 0;
    for (int j = 0; j < n_raters; ++j) {
      const int code = ratings(i, j);
      if (code != NA_INTEGER && code >= 1 && code <= n_levels) {
        observed += 1;
      }
    }
    if (na_code == 3) {
      if (observed >= min_raters) {
        keep_rows.push_back(i);
      }
    } else {
      keep_rows.push_back(i);
    }
  }

  Rcpp::IntegerMatrix counts(static_cast<int>(keep_rows.size()), n_levels);
  for (int out_i = 0; out_i < static_cast<int>(keep_rows.size()); ++out_i) {
    const int src_i = keep_rows[static_cast<std::size_t>(out_i)];
    for (int j = 0; j < n_raters; ++j) {
      const int code = ratings(src_i, j);
      if (code != NA_INTEGER && code >= 1 && code <= n_levels) {
        counts(out_i, code - 1) += 1;
      }
    }
  }
  return counts;
}

Rcpp::List build_counts_result(const Rcpp::IntegerMatrix& counts,
                               const Rcpp::NumericMatrix& weights,
                               const bool by_category,
                               const int inference_code,
                               const double conf_level,
                               const int n_threads) {
  const GwetCore full = compute_gwet_from_counts(counts, weights);
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("estimate") = full.estimate,
    Rcpp::Named("observed_agreement") = full.observed_agreement,
    Rcpp::Named("expected_agreement") = full.expected_agreement,
    Rcpp::Named("n_items") = full.n_items,
    Rcpp::Named("n_categories") = full.n_categories,
    Rcpp::Named("n_ratings_total") = full.n_ratings_total,
    Rcpp::Named("n_raters_min") = full.n_raters_min,
    Rcpp::Named("n_raters_max") = full.n_raters_max,
    Rcpp::Named("n_raters") = full.n_raters,
    Rcpp::Named("balanced") = full.balanced,
    Rcpp::Named("category_counts") = Rcpp::wrap(full.category_counts),
    Rcpp::Named("category_proportions") = Rcpp::wrap(full.category_proportions),
    Rcpp::Named("item_category_counts") = counts,
    Rcpp::Named("item_agreement") = Rcpp::wrap(full.item_agreement),
    Rcpp::Named("n_raters_by_item") = Rcpp::wrap(full.n_raters_by_item)
  );

  if (inference_code != 0) {
    const InferenceStats inf = (inference_code == 1)
      ? compute_asymptotic_se_counts(full, weights, conf_level)
      : compute_jackknife_se_counts(counts, weights, conf_level, n_threads);
    out["se"] = inf.se;
    out["lwr"] = inf.lwr;
    out["upr"] = inf.upr;
    out["statistic"] = inf.statistic;
    out["p_value"] = inf.p_value;
    out["conf_level"] = conf_level;
    out["se_method"] = (inference_code == 1) ? "asymptotic" : "jackknife";
  }

  if (by_category) {
    const int q = counts.ncol();
    Rcpp::NumericVector by_estimate(q, NA_REAL);
    Rcpp::NumericVector by_observed(q, NA_REAL);
    Rcpp::NumericVector by_expected(q, NA_REAL);
    Rcpp::NumericVector by_se(q, NA_REAL);
    Rcpp::NumericVector by_lwr(q, NA_REAL);
    Rcpp::NumericVector by_upr(q, NA_REAL);
    Rcpp::NumericVector by_statistic(q, NA_REAL);
    Rcpp::NumericVector by_p_value(q, NA_REAL);

    Rcpp::NumericMatrix binary_weights(2, 2);
    binary_weights(0, 0) = 1.0;
    binary_weights(1, 1) = 1.0;
    binary_weights(0, 1) = 0.0;
    binary_weights(1, 0) = 0.0;

    for (int k = 0; k < q; ++k) {
      const Rcpp::IntegerMatrix collapsed = collapse_counts_for_category(counts, k);
      const GwetCore cat_full = compute_gwet_from_counts(collapsed, binary_weights);
      by_estimate[k] = cat_full.estimate;
      by_observed[k] = cat_full.observed_agreement;
      by_expected[k] = cat_full.expected_agreement;
      if (inference_code != 0) {
        const InferenceStats cat_inf = (inference_code == 1)
          ? compute_asymptotic_se_counts(cat_full, binary_weights, conf_level)
          : compute_jackknife_se_counts(collapsed, binary_weights, conf_level, n_threads);
        by_se[k] = cat_inf.se;
        by_lwr[k] = cat_inf.lwr;
        by_upr[k] = cat_inf.upr;
        by_statistic[k] = cat_inf.statistic;
        by_p_value[k] = cat_inf.p_value;
      }
    }

    out["by_category_estimate"] = by_estimate;
    out["by_category_observed_agreement"] = by_observed;
    out["by_category_expected_agreement"] = by_expected;
    out["by_category_se"] = by_se;
    out["by_category_lwr"] = by_lwr;
    out["by_category_upr"] = by_upr;
    out["by_category_statistic"] = by_statistic;
    out["by_category_p_value"] = by_p_value;
  }

  return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List gwet_ac_pair_cpp(Rcpp::IntegerVector x,
                            Rcpp::IntegerVector y,
                            Rcpp::NumericMatrix weights,
                            bool drop_unused_levels = false,
                            bool return_inference = false,
                            double conf_level = 0.95) {
  GwetPairStats stats = compute_gwet_pair(
    x,
    y,
    weights,
    drop_unused_levels,
    return_inference,
    conf_level
  );
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["est"] = stats.est,
    Rcpp::_["n_complete"] = stats.n_complete,
    Rcpp::_["po"] = stats.po,
    Rcpp::_["pe"] = stats.pe
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
Rcpp::List gwet_ac_matrix_cpp(Rcpp::IntegerMatrix X,
                              Rcpp::NumericMatrix weights,
                              bool drop_unused_levels = false,
                              bool pairwise_complete = false,
                              bool return_inference = false,
                              double conf_level = 0.95,
                              int n_threads = 1) {
  const int p = X.ncol();
  Rcpp::NumericMatrix est(p, p);
  Rcpp::IntegerMatrix n_complete(p, p);
  Rcpp::NumericMatrix po(p, p);
  Rcpp::NumericMatrix pe(p, p);
  std::fill(est.begin(), est.end(), NA_REAL);
  std::fill(n_complete.begin(), n_complete.end(), NA_INTEGER);
  std::fill(po.begin(), po.end(), NA_REAL);
  std::fill(pe.begin(), pe.end(), NA_REAL);

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
    GwetPairStats diag_stats = compute_gwet_diag(X, j, weights, return_inference, conf_level);
    est(j, j) = diag_stats.est;
    n_complete(j, j) = diag_stats.n_complete;
    po(j, j) = diag_stats.po;
    pe(j, j) = diag_stats.pe;
    if (return_inference) {
      se(j, j) = diag_stats.se;
      lwr(j, j) = diag_stats.lwr;
      upr(j, j) = diag_stats.upr;
      statistic(j, j) = diag_stats.statistic;
      p_value(j, j) = diag_stats.p_value;
    }

    for (int k = j + 1; k < p; ++k) {
      GwetPairStats stats = compute_gwet_pair_matrix(
        X,
        j,
        k,
        weights,
        drop_unused_levels,
        return_inference,
        conf_level
      );
      store_sym_double(est, j, k, stats.est);
      store_sym_int(n_complete, j, k, stats.n_complete);
      store_sym_double(po, j, k, stats.po);
      store_sym_double(pe, j, k, stats.pe);
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
    Rcpp::_["po"] = po,
    Rcpp::_["pe"] = pe
  );
  if (return_inference) {
    out["se"] = se;
    out["lwr"] = lwr;
    out["upr"] = upr;
    out["statistic"] = statistic;
    out["p_value"] = p_value;
    out["conf_level"] = conf_level;
  }
  out["pairwise_complete"] = pairwise_complete;
  return out;
}

// [[Rcpp::export]]
Rcpp::List gwet_ac_threshold_triplets_cpp(Rcpp::IntegerMatrix X,
                                          Rcpp::NumericMatrix weights,
                                          bool drop_unused_levels = false,
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

  const int p = X.ncol();
  const int n = X.nrow();
  if (p < 1 || n < 1) {
    return Rcpp::List::create(
      Rcpp::_["i"] = Rcpp::IntegerVector(),
      Rcpp::_["j"] = Rcpp::IntegerVector(),
      Rcpp::_["x"] = Rcpp::NumericVector()
    );
  }

  const int max_threads =
#ifdef _OPENMP
    std::max(1, n_threads);
#else
    1;
#endif
  std::vector< std::vector<int> > i_buf(static_cast<std::size_t>(max_threads));
  std::vector< std::vector<int> > j_buf(static_cast<std::size_t>(max_threads));
  std::vector< std::vector<double> > x_buf(static_cast<std::size_t>(max_threads));

#ifdef _OPENMP
  omp_set_num_threads(max_threads);
#pragma omp parallel
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    std::vector<int>& i_local = i_buf[static_cast<std::size_t>(tid)];
    std::vector<int>& j_local = j_buf[static_cast<std::size_t>(tid)];
    std::vector<double>& x_local = x_buf[static_cast<std::size_t>(tid)];

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int j = 0; j < p; ++j) {
      if (diag) {
        GwetPairStats d = compute_gwet_diag(X, j, weights, false, 0.95);
        if (std::isfinite(d.est) && std::fabs(d.est) >= threshold) {
          i_local.push_back(j + 1);
          j_local.push_back(j + 1);
          x_local.push_back(d.est);
        }
      }
      for (int k = j + 1; k < p; ++k) {
        GwetPairStats stats = compute_gwet_pair_matrix(
          X,
          j,
          k,
          weights,
          drop_unused_levels,
          false,
          0.95
        );
        if (std::isfinite(stats.est) && std::fabs(stats.est) >= threshold) {
          i_local.push_back(j + 1);
          j_local.push_back(k + 1);
          x_local.push_back(stats.est);
        }
      }
    }
  }

  std::size_t total = 0u;
  for (int t = 0; t < max_threads; ++t) {
    total += i_buf[static_cast<std::size_t>(t)].size();
  }

  std::vector<int> i_out;
  std::vector<int> j_out;
  std::vector<double> x_out;
  i_out.reserve(total);
  j_out.reserve(total);
  x_out.reserve(total);
  for (int t = 0; t < max_threads; ++t) {
    i_out.insert(i_out.end(), i_buf[static_cast<std::size_t>(t)].begin(), i_buf[static_cast<std::size_t>(t)].end());
    j_out.insert(j_out.end(), j_buf[static_cast<std::size_t>(t)].begin(), j_buf[static_cast<std::size_t>(t)].end());
    x_out.insert(x_out.end(), x_buf[static_cast<std::size_t>(t)].begin(), x_buf[static_cast<std::size_t>(t)].end());
  }

  return Rcpp::List::create(
    Rcpp::_["i"] = Rcpp::wrap(i_out),
    Rcpp::_["j"] = Rcpp::wrap(j_out),
    Rcpp::_["x"] = Rcpp::wrap(x_out)
  );
}

// [[Rcpp::export]]
Rcpp::List gwet_ac_counts_cpp(Rcpp::IntegerMatrix counts,
                              Rcpp::NumericMatrix weights,
                              bool by_category = false,
                              int inference_code = 0,
                              double conf_level = 0.95,
                              int n_threads = 1) {
  return build_counts_result(
    counts,
    weights,
    by_category,
    inference_code,
    conf_level,
    n_threads
  );
}

// [[Rcpp::export]]
Rcpp::List gwet_ac_ratings_cpp(Rcpp::IntegerMatrix ratings,
                               int n_levels,
                               Rcpp::NumericMatrix weights,
                               int na_code = 1,
                               int min_raters = 2,
                               bool by_category = false,
                               int inference_code = 0,
                               double conf_level = 0.95,
                               int n_threads = 1) {
  const Rcpp::IntegerMatrix counts = ratings_to_counts(
    ratings,
    n_levels,
    na_code,
    min_raters
  );
  return build_counts_result(
    counts,
    weights,
    by_category,
    inference_code,
    conf_level,
    n_threads
  );
}
