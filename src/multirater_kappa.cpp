#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "matrixCorr_omp.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

namespace {

struct KappaScalar {
  double estimate = NA_REAL;
  double observed_agreement = NA_REAL;
  double expected_agreement = NA_REAL;
};

struct JackknifeStats {
  double se = NA_REAL;
  double lwr = NA_REAL;
  double upr = NA_REAL;
  double statistic = NA_REAL;
  double p_value = NA_REAL;
};

inline JackknifeStats finalize_inference_from_se(const double estimate,
                                                 const double se,
                                                 const double conf_level,
                                                 const bool irr_style_p_value = false) {
  JackknifeStats out;
  if (!std::isfinite(se) || se <= 0.0 || !std::isfinite(estimate)) {
    return out;
  }

  out.se = se;
  out.statistic = estimate / se;
  if (irr_style_p_value) {
    out.p_value = 2.0 * (1.0 - R::pnorm5(std::fabs(out.statistic), 0.0, 1.0, 1, 0));
  } else {
    out.p_value = 2.0 * R::pnorm5(std::fabs(out.statistic), 0.0, 1.0, 0, 0);
  }
  const double z = R::qnorm5(0.5 * (1.0 + conf_level), 0.0, 1.0, 1, 0);
  out.lwr = std::max(-1.0, estimate - z * se);
  out.upr = std::min(1.0, estimate + z * se);
  return out;
}

struct MultiKappaCore {
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
  bool exact = false;
  std::vector<int> category_counts;
  std::vector<int> rater_category_counts;
  std::vector<double> category_proportions;
  std::vector<double> item_agreement;
  std::vector<int> n_raters_by_item;
};

inline double compute_fleiss_expected(const std::vector<int>& category_counts,
                                      const int total_ratings,
                                      const int n_categories) {
  if (n_categories < 2 || total_ratings <= 0) {
    return NA_REAL;
  }
  const double denom = static_cast<double>(total_ratings);
  double expected = 0.0;
  for (int j = 0; j < n_categories; ++j) {
    const double pj = static_cast<double>(category_counts[static_cast<std::size_t>(j)]) / denom;
    expected += pj * pj;
  }
  return expected;
}

inline double compute_randolph_expected(const int n_categories) {
  if (n_categories < 2) {
    return NA_REAL;
  }
  return 1.0 / static_cast<double>(n_categories);
}

double compute_fleiss_exact_expected(const std::vector<int>& category_counts,
                                     const std::vector<int>& rater_category_counts,
                                     const int n_items,
                                     const int n_raters,
                                     const int n_categories) {
  if (n_items < 2 || n_raters < 2 || n_categories < 2) {
    return NA_REAL;
  }

  const double ns = static_cast<double>(n_items);
  const double nr = static_cast<double>(n_raters);
  const double total = ns * nr;

  double chance = 0.0;
  for (int j = 0; j < n_categories; ++j) {
    const double col_sum = static_cast<double>(category_counts[static_cast<std::size_t>(j)]);
    chance += (col_sum * col_sum) / (total * total);
  }

  double variance_term = 0.0;
  for (int j = 0; j < n_categories; ++j) {
    double mean_prop = 0.0;
    for (int r = 0; r < n_raters; ++r) {
      mean_prop += static_cast<double>(
        rater_category_counts[static_cast<std::size_t>(r * n_categories + j)]
      ) / ns;
    }
    mean_prop /= nr;

    double ss = 0.0;
    for (int r = 0; r < n_raters; ++r) {
      const double prop = static_cast<double>(
        rater_category_counts[static_cast<std::size_t>(r * n_categories + j)]
      ) / ns;
      const double diff = prop - mean_prop;
      ss += diff * diff;
    }
    const double sample_var = ss / static_cast<double>(n_raters - 1);
    variance_term += sample_var / nr;
  }

  return chance - variance_term;
}

inline KappaScalar finalize_kappa_scalar(const double observed,
                                         const double expected) {
  KappaScalar out;
  out.observed_agreement = observed;
  out.expected_agreement = expected;

  const double denom = 1.0 - expected;
  if (!std::isfinite(observed) || !std::isfinite(expected) ||
      !std::isfinite(denom) || denom <= 1e-15) {
    return out;
  }

  const double estimate = (observed - expected) / denom;
  if (!std::isfinite(estimate)) {
    return out;
  }
  out.estimate = estimate;
  return out;
}

inline KappaScalar compute_multirater_kappa_scalar(const double observed,
                                                   const std::vector<int>& category_counts,
                                                   const int total_ratings,
                                                   const int n_categories,
                                                   const int method_code) {
  const double expected =
    (method_code == 2)
      ? compute_randolph_expected(n_categories)
      : compute_fleiss_expected(category_counts, total_ratings, n_categories);
  return finalize_kappa_scalar(observed, expected);
}

inline KappaScalar compute_multirater_kappa_scalar_exact(const double observed,
                                                         const std::vector<int>& category_counts,
                                                         const std::vector<int>& rater_category_counts,
                                                         const int n_items,
                                                         const int n_raters,
                                                         const int n_categories) {
  const double expected = compute_fleiss_exact_expected(
    category_counts,
    rater_category_counts,
    n_items,
    n_raters,
    n_categories
  );
  return finalize_kappa_scalar(observed, expected);
}

MultiKappaCore compute_multirater_kappa_from_counts(const Rcpp::IntegerMatrix& counts,
                                                    const int method_code) {
  MultiKappaCore out;
  const int n_items = counts.nrow();
  const int n_categories = counts.ncol();
  out.n_items = n_items;
  out.n_categories = n_categories;
  out.category_counts.assign(static_cast<std::size_t>(std::max(0, n_categories)), 0);
  out.category_proportions.assign(static_cast<std::size_t>(std::max(0, n_categories)), NA_REAL);
  out.item_agreement.assign(static_cast<std::size_t>(std::max(0, n_items)), NA_REAL);
  out.n_raters_by_item.assign(static_cast<std::size_t>(std::max(0, n_items)), 0);

  if (n_items <= 0 || n_categories <= 0) {
    return out;
  }

  double sum_item_agreement = 0.0;
  bool first_row = true;

  for (int i = 0; i < n_items; ++i) {
    int row_sum = 0;
    double numer = 0.0;
    for (int j = 0; j < n_categories; ++j) {
      const int nij = counts(i, j);
      row_sum += nij;
      out.category_counts[static_cast<std::size_t>(j)] += nij;
      numer += static_cast<double>(nij) * static_cast<double>(nij - 1);
    }

    out.n_raters_by_item[static_cast<std::size_t>(i)] = row_sum;
    out.n_ratings_total += row_sum;
    if (first_row) {
      out.n_raters_min = row_sum;
      out.n_raters_max = row_sum;
      first_row = false;
    } else {
      out.n_raters_min = std::min(out.n_raters_min, row_sum);
      out.n_raters_max = std::max(out.n_raters_max, row_sum);
    }

    if (row_sum >= 2) {
      const double denom = static_cast<double>(row_sum) * static_cast<double>(row_sum - 1);
      const double item_agreement = numer / denom;
      out.item_agreement[static_cast<std::size_t>(i)] = item_agreement;
      sum_item_agreement += item_agreement;
    }
  }

  out.balanced = (out.n_raters_min == out.n_raters_max);
  if (out.balanced) {
    out.n_raters = out.n_raters_max;
  }
  if (out.n_items < 2 || out.n_categories < 2 || out.n_ratings_total <= 0) {
    return out;
  }

  for (int j = 0; j < n_categories; ++j) {
    out.category_proportions[static_cast<std::size_t>(j)] =
      static_cast<double>(out.category_counts[static_cast<std::size_t>(j)]) /
      static_cast<double>(out.n_ratings_total);
  }

  const double observed = sum_item_agreement / static_cast<double>(out.n_items);
  const KappaScalar scalar = compute_multirater_kappa_scalar(
    observed,
    out.category_counts,
    out.n_ratings_total,
    out.n_categories,
    method_code
  );
  out.estimate = scalar.estimate;
  out.observed_agreement = scalar.observed_agreement;
  out.expected_agreement = scalar.expected_agreement;
  return out;
}

MultiKappaCore compute_multirater_kappa_from_ratings_complete(const Rcpp::IntegerMatrix& ratings,
                                                              const int n_levels,
                                                              const bool exact) {
  MultiKappaCore out;
  const int n_items = ratings.nrow();
  const int n_raters = ratings.ncol();
  out.n_items = n_items;
  out.n_categories = n_levels;
  out.n_ratings_total = n_items * n_raters;
  out.n_raters = n_raters;
  out.n_raters_min = n_raters;
  out.n_raters_max = n_raters;
  out.balanced = true;
  out.exact = exact;
  out.category_counts.assign(static_cast<std::size_t>(std::max(0, n_levels)), 0);
  out.category_proportions.assign(static_cast<std::size_t>(std::max(0, n_levels)), NA_REAL);
  out.item_agreement.assign(static_cast<std::size_t>(std::max(0, n_items)), NA_REAL);
  out.n_raters_by_item.assign(static_cast<std::size_t>(std::max(0, n_items)), n_raters);
  if (exact) {
    out.rater_category_counts.assign(
      static_cast<std::size_t>(std::max(0, n_raters * n_levels)),
      0
    );
  }

  if (n_items <= 0 || n_levels <= 0 || n_raters <= 0) {
    return out;
  }

  Rcpp::IntegerMatrix item_counts(n_items, n_levels);
  double sum_item_agreement = 0.0;

  for (int i = 0; i < n_items; ++i) {
    for (int r = 0; r < n_raters; ++r) {
      const int code = ratings(i, r);
      if (code == NA_INTEGER || code < 1 || code > n_levels) {
        continue;
      }
      item_counts(i, code - 1) += 1;
      out.category_counts[static_cast<std::size_t>(code - 1)] += 1;
      if (exact) {
        out.rater_category_counts[static_cast<std::size_t>(r * n_levels + (code - 1))] += 1;
      }
    }

    double numer = 0.0;
    for (int j = 0; j < n_levels; ++j) {
      const int nij = item_counts(i, j);
      numer += static_cast<double>(nij) * static_cast<double>(nij - 1);
    }
    const double denom = static_cast<double>(n_raters) * static_cast<double>(n_raters - 1);
    const double agreement = numer / denom;
    out.item_agreement[static_cast<std::size_t>(i)] = agreement;
    sum_item_agreement += agreement;
  }

  if (out.n_items < 2 || out.n_categories < 2 || out.n_ratings_total <= 0) {
    return out;
  }

  for (int j = 0; j < n_levels; ++j) {
    out.category_proportions[static_cast<std::size_t>(j)] =
      static_cast<double>(out.category_counts[static_cast<std::size_t>(j)]) /
      static_cast<double>(out.n_ratings_total);
  }

  const double observed = sum_item_agreement / static_cast<double>(out.n_items);
  const KappaScalar scalar = exact
    ? compute_multirater_kappa_scalar_exact(
        observed,
        out.category_counts,
        out.rater_category_counts,
        out.n_items,
        out.n_raters,
        out.n_categories
      )
    : compute_multirater_kappa_scalar(
        observed,
        out.category_counts,
        out.n_ratings_total,
        out.n_categories,
        1
      );
  out.estimate = scalar.estimate;
  out.observed_agreement = scalar.observed_agreement;
  out.expected_agreement = scalar.expected_agreement;
  return out;
}

JackknifeStats compute_jackknife_se_counts(const Rcpp::IntegerMatrix& counts,
                                           const MultiKappaCore& full,
                                           const int method_code,
                                           const double conf_level,
                                           const int n_threads) {
  JackknifeStats out;
  const int n_items = counts.nrow();
  const int n_categories = counts.ncol();
  if (n_items < 3 || n_categories < 2 || !std::isfinite(full.estimate)) {
    return out;
  }

  const double sum_item_agreement = std::accumulate(
    full.item_agreement.begin(),
    full.item_agreement.end(),
    0.0
  );
  std::vector<double> theta_minus(static_cast<std::size_t>(n_items), NA_REAL);

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) if (n_items >= 32)
#endif
  for (int i = 0; i < n_items; ++i) {
    const int n_items_minus = n_items - 1;
    const int row_raters = full.n_raters_by_item[static_cast<std::size_t>(i)];
    const int total_ratings_minus = full.n_ratings_total - row_raters;
    if (n_items_minus < 2 || total_ratings_minus <= 0) {
      continue;
    }

    const double observed_minus =
      (sum_item_agreement - full.item_agreement[static_cast<std::size_t>(i)]) /
      static_cast<double>(n_items_minus);

    std::vector<int> category_counts_minus = full.category_counts;
    for (int j = 0; j < n_categories; ++j) {
      category_counts_minus[static_cast<std::size_t>(j)] -= counts(i, j);
    }

    const KappaScalar scalar = compute_multirater_kappa_scalar(
      observed_minus,
      category_counts_minus,
      total_ratings_minus,
      n_categories,
      method_code
    );
    theta_minus[static_cast<std::size_t>(i)] = scalar.estimate;
  }

  for (int i = 0; i < n_items; ++i) {
    if (!std::isfinite(theta_minus[static_cast<std::size_t>(i)])) {
      return out;
    }
  }

  const double theta_bar = std::accumulate(
    theta_minus.begin(),
    theta_minus.end(),
    0.0
  ) / static_cast<double>(n_items);

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

  return finalize_inference_from_se(
    full.estimate,
    std::sqrt(std::max(variance, 0.0)),
    conf_level
  );
}

JackknifeStats compute_jackknife_se_exact_ratings(const Rcpp::IntegerMatrix& ratings,
                                                  const Rcpp::IntegerMatrix& counts,
                                                  const MultiKappaCore& full,
                                                  const double conf_level,
                                                  const int n_threads) {
  JackknifeStats out;
  const int n_items = ratings.nrow();
  const int n_raters = ratings.ncol();
  const int n_categories = full.n_categories;
  if (n_items < 3 || n_raters < 2 || n_categories < 2 || !std::isfinite(full.estimate)) {
    return out;
  }

  const double sum_item_agreement = std::accumulate(
    full.item_agreement.begin(),
    full.item_agreement.end(),
    0.0
  );
  std::vector<double> theta_minus(static_cast<std::size_t>(n_items), NA_REAL);

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) if (n_items >= 32)
#endif
  for (int i = 0; i < n_items; ++i) {
    const int n_items_minus = n_items - 1;
    if (n_items_minus < 2) {
      continue;
    }

    const double observed_minus =
      (sum_item_agreement - full.item_agreement[static_cast<std::size_t>(i)]) /
      static_cast<double>(n_items_minus);

    std::vector<int> category_counts_minus = full.category_counts;
    for (int j = 0; j < n_categories; ++j) {
      category_counts_minus[static_cast<std::size_t>(j)] -= counts(i, j);
    }

    std::vector<int> rater_category_counts_minus = full.rater_category_counts;
    for (int r = 0; r < n_raters; ++r) {
      const int code = ratings(i, r);
      if (code == NA_INTEGER || code < 1 || code > n_categories) {
        continue;
      }
      rater_category_counts_minus[static_cast<std::size_t>(r * n_categories + (code - 1))] -= 1;
    }

    const KappaScalar scalar = compute_multirater_kappa_scalar_exact(
      observed_minus,
      category_counts_minus,
      rater_category_counts_minus,
      n_items_minus,
      n_raters,
      n_categories
    );
    theta_minus[static_cast<std::size_t>(i)] = scalar.estimate;
  }

  for (int i = 0; i < n_items; ++i) {
    if (!std::isfinite(theta_minus[static_cast<std::size_t>(i)])) {
      return out;
    }
  }

  const double theta_bar = std::accumulate(
    theta_minus.begin(),
    theta_minus.end(),
    0.0
  ) / static_cast<double>(n_items);

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

  return finalize_inference_from_se(
    full.estimate,
    std::sqrt(std::max(variance, 0.0)),
    conf_level
  );
}

JackknifeStats compute_asymptotic_se_fleiss_counts(const MultiKappaCore& full,
                                                   const double conf_level) {
  JackknifeStats out;
  if (!full.balanced || full.n_raters < 2 || full.n_items < 1 ||
      full.n_categories < 2 || !std::isfinite(full.estimate)) {
    return out;
  }

  const double ns = static_cast<double>(full.n_items);
  const double nr = static_cast<double>(full.n_raters);
  double sum_pjqj = 0.0;
  double correction = 0.0;
  for (int j = 0; j < full.n_categories; ++j) {
    const double pj = full.category_proportions[static_cast<std::size_t>(j)];
    const double qj = 1.0 - pj;
    const double term = pj * qj;
    sum_pjqj += term;
    correction += term * (qj - pj);
  }

  const double denom = sum_pjqj * sum_pjqj * (ns * nr * (nr - 1.0));
  if (!std::isfinite(denom) || denom <= 0.0) {
    return out;
  }
  const double variance = (2.0 / denom) * ((sum_pjqj * sum_pjqj) - correction);
  if (!std::isfinite(variance)) {
    return out;
  }

  return finalize_inference_from_se(
    full.estimate,
    std::sqrt(std::max(variance, 0.0)),
    conf_level,
    true
  );
}

Rcpp::IntegerMatrix collapse_counts_for_category(const Rcpp::IntegerMatrix& counts,
                                                 const int category_index) {
  const int n_items = counts.nrow();
  const int n_categories = counts.ncol();
  Rcpp::IntegerMatrix out(n_items, 2);
  for (int i = 0; i < n_items; ++i) {
    int target = counts(i, category_index);
    int total = 0;
    for (int j = 0; j < n_categories; ++j) {
      total += counts(i, j);
    }
    out(i, 0) = target;
    out(i, 1) = total - target;
  }
  return out;
}

Rcpp::IntegerMatrix collapse_ratings_for_category(const Rcpp::IntegerMatrix& ratings,
                                                  const int category_index) {
  const int n_items = ratings.nrow();
  const int n_raters = ratings.ncol();
  Rcpp::IntegerMatrix out(n_items, n_raters);
  const int target_code = category_index + 1;
  for (int i = 0; i < n_items; ++i) {
    for (int r = 0; r < n_raters; ++r) {
      const int code = ratings(i, r);
      if (code == NA_INTEGER) {
        out(i, r) = NA_INTEGER;
      } else {
        out(i, r) = (code == target_code) ? 1 : 2;
      }
    }
  }
  return out;
}

Rcpp::IntegerMatrix ratings_to_counts_complete(const Rcpp::IntegerMatrix& ratings,
                                               const int n_levels) {
  const int n_items = ratings.nrow();
  const int n_raters = ratings.ncol();
  Rcpp::IntegerMatrix counts(n_items, n_levels);
  for (int i = 0; i < n_items; ++i) {
    for (int r = 0; r < n_raters; ++r) {
      const int code = ratings(i, r);
      if (code != NA_INTEGER && code >= 1 && code <= n_levels) {
        counts(i, code - 1) += 1;
      }
    }
  }
  return counts;
}

Rcpp::List build_multirater_result_from_counts(const Rcpp::IntegerMatrix& counts,
                                               const int method_code,
                                               const bool by_category,
                                               const int inference_code,
                                               const double conf_level,
                                               const int n_threads) {
  const MultiKappaCore full = compute_multirater_kappa_from_counts(counts, method_code);
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("estimate") = full.estimate,
    Rcpp::Named("observed_agreement") = full.observed_agreement,
    Rcpp::Named("expected_agreement") = full.expected_agreement,
    Rcpp::Named("n_items") = full.n_items,
    Rcpp::Named("n_categories") = full.n_categories,
    Rcpp::Named("n_ratings_total") = full.n_ratings_total,
    Rcpp::Named("n_raters_min") = full.n_raters_min,
    Rcpp::Named("n_raters_max") = full.n_raters_max,
    Rcpp::Named("balanced") = full.balanced,
    Rcpp::Named("exact") = false,
    Rcpp::Named("asymptotic_available") = (method_code == 1 && full.balanced && full.n_raters >= 2),
    Rcpp::Named("category_counts") = Rcpp::wrap(full.category_counts),
    Rcpp::Named("category_proportions") = Rcpp::wrap(full.category_proportions),
    Rcpp::Named("item_category_counts") = counts,
    Rcpp::Named("item_agreement") = Rcpp::wrap(full.item_agreement),
    Rcpp::Named("n_raters_by_item") = Rcpp::wrap(full.n_raters_by_item)
  );

  if (inference_code != 0) {
    const JackknifeStats inf = (inference_code == 1)
      ? compute_asymptotic_se_fleiss_counts(full, conf_level)
      : compute_jackknife_se_counts(
          counts,
          full,
          method_code,
          conf_level,
          n_threads
        );
    out["se"] = inf.se;
    out["lwr"] = inf.lwr;
    out["upr"] = inf.upr;
    out["statistic"] = inf.statistic;
    out["p_value"] = inf.p_value;
    out["conf_level"] = conf_level;
    out["se_method"] = (inference_code == 1) ? "asymptotic" : "jackknife";
  }

  if (by_category) {
    const int n_categories = counts.ncol();
    Rcpp::NumericVector by_estimate(n_categories, NA_REAL);
    Rcpp::NumericVector by_observed(n_categories, NA_REAL);
    Rcpp::NumericVector by_expected(n_categories, NA_REAL);
    Rcpp::NumericVector by_se(n_categories, NA_REAL);
    Rcpp::NumericVector by_lwr(n_categories, NA_REAL);
    Rcpp::NumericVector by_upr(n_categories, NA_REAL);
    Rcpp::NumericVector by_statistic(n_categories, NA_REAL);
    Rcpp::NumericVector by_p_value(n_categories, NA_REAL);

    for (int j = 0; j < n_categories; ++j) {
      const Rcpp::IntegerMatrix collapsed = collapse_counts_for_category(counts, j);
      const MultiKappaCore cat_full = compute_multirater_kappa_from_counts(collapsed, method_code);
      by_estimate[j] = cat_full.estimate;
      by_observed[j] = cat_full.observed_agreement;
      by_expected[j] = cat_full.expected_agreement;

      if (inference_code != 0) {
        const JackknifeStats cat_inf = (inference_code == 1)
          ? compute_asymptotic_se_fleiss_counts(cat_full, conf_level)
          : compute_jackknife_se_counts(
              collapsed,
              cat_full,
              method_code,
              conf_level,
              n_threads
            );
        by_se[j] = cat_inf.se;
        by_lwr[j] = cat_inf.lwr;
        by_upr[j] = cat_inf.upr;
        by_statistic[j] = cat_inf.statistic;
        by_p_value[j] = cat_inf.p_value;
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

Rcpp::List build_multirater_result_from_ratings_complete(const Rcpp::IntegerMatrix& ratings,
                                                         const int n_levels,
                                                         const bool exact,
                                                         const bool by_category,
                                                         const int inference_code,
                                                         const double conf_level,
                                                         const int n_threads) {
  const MultiKappaCore full = compute_multirater_kappa_from_ratings_complete(
    ratings,
    n_levels,
    exact
  );
  const Rcpp::IntegerMatrix counts = ratings_to_counts_complete(ratings, n_levels);

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("estimate") = full.estimate,
    Rcpp::Named("observed_agreement") = full.observed_agreement,
    Rcpp::Named("expected_agreement") = full.expected_agreement,
    Rcpp::Named("n_items") = full.n_items,
    Rcpp::Named("n_categories") = full.n_categories,
    Rcpp::Named("n_ratings_total") = full.n_ratings_total,
    Rcpp::Named("n_raters_min") = full.n_raters_min,
    Rcpp::Named("n_raters_max") = full.n_raters_max,
    Rcpp::Named("balanced") = full.balanced,
    Rcpp::Named("exact") = true,
    Rcpp::Named("asymptotic_available") = false,
    Rcpp::Named("category_counts") = Rcpp::wrap(full.category_counts),
    Rcpp::Named("category_proportions") = Rcpp::wrap(full.category_proportions),
    Rcpp::Named("item_category_counts") = counts,
    Rcpp::Named("item_agreement") = Rcpp::wrap(full.item_agreement),
    Rcpp::Named("n_raters_by_item") = Rcpp::wrap(full.n_raters_by_item)
  );

  if (inference_code == 2) {
    const JackknifeStats jk = compute_jackknife_se_exact_ratings(
      ratings,
      counts,
      full,
      conf_level,
      n_threads
    );
    out["se"] = jk.se;
    out["lwr"] = jk.lwr;
    out["upr"] = jk.upr;
    out["statistic"] = jk.statistic;
    out["p_value"] = jk.p_value;
    out["conf_level"] = conf_level;
    out["se_method"] = "jackknife";
  }

  if (by_category) {
    const int n_categories = n_levels;
    Rcpp::NumericVector by_estimate(n_categories, NA_REAL);
    Rcpp::NumericVector by_observed(n_categories, NA_REAL);
    Rcpp::NumericVector by_expected(n_categories, NA_REAL);
    Rcpp::NumericVector by_se(n_categories, NA_REAL);
    Rcpp::NumericVector by_lwr(n_categories, NA_REAL);
    Rcpp::NumericVector by_upr(n_categories, NA_REAL);
    Rcpp::NumericVector by_statistic(n_categories, NA_REAL);
    Rcpp::NumericVector by_p_value(n_categories, NA_REAL);

    for (int j = 0; j < n_categories; ++j) {
      const Rcpp::IntegerMatrix collapsed_ratings = collapse_ratings_for_category(ratings, j);
      const MultiKappaCore cat_full = compute_multirater_kappa_from_ratings_complete(
        collapsed_ratings,
        2,
        true
      );
      by_estimate[j] = cat_full.estimate;
      by_observed[j] = cat_full.observed_agreement;
      by_expected[j] = cat_full.expected_agreement;

      if (inference_code == 2) {
        const Rcpp::IntegerMatrix collapsed_counts = ratings_to_counts_complete(collapsed_ratings, 2);
        const JackknifeStats cat_jk = compute_jackknife_se_exact_ratings(
          collapsed_ratings,
          collapsed_counts,
          cat_full,
          conf_level,
          n_threads
        );
        by_se[j] = cat_jk.se;
        by_lwr[j] = cat_jk.lwr;
        by_upr[j] = cat_jk.upr;
        by_statistic[j] = cat_jk.statistic;
        by_p_value[j] = cat_jk.p_value;
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
Rcpp::List multirater_kappa_counts_cpp(Rcpp::IntegerMatrix counts,
                                       int method_code = 1,
                                       bool by_category = false,
                                       int inference_code = 0,
                                       double conf_level = 0.95,
                                       int n_threads = 1) {
  return build_multirater_result_from_counts(
    counts,
    method_code,
    by_category,
    inference_code,
    conf_level,
    n_threads
  );
}

// [[Rcpp::export]]
Rcpp::List multirater_kappa_ratings_cpp(Rcpp::IntegerMatrix ratings,
                                        int n_levels,
                                        int method_code = 1,
                                        int na_code = 1,
                                        int min_raters = 2,
                                        bool by_category = false,
                                        int inference_code = 0,
                                        double conf_level = 0.95,
                                        int n_threads = 1,
                                        bool exact = false) {
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

  Rcpp::IntegerMatrix retained(static_cast<int>(keep_rows.size()), n_raters);
  for (int out_i = 0; out_i < static_cast<int>(keep_rows.size()); ++out_i) {
    const int src_i = keep_rows[static_cast<std::size_t>(out_i)];
    for (int j = 0; j < n_raters; ++j) {
      retained(out_i, j) = ratings(src_i, j);
    }
  }

  if (exact) {
    return build_multirater_result_from_ratings_complete(
      retained,
      n_levels,
      true,
      by_category,
      inference_code,
      conf_level,
      n_threads
    );
  }

  Rcpp::IntegerMatrix counts(static_cast<int>(keep_rows.size()), n_levels);
  for (int out_i = 0; out_i < static_cast<int>(keep_rows.size()); ++out_i) {
    for (int j = 0; j < n_raters; ++j) {
      const int code = retained(out_i, j);
      if (code != NA_INTEGER && code >= 1 && code <= n_levels) {
        counts(out_i, code - 1) += 1;
      }
    }
  }

  return build_multirater_result_from_counts(
    counts,
    method_code,
    by_category,
    inference_code,
    conf_level,
    n_threads
  );
}
