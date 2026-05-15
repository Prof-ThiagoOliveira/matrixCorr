#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "matrixCorr_omp.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

namespace {

inline double clamp_unit(const double x, const double tol = 1e-12) {
  if (!std::isfinite(x)) {
    return NA_REAL;
  }
  if (x > 1.0 && x <= 1.0 + tol) {
    return 1.0;
  }
  if (x < -1.0 && x >= -1.0 - tol) {
    return -1.0;
  }
  return x;
}

struct AlphaCore {
  double customary_alpha = NA_REAL;
  double analytical_alpha = NA_REAL;
  double observed_disagreement = NA_REAL;
  double expected_disagreement = NA_REAL;
  double mse = NA_REAL;
  double sst = NA_REAL;
  double ssa = NA_REAL;
  double msa = NA_REAL;
  double theta = NA_REAL;
  double eta = NA_REAL;
  double n_bar = NA_REAL;
  int n_items = 0;
  int n_categories = 0;
  int n_ratings_total = 0;
  int n_raters_min = 0;
  int n_raters_max = 0;
  bool balanced = true;
  std::vector<int> category_counts;
  std::vector<double> category_proportions;
  std::vector<int> n_raters_by_item;
  std::vector<double> item_disagreement;
  std::vector<double> coincidence;
  std::vector<double> expected_coincidence;
};

AlphaCore compute_alpha_core(const Rcpp::IntegerMatrix& counts,
                             const Rcpp::NumericMatrix& delta2,
                             const bool need_matrices,
                             const int n_threads) {
  AlphaCore out;
  const int n_items = counts.nrow();
  const int n_categories = counts.ncol();
  out.n_items = n_items;
  out.n_categories = n_categories;
  out.category_counts.assign(static_cast<std::size_t>(std::max(0, n_categories)), 0);
  out.category_proportions.assign(static_cast<std::size_t>(std::max(0, n_categories)), NA_REAL);
  out.n_raters_by_item.assign(static_cast<std::size_t>(std::max(0, n_items)), 0);
  out.item_disagreement.assign(static_cast<std::size_t>(std::max(0, n_items)), NA_REAL);
  out.coincidence.assign(static_cast<std::size_t>(std::max(0, n_categories * n_categories)), 0.0);

  if (n_items <= 0 || n_categories <= 0) {
    return out;
  }

  bool first_row = true;
  for (int i = 0; i < n_items; ++i) {
    int m_i = 0;
    for (int j = 0; j < n_categories; ++j) {
      const int nij = counts(i, j);
      m_i += nij;
      out.category_counts[static_cast<std::size_t>(j)] += nij;
    }
    out.n_raters_by_item[static_cast<std::size_t>(i)] = m_i;
    out.n_ratings_total += m_i;
    if (first_row) {
      out.n_raters_min = m_i;
      out.n_raters_max = m_i;
      first_row = false;
    } else {
      out.n_raters_min = std::min(out.n_raters_min, m_i);
      out.n_raters_max = std::max(out.n_raters_max, m_i);
    }
  }
  out.balanced = (out.n_raters_min == out.n_raters_max);

  if (out.n_ratings_total <= 0) {
    return out;
  }

  for (int j = 0; j < n_categories; ++j) {
    out.category_proportions[static_cast<std::size_t>(j)] =
      static_cast<double>(out.category_counts[static_cast<std::size_t>(j)]) /
      static_cast<double>(out.n_ratings_total);
  }

  std::vector<double> within_sum(static_cast<std::size_t>(n_items), 0.0);

#ifdef _OPENMP
  const int nthreads_use = std::max(1, n_threads);
  std::vector< std::vector<double> > coincidence_tls(
    static_cast<std::size_t>(nthreads_use),
    std::vector<double>(static_cast<std::size_t>(n_categories * n_categories), 0.0)
  );
#pragma omp parallel for num_threads(nthreads_use) if(n_items >= 32)
  for (int i = 0; i < n_items; ++i) {
    const int tid = omp_get_thread_num();
    std::vector<double>& coinc_local = coincidence_tls[static_cast<std::size_t>(tid)];
    const int m_i = out.n_raters_by_item[static_cast<std::size_t>(i)];
    if (m_i <= 1) {
      continue;
    }

    double item_within = 0.0;
    const double denom = static_cast<double>(m_i - 1);
    for (int c = 0; c < n_categories; ++c) {
      const int n_ic = counts(i, c);
      if (n_ic <= 0) {
        continue;
      }
      for (int k = 0; k < n_categories; ++k) {
        const int n_ik = counts(i, k);
        if (n_ik <= 0) {
          continue;
        }
        coinc_local[static_cast<std::size_t>(c * n_categories + k)] +=
          static_cast<double>(n_ic) *
          (static_cast<double>(n_ik) - (c == k ? 1.0 : 0.0)) / denom;
      }
      for (int k = c + 1; k < n_categories; ++k) {
        const int n_ik = counts(i, k);
        if (n_ik <= 0) {
          continue;
        }
        item_within += static_cast<double>(n_ic) *
          static_cast<double>(n_ik) * delta2(c, k);
      }
    }
    within_sum[static_cast<std::size_t>(i)] = item_within;
    out.item_disagreement[static_cast<std::size_t>(i)] =
      (2.0 * item_within) /
      (static_cast<double>(m_i) * static_cast<double>(m_i - 1));
  }

  for (int t = 0; t < nthreads_use; ++t) {
    const std::vector<double>& src = coincidence_tls[static_cast<std::size_t>(t)];
    for (int idx = 0; idx < n_categories * n_categories; ++idx) {
      out.coincidence[static_cast<std::size_t>(idx)] += src[static_cast<std::size_t>(idx)];
    }
  }
#else
  (void)n_threads;
  for (int i = 0; i < n_items; ++i) {
    const int m_i = out.n_raters_by_item[static_cast<std::size_t>(i)];
    if (m_i <= 1) {
      continue;
    }

    double item_within = 0.0;
    const double denom = static_cast<double>(m_i - 1);
    for (int c = 0; c < n_categories; ++c) {
      const int n_ic = counts(i, c);
      if (n_ic <= 0) {
        continue;
      }
      for (int k = 0; k < n_categories; ++k) {
        const int n_ik = counts(i, k);
        if (n_ik <= 0) {
          continue;
        }
        out.coincidence[static_cast<std::size_t>(c * n_categories + k)] +=
          static_cast<double>(n_ic) *
          (static_cast<double>(n_ik) - (c == k ? 1.0 : 0.0)) / denom;
      }
      for (int k = c + 1; k < n_categories; ++k) {
        const int n_ik = counts(i, k);
        if (n_ik <= 0) {
          continue;
        }
        item_within += static_cast<double>(n_ic) *
          static_cast<double>(n_ik) * delta2(c, k);
      }
    }
    within_sum[static_cast<std::size_t>(i)] = item_within;
    out.item_disagreement[static_cast<std::size_t>(i)] =
      (2.0 * item_within) /
      (static_cast<double>(m_i) * static_cast<double>(m_i - 1));
  }
#endif

  double observed_num = 0.0;
  for (int c = 0; c < n_categories; ++c) {
    for (int k = 0; k < n_categories; ++k) {
      observed_num += out.coincidence[static_cast<std::size_t>(c * n_categories + k)] * delta2(c, k);
    }
  }
  if (out.n_ratings_total > 0) {
    out.observed_disagreement = observed_num / static_cast<double>(out.n_ratings_total);
  }

  double expected_num = 0.0;
  double total_unordered = 0.0;
  for (int c = 0; c < n_categories; ++c) {
    const double n_c = static_cast<double>(out.category_counts[static_cast<std::size_t>(c)]);
    for (int k = 0; k < n_categories; ++k) {
      const double n_k = static_cast<double>(out.category_counts[static_cast<std::size_t>(k)]);
      expected_num += n_c * n_k * delta2(c, k);
    }
    for (int k = c + 1; k < n_categories; ++k) {
      total_unordered += n_c *
        static_cast<double>(out.category_counts[static_cast<std::size_t>(k)]) * delta2(c, k);
    }
  }

  if (out.n_ratings_total > 1) {
    out.expected_disagreement =
      expected_num /
      (static_cast<double>(out.n_ratings_total) * static_cast<double>(out.n_ratings_total - 1));
  }

  if (std::isfinite(out.observed_disagreement) &&
      std::isfinite(out.expected_disagreement) &&
      out.expected_disagreement > 0.0) {
    out.customary_alpha = clamp_unit(
      1.0 - (out.observed_disagreement / out.expected_disagreement)
    );
  }

  const int a = n_items;
  const int N = out.n_ratings_total;
  if (a >= 2 && N > a) {
    double sum_m_sq = 0.0;
    double mse_num = 0.0;
    for (int i = 0; i < n_items; ++i) {
      const int m_i = out.n_raters_by_item[static_cast<std::size_t>(i)];
      sum_m_sq += static_cast<double>(m_i) * static_cast<double>(m_i);
      if (m_i > 1) {
        mse_num += within_sum[static_cast<std::size_t>(i)] /
          static_cast<double>(m_i - 1);
      }
    }

    out.n_bar =
      (static_cast<double>(N) - sum_m_sq / static_cast<double>(N)) /
      static_cast<double>(a - 1);
    out.mse = mse_num / static_cast<double>(N);
    out.sst = total_unordered / static_cast<double>(N);
    out.ssa = out.sst - static_cast<double>(N - a) * out.mse;
    out.msa = out.ssa / static_cast<double>(a - 1);

    if (std::isfinite(out.mse) && out.mse > 0.0 &&
        std::isfinite(out.msa) && out.msa > 0.0) {
      out.theta = out.msa / out.mse;
      if (std::isfinite(out.theta) && out.theta > 0.0) {
        out.eta = std::log(out.theta);
        const double denom = out.theta + out.n_bar - 1.0;
        if (std::isfinite(denom) && denom != 0.0) {
          out.analytical_alpha = clamp_unit(
            (out.theta - 1.0) / denom
          );
        }
      }
    }
  }

  if (need_matrices && out.n_ratings_total > 1) {
    out.expected_coincidence.assign(
      static_cast<std::size_t>(std::max(0, n_categories * n_categories)),
      0.0
    );
    const double denom = static_cast<double>(out.n_ratings_total - 1);
    for (int c = 0; c < n_categories; ++c) {
      const double n_c = static_cast<double>(out.category_counts[static_cast<std::size_t>(c)]);
      for (int k = 0; k < n_categories; ++k) {
        const double n_k = static_cast<double>(out.category_counts[static_cast<std::size_t>(k)]);
        out.expected_coincidence[static_cast<std::size_t>(c * n_categories + k)] =
          n_c * (n_k - (c == k ? 1.0 : 0.0)) / denom;
      }
    }
  }

  return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List krippendorff_alpha_core_cpp(Rcpp::IntegerMatrix counts,
                                       Rcpp::NumericMatrix delta2,
                                       int method_code = 1,
                                       bool return_matrices = false,
                                       int n_threads = 1) {
  const AlphaCore core = compute_alpha_core(counts, delta2, return_matrices, n_threads);
  const double estimate = (method_code == 2)
    ? core.analytical_alpha
    : core.customary_alpha;

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("estimate") = estimate,
    Rcpp::Named("customary_alpha") = core.customary_alpha,
    Rcpp::Named("analytical_alpha") = core.analytical_alpha,
    Rcpp::Named("observed_disagreement") = core.observed_disagreement,
    Rcpp::Named("expected_disagreement") = core.expected_disagreement,
    Rcpp::Named("mse") = core.mse,
    Rcpp::Named("sst") = core.sst,
    Rcpp::Named("ssa") = core.ssa,
    Rcpp::Named("msa") = core.msa,
    Rcpp::Named("theta") = core.theta,
    Rcpp::Named("eta") = core.eta,
    Rcpp::Named("n_bar") = core.n_bar,
    Rcpp::Named("n_items") = core.n_items,
    Rcpp::Named("n_categories") = core.n_categories,
    Rcpp::Named("n_ratings_total") = core.n_ratings_total,
    Rcpp::Named("n_raters_min") = core.n_raters_min,
    Rcpp::Named("n_raters_max") = core.n_raters_max,
    Rcpp::Named("balanced") = core.balanced,
    Rcpp::Named("category_counts") = Rcpp::wrap(core.category_counts),
    Rcpp::Named("category_proportions") = Rcpp::wrap(core.category_proportions),
    Rcpp::Named("n_raters_by_item") = Rcpp::wrap(core.n_raters_by_item),
    Rcpp::Named("item_disagreement") = Rcpp::wrap(core.item_disagreement)
  );

  if (return_matrices) {
    Rcpp::NumericMatrix coincidence(core.n_categories, core.n_categories);
    Rcpp::NumericMatrix expected(core.n_categories, core.n_categories);
    for (int c = 0; c < core.n_categories; ++c) {
      for (int k = 0; k < core.n_categories; ++k) {
        coincidence(c, k) = core.coincidence[static_cast<std::size_t>(c * core.n_categories + k)];
        if (!core.expected_coincidence.empty()) {
          expected(c, k) = core.expected_coincidence[static_cast<std::size_t>(c * core.n_categories + k)];
        }
      }
    }
    out["coincidence"] = coincidence;
    out["expected_coincidence"] = expected;
    out["delta2"] = delta2;
    out["item_category_counts"] = counts;
  }

  return out;
}
