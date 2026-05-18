// Thiago de Paula Oliveira
// cia.cpp
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

const char* kOverallBalancedError =
  "Overall CIA currently implements the balanced-replication moment estimator from the cited papers; use pairwise CIA for unequal replication.";

double within_msd_subject(const std::vector<double>& x) {
  const std::size_t n = x.size();
  if (n < 2) {
    return NA_REAL;
  }
  double sum = 0.0;
  std::size_t count = 0;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const double d = x[i] - x[j];
      sum += d * d;
      count += 1;
    }
  }
  return count > 0 ? sum / static_cast<double>(count) : NA_REAL;
}

double between_msd_subject(const std::vector<double>& x, const std::vector<double>& y) {
  if (x.empty() || y.empty()) {
    return NA_REAL;
  }
  double sum = 0.0;
  std::size_t count = 0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    for (std::size_t j = 0; j < y.size(); ++j) {
      const double d = x[i] - y[j];
      sum += d * d;
      count += 1;
    }
  }
  return count > 0 ? sum / static_cast<double>(count) : NA_REAL;
}

double sample_var_from_sums(const double sum, const double sum_sq, const int n) {
  if (n < 2) {
    return NA_REAL;
  }
  const double nn = static_cast<double>(n);
  const double centered = sum_sq - (sum * sum) / nn;
  return centered / static_cast<double>(n - 1);
}

double sample_cov_from_sums(const double sum_xy,
                            const double sum_x,
                            const double sum_y,
                            const int n) {
  if (n < 2) {
    return NA_REAL;
  }
  const double nn = static_cast<double>(n);
  const double centered = sum_xy - (sum_x * sum_y) / nn;
  return centered / static_cast<double>(n - 1);
}

std::vector< std::vector< std::vector<double> > > build_subject_method_values(
    const Rcpp::NumericVector& y,
    const Rcpp::IntegerVector& subject,
    const Rcpp::IntegerVector& method,
    const int n_methods,
    int* n_subjects_out) {
  const R_xlen_t n_obs = y.size();
  if (subject.size() != n_obs || method.size() != n_obs) {
    Rcpp::stop("y, subject, and method must have the same length.");
  }
  if (n_methods <= 0) {
    Rcpp::stop("n_methods must be positive.");
  }

  int n_subjects = 0;
  for (R_xlen_t i = 0; i < n_obs; ++i) {
    if (subject[i] != NA_INTEGER && subject[i] > n_subjects) {
      n_subjects = subject[i];
    }
  }
  *n_subjects_out = n_subjects;

  std::vector< std::vector< std::vector<double> > > by_subject(
    static_cast<std::size_t>(std::max(n_subjects, 0)),
    std::vector< std::vector<double> >(static_cast<std::size_t>(n_methods))
  );

  for (R_xlen_t i = 0; i < n_obs; ++i) {
    if (Rcpp::NumericVector::is_na(y[i]) ||
        subject[i] == NA_INTEGER ||
        method[i] == NA_INTEGER) {
      continue;
    }
    const int s = subject[i];
    const int m = method[i];
    if (s < 1 || s > n_subjects || m < 1 || m > n_methods) {
      continue;
    }
    by_subject[static_cast<std::size_t>(s - 1)][static_cast<std::size_t>(m - 1)].push_back(y[i]);
  }

  return by_subject;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cia_moments_cpp(Rcpp::NumericVector y,
                           Rcpp::IntegerVector subject,
                           Rcpp::IntegerVector method,
                           Rcpp::IntegerVector replicate,
                           int n_methods,
                           int reference_method,
                           bool has_reference,
                           bool pairwise,
                           int n_threads = 1) {
  (void)replicate;
  (void)reference_method;
  (void)has_reference;
  (void)pairwise;
  (void)n_threads;

  int n_subjects = 0;
  const std::vector< std::vector< std::vector<double> > > by_subject =
    build_subject_method_values(y, subject, method, n_methods, &n_subjects);

  std::vector<double> within_sum(static_cast<std::size_t>(n_methods), 0.0);
  std::vector<int> n_replicate_pairs(static_cast<std::size_t>(n_methods), 0);
  std::vector< std::vector<double> > between_sum(
    static_cast<std::size_t>(n_methods),
    std::vector<double>(static_cast<std::size_t>(n_methods), 0.0)
  );
  std::vector< std::vector<int> > n_between_pairs(
    static_cast<std::size_t>(n_methods),
    std::vector<int>(static_cast<std::size_t>(n_methods), 0)
  );

  for (int s = 0; s < n_subjects; ++s) {
    const std::vector< std::vector<double> >& subj = by_subject[static_cast<std::size_t>(s)];

    for (int j = 0; j < n_methods; ++j) {
      const std::vector<double>& vals = subj[static_cast<std::size_t>(j)];
      const int nv = static_cast<int>(vals.size());
      if (nv < 2) {
        continue;
      }
      for (int a = 0; a < nv - 1; ++a) {
        for (int b = a + 1; b < nv; ++b) {
          const double d = vals[static_cast<std::size_t>(a)] - vals[static_cast<std::size_t>(b)];
          within_sum[static_cast<std::size_t>(j)] += d * d;
          n_replicate_pairs[static_cast<std::size_t>(j)] += 1;
        }
      }
    }

    for (int j = 0; j < n_methods - 1; ++j) {
      const std::vector<double>& lhs = subj[static_cast<std::size_t>(j)];
      if (lhs.empty()) {
        continue;
      }
      for (int k = j + 1; k < n_methods; ++k) {
        const std::vector<double>& rhs = subj[static_cast<std::size_t>(k)];
        if (rhs.empty()) {
          continue;
        }
        for (std::size_t a = 0; a < lhs.size(); ++a) {
          for (std::size_t b = 0; b < rhs.size(); ++b) {
            const double d = lhs[a] - rhs[b];
            between_sum[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] += d * d;
            n_between_pairs[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] += 1;
          }
        }
      }
    }
  }

  Rcpp::NumericVector within_msd(n_methods, NA_REAL);
  Rcpp::NumericVector sigma_within(n_methods, NA_REAL);
  Rcpp::IntegerVector replicate_pairs(n_methods);
  for (int j = 0; j < n_methods; ++j) {
    const int count = n_replicate_pairs[static_cast<std::size_t>(j)];
    replicate_pairs[j] = count;
    if (count > 0) {
      const double msd = within_sum[static_cast<std::size_t>(j)] / static_cast<double>(count);
      within_msd[j] = msd;
      sigma_within[j] = msd / 2.0;
    }
  }

  Rcpp::NumericMatrix between_msd(n_methods, n_methods);
  Rcpp::IntegerMatrix between_pairs(n_methods, n_methods);
  std::fill(between_msd.begin(), between_msd.end(), NA_REAL);
  std::fill(between_pairs.begin(), between_pairs.end(), 0);
  for (int j = 0; j < n_methods; ++j) {
    between_msd(j, j) = 0.0;
  }
  for (int j = 0; j < n_methods - 1; ++j) {
    for (int k = j + 1; k < n_methods; ++k) {
      const int count = n_between_pairs[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
      between_pairs(j, k) = count;
      between_pairs(k, j) = count;
      if (count > 0) {
        const double msd =
          between_sum[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] /
          static_cast<double>(count);
        between_msd(j, k) = msd;
        between_msd(k, j) = msd;
      }
    }
  }

  // Diagnostics-only backend: authoritative CIA estimates are computed in
  // cia_pairwise_stats_cpp() and cia_overall_balanced_cpp().
  return Rcpp::List::create(
    Rcpp::_["estimate"] = R_NilValue,
    Rcpp::_["within_msd"] = within_msd,
    Rcpp::_["between_msd"] = between_msd,
    Rcpp::_["sigma_within"] = sigma_within,
    Rcpp::_["n_subjects"] = n_subjects,
    Rcpp::_["n_obs"] = static_cast<int>(y.size()),
    Rcpp::_["n_replicate_pairs"] = replicate_pairs,
    Rcpp::_["n_between_pairs"] = between_pairs
  );
}

// [[Rcpp::export]]
Rcpp::List cia_pairwise_stats_cpp(Rcpp::NumericVector y,
                                  Rcpp::IntegerVector subject,
                                  Rcpp::IntegerVector method,
                                  Rcpp::IntegerVector replicate,
                                  int n_methods,
                                  int reference_method,
                                  bool has_reference,
                                  int n_threads = 1) {
  (void)replicate;
  (void)n_threads;

  int n_subjects = 0;
  const std::vector< std::vector< std::vector<double> > > by_subject =
    build_subject_method_values(y, subject, method, n_methods, &n_subjects);

  Rcpp::NumericMatrix estimate_raw(n_methods, n_methods);
  Rcpp::NumericMatrix A_bar(n_methods, n_methods);
  Rcpp::NumericMatrix B_bar(n_methods, n_methods);
  Rcpp::NumericMatrix S11(n_methods, n_methods);
  Rcpp::NumericMatrix S22(n_methods, n_methods);
  Rcpp::NumericMatrix S33(n_methods, n_methods);
  Rcpp::NumericMatrix S12(n_methods, n_methods);
  Rcpp::NumericMatrix S13(n_methods, n_methods);
  Rcpp::NumericMatrix S23(n_methods, n_methods);
  Rcpp::IntegerMatrix n_eligible(n_methods, n_methods);

  std::fill(estimate_raw.begin(), estimate_raw.end(), NA_REAL);
  std::fill(A_bar.begin(), A_bar.end(), NA_REAL);
  std::fill(B_bar.begin(), B_bar.end(), NA_REAL);
  std::fill(S11.begin(), S11.end(), NA_REAL);
  std::fill(S22.begin(), S22.end(), NA_REAL);
  std::fill(S33.begin(), S33.end(), NA_REAL);
  std::fill(S12.begin(), S12.end(), NA_REAL);
  std::fill(S13.begin(), S13.end(), NA_REAL);
  std::fill(S23.begin(), S23.end(), NA_REAL);
  std::fill(n_eligible.begin(), n_eligible.end(), 0);

  for (int j = 0; j < n_methods; ++j) {
    estimate_raw(j, j) = 1.0;
  }

  if (has_reference) {
    const int ref = reference_method - 1;
    if (ref < 0 || ref >= n_methods) {
      Rcpp::stop("reference_method is out of bounds.");
    }
    for (int j = 0; j < n_methods; ++j) {
      if (j == ref) {
        continue;
      }
      double sum_g1 = 0.0;
      double sum_g3 = 0.0;
      double sum_g1_sq = 0.0;
      double sum_g3_sq = 0.0;
      double sum_g1g3 = 0.0;
      int n = 0;

      for (int s = 0; s < n_subjects; ++s) {
        const std::vector<double>& ref_vals =
          by_subject[static_cast<std::size_t>(s)][static_cast<std::size_t>(ref)];
        const std::vector<double>& cmp_vals =
          by_subject[static_cast<std::size_t>(s)][static_cast<std::size_t>(j)];
        if (ref_vals.size() < 2 || cmp_vals.empty()) {
          continue;
        }
        const double g1 = within_msd_subject(ref_vals);
        const double g3 = between_msd_subject(ref_vals, cmp_vals);
        if (!R_finite(g1) || !R_finite(g3)) {
          continue;
        }
        sum_g1 += g1;
        sum_g3 += g3;
        sum_g1_sq += g1 * g1;
        sum_g3_sq += g3 * g3;
        sum_g1g3 += g1 * g3;
        n += 1;
      }

      n_eligible(j, ref) = n;
      n_eligible(ref, j) = n;
      if (n <= 0) {
        continue;
      }

      const double mean_g1 = sum_g1 / static_cast<double>(n);
      const double mean_g3 = sum_g3 / static_cast<double>(n);
      const double est = (mean_g3 != 0.0) ? mean_g1 / mean_g3 : NA_REAL;
      const double s11 = sample_var_from_sums(sum_g1, sum_g1_sq, n);
      const double s33 = sample_var_from_sums(sum_g3, sum_g3_sq, n);
      const double s13 = sample_cov_from_sums(sum_g1g3, sum_g1, sum_g3, n);

      estimate_raw(j, ref) = est;
      estimate_raw(ref, j) = est;
      A_bar(j, ref) = mean_g1;
      A_bar(ref, j) = mean_g1;
      B_bar(j, ref) = mean_g3;
      B_bar(ref, j) = mean_g3;
      S11(j, ref) = s11;
      S11(ref, j) = s11;
      S33(j, ref) = s33;
      S33(ref, j) = s33;
      S13(j, ref) = s13;
      S13(ref, j) = s13;
    }
  } else {
    for (int j = 0; j < n_methods - 1; ++j) {
      for (int k = j + 1; k < n_methods; ++k) {
        double sum_g1 = 0.0;
        double sum_g2 = 0.0;
        double sum_g3 = 0.0;
        double sum_g1_sq = 0.0;
        double sum_g2_sq = 0.0;
        double sum_g3_sq = 0.0;
        double sum_g1g2 = 0.0;
        double sum_g1g3 = 0.0;
        double sum_g2g3 = 0.0;
        int n = 0;

        for (int s = 0; s < n_subjects; ++s) {
          const std::vector<double>& lhs =
            by_subject[static_cast<std::size_t>(s)][static_cast<std::size_t>(j)];
          const std::vector<double>& rhs =
            by_subject[static_cast<std::size_t>(s)][static_cast<std::size_t>(k)];
          if (lhs.size() < 2 || rhs.size() < 2) {
            continue;
          }
          const double g1 = within_msd_subject(lhs);
          const double g2 = within_msd_subject(rhs);
          const double g3 = between_msd_subject(lhs, rhs);
          if (!R_finite(g1) || !R_finite(g2) || !R_finite(g3)) {
            continue;
          }
          sum_g1 += g1;
          sum_g2 += g2;
          sum_g3 += g3;
          sum_g1_sq += g1 * g1;
          sum_g2_sq += g2 * g2;
          sum_g3_sq += g3 * g3;
          sum_g1g2 += g1 * g2;
          sum_g1g3 += g1 * g3;
          sum_g2g3 += g2 * g3;
          n += 1;
        }

        n_eligible(j, k) = n;
        n_eligible(k, j) = n;
        if (n <= 0) {
          continue;
        }

        const double mean_g1 = sum_g1 / static_cast<double>(n);
        const double mean_g2 = sum_g2 / static_cast<double>(n);
        const double mean_g3 = sum_g3 / static_cast<double>(n);
        const double mean_a = (mean_g1 + mean_g2) / 2.0;
        const double est = (mean_g3 != 0.0) ? mean_a / mean_g3 : NA_REAL;
        const double s11 = sample_var_from_sums(sum_g1, sum_g1_sq, n);
        const double s22 = sample_var_from_sums(sum_g2, sum_g2_sq, n);
        const double s33 = sample_var_from_sums(sum_g3, sum_g3_sq, n);
        const double s12 = sample_cov_from_sums(sum_g1g2, sum_g1, sum_g2, n);
        const double s13 = sample_cov_from_sums(sum_g1g3, sum_g1, sum_g3, n);
        const double s23 = sample_cov_from_sums(sum_g2g3, sum_g2, sum_g3, n);

        estimate_raw(j, k) = est;
        estimate_raw(k, j) = est;
        A_bar(j, k) = mean_a;
        A_bar(k, j) = mean_a;
        B_bar(j, k) = mean_g3;
        B_bar(k, j) = mean_g3;
        S11(j, k) = s11;
        S11(k, j) = s11;
        S22(j, k) = s22;
        S22(k, j) = s22;
        S33(j, k) = s33;
        S33(k, j) = s33;
        S12(j, k) = s12;
        S12(k, j) = s12;
        S13(j, k) = s13;
        S13(k, j) = s13;
        S23(j, k) = s23;
        S23(k, j) = s23;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::_["estimate_raw"] = estimate_raw,
    Rcpp::_["A_bar"] = A_bar,
    Rcpp::_["B_bar"] = B_bar,
    Rcpp::_["S11"] = S11,
    Rcpp::_["S22"] = S22,
    Rcpp::_["S33"] = S33,
    Rcpp::_["S12"] = S12,
    Rcpp::_["S13"] = S13,
    Rcpp::_["S23"] = S23,
    Rcpp::_["n_eligible"] = n_eligible,
    Rcpp::_["n_subjects"] = n_subjects
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cia_pairwise_bootstrap_est_cpp(Rcpp::NumericVector y,
                                                   Rcpp::IntegerVector subject,
                                                   Rcpp::IntegerVector method,
                                                   Rcpp::IntegerVector replicate,
                                                   int n_methods,
                                                   int reference_method,
                                                   bool has_reference,
                                                   Rcpp::IntegerVector sampled_subjects,
                                                   std::string estimator,
                                                   int n_threads = 1) {
  const R_xlen_t n_obs = y.size();
  if (subject.size() != n_obs || method.size() != n_obs || replicate.size() != n_obs) {
    Rcpp::stop("y, subject, method, and replicate must have the same length.");
  }

  int n_subjects = 0;
  for (R_xlen_t i = 0; i < n_obs; ++i) {
    if (subject[i] != NA_INTEGER && subject[i] > n_subjects) {
      n_subjects = subject[i];
    }
  }
  std::vector< std::vector<R_xlen_t> > subject_rows(static_cast<std::size_t>(n_subjects));
  for (R_xlen_t i = 0; i < n_obs; ++i) {
    const int s = subject[i];
    if (s != NA_INTEGER && s >= 1 && s <= n_subjects) {
      subject_rows[static_cast<std::size_t>(s - 1)].push_back(i);
    }
  }

  R_xlen_t total = 0;
  for (R_xlen_t draw = 0; draw < sampled_subjects.size(); ++draw) {
    const int sampled = sampled_subjects[draw];
    if (sampled != NA_INTEGER && sampled >= 1 && sampled <= n_subjects) {
      total += static_cast<R_xlen_t>(subject_rows[static_cast<std::size_t>(sampled - 1)].size());
    }
  }

  Rcpp::NumericVector boot_y(total);
  Rcpp::IntegerVector boot_subject(total);
  Rcpp::IntegerVector boot_method(total);
  Rcpp::IntegerVector boot_replicate(total);
  R_xlen_t pos = 0;
  for (R_xlen_t draw = 0; draw < sampled_subjects.size(); ++draw) {
    const int sampled = sampled_subjects[draw];
    if (sampled == NA_INTEGER || sampled < 1 || sampled > n_subjects) {
      continue;
    }
    const std::vector<R_xlen_t>& rows = subject_rows[static_cast<std::size_t>(sampled - 1)];
    for (std::size_t row_idx = 0; row_idx < rows.size(); ++row_idx) {
      const R_xlen_t i = rows[row_idx];
      boot_y[pos] = y[i];
      boot_subject[pos] = static_cast<int>(draw + 1);
      boot_method[pos] = method[i];
      boot_replicate[pos] = replicate[i];
      pos += 1;
    }
  }

  Rcpp::List raw = cia_pairwise_stats_cpp(
    boot_y,
    boot_subject,
    boot_method,
    boot_replicate,
    n_methods,
    reference_method,
    has_reference,
    n_threads
  );
  Rcpp::NumericMatrix A_bar = raw["A_bar"];
  Rcpp::NumericMatrix B_bar = raw["B_bar"];
  Rcpp::IntegerMatrix n_eligible = raw["n_eligible"];

  Rcpp::NumericMatrix est(n_methods, n_methods);
  std::fill(est.begin(), est.end(), NA_REAL);
  for (int j = 0; j < n_methods; ++j) {
    est(j, j) = 1.0;
  }

  const bool constrained = estimator == "vc_constrained";
  for (int j = 0; j < n_methods - 1; ++j) {
    for (int k = j + 1; k < n_methods; ++k) {
      if (n_eligible(j, k) <= 0) {
        continue;
      }
      const double A = A_bar(j, k);
      const double B = B_bar(j, k);
      double value = NA_REAL;
      if (constrained) {
        const double tau2_raw = B - A;
        const double tau2 = R_finite(tau2_raw) ? std::max(tau2_raw, 0.0) : NA_REAL;
        const double denom = A + tau2;
        if (R_finite(denom) && denom != 0.0) {
          value = A / denom;
        }
      } else if (B != 0.0) {
        value = A / B;
      }
      est(j, k) = value;
      est(k, j) = value;
    }
  }

  return est;
}

// [[Rcpp::export]]
Rcpp::List cia_overall_balanced_cpp(Rcpp::NumericVector y,
                                    Rcpp::IntegerVector subject,
                                    Rcpp::IntegerVector method,
                                    Rcpp::IntegerVector replicate,
                                    int n_methods,
                                    int reference_method,
                                    bool has_reference,
                                    int n_threads = 1) {
  (void)replicate;
  (void)n_threads;

  int n_subjects = 0;
  const std::vector< std::vector< std::vector<double> > > by_subject =
    build_subject_method_values(y, subject, method, n_methods, &n_subjects);

  if (n_subjects <= 0) {
    Rcpp::stop(kOverallBalancedError);
  }

  int K = -1;
  Rcpp::NumericMatrix Ybar(n_subjects, n_methods);
  Rcpp::NumericMatrix A(n_subjects, n_methods);

  for (int s = 0; s < n_subjects; ++s) {
    for (int j = 0; j < n_methods; ++j) {
      const std::vector<double>& vals =
        by_subject[static_cast<std::size_t>(s)][static_cast<std::size_t>(j)];
      const int nv = static_cast<int>(vals.size());
      if (nv < 2) {
        Rcpp::stop(kOverallBalancedError);
      }
      if (K < 0) {
        K = nv;
      } else if (nv != K) {
        Rcpp::stop(kOverallBalancedError);
      }
      double mean_ij = 0.0;
      for (int k = 0; k < nv; ++k) {
        mean_ij += vals[static_cast<std::size_t>(k)];
      }
      mean_ij /= static_cast<double>(nv);
      Ybar(s, j) = mean_ij;
      double ss = 0.0;
      for (int k = 0; k < nv; ++k) {
        const double d = vals[static_cast<std::size_t>(k)] - mean_ij;
        ss += d * d;
      }
      A(s, j) = ss / static_cast<double>(K - 1);
    }
  }

  if (K < 2) {
    Rcpp::stop(kOverallBalancedError);
  }

  Rcpp::NumericVector A_i(n_subjects, NA_REAL);
  Rcpp::NumericVector B_i(n_subjects, NA_REAL);
  const double Jm1 = static_cast<double>(n_methods - 1);
  const double one_minus_invK = 1.0 - 1.0 / static_cast<double>(K);

  if (has_reference) {
    const int ref = reference_method - 1;
    if (ref < 0 || ref >= n_methods) {
      Rcpp::stop("reference_method is out of bounds.");
    }
    for (int s = 0; s < n_subjects; ++s) {
      A_i[s] = A(s, ref);
      double mean_term = 0.0;
      double within_nonref = 0.0;
      for (int j = 0; j < n_methods; ++j) {
        if (j == ref) {
          continue;
        }
        const double d = Ybar(s, j) - Ybar(s, ref);
        mean_term += d * d;
        within_nonref += A(s, j);
      }
      B_i[s] =
        mean_term / Jm1 +
        one_minus_invK * within_nonref / Jm1 +
        one_minus_invK * A(s, ref);
    }
  } else {
    for (int s = 0; s < n_subjects; ++s) {
      double A_mean = 0.0;
      double Y_mean = 0.0;
      for (int j = 0; j < n_methods; ++j) {
        A_mean += A(s, j);
        Y_mean += Ybar(s, j);
      }
      A_mean /= static_cast<double>(n_methods);
      Y_mean /= static_cast<double>(n_methods);
      A_i[s] = A_mean;
      double mean_term = 0.0;
      for (int j = 0; j < n_methods; ++j) {
        const double d = Ybar(s, j) - Y_mean;
        mean_term += d * d;
      }
      B_i[s] = mean_term / Jm1 + one_minus_invK * A_mean;
    }
  }

  double A_bar_val = 0.0;
  double B_bar_val = 0.0;
  for (int s = 0; s < n_subjects; ++s) {
    A_bar_val += A_i[s];
    B_bar_val += B_i[s];
  }
  A_bar_val /= static_cast<double>(n_subjects);
  B_bar_val /= static_cast<double>(n_subjects);

  const double estimate_raw =
    has_reference ? (2.0 * A_bar_val) / B_bar_val : A_bar_val / B_bar_val;

  return Rcpp::List::create(
    Rcpp::_["estimate_raw"] = estimate_raw,
    Rcpp::_["A_i"] = A_i,
    Rcpp::_["B_i"] = B_i,
    Rcpp::_["A_bar"] = A_bar_val,
    Rcpp::_["B_bar"] = B_bar_val,
    Rcpp::_["Ybar"] = Ybar,
    Rcpp::_["A"] = A,
    Rcpp::_["n_subjects"] = n_subjects,
    Rcpp::_["n_methods"] = n_methods,
    Rcpp::_["K"] = K,
    Rcpp::_["ref_idx"] = has_reference ? reference_method : 0
  );
}
