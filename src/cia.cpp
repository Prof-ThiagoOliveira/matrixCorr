// Thiago de Paula Oliveira
// cia.cpp
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {

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
  (void) replicate;
  (void) n_threads;

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
    if (s < 1 || m < 1 || m > n_methods || s > n_subjects) {
      continue;
    }
    by_subject[static_cast<std::size_t>(s - 1)][static_cast<std::size_t>(m - 1)].push_back(y[i]);
  }

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

  Rcpp::RObject estimate;
  if (pairwise) {
    Rcpp::NumericMatrix est(n_methods, n_methods);
    std::fill(est.begin(), est.end(), NA_REAL);
    for (int j = 0; j < n_methods; ++j) {
      est(j, j) = 1.0;
    }

    if (has_reference) {
      const int ref = reference_method - 1;
      if (ref >= 0 && ref < n_methods) {
        for (int j = 0; j < n_methods; ++j) {
          if (j == ref) {
            continue;
          }
          double g1_sum = 0.0;
          double g3_sum = 0.0;
          int eligible_n = 0;
          for (int s = 0; s < n_subjects; ++s) {
            const std::vector<double>& ref_vals = by_subject[static_cast<std::size_t>(s)][static_cast<std::size_t>(ref)];
            const std::vector<double>& cmp_vals = by_subject[static_cast<std::size_t>(s)][static_cast<std::size_t>(j)];
            if (ref_vals.size() < 2 || cmp_vals.empty()) {
              continue;
            }
            const double g1 = within_msd_subject(ref_vals);
            const double g3 = between_msd_subject(ref_vals, cmp_vals);
            if (!R_finite(g1) || !R_finite(g3)) {
              continue;
            }
            g1_sum += g1;
            g3_sum += g3;
            eligible_n += 1;
          }
          if (eligible_n > 0) {
            const double denom = g3_sum / static_cast<double>(eligible_n);
            if (denom != 0.0) {
              est(j, ref) = (g1_sum / static_cast<double>(eligible_n)) / denom;
            }
            est(ref, j) = est(j, ref);
          }
        }
      }
    } else {
      for (int j = 0; j < n_methods - 1; ++j) {
        for (int k = j + 1; k < n_methods; ++k) {
          double g1_sum = 0.0;
          double g2_sum = 0.0;
          double g3_sum = 0.0;
          int eligible_n = 0;
          for (int s = 0; s < n_subjects; ++s) {
            const std::vector<double>& lhs = by_subject[static_cast<std::size_t>(s)][static_cast<std::size_t>(j)];
            const std::vector<double>& rhs = by_subject[static_cast<std::size_t>(s)][static_cast<std::size_t>(k)];
            if (lhs.size() < 2 || rhs.size() < 2) {
              continue;
            }
            const double g1 = within_msd_subject(lhs);
            const double g2 = within_msd_subject(rhs);
            const double g3 = between_msd_subject(lhs, rhs);
            if (!R_finite(g1) || !R_finite(g2) || !R_finite(g3)) {
              continue;
            }
            g1_sum += g1;
            g2_sum += g2;
            g3_sum += g3;
            eligible_n += 1;
          }
          if (eligible_n <= 0) {
            continue;
          }
          const double denom = g3_sum / static_cast<double>(eligible_n);
          if (denom == 0.0) {
            continue;
          }
          const double numer =
            ((g1_sum / static_cast<double>(eligible_n)) +
             (g2_sum / static_cast<double>(eligible_n))) / 2.0;
          est(j, k) = numer / denom;
          est(k, j) = est(j, k);
        }
      }
    }
    estimate = est;
  } else {
    double est = NA_REAL;
    if (has_reference) {
      const int ref = reference_method - 1;
      if (ref >= 0 && ref < n_methods && R_finite(sigma_within[ref])) {
        double denom_sum = 0.0;
        int denom_n = 0;
        for (int j = 0; j < n_methods; ++j) {
          if (j == ref || !R_finite(between_msd(j, ref))) {
            continue;
          }
          denom_sum += between_msd(j, ref);
          denom_n += 1;
        }
        if (denom_n > 0) {
          const double denom = denom_sum / static_cast<double>(denom_n);
          if (denom != 0.0) {
            est = sigma_within[ref] / denom;
          }
        }
      }
    } else if (n_methods > 1) {
      double within_term = 0.0;
      bool valid_within = true;
      for (int j = 0; j < n_methods; ++j) {
        if (!R_finite(within_msd[j])) {
          valid_within = false;
          break;
        }
        within_term += within_msd[j] / 2.0;
      }

      double between_sum_all = 0.0;
      bool valid_between = true;
      for (int j = 0; j < n_methods - 1; ++j) {
        for (int k = j + 1; k < n_methods; ++k) {
          if (!R_finite(between_msd(j, k))) {
            valid_between = false;
            break;
          }
          between_sum_all += between_msd(j, k);
        }
        if (!valid_between) {
          break;
        }
      }

      if (valid_within && valid_between) {
        const double denom = between_sum_all / static_cast<double>(n_methods - 1);
        if (denom != 0.0) {
          est = within_term / denom;
        }
      }
    }
    estimate = Rcpp::wrap(est);
  }

  return Rcpp::List::create(
    Rcpp::_["estimate"] = estimate,
    Rcpp::_["within_msd"] = within_msd,
    Rcpp::_["between_msd"] = between_msd,
    Rcpp::_["sigma_within"] = sigma_within,
    Rcpp::_["n_subjects"] = n_subjects,
    Rcpp::_["n_obs"] = static_cast<int>(n_obs),
    Rcpp::_["n_replicate_pairs"] = replicate_pairs,
    Rcpp::_["n_between_pairs"] = between_pairs
  );
}
