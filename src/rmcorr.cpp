// Thiago de Paula Oliveira
// repeated-measures correlation kernels
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

#include "matrixCorr_omp.h"

using namespace Rcpp;

namespace {

struct RmStats {
  double estimate{NA_REAL};
  double slope{NA_REAL};
  double p_value{NA_REAL};
  double df{NA_REAL};
  double conf_low{NA_REAL};
  double conf_high{NA_REAL};
  double weighted_ss_x{NA_REAL};
  double weighted_ss_e{NA_REAL};
  int n_complete{0};
  int n_subjects{0};
  int obs_per_subject_min{NA_INTEGER};
  int obs_per_subject_max{NA_INTEGER};
  bool valid{false};
};

struct SubjectIndex {
  std::vector<int> row_index;
  std::vector<int> start;
};

struct CompleteSubjectIndex {
  std::vector<int> row_index;
  std::vector<int> start;
  int n_complete{0};
  int n_subjects{0};
  double df{NA_REAL};
};

inline double clamp_corr(double x) {
  if (!std::isfinite(x)) return NA_REAL;
  if (x > 1.0) return 1.0;
  if (x < -1.0) return -1.0;
  return x;
}

inline bool all_finite_ptr(const double* x, const R_xlen_t n) {
  for (R_xlen_t i = 0; i < n; ++i) {
    if (!std::isfinite(x[i])) return false;
  }
  return true;
}

inline SubjectIndex build_subject_index(const IntegerVector& subject) {
  SubjectIndex out;
  int max_id = 0;
  for (int i = 0; i < subject.size(); ++i) {
    const int s = subject[i];
    if (s != NA_INTEGER && s > max_id) max_id = s;
  }

  if (max_id <= 0) return out;

  std::vector<int> counts(static_cast<std::size_t>(max_id), 0);
  for (int i = 0; i < subject.size(); ++i) {
    const int s = subject[i];
    if (s != NA_INTEGER && s >= 1) {
      counts[static_cast<std::size_t>(s - 1)] += 1;
    }
  }

  out.start.resize(static_cast<std::size_t>(max_id) + 1u, 0);
  for (int s = 0; s < max_id; ++s) {
    out.start[static_cast<std::size_t>(s) + 1u] =
      out.start[static_cast<std::size_t>(s)] + counts[static_cast<std::size_t>(s)];
  }

  out.row_index.resize(static_cast<std::size_t>(out.start.back()));
  std::vector<int> cursor(out.start.begin(), out.start.end() - 1);
  for (int i = 0; i < subject.size(); ++i) {
    const int s = subject[i];
    if (s == NA_INTEGER || s < 1) continue;
    const std::size_t pos = static_cast<std::size_t>(cursor[static_cast<std::size_t>(s - 1)]++);
    out.row_index[pos] = i;
  }

  return out;
}

inline CompleteSubjectIndex build_complete_subject_index(const IntegerVector& subject) {
  CompleteSubjectIndex out;
  int max_id = 0;
  for (int i = 0; i < subject.size(); ++i) {
    const int s = subject[i];
    if (s != NA_INTEGER && s > max_id) max_id = s;
  }

  if (max_id <= 0) return out;

  std::vector<int> counts(static_cast<std::size_t>(max_id), 0);
  for (int i = 0; i < subject.size(); ++i) {
    const int s = subject[i];
    if (s != NA_INTEGER && s >= 1) {
      counts[static_cast<std::size_t>(s - 1)] += 1;
    }
  }

  out.start.reserve(static_cast<std::size_t>(max_id) + 1u);
  out.start.push_back(0);
  for (int s = 0; s < max_id; ++s) {
    const int count = counts[static_cast<std::size_t>(s)];
    if (count >= 2) {
      out.n_subjects += 1;
      out.n_complete += count;
      out.start.push_back(out.n_complete);
    }
  }

  if (out.n_subjects == 0) {
    out.start.resize(1u, 0);
    return out;
  }

  out.df = static_cast<double>(out.n_complete - out.n_subjects - 1);
  out.row_index.resize(static_cast<std::size_t>(out.n_complete));

  std::vector<int> valid_pos(static_cast<std::size_t>(max_id), -1);
  int idx = 0;
  for (int s = 0; s < max_id; ++s) {
    if (counts[static_cast<std::size_t>(s)] >= 2) {
      valid_pos[static_cast<std::size_t>(s)] = idx++;
    }
  }

  std::vector<int> cursor(out.start.begin(), out.start.end() - 1);
  for (int i = 0; i < subject.size(); ++i) {
    const int s = subject[i];
    if (s == NA_INTEGER || s < 1) continue;
    const int pos = valid_pos[static_cast<std::size_t>(s - 1)];
    if (pos < 0) continue;
    out.row_index[static_cast<std::size_t>(cursor[static_cast<std::size_t>(pos)]++)] = i;
  }

  return out;
}

inline RmStats compute_rmcorr_stats(
  const double* x,
  const double* y,
  const SubjectIndex& groups,
  const bool use_confint,
  const double zcrit
) {
  RmStats out;
  const std::size_t m = groups.start.empty() ? 0u : groups.start.size() - 1u;
  if (m == 0) return out;

  static thread_local std::vector<int> count_buf;
  static thread_local std::vector<double> sumx_buf;
  static thread_local std::vector<double> sumy_buf;

  if (count_buf.size() < m) count_buf.resize(m);
  if (sumx_buf.size() < m) sumx_buf.resize(m);
  if (sumy_buf.size() < m) sumy_buf.resize(m);
  std::fill(count_buf.begin(), count_buf.begin() + static_cast<std::ptrdiff_t>(m), 0);
  std::fill(sumx_buf.begin(), sumx_buf.begin() + static_cast<std::ptrdiff_t>(m), 0.0);
  std::fill(sumy_buf.begin(), sumy_buf.begin() + static_cast<std::ptrdiff_t>(m), 0.0);

  for (std::size_t s = 0; s < m; ++s) {
    const int begin = groups.start[s];
    const int end = groups.start[s + 1u];
    for (int pos = begin; pos < end; ++pos) {
      const int row = groups.row_index[static_cast<std::size_t>(pos)];
      const double xv = x[row];
      const double yv = y[row];
      if (!std::isfinite(xv) || !std::isfinite(yv)) continue;
      count_buf[s] += 1;
      sumx_buf[s] += xv;
      sumy_buf[s] += yv;
    }
    if (count_buf[s] >= 2) {
      out.n_complete += count_buf[s];
      out.n_subjects += 1;
      if (out.obs_per_subject_min == NA_INTEGER || count_buf[s] < out.obs_per_subject_min) {
        out.obs_per_subject_min = count_buf[s];
      }
      if (out.obs_per_subject_max == NA_INTEGER || count_buf[s] > out.obs_per_subject_max) {
        out.obs_per_subject_max = count_buf[s];
      }
    }
  }

  if (out.n_subjects < 2) return out;

  out.df = static_cast<double>(out.n_complete - out.n_subjects - 1);
  if (!(out.df > 0.0)) return out;

  double sxx = 0.0;
  double syy = 0.0;
  double sxy = 0.0;

  for (std::size_t s = 0; s < m; ++s) {
    const int count = count_buf[s];
    if (count < 2) continue;

    const double mean_x = sumx_buf[s] / static_cast<double>(count);
    const double mean_y = sumy_buf[s] / static_cast<double>(count);
    const int begin = groups.start[s];
    const int end = groups.start[s + 1u];
    for (int pos = begin; pos < end; ++pos) {
      const int row = groups.row_index[static_cast<std::size_t>(pos)];
      const double xv = x[row];
      const double yv = y[row];
      if (!std::isfinite(xv) || !std::isfinite(yv)) continue;
      const double dx = xv - mean_x;
      const double dy = yv - mean_y;
      sxx += dx * dx;
      syy += dy * dy;
      sxy += dx * dy;
    }
  }

  if (!(sxx > 0.0) || !(syy > 0.0)) return out;

  out.estimate = clamp_corr(sxy / std::sqrt(sxx * syy));
  out.slope = sxy / sxx;
  out.valid = std::isfinite(out.estimate) && std::isfinite(out.slope);
  if (!out.valid) return out;

  double sse = syy - (sxy * sxy) / sxx;
  if (sse < 0.0 && sse > -1e-12) sse = 0.0;
  if (sse < 0.0) return out;

  const double mse = sse / out.df;
  if (!(mse >= 0.0) || !std::isfinite(mse)) return out;

  const double se_slope = std::sqrt(mse / sxx);
  double t_value = NA_REAL;
  if (se_slope == 0.0) {
    t_value = (out.slope == 0.0) ? 0.0 : R_PosInf;
  } else {
    t_value = out.slope / se_slope;
  }
  if (std::isfinite(t_value) || t_value == R_PosInf || t_value == R_NegInf) {
    out.p_value = 2.0 * R::pt(-std::abs(t_value), out.df, 1, 0);
  }

  if (use_confint) {
    const double eff_n = out.df + 2.0;
    if (std::abs(out.estimate) >= 1.0) {
      out.conf_low = out.estimate;
      out.conf_high = out.estimate;
    } else if (eff_n > 3.0) {
      const double se_z = 1.0 / std::sqrt(eff_n - 3.0);
      if (std::isfinite(zcrit) && std::isfinite(se_z) && se_z > 0.0) {
        const double zr = std::atanh(out.estimate);
        out.conf_low = clamp_corr(std::tanh(zr - zcrit * se_z));
        out.conf_high = clamp_corr(std::tanh(zr + zcrit * se_z));
      }
    }
  }

  return out;
}

inline RmStats compute_rmcorr_stats_complete(
  const double* x,
  const double* y,
  const CompleteSubjectIndex& groups,
  const bool use_confint,
  const double zcrit
) {
  RmStats out;
  const std::size_t m = groups.start.empty() ? 0u : groups.start.size() - 1u;
  if (m == 0) return out;

  out.n_complete = groups.n_complete;
  out.n_subjects = groups.n_subjects;
  out.df = groups.df;
  if (out.n_subjects < 2 || !(out.df > 0.0)) return out;

  static thread_local std::vector<double> sumx_buf;
  static thread_local std::vector<double> sumy_buf;

  if (sumx_buf.size() < m) sumx_buf.resize(m);
  if (sumy_buf.size() < m) sumy_buf.resize(m);
  std::fill(sumx_buf.begin(), sumx_buf.begin() + static_cast<std::ptrdiff_t>(m), 0.0);
  std::fill(sumy_buf.begin(), sumy_buf.begin() + static_cast<std::ptrdiff_t>(m), 0.0);

  for (std::size_t s = 0; s < m; ++s) {
    const int begin = groups.start[s];
    const int end = groups.start[s + 1u];
    for (int pos = begin; pos < end; ++pos) {
      const int row = groups.row_index[static_cast<std::size_t>(pos)];
      sumx_buf[s] += x[row];
      sumy_buf[s] += y[row];
    }
  }

  double sxx = 0.0;
  double syy = 0.0;
  double sxy = 0.0;

  for (std::size_t s = 0; s < m; ++s) {
    const int begin = groups.start[s];
    const int end = groups.start[s + 1u];
    const int count = end - begin;
    if (out.obs_per_subject_min == NA_INTEGER || count < out.obs_per_subject_min) {
      out.obs_per_subject_min = count;
    }
    if (out.obs_per_subject_max == NA_INTEGER || count > out.obs_per_subject_max) {
      out.obs_per_subject_max = count;
    }
    const double mean_x = sumx_buf[s] / static_cast<double>(count);
    const double mean_y = sumy_buf[s] / static_cast<double>(count);
    for (int pos = begin; pos < end; ++pos) {
      const int row = groups.row_index[static_cast<std::size_t>(pos)];
      const double dx = x[row] - mean_x;
      const double dy = y[row] - mean_y;
      sxx += dx * dx;
      syy += dy * dy;
      sxy += dx * dy;
    }
  }

  if (!(sxx > 0.0) || !(syy > 0.0)) return out;

  out.estimate = clamp_corr(sxy / std::sqrt(sxx * syy));
  out.slope = sxy / sxx;
  out.valid = std::isfinite(out.estimate) && std::isfinite(out.slope);
  if (!out.valid) return out;

  double sse = syy - (sxy * sxy) / sxx;
  if (sse < 0.0 && sse > -1e-12) sse = 0.0;
  if (sse < 0.0) return out;

  const double mse = sse / out.df;
  if (!(mse >= 0.0) || !std::isfinite(mse)) return out;

  const double se_slope = std::sqrt(mse / sxx);
  double t_value = NA_REAL;
  if (se_slope == 0.0) {
    t_value = (out.slope == 0.0) ? 0.0 : R_PosInf;
  } else {
    t_value = out.slope / se_slope;
  }
  if (std::isfinite(t_value) || t_value == R_PosInf || t_value == R_NegInf) {
    out.p_value = 2.0 * R::pt(-std::abs(t_value), out.df, 1, 0);
  }

  if (use_confint) {
    const double eff_n = out.df + 2.0;
    if (std::abs(out.estimate) >= 1.0) {
      out.conf_low = out.estimate;
      out.conf_high = out.estimate;
    } else if (eff_n > 3.0) {
      const double se_z = 1.0 / std::sqrt(eff_n - 3.0);
      if (std::isfinite(zcrit) && std::isfinite(se_z) && se_z > 0.0) {
        const double zr = std::atanh(out.estimate);
        out.conf_low = clamp_corr(std::tanh(zr - zcrit * se_z));
        out.conf_high = clamp_corr(std::tanh(zr + zcrit * se_z));
      }
    }
  }

  return out;
}

inline RmStats compute_rmcorr_weighted_stats(
  const double* x,
  const double* y,
  const SubjectIndex& groups,
  const bool use_confint,
  const double zcrit
) {
  RmStats out;
  const std::size_t m = groups.start.empty() ? 0u : groups.start.size() - 1u;
  if (m == 0) return out;

  static thread_local std::vector<int> count_buf;
  static thread_local std::vector<double> sumx_buf;
  static thread_local std::vector<double> sumy_buf;

  if (count_buf.size() < m) count_buf.resize(m);
  if (sumx_buf.size() < m) sumx_buf.resize(m);
  if (sumy_buf.size() < m) sumy_buf.resize(m);
  std::fill(count_buf.begin(), count_buf.begin() + static_cast<std::ptrdiff_t>(m), 0);
  std::fill(sumx_buf.begin(), sumx_buf.begin() + static_cast<std::ptrdiff_t>(m), 0.0);
  std::fill(sumy_buf.begin(), sumy_buf.begin() + static_cast<std::ptrdiff_t>(m), 0.0);

  for (std::size_t s = 0; s < m; ++s) {
    const int begin = groups.start[s];
    const int end = groups.start[s + 1u];
    for (int pos = begin; pos < end; ++pos) {
      const int row = groups.row_index[static_cast<std::size_t>(pos)];
      const double xv = x[row];
      const double yv = y[row];
      if (!std::isfinite(xv) || !std::isfinite(yv)) continue;
      count_buf[s] += 1;
      sumx_buf[s] += xv;
      sumy_buf[s] += yv;
    }
    if (count_buf[s] >= 2) {
      out.n_complete += count_buf[s];
      out.n_subjects += 1;
      if (out.obs_per_subject_min == NA_INTEGER || count_buf[s] < out.obs_per_subject_min) {
        out.obs_per_subject_min = count_buf[s];
      }
      if (out.obs_per_subject_max == NA_INTEGER || count_buf[s] > out.obs_per_subject_max) {
        out.obs_per_subject_max = count_buf[s];
      }
    }
  }

  if (out.n_subjects < 2) return out;

  out.df = static_cast<double>(out.n_complete - out.n_subjects - 1);
  if (!(out.df > 0.0)) return out;

  double sxx = 0.0;
  double syy = 0.0;
  double sxy = 0.0;

  for (std::size_t s = 0; s < m; ++s) {
    const int count = count_buf[s];
    if (count < 2) continue;

    const double mean_x = sumx_buf[s] / static_cast<double>(count);
    const double mean_y = sumy_buf[s] / static_cast<double>(count);
    const int begin = groups.start[s];
    const int end = groups.start[s + 1u];
    for (int pos = begin; pos < end; ++pos) {
      const int row = groups.row_index[static_cast<std::size_t>(pos)];
      const double xv = x[row];
      const double yv = y[row];
      if (!std::isfinite(xv) || !std::isfinite(yv)) continue;
      const double dx = xv - mean_x;
      const double dy = yv - mean_y;
      sxx += dx * dx;
      syy += dy * dy;
      sxy += dx * dy;
    }
  }

  if (!(sxx > 0.0)) return out;

  out.slope = sxy / sxx;
  if (!std::isfinite(out.slope)) return out;

  double weighted_ss_x = 0.0;
  double weighted_ss_e = 0.0;
  double sse = 0.0;

  for (std::size_t s = 0; s < m; ++s) {
    const int count = count_buf[s];
    if (count < 2) continue;

    const double mean_x = sumx_buf[s] / static_cast<double>(count);
    const double mean_y = sumy_buf[s] / static_cast<double>(count);
    const double inv_count = 1.0 / static_cast<double>(count);
    const int begin = groups.start[s];
    const int end = groups.start[s + 1u];
    double sum_e2 = 0.0;
    double sum_y2 = 0.0;

    for (int pos = begin; pos < end; ++pos) {
      const int row = groups.row_index[static_cast<std::size_t>(pos)];
      const double xv = x[row];
      const double yv = y[row];
      if (!std::isfinite(xv) || !std::isfinite(yv)) continue;
      const double dx = xv - mean_x;
      const double dy = yv - mean_y;
      const double e = dy - out.slope * dx;
      const double e2 = e * e;
      sum_e2 += e2;
      sum_y2 += dy * dy;
    }

    sse += sum_e2;
    weighted_ss_e += inv_count * sum_e2;
    weighted_ss_x += inv_count * (sum_y2 - sum_e2);
  }

  const double neg_tol = std::sqrt(std::numeric_limits<double>::epsilon());
  if (weighted_ss_x < 0.0 && weighted_ss_x > -neg_tol) {
    weighted_ss_x = 0.0;
  }
  if (weighted_ss_x < 0.0) return out;

  const double denom = weighted_ss_x + weighted_ss_e;
  if (!(denom > 0.0) || !std::isfinite(denom)) return out;

  const double mag = std::sqrt(weighted_ss_x / denom);
  if (!std::isfinite(mag)) return out;
  const double sign_beta = (out.slope > 0.0) ? 1.0 : ((out.slope < 0.0) ? -1.0 : 0.0);
  out.estimate = clamp_corr(sign_beta * mag);
  if (!std::isfinite(out.estimate)) return out;

  if (sse < 0.0 && sse > -1e-12) sse = 0.0;
  if (sse < 0.0) return out;
  const double mse = sse / out.df;
  if (!(mse >= 0.0) || !std::isfinite(mse)) return out;

  const double se_slope = std::sqrt(mse / sxx);
  double t_value = NA_REAL;
  if (se_slope == 0.0) {
    t_value = (out.slope == 0.0) ? 0.0 : R_PosInf;
  } else {
    t_value = out.slope / se_slope;
  }
  if (std::isfinite(t_value) || t_value == R_PosInf || t_value == R_NegInf) {
    out.p_value = 2.0 * R::pt(-std::abs(t_value), out.df, 1, 0);
  }

  if (use_confint) {
    const double eff_n = out.df + 2.0;
    if (std::abs(out.estimate) >= 1.0) {
      out.conf_low = out.estimate;
      out.conf_high = out.estimate;
    } else if (eff_n > 3.0) {
      const double se_z = 1.0 / std::sqrt(eff_n - 3.0);
      if (std::isfinite(zcrit) && std::isfinite(se_z) && se_z > 0.0) {
        const double zr = std::atanh(out.estimate);
        out.conf_low = clamp_corr(std::tanh(zr - zcrit * se_z));
        out.conf_high = clamp_corr(std::tanh(zr + zcrit * se_z));
      }
    }
  }

  out.weighted_ss_x = weighted_ss_x;
  out.weighted_ss_e = weighted_ss_e;
  out.valid = true;
  return out;
}

inline List stats_to_list(const RmStats& stats) {
  return List::create(
    _["estimate"] = stats.valid ? stats.estimate : NA_REAL,
    _["slope"] = stats.valid ? stats.slope : NA_REAL,
    _["p_value"] = stats.valid ? stats.p_value : NA_REAL,
    _["df"] = (stats.n_subjects >= 2) ? stats.df : NA_REAL,
    _["conf_low"] = stats.valid ? stats.conf_low : NA_REAL,
    _["conf_high"] = stats.valid ? stats.conf_high : NA_REAL,
    _["weighted_ss_x"] = stats.valid ? stats.weighted_ss_x : NA_REAL,
    _["weighted_ss_e"] = stats.valid ? stats.weighted_ss_e : NA_REAL,
    _["n_complete"] = stats.n_complete,
    _["n_subjects"] = stats.n_subjects,
    _["obs_per_subject_min"] = (stats.n_subjects >= 1) ? stats.obs_per_subject_min : NA_INTEGER,
    _["obs_per_subject_max"] = (stats.n_subjects >= 1) ? stats.obs_per_subject_max : NA_INTEGER,
    _["valid"] = stats.valid
  );
}

} // namespace

// [[Rcpp::export]]
List rmcorr_pair_cpp(NumericVector x, NumericVector y, IntegerVector subject,
                     double conf_level = 0.95) {
  if (x.size() != y.size() || x.size() != subject.size()) {
    stop("x, y, and subject must have the same length.");
  }
  const bool use_confint = conf_level > 0.0 && conf_level < 1.0;
  const double zcrit = use_confint ?
    R::qnorm(0.5 * (1.0 + conf_level), 0.0, 1.0, 1, 0) : NA_REAL;
  const double* x_ptr = REAL(x);
  const double* y_ptr = REAL(y);
  const bool complete = all_finite_ptr(x_ptr, x.size()) && all_finite_ptr(y_ptr, y.size());
  const RmStats stats = complete ?
    compute_rmcorr_stats_complete(
      x_ptr, y_ptr, build_complete_subject_index(subject), use_confint, zcrit
    ) :
    compute_rmcorr_stats(
      x_ptr, y_ptr, build_subject_index(subject), use_confint, zcrit
    );
  return stats_to_list(stats);
}

// [[Rcpp::export]]
List rmcorr_pair_weighted_cpp(NumericVector x, NumericVector y, IntegerVector subject,
                              double conf_level = 0.95) {
  if (x.size() != y.size() || x.size() != subject.size()) {
    stop("x, y, and subject must have the same length.");
  }
  const bool use_confint = conf_level > 0.0 && conf_level < 1.0;
  const double zcrit = use_confint ?
    R::qnorm(0.5 * (1.0 + conf_level), 0.0, 1.0, 1, 0) : NA_REAL;
  const double* x_ptr = REAL(x);
  const double* y_ptr = REAL(y);
  const SubjectIndex groups = build_subject_index(subject);
  const RmStats stats = compute_rmcorr_weighted_stats(
    x_ptr, y_ptr, groups, use_confint, zcrit
  );
  return stats_to_list(stats);
}

// [[Rcpp::export]]
List rmcorr_matrix_cpp(NumericMatrix x, NumericMatrix y, IntegerVector subject,
                       bool symmetric = false, double conf_level = 0.95,
                       int n_threads = 1) {
  const int n = x.nrow();
  if (y.nrow() != n || subject.size() != n) {
    stop("x, y, and subject must have the same number of rows.");
  }

  const int px = x.ncol();
  const int py = y.ncol();
  if (px < 1 || py < 1) {
    stop("x and y must each have at least one column.");
  }
  if (symmetric && px != py) {
    stop("symmetric = TRUE requires x and y to have the same number of columns.");
  }

  NumericMatrix estimate(px, py);
  NumericMatrix slope(px, py);
  NumericMatrix p_value(px, py);
  NumericMatrix df(px, py);
  NumericMatrix conf_low(px, py);
  NumericMatrix conf_high(px, py);
  IntegerMatrix n_complete(px, py);
  IntegerMatrix n_subjects(px, py);

  const bool use_confint = conf_level > 0.0 && conf_level < 1.0;
  const double zcrit = use_confint ?
    R::qnorm(0.5 * (1.0 + conf_level), 0.0, 1.0, 1, 0) : NA_REAL;
  const double* x_ptr = REAL(x);
  const double* y_ptr = REAL(y);
  const bool complete =
    all_finite_ptr(x_ptr, static_cast<R_xlen_t>(n) * px) &&
    all_finite_ptr(y_ptr, static_cast<R_xlen_t>(n) * py);
  const SubjectIndex groups = complete ? SubjectIndex() : build_subject_index(subject);
  const CompleteSubjectIndex groups_complete = complete ? build_complete_subject_index(subject) : CompleteSubjectIndex();

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
  for (int j = 0; j < px; ++j) {
    const int k_start = symmetric ? j : 0;
    const double* col_x = x_ptr + static_cast<R_xlen_t>(j) * n;
    for (int k = k_start; k < py; ++k) {
      const double* col_y = y_ptr + static_cast<R_xlen_t>(k) * n;
      const RmStats stats = complete ?
        compute_rmcorr_stats_complete(col_x, col_y, groups_complete, use_confint, zcrit) :
        compute_rmcorr_stats(col_x, col_y, groups, use_confint, zcrit);

      estimate(j, k) = stats.valid ? stats.estimate : NA_REAL;
      slope(j, k) = stats.valid ? stats.slope : NA_REAL;
      p_value(j, k) = stats.valid ? stats.p_value : NA_REAL;
      df(j, k) = (stats.n_subjects >= 2) ? stats.df : NA_REAL;
      conf_low(j, k) = stats.valid ? stats.conf_low : NA_REAL;
      conf_high(j, k) = stats.valid ? stats.conf_high : NA_REAL;
      n_complete(j, k) = stats.n_complete;
      n_subjects(j, k) = stats.n_subjects;

      if (symmetric && k != j) {
        estimate(k, j) = estimate(j, k);
        slope(k, j) = slope(j, k);
        p_value(k, j) = p_value(j, k);
        df(k, j) = df(j, k);
        conf_low(k, j) = conf_low(j, k);
        conf_high(k, j) = conf_high(j, k);
        n_complete(k, j) = n_complete(j, k);
        n_subjects(k, j) = n_subjects(j, k);
      }
    }
  }

  if (symmetric) {
    for (int j = 0; j < px; ++j) {
      if (std::isfinite(estimate(j, j))) {
        estimate(j, j) = 1.0;
      }
    }
  }

  return List::create(
    _["estimate"] = estimate,
    _["slope"] = slope,
    _["p_value"] = p_value,
    _["df"] = df,
    _["conf_low"] = conf_low,
    _["conf_high"] = conf_high,
    _["n_complete"] = n_complete,
    _["n_subjects"] = n_subjects
  );
}

// [[Rcpp::export]]
List rmcorr_matrix_weighted_cpp(NumericMatrix x, NumericMatrix y, IntegerVector subject,
                                bool symmetric = false, double conf_level = 0.95,
                                int n_threads = 1) {
  const int n = x.nrow();
  if (y.nrow() != n || subject.size() != n) {
    stop("x, y, and subject must have the same number of rows.");
  }

  const int px = x.ncol();
  const int py = y.ncol();
  if (px < 1 || py < 1) {
    stop("x and y must each have at least one column.");
  }
  if (symmetric && px != py) {
    stop("symmetric = TRUE requires x and y to have the same number of columns.");
  }

  NumericMatrix estimate(px, py);
  NumericMatrix slope(px, py);
  NumericMatrix p_value(px, py);
  NumericMatrix df(px, py);
  NumericMatrix conf_low(px, py);
  NumericMatrix conf_high(px, py);
  NumericMatrix weighted_ss_x(px, py);
  NumericMatrix weighted_ss_e(px, py);
  IntegerMatrix n_complete(px, py);
  IntegerMatrix n_subjects(px, py);
  IntegerMatrix obs_per_subject_min(px, py);
  IntegerMatrix obs_per_subject_max(px, py);
  LogicalMatrix valid(px, py);

  const bool use_confint = conf_level > 0.0 && conf_level < 1.0;
  const double zcrit = use_confint ?
    R::qnorm(0.5 * (1.0 + conf_level), 0.0, 1.0, 1, 0) : NA_REAL;
  const double* x_ptr = REAL(x);
  const double* y_ptr = REAL(y);
  const SubjectIndex groups = build_subject_index(subject);

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
  for (int j = 0; j < px; ++j) {
    const int k_start = symmetric ? j : 0;
    const double* col_x = x_ptr + static_cast<R_xlen_t>(j) * n;
    for (int k = k_start; k < py; ++k) {
      const double* col_y = y_ptr + static_cast<R_xlen_t>(k) * n;
      const RmStats stats = compute_rmcorr_weighted_stats(col_x, col_y, groups, use_confint, zcrit);

      estimate(j, k) = stats.valid ? stats.estimate : NA_REAL;
      slope(j, k) = stats.valid ? stats.slope : NA_REAL;
      p_value(j, k) = stats.valid ? stats.p_value : NA_REAL;
      df(j, k) = (stats.n_subjects >= 2) ? stats.df : NA_REAL;
      conf_low(j, k) = stats.valid ? stats.conf_low : NA_REAL;
      conf_high(j, k) = stats.valid ? stats.conf_high : NA_REAL;
      weighted_ss_x(j, k) = stats.valid ? stats.weighted_ss_x : NA_REAL;
      weighted_ss_e(j, k) = stats.valid ? stats.weighted_ss_e : NA_REAL;
      n_complete(j, k) = stats.n_complete;
      n_subjects(j, k) = stats.n_subjects;
      obs_per_subject_min(j, k) = (stats.n_subjects >= 1) ? stats.obs_per_subject_min : NA_INTEGER;
      obs_per_subject_max(j, k) = (stats.n_subjects >= 1) ? stats.obs_per_subject_max : NA_INTEGER;
      valid(j, k) = stats.valid;

      if (symmetric && k != j) {
        estimate(k, j) = estimate(j, k);
        slope(k, j) = slope(j, k);
        p_value(k, j) = p_value(j, k);
        df(k, j) = df(j, k);
        conf_low(k, j) = conf_low(j, k);
        conf_high(k, j) = conf_high(j, k);
        weighted_ss_x(k, j) = weighted_ss_x(j, k);
        weighted_ss_e(k, j) = weighted_ss_e(j, k);
        n_complete(k, j) = n_complete(j, k);
        n_subjects(k, j) = n_subjects(j, k);
        obs_per_subject_min(k, j) = obs_per_subject_min(j, k);
        obs_per_subject_max(k, j) = obs_per_subject_max(j, k);
        valid(k, j) = valid(j, k);
      }
    }
  }

  if (symmetric) {
    for (int j = 0; j < px; ++j) {
      if (valid(j, j) == TRUE && std::isfinite(estimate(j, j))) {
        estimate(j, j) = 1.0;
      }
    }
  }

  return List::create(
    _["estimate"] = estimate,
    _["slope"] = slope,
    _["p_value"] = p_value,
    _["df"] = df,
    _["conf_low"] = conf_low,
    _["conf_high"] = conf_high,
    _["weighted_ss_x"] = weighted_ss_x,
    _["weighted_ss_e"] = weighted_ss_e,
    _["n_complete"] = n_complete,
    _["n_subjects"] = n_subjects,
    _["obs_per_subject_min"] = obs_per_subject_min,
    _["obs_per_subject_max"] = obs_per_subject_max,
    _["valid"] = valid
  );
}
