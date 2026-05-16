// Thiago de Paula Oliveira
// cia_repeated.cpp
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

namespace {

const char* kCiaRmBalancedError =
  "cia_rm() currently implements the balanced categorical repeated-measures ANOVA estimator; each subject-method-time cell must contain exactly one observation.";

inline int cube_index(const int s,
                      const int m,
                      const int t,
                      const int n_methods,
                      const int n_times) {
  return (s * n_methods + m) * n_times + t;
}

inline int est_index(const int i,
                     const int j,
                     const int k,
                     const int n_methods,
                     const int n_times) {
  return i + n_methods * (j + n_methods * k);
}

inline int matrix_index(const int row, const int col, const int nrow) {
  return row + nrow * col;
}

std::vector<double> build_cube(const Rcpp::NumericVector& y,
                               const Rcpp::IntegerVector& subject,
                               const Rcpp::IntegerVector& method,
                               const Rcpp::IntegerVector& time,
                               const int n_subjects,
                               const int n_methods,
                               const int n_times) {
  const R_xlen_t n_obs = y.size();
  if (subject.size() != n_obs || method.size() != n_obs || time.size() != n_obs) {
    Rcpp::stop("y, subject, method, and time must have the same length.");
  }

  const int total = n_subjects * n_methods * n_times;
  std::vector<double> cube(static_cast<std::size_t>(total), NA_REAL);
  std::vector<int> counts(static_cast<std::size_t>(total), 0);

  for (R_xlen_t i = 0; i < n_obs; ++i) {
    if (!R_finite(y[i])) {
      Rcpp::stop("response must contain only finite numeric values.");
    }
    if (subject[i] == NA_INTEGER || method[i] == NA_INTEGER || time[i] == NA_INTEGER) {
      Rcpp::stop(kCiaRmBalancedError);
    }

    const int s = subject[i] - 1;
    const int m = method[i] - 1;
    const int t = time[i] - 1;
    if (s < 0 || s >= n_subjects ||
        m < 0 || m >= n_methods ||
        t < 0 || t >= n_times) {
      Rcpp::stop(kCiaRmBalancedError);
    }

    const int idx = cube_index(s, m, t, n_methods, n_times);
    counts[static_cast<std::size_t>(idx)] += 1;
    if (counts[static_cast<std::size_t>(idx)] > 1) {
      Rcpp::stop(kCiaRmBalancedError);
    }
    cube[static_cast<std::size_t>(idx)] = y[i];
  }

  for (int idx = 0; idx < total; ++idx) {
    if (counts[static_cast<std::size_t>(idx)] != 1 || !R_finite(cube[static_cast<std::size_t>(idx)])) {
      Rcpp::stop(kCiaRmBalancedError);
    }
  }

  return cube;
}

struct OverallRmStats {
  Rcpp::NumericVector overall;
  double overall_common;
  double sigma2_error;
  double sigma2_subject_method;
  double repeatability;
  double ms_subject_method;
  double ms_method_time;
  double ms_error;
  double ss_subject_method;
  double ss_method_time;
  double ss_error;
  double df_subject_method;
  double df_method_time;
  double df_error;
  double homogeneity_F;
  double homogeneity_p;
  double common_sigma2_error;
  double common_sigma2_subject_method;
  double common_repeatability;

  OverallRmStats(const int n_times)
    : overall(n_times, NA_REAL),
      overall_common(NA_REAL),
      sigma2_error(NA_REAL),
      sigma2_subject_method(NA_REAL),
      repeatability(NA_REAL),
      ms_subject_method(NA_REAL),
      ms_method_time(NA_REAL),
      ms_error(NA_REAL),
      ss_subject_method(NA_REAL),
      ss_method_time(NA_REAL),
      ss_error(NA_REAL),
      df_subject_method(NA_REAL),
      df_method_time(NA_REAL),
      df_error(NA_REAL),
      homogeneity_F(NA_REAL),
      homogeneity_p(NA_REAL),
      common_sigma2_error(NA_REAL),
      common_sigma2_subject_method(NA_REAL),
      common_repeatability(NA_REAL) {}
};

OverallRmStats compute_overall_rm_stats(const std::vector<double>& cube,
                                        const int n_subjects,
                                        const int n_methods,
                                        const int n_times,
                                        const bool homogeneous,
                                        const bool constrain_vc) {
  OverallRmStats out(n_times);

  const double n_subj_d = static_cast<double>(n_subjects);
  const double n_methods_d = static_cast<double>(n_methods);
  const double n_times_d = static_cast<double>(n_times);
  const double n_pairs_d =
    static_cast<double>(n_methods * (n_methods - 1)) / 2.0;

  std::vector<double> subject_mean(static_cast<std::size_t>(n_subjects), 0.0);
  std::vector<double> method_mean(static_cast<std::size_t>(n_methods), 0.0);
  std::vector<double> time_mean(static_cast<std::size_t>(n_times), 0.0);
  std::vector<double> subject_method_mean(
    static_cast<std::size_t>(n_subjects) * n_methods, 0.0
  );
  std::vector<double> subject_time_mean(
    static_cast<std::size_t>(n_subjects) * n_times, 0.0
  );
  std::vector<double> method_time_mean(
    static_cast<std::size_t>(n_methods) * n_times, 0.0
  );
  double grand_sum = 0.0;

  for (int s = 0; s < n_subjects; ++s) {
    for (int j = 0; j < n_methods; ++j) {
      for (int k = 0; k < n_times; ++k) {
        const double value =
          cube[static_cast<std::size_t>(cube_index(s, j, k, n_methods, n_times))];
        subject_mean[static_cast<std::size_t>(s)] += value;
        method_mean[static_cast<std::size_t>(j)] += value;
        time_mean[static_cast<std::size_t>(k)] += value;
        subject_method_mean[static_cast<std::size_t>(matrix_index(s, j, n_subjects))] += value;
        subject_time_mean[static_cast<std::size_t>(matrix_index(s, k, n_subjects))] += value;
        method_time_mean[static_cast<std::size_t>(matrix_index(j, k, n_methods))] += value;
        grand_sum += value;
      }
    }
  }

  for (int s = 0; s < n_subjects; ++s) {
    subject_mean[static_cast<std::size_t>(s)] /= (n_methods_d * n_times_d);
  }
  for (int j = 0; j < n_methods; ++j) {
    method_mean[static_cast<std::size_t>(j)] /= (n_subj_d * n_times_d);
    for (int s = 0; s < n_subjects; ++s) {
      subject_method_mean[static_cast<std::size_t>(matrix_index(s, j, n_subjects))] /= n_times_d;
    }
  }
  for (int k = 0; k < n_times; ++k) {
    time_mean[static_cast<std::size_t>(k)] /= (n_subj_d * n_methods_d);
    for (int s = 0; s < n_subjects; ++s) {
      subject_time_mean[static_cast<std::size_t>(matrix_index(s, k, n_subjects))] /= n_methods_d;
    }
    for (int j = 0; j < n_methods; ++j) {
      method_time_mean[static_cast<std::size_t>(matrix_index(j, k, n_methods))] /= n_subj_d;
    }
  }

  const double grand_mean = grand_sum / (n_subj_d * n_methods_d * n_times_d);

  double ss_sm = 0.0;
  double ss_mt = 0.0;
  double ss_e = 0.0;

  for (int s = 0; s < n_subjects; ++s) {
    const double zi = subject_mean[static_cast<std::size_t>(s)];
    for (int j = 0; j < n_methods; ++j) {
      const double zij =
        subject_method_mean[static_cast<std::size_t>(matrix_index(s, j, n_subjects))];
      const double zj = method_mean[static_cast<std::size_t>(j)];
      ss_sm += n_times_d * std::pow(zij - zi - zj + grand_mean, 2.0);
    }
  }

  for (int j = 0; j < n_methods; ++j) {
    const double zj = method_mean[static_cast<std::size_t>(j)];
    for (int k = 0; k < n_times; ++k) {
      const double zjk =
        method_time_mean[static_cast<std::size_t>(matrix_index(j, k, n_methods))];
      const double zk = time_mean[static_cast<std::size_t>(k)];
      ss_mt += n_subj_d * std::pow(zjk - zj - zk + grand_mean, 2.0);
    }
  }

  for (int s = 0; s < n_subjects; ++s) {
    const double zi = subject_mean[static_cast<std::size_t>(s)];
    for (int j = 0; j < n_methods; ++j) {
      const double zij =
        subject_method_mean[static_cast<std::size_t>(matrix_index(s, j, n_subjects))];
      const double zj = method_mean[static_cast<std::size_t>(j)];
      for (int k = 0; k < n_times; ++k) {
        const double zik =
          subject_time_mean[static_cast<std::size_t>(matrix_index(s, k, n_subjects))];
        const double zjk =
          method_time_mean[static_cast<std::size_t>(matrix_index(j, k, n_methods))];
        const double zk = time_mean[static_cast<std::size_t>(k)];
        const double value =
          cube[static_cast<std::size_t>(cube_index(s, j, k, n_methods, n_times))];
        const double resid =
          value - zij - zik - zjk + zi + zj + zk - grand_mean;
        ss_e += resid * resid;
      }
    }
  }

  const double df_sm =
    static_cast<double>((n_subjects - 1) * (n_methods - 1));
  const double df_mt =
    static_cast<double>((n_methods - 1) * (n_times - 1));
  const double df_err =
    static_cast<double>((n_subjects - 1) * (n_methods - 1) * (n_times - 1));

  const double ms_sm = ss_sm / df_sm;
  const double ms_mt = ss_mt / df_mt;
  const double ms_err = ss_e / df_err;

  double sigma_e2 = ms_err;
  double sigma_sm2 = (ms_sm - ms_err) / n_times_d;
  if (constrain_vc) {
    sigma_e2 = std::max(sigma_e2, 0.0);
    sigma_sm2 = std::max(sigma_sm2, 0.0);
  }

  double f_stat = NA_REAL;
  double p_val = NA_REAL;
  if (R_finite(ms_mt) && R_finite(ms_err)) {
    if (ms_err > 0.0) {
      f_stat = ms_mt / ms_err;
      p_val = R::pf(f_stat, df_mt, df_err, false, false);
    } else if (ms_mt > 0.0) {
      f_stat = R_PosInf;
      p_val = 0.0;
    }
  }

  for (int k = 0; k < n_times; ++k) {
    double mean_sq_diff = 0.0;
    for (int a = 0; a < n_methods - 1; ++a) {
      const double za =
        method_time_mean[static_cast<std::size_t>(matrix_index(a, k, n_methods))];
      for (int b = a + 1; b < n_methods; ++b) {
        const double zb =
          method_time_mean[static_cast<std::size_t>(matrix_index(b, k, n_methods))];
        const double diff = za - zb;
        mean_sq_diff += diff * diff;
      }
    }
    mean_sq_diff /= n_pairs_d;
    const double denom = mean_sq_diff + 2.0 * sigma_sm2 + 2.0 * sigma_e2;
    if (R_finite(sigma_e2) && R_finite(denom) && denom > 0.0) {
      out.overall[k] = 2.0 * sigma_e2 / denom;
    }
  }

  if (R_finite(sigma_e2) && sigma_e2 >= 0.0) {
    out.repeatability = 1.96 * std::sqrt(2.0 * sigma_e2);
  }

  out.sigma2_error = sigma_e2;
  out.sigma2_subject_method = sigma_sm2;
  out.ms_subject_method = ms_sm;
  out.ms_method_time = ms_mt;
  out.ms_error = ms_err;
  out.ss_subject_method = ss_sm;
  out.ss_method_time = ss_mt;
  out.ss_error = ss_e;
  out.df_subject_method = df_sm;
  out.df_method_time = df_mt;
  out.df_error = df_err;
  out.homogeneity_F = f_stat;
  out.homogeneity_p = p_val;

  if (homogeneous) {
    const double ss_err_reduced = ss_mt + ss_e;
    const double df_err_reduced = df_mt + df_err;
    const double ms_err_reduced = ss_err_reduced / df_err_reduced;
    double common_sigma_e2 = ms_err_reduced;
    double common_sigma_sm2 = (ms_sm - ms_err_reduced) / n_times_d;
    if (constrain_vc) {
      common_sigma_e2 = std::max(common_sigma_e2, 0.0);
      common_sigma_sm2 = std::max(common_sigma_sm2, 0.0);
    }

    double mean_sq_diff_common = 0.0;
    for (int a = 0; a < n_methods - 1; ++a) {
      const double za = method_mean[static_cast<std::size_t>(a)];
      for (int b = a + 1; b < n_methods; ++b) {
        const double zb = method_mean[static_cast<std::size_t>(b)];
        const double diff = za - zb;
        mean_sq_diff_common += diff * diff;
      }
    }
    mean_sq_diff_common /= n_pairs_d;

    const double common_denom =
      mean_sq_diff_common + 2.0 * common_sigma_sm2 + 2.0 * common_sigma_e2;
    if (R_finite(common_sigma_e2) &&
        R_finite(common_denom) &&
        common_denom > 0.0) {
      out.overall_common = 2.0 * common_sigma_e2 / common_denom;
    }
    if (R_finite(common_sigma_e2) && common_sigma_e2 >= 0.0) {
      out.common_repeatability = 1.96 * std::sqrt(2.0 * common_sigma_e2);
    }
    out.common_sigma2_error = common_sigma_e2;
    out.common_sigma2_subject_method = common_sigma_sm2;
  }

  return out;
}

struct PairMomentPayload {
  std::vector<double> est;
  double common;

  explicit PairMomentPayload(const int n_times)
    : est(static_cast<std::size_t>(n_times), NA_REAL),
      common(NA_REAL) {}
};

inline void apply_bounds_scalar(double& lwr,
                                double& upr,
                                const bool constrain_vc) {
  if (constrain_vc) {
    lwr = std::max(0.0, lwr);
    upr = std::min(1.0, upr);
  }
}

double mean_pair_sqdiff_cpp(const std::vector<double>& x) {
  const int J = static_cast<int>(x.size());
  if (J < 2) {
    return NA_REAL;
  }

  double acc = 0.0;
  int count = 0;
  for (int i = 0; i < J - 1; ++i) {
    for (int j = i + 1; j < J; ++j) {
      const double diff = x[static_cast<std::size_t>(i)] - x[static_cast<std::size_t>(j)];
      acc += diff * diff;
      ++count;
    }
  }
  return acc / static_cast<double>(count);
}

PairMomentPayload pair_payload_from_moments(const std::vector<double>& m,
                                            const int K,
                                            const int n,
                                            const bool constrain_vc,
                                            const bool homogeneous) {
  PairMomentPayload out(K);
  const double fac = static_cast<double>(n) / static_cast<double>(n - 1);

  std::vector<double> d(static_cast<std::size_t>(K), 0.0);
  std::copy(m.begin(), m.begin() + K, d.begin());

  double off_sum = 0.0;
  double diag_sum = 0.0;
  double mean_a2 = 0.0;
  for (int r = 0; r < K; ++r) {
    for (int c = 0; c < K; ++c) {
      const double Mrc = m[static_cast<std::size_t>(K + r + K * c)];
      const double Crc =
        fac * (Mrc - d[static_cast<std::size_t>(r)] * d[static_cast<std::size_t>(c)]);
      mean_a2 += Mrc;
      if (r == c) {
        diag_sum += Crc;
      } else {
        off_sum += Crc;
      }
    }
  }
  mean_a2 /= static_cast<double>(K * K);

  const double off_mean = off_sum / static_cast<double>(K * (K - 1));
  const double diag_mean = diag_sum / static_cast<double>(K);
  double sigma2_subject_method = off_mean / 2.0;
  const double sigma2_error = (diag_mean - off_mean) / 2.0;
  if (constrain_vc) {
    sigma2_subject_method = std::max(sigma2_subject_method, 0.0);
  }

  for (int k = 0; k < K; ++k) {
    const double denom =
      d[static_cast<std::size_t>(k)] * d[static_cast<std::size_t>(k)] +
      2.0 * sigma2_subject_method +
      2.0 * sigma2_error;
    if (R_finite(denom) && denom > 0.0) {
      out.est[static_cast<std::size_t>(k)] = 2.0 * sigma2_error / denom;
    }
  }

  if (homogeneous) {
    const double d_bar =
      std::accumulate(d.begin(), d.end(), 0.0) / static_cast<double>(K);
    double ss_time = 0.0;
    double ss_error_reduced = 0.0;
    for (int r = 0; r < K; ++r) {
      const double dr = d[static_cast<std::size_t>(r)] - d_bar;
      ss_time += dr * dr;
      for (int c = 0; c < K; ++c) {
        const double Mrc = m[static_cast<std::size_t>(K + r + K * c)];
        const double Qrc =
          Mrc - d[static_cast<std::size_t>(r)] * d[static_cast<std::size_t>(c)];
        const double Prc =
          (r == c ? 1.0 : 0.0) - 1.0 / static_cast<double>(K);
        ss_error_reduced += Prc * Qrc;
      }
    }
    ss_time *= static_cast<double>(n);
    ss_error_reduced *= static_cast<double>(n);
    const double ms_error_reduced =
      (ss_time + ss_error_reduced) /
      (static_cast<double>(n) * static_cast<double>(K - 1));
    const double ms_subject_method =
      static_cast<double>(K) * fac * (mean_a2 - d_bar * d_bar);

    double common_sigma2_subject_method =
      (ms_subject_method - ms_error_reduced) / (2.0 * static_cast<double>(K));
    if (constrain_vc) {
      common_sigma2_subject_method = std::max(common_sigma2_subject_method, 0.0);
    }
    const double common_sigma2_error = ms_error_reduced / 2.0;
    const double common_denom =
      d_bar * d_bar +
      2.0 * common_sigma2_subject_method +
      2.0 * common_sigma2_error;
    if (R_finite(common_denom) && common_denom > 0.0) {
      out.common = 2.0 * common_sigma2_error / common_denom;
    }
  }

  return out;
}

struct OverallMomentPayload {
  std::vector<double> overall;
  double overall_common;

  explicit OverallMomentPayload(const int n_times)
    : overall(static_cast<std::size_t>(n_times), NA_REAL),
      overall_common(NA_REAL) {}
};

OverallMomentPayload overall_payload_from_moments(const std::vector<double>& m,
                                                  const int J,
                                                  const int K,
                                                  const int n,
                                                  const bool homogeneous,
                                                  const bool constrain_vc) {
  const int P = J * K;
  const double fac = static_cast<double>(n) / static_cast<double>(n - 1);
  OverallMomentPayload out(K);

  std::vector<double> mu(static_cast<std::size_t>(P), 0.0);
  std::copy(m.begin(), m.begin() + P, mu.begin());

  std::vector<double> Sx(static_cast<std::size_t>(P * P), 0.0);
  for (int r = 0; r < P; ++r) {
    for (int c = 0; c < P; ++c) {
      const double Mrc = m[static_cast<std::size_t>(P + r + P * c)];
      Sx[static_cast<std::size_t>(r + P * c)] =
        fac * (Mrc - mu[static_cast<std::size_t>(r)] * mu[static_cast<std::size_t>(c)]);
    }
  }

  std::vector<double> method_mean(static_cast<std::size_t>(J), 0.0);
  std::vector<double> time_mean(static_cast<std::size_t>(K), 0.0);
  double grand_mean = 0.0;
  for (int k = 0; k < K; ++k) {
    for (int j = 0; j < J; ++j) {
      const double value = mu[static_cast<std::size_t>(j + J * k)];
      method_mean[static_cast<std::size_t>(j)] += value;
      time_mean[static_cast<std::size_t>(k)] += value;
      grand_mean += value;
    }
  }
  for (int j = 0; j < J; ++j) {
    method_mean[static_cast<std::size_t>(j)] /= static_cast<double>(K);
  }
  for (int k = 0; k < K; ++k) {
    time_mean[static_cast<std::size_t>(k)] /= static_cast<double>(J);
  }
  grand_mean /= static_cast<double>(P);

  std::vector<double> S_method(static_cast<std::size_t>(J * J), 0.0);
  for (int a = 0; a < J; ++a) {
    for (int b = 0; b < J; ++b) {
      double acc = 0.0;
      for (int k1 = 0; k1 < K; ++k1) {
        const int idx1 = a + J * k1;
        for (int k2 = 0; k2 < K; ++k2) {
          const int idx2 = b + J * k2;
          acc += Sx[static_cast<std::size_t>(idx1 + P * idx2)];
        }
      }
      S_method[static_cast<std::size_t>(a + J * b)] =
        acc / static_cast<double>(K * K);
    }
  }

  double ss_subject_method = 0.0;
  for (int a = 0; a < J; ++a) {
    for (int b = 0; b < J; ++b) {
      const double Pjab =
        (a == b ? 1.0 : 0.0) - 1.0 / static_cast<double>(J);
      ss_subject_method +=
        Pjab * S_method[static_cast<std::size_t>(a + J * b)];
    }
  }
  ss_subject_method *= static_cast<double>(K * (n - 1));

  double ss_error = 0.0;
  for (int p = 0; p < P; ++p) {
    const int jp = p % J;
    const int kp = p / J;
    for (int q = 0; q < P; ++q) {
      const int jq = q % J;
      const int kq = q / J;
      const double Pjab =
        (jp == jq ? 1.0 : 0.0) - 1.0 / static_cast<double>(J);
      const double Pkab =
        (kp == kq ? 1.0 : 0.0) - 1.0 / static_cast<double>(K);
      ss_error +=
        Pjab * Pkab * Sx[static_cast<std::size_t>(p + P * q)];
    }
  }
  ss_error *= static_cast<double>(n - 1);

  double ss_method_time = 0.0;
  for (int j = 0; j < J; ++j) {
    for (int k = 0; k < K; ++k) {
      const double value = mu[static_cast<std::size_t>(j + J * k)];
      const double interaction =
        value -
        method_mean[static_cast<std::size_t>(j)] -
        time_mean[static_cast<std::size_t>(k)] +
        grand_mean;
      ss_method_time += interaction * interaction;
    }
  }
  ss_method_time *= static_cast<double>(n);

  const double df_subject_method =
    static_cast<double>((n - 1) * (J - 1));
  const double df_method_time =
    static_cast<double>((J - 1) * (K - 1));
  const double df_error =
    static_cast<double>((n - 1) * (J - 1) * (K - 1));

  const double ms_subject_method = ss_subject_method / df_subject_method;
  const double ms_error = ss_error / df_error;
  double sigma2_error = ms_error;
  double sigma2_subject_method =
    (ms_subject_method - ms_error) / static_cast<double>(K);
  if (constrain_vc) {
    sigma2_error = std::max(sigma2_error, 0.0);
    sigma2_subject_method = std::max(sigma2_subject_method, 0.0);
  }

  for (int k = 0; k < K; ++k) {
    std::vector<double> method_slice(static_cast<std::size_t>(J), 0.0);
    for (int j = 0; j < J; ++j) {
      method_slice[static_cast<std::size_t>(j)] =
        mu[static_cast<std::size_t>(j + J * k)];
    }
    const double mean_sq_diff = mean_pair_sqdiff_cpp(method_slice);
    const double denom =
      mean_sq_diff + 2.0 * sigma2_subject_method + 2.0 * sigma2_error;
    if (R_finite(denom) && denom > 0.0) {
      out.overall[static_cast<std::size_t>(k)] = 2.0 * sigma2_error / denom;
    }
  }

  if (homogeneous) {
    const double ss_error_reduced = ss_method_time + ss_error;
    const double df_error_reduced = df_method_time + df_error;
    const double ms_error_reduced = ss_error_reduced / df_error_reduced;
    double common_sigma2_error = ms_error_reduced;
    double common_sigma2_subject_method =
      (ms_subject_method - ms_error_reduced) / static_cast<double>(K);
    if (constrain_vc) {
      common_sigma2_error = std::max(common_sigma2_error, 0.0);
      common_sigma2_subject_method = std::max(common_sigma2_subject_method, 0.0);
    }
    const double mean_sq_diff_common = mean_pair_sqdiff_cpp(method_mean);
    const double denom_common =
      mean_sq_diff_common +
      2.0 * common_sigma2_subject_method +
      2.0 * common_sigma2_error;
    if (R_finite(denom_common) && denom_common > 0.0) {
      out.overall_common = 2.0 * common_sigma2_error / denom_common;
    }
  }

  return out;
}

std::vector<double> numerical_gradient_matrix(
    const std::function<std::vector<double>(const std::vector<double>&)>& fun,
    const std::vector<double>& m,
    const int out_dim) {
  const int p = static_cast<int>(m.size());
  std::vector<double> grad(static_cast<std::size_t>(p * out_dim), 0.0);
  std::vector<double> plus = m;
  std::vector<double> minus = m;

  for (int idx = 0; idx < p; ++idx) {
    const double step = 1e-6 * std::max(1.0, std::abs(m[static_cast<std::size_t>(idx)]));
    plus[static_cast<std::size_t>(idx)] = m[static_cast<std::size_t>(idx)] + step;
    minus[static_cast<std::size_t>(idx)] = m[static_cast<std::size_t>(idx)] - step;
    const std::vector<double> f_plus = fun(plus);
    const std::vector<double> f_minus = fun(minus);
    for (int out_idx = 0; out_idx < out_dim; ++out_idx) {
      grad[static_cast<std::size_t>(idx + p * out_idx)] =
        (f_plus[static_cast<std::size_t>(out_idx)] -
         f_minus[static_cast<std::size_t>(out_idx)]) / (2.0 * step);
    }
    plus[static_cast<std::size_t>(idx)] = m[static_cast<std::size_t>(idx)];
    minus[static_cast<std::size_t>(idx)] = m[static_cast<std::size_t>(idx)];
  }

  return grad;
}

double projection_variance_quadratic(const std::vector<double>& rows,
                                     const int n_rows,
                                     const int dim,
                                     const std::vector<double>& grad) {
  if (n_rows < 2) {
    return NA_REAL;
  }

  const int offset = dim;
  std::vector<double> h(static_cast<std::size_t>(n_rows), 0.0);
  for (int i = 0; i < n_rows; ++i) {
    const double* row_ptr = rows.data() + static_cast<std::size_t>(i * dim);
    double value = 0.0;
    for (int a = 0; a < dim; ++a) {
      value += grad[static_cast<std::size_t>(a)] * row_ptr[a];
    }
    for (int a = 0; a < dim; ++a) {
      for (int b = 0; b < dim; ++b) {
        value +=
          grad[static_cast<std::size_t>(offset + a + dim * b)] *
          row_ptr[a] * row_ptr[b];
      }
    }
    h[static_cast<std::size_t>(i)] = value;
  }

  const double mean_h =
    std::accumulate(h.begin(), h.end(), 0.0) / static_cast<double>(n_rows);
  double ss = 0.0;
  for (int i = 0; i < n_rows; ++i) {
    const double centered = h[static_cast<std::size_t>(i)] - mean_h;
    ss += centered * centered;
  }
  return (ss / static_cast<double>(n_rows - 1)) / static_cast<double>(n_rows);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cia_rm_anova_cpp(Rcpp::NumericVector y,
                            Rcpp::IntegerVector subject,
                            Rcpp::IntegerVector method,
                            Rcpp::IntegerVector time,
                            int n_subjects,
                            int n_methods,
                            int n_times,
                            bool homogeneous,
                            bool constrain_vc,
                            int n_threads = 1) {
  (void)n_threads;

  if (n_subjects < 2 || n_methods < 2 || n_times < 2) {
    Rcpp::stop(kCiaRmBalancedError);
  }

  const std::vector<double> cube =
    build_cube(y, subject, method, time, n_subjects, n_methods, n_times);
  const OverallRmStats overall_stats =
    compute_overall_rm_stats(
      cube,
      n_subjects,
      n_methods,
      n_times,
      homogeneous,
      constrain_vc
    );

  Rcpp::NumericVector est(
    static_cast<R_xlen_t>(n_methods) * n_methods * n_times,
    NA_REAL
  );
  est.attr("dim") = Rcpp::IntegerVector::create(n_methods, n_methods, n_times);

  Rcpp::NumericVector method_time_diff(
    static_cast<R_xlen_t>(n_methods) * n_methods * n_times,
    NA_REAL
  );
  method_time_diff.attr("dim") = Rcpp::IntegerVector::create(n_methods, n_methods, n_times);

  Rcpp::NumericMatrix sigma2_error(n_methods, n_methods);
  Rcpp::NumericMatrix sigma2_subject_method(n_methods, n_methods);
  Rcpp::NumericMatrix repeatability(n_methods, n_methods);
  Rcpp::NumericMatrix ms_subject_method(n_methods, n_methods);
  Rcpp::NumericMatrix ms_method_time(n_methods, n_methods);
  Rcpp::NumericMatrix ms_error(n_methods, n_methods);
  Rcpp::NumericMatrix ss_subject_method(n_methods, n_methods);
  Rcpp::NumericMatrix ss_method_time(n_methods, n_methods);
  Rcpp::NumericMatrix ss_error(n_methods, n_methods);
  Rcpp::NumericMatrix df_subject_method(n_methods, n_methods);
  Rcpp::NumericMatrix df_method_time(n_methods, n_methods);
  Rcpp::NumericMatrix df_error(n_methods, n_methods);
  Rcpp::NumericMatrix homogeneity_F(n_methods, n_methods);
  Rcpp::NumericMatrix homogeneity_p(n_methods, n_methods);
  Rcpp::NumericMatrix common(n_methods, n_methods);
  Rcpp::NumericMatrix common_sigma2_error(n_methods, n_methods);
  Rcpp::NumericMatrix common_sigma2_subject_method(n_methods, n_methods);
  Rcpp::NumericMatrix common_repeatability(n_methods, n_methods);

  std::fill(sigma2_error.begin(), sigma2_error.end(), NA_REAL);
  std::fill(sigma2_subject_method.begin(), sigma2_subject_method.end(), NA_REAL);
  std::fill(repeatability.begin(), repeatability.end(), NA_REAL);
  std::fill(ms_subject_method.begin(), ms_subject_method.end(), NA_REAL);
  std::fill(ms_method_time.begin(), ms_method_time.end(), NA_REAL);
  std::fill(ms_error.begin(), ms_error.end(), NA_REAL);
  std::fill(ss_subject_method.begin(), ss_subject_method.end(), NA_REAL);
  std::fill(ss_method_time.begin(), ss_method_time.end(), NA_REAL);
  std::fill(ss_error.begin(), ss_error.end(), NA_REAL);
  std::fill(df_subject_method.begin(), df_subject_method.end(), NA_REAL);
  std::fill(df_method_time.begin(), df_method_time.end(), NA_REAL);
  std::fill(df_error.begin(), df_error.end(), NA_REAL);
  std::fill(homogeneity_F.begin(), homogeneity_F.end(), NA_REAL);
  std::fill(homogeneity_p.begin(), homogeneity_p.end(), NA_REAL);
  std::fill(common.begin(), common.end(), NA_REAL);
  std::fill(common_sigma2_error.begin(), common_sigma2_error.end(), NA_REAL);
  std::fill(common_sigma2_subject_method.begin(), common_sigma2_subject_method.end(), NA_REAL);
  std::fill(common_repeatability.begin(), common_repeatability.end(), NA_REAL);

  for (int j = 0; j < n_methods; ++j) {
    for (int k = 0; k < n_times; ++k) {
      est[est_index(j, j, k, n_methods, n_times)] = 1.0;
      method_time_diff[est_index(j, j, k, n_methods, n_times)] = 0.0;
    }
    common(j, j) = 1.0;
  }

  const double n_subj_d = static_cast<double>(n_subjects);
  const double n_times_d = static_cast<double>(n_times);
  const double two_times_d = 2.0 * n_times_d;

  for (int a = 0; a < n_methods - 1; ++a) {
    for (int b = a + 1; b < n_methods; ++b) {
      std::vector<double> subject_mean(static_cast<std::size_t>(n_subjects), 0.0);
      std::vector<double> subject_method_mean(static_cast<std::size_t>(n_subjects) * 2, 0.0);
      std::vector<double> subject_time_mean(static_cast<std::size_t>(n_subjects) * n_times, 0.0);
      std::vector<double> method_time_mean(static_cast<std::size_t>(2) * n_times, 0.0);
      std::vector<double> time_mean(static_cast<std::size_t>(n_times), 0.0);
      double grand_sum = 0.0;
      double method_sum[2] = {0.0, 0.0};

      for (int s = 0; s < n_subjects; ++s) {
        double subj_sum = 0.0;
        double subj_m0 = 0.0;
        double subj_m1 = 0.0;
        for (int k = 0; k < n_times; ++k) {
          const double z0 = cube[static_cast<std::size_t>(cube_index(s, a, k, n_methods, n_times))];
          const double z1 = cube[static_cast<std::size_t>(cube_index(s, b, k, n_methods, n_times))];
          subj_sum += z0 + z1;
          subj_m0 += z0;
          subj_m1 += z1;
          subject_time_mean[static_cast<std::size_t>(s * n_times + k)] = (z0 + z1) / 2.0;
          method_time_mean[static_cast<std::size_t>(k)] += z0;
          method_time_mean[static_cast<std::size_t>(n_times + k)] += z1;
          time_mean[static_cast<std::size_t>(k)] += z0 + z1;
          method_sum[0] += z0;
          method_sum[1] += z1;
          grand_sum += z0 + z1;
        }
        subject_mean[static_cast<std::size_t>(s)] = subj_sum / two_times_d;
        subject_method_mean[static_cast<std::size_t>(s * 2)] = subj_m0 / n_times_d;
        subject_method_mean[static_cast<std::size_t>(s * 2 + 1)] = subj_m1 / n_times_d;
      }

      const double grand_mean = grand_sum / (n_subj_d * 2.0 * n_times_d);
      const double method_mean0 = method_sum[0] / (n_subj_d * n_times_d);
      const double method_mean1 = method_sum[1] / (n_subj_d * n_times_d);

      for (int k = 0; k < n_times; ++k) {
        method_time_mean[static_cast<std::size_t>(k)] /= n_subj_d;
        method_time_mean[static_cast<std::size_t>(n_times + k)] /= n_subj_d;
        time_mean[static_cast<std::size_t>(k)] /= (2.0 * n_subj_d);
      }

      double ss_sm = 0.0;
      double ss_mt = 0.0;
      double ss_e = 0.0;

      for (int s = 0; s < n_subjects; ++s) {
        const double zi = subject_mean[static_cast<std::size_t>(s)];
        const double zir0 = subject_method_mean[static_cast<std::size_t>(s * 2)];
        const double zir1 = subject_method_mean[static_cast<std::size_t>(s * 2 + 1)];
        ss_sm += n_times_d * std::pow(zir0 - zi - method_mean0 + grand_mean, 2.0);
        ss_sm += n_times_d * std::pow(zir1 - zi - method_mean1 + grand_mean, 2.0);

        for (int k = 0; k < n_times; ++k) {
          const double zi_k = subject_time_mean[static_cast<std::size_t>(s * n_times + k)];
          const double z0 = cube[static_cast<std::size_t>(cube_index(s, a, k, n_methods, n_times))];
          const double z1 = cube[static_cast<std::size_t>(cube_index(s, b, k, n_methods, n_times))];
          const double z0k = method_time_mean[static_cast<std::size_t>(k)];
          const double z1k = method_time_mean[static_cast<std::size_t>(n_times + k)];
          const double z_k = time_mean[static_cast<std::size_t>(k)];

          const double resid0 = z0 - zir0 - zi_k - z0k + zi + method_mean0 + z_k - grand_mean;
          const double resid1 = z1 - zir1 - zi_k - z1k + zi + method_mean1 + z_k - grand_mean;
          ss_e += resid0 * resid0 + resid1 * resid1;
        }
      }

      for (int k = 0; k < n_times; ++k) {
        const double z0k = method_time_mean[static_cast<std::size_t>(k)];
        const double z1k = method_time_mean[static_cast<std::size_t>(n_times + k)];
        const double z_k = time_mean[static_cast<std::size_t>(k)];
        ss_mt += n_subj_d * std::pow(z0k - method_mean0 - z_k + grand_mean, 2.0);
        ss_mt += n_subj_d * std::pow(z1k - method_mean1 - z_k + grand_mean, 2.0);
      }

      const double df_sm = static_cast<double>(n_subjects - 1);
      const double df_mt = static_cast<double>(n_times - 1);
      const double df_err = static_cast<double>((n_subjects - 1) * (n_times - 1));
      const double ms_sm = ss_sm / df_sm;
      const double ms_mt = ss_mt / df_mt;
      const double ms_err = ss_e / df_err;

      double sigma_e2 = ms_err;
      double sigma_ab2 = (ms_sm - ms_err) / n_times_d;
      if (constrain_vc) {
        sigma_e2 = std::max(sigma_e2, 0.0);
        sigma_ab2 = std::max(sigma_ab2, 0.0);
      }

      double f_stat = NA_REAL;
      double p_val = NA_REAL;
      if (R_finite(ms_mt) && R_finite(ms_err)) {
        if (ms_err > 0.0) {
          f_stat = ms_mt / ms_err;
          p_val = R::pf(f_stat, df_mt, df_err, false, false);
        } else if (ms_mt > 0.0) {
          f_stat = R_PosInf;
          p_val = 0.0;
        }
      }

      double common_sigma_e2 = NA_REAL;
      double common_sigma_ab2 = NA_REAL;
      double common_repeat = NA_REAL;
      double common_est = NA_REAL;
      if (homogeneous) {
        const double ss_err_reduced = ss_mt + ss_e;
        const double df_err_reduced = df_mt + df_err;
        const double ms_err_reduced = ss_err_reduced / df_err_reduced;
        common_sigma_e2 = ms_err_reduced;
        common_sigma_ab2 = (ms_sm - ms_err_reduced) / n_times_d;
        if (constrain_vc) {
          common_sigma_e2 = std::max(common_sigma_e2, 0.0);
          common_sigma_ab2 = std::max(common_sigma_ab2, 0.0);
        }
        double d_common = 0.0;
        for (int s = 0; s < n_subjects; ++s) {
          for (int k = 0; k < n_times; ++k) {
            const double z0 = cube[static_cast<std::size_t>(cube_index(s, a, k, n_methods, n_times))];
            const double z1 = cube[static_cast<std::size_t>(cube_index(s, b, k, n_methods, n_times))];
            d_common += (z0 - z1);
          }
        }
        d_common /= (n_subj_d * n_times_d);
        const double common_denom = d_common * d_common + 2.0 * common_sigma_ab2 + 2.0 * common_sigma_e2;
        if (R_finite(common_sigma_e2) && R_finite(common_denom) && common_denom > 0.0) {
          common_est = 2.0 * common_sigma_e2 / common_denom;
        }
        if (R_finite(common_sigma_e2) && common_sigma_e2 >= 0.0) {
          common_repeat = 1.96 * std::sqrt(2.0 * common_sigma_e2);
        }
      }

      for (int k = 0; k < n_times; ++k) {
        double d_k = 0.0;
        for (int s = 0; s < n_subjects; ++s) {
          const double z0 = cube[static_cast<std::size_t>(cube_index(s, a, k, n_methods, n_times))];
          const double z1 = cube[static_cast<std::size_t>(cube_index(s, b, k, n_methods, n_times))];
          d_k += (z0 - z1);
        }
        d_k /= n_subj_d;
        const double denom = d_k * d_k + 2.0 * sigma_ab2 + 2.0 * sigma_e2;
        double psi_k = NA_REAL;
        if (R_finite(sigma_e2) && R_finite(denom) && denom > 0.0) {
          psi_k = 2.0 * sigma_e2 / denom;
        }
        est[est_index(a, b, k, n_methods, n_times)] = psi_k;
        est[est_index(b, a, k, n_methods, n_times)] = psi_k;
        method_time_diff[est_index(a, b, k, n_methods, n_times)] = d_k;
        method_time_diff[est_index(b, a, k, n_methods, n_times)] = -d_k;
      }

      double repeat = NA_REAL;
      if (R_finite(sigma_e2) && sigma_e2 >= 0.0) {
        repeat = 1.96 * std::sqrt(2.0 * sigma_e2);
      }

      sigma2_error(a, b) = sigma2_error(b, a) = sigma_e2;
      sigma2_subject_method(a, b) = sigma2_subject_method(b, a) = sigma_ab2;
      repeatability(a, b) = repeatability(b, a) = repeat;
      ms_subject_method(a, b) = ms_subject_method(b, a) = ms_sm;
      ms_method_time(a, b) = ms_method_time(b, a) = ms_mt;
      ms_error(a, b) = ms_error(b, a) = ms_err;
      ss_subject_method(a, b) = ss_subject_method(b, a) = ss_sm;
      ss_method_time(a, b) = ss_method_time(b, a) = ss_mt;
      ss_error(a, b) = ss_error(b, a) = ss_e;
      df_subject_method(a, b) = df_subject_method(b, a) = df_sm;
      df_method_time(a, b) = df_method_time(b, a) = df_mt;
      df_error(a, b) = df_error(b, a) = df_err;
      homogeneity_F(a, b) = homogeneity_F(b, a) = f_stat;
      homogeneity_p(a, b) = homogeneity_p(b, a) = p_val;
      common(a, b) = common(b, a) = common_est;
      common_sigma2_error(a, b) = common_sigma2_error(b, a) = common_sigma_e2;
      common_sigma2_subject_method(a, b) = common_sigma2_subject_method(b, a) = common_sigma_ab2;
      common_repeatability(a, b) = common_repeatability(b, a) = common_repeat;
    }
  }

  return Rcpp::List::create(
    Rcpp::_["est"] = est,
    Rcpp::_["sigma2_error"] = sigma2_error,
    Rcpp::_["sigma2_subject_method"] = sigma2_subject_method,
    Rcpp::_["repeatability"] = repeatability,
    Rcpp::_["method_time_diff"] = method_time_diff,
    Rcpp::_["ms_subject_method"] = ms_subject_method,
    Rcpp::_["ms_method_time"] = ms_method_time,
    Rcpp::_["ms_error"] = ms_error,
    Rcpp::_["ss_subject_method"] = ss_subject_method,
    Rcpp::_["ss_method_time"] = ss_method_time,
    Rcpp::_["ss_error"] = ss_error,
    Rcpp::_["df_subject_method"] = df_subject_method,
    Rcpp::_["df_method_time"] = df_method_time,
    Rcpp::_["df_error"] = df_error,
    Rcpp::_["homogeneity_F"] = homogeneity_F,
    Rcpp::_["homogeneity_p"] = homogeneity_p,
    Rcpp::_["common"] = common,
    Rcpp::_["common_sigma2_error"] = common_sigma2_error,
    Rcpp::_["common_sigma2_subject_method"] = common_sigma2_subject_method,
    Rcpp::_["common_repeatability"] = common_repeatability,
    Rcpp::_["overall"] = overall_stats.overall,
    Rcpp::_["overall_common"] = overall_stats.overall_common,
    Rcpp::_["overall_sigma2_error"] = overall_stats.sigma2_error,
    Rcpp::_["overall_sigma2_subject_method"] = overall_stats.sigma2_subject_method,
    Rcpp::_["overall_repeatability"] = overall_stats.repeatability,
    Rcpp::_["overall_ms_subject_method"] = overall_stats.ms_subject_method,
    Rcpp::_["overall_ms_method_time"] = overall_stats.ms_method_time,
    Rcpp::_["overall_ms_error"] = overall_stats.ms_error,
    Rcpp::_["overall_ss_subject_method"] = overall_stats.ss_subject_method,
    Rcpp::_["overall_ss_method_time"] = overall_stats.ss_method_time,
    Rcpp::_["overall_ss_error"] = overall_stats.ss_error,
    Rcpp::_["overall_df_subject_method"] = overall_stats.df_subject_method,
    Rcpp::_["overall_df_method_time"] = overall_stats.df_method_time,
    Rcpp::_["overall_df_error"] = overall_stats.df_error,
    Rcpp::_["overall_homogeneity_F"] = overall_stats.homogeneity_F,
    Rcpp::_["overall_homogeneity_p"] = overall_stats.homogeneity_p,
    Rcpp::_["overall_common_sigma2_error"] = overall_stats.common_sigma2_error,
    Rcpp::_["overall_common_sigma2_subject_method"] = overall_stats.common_sigma2_subject_method,
    Rcpp::_["overall_common_repeatability"] = overall_stats.common_repeatability
  );
}

// [[Rcpp::export]]
Rcpp::List cia_rm_delta_cpp(Rcpp::NumericVector y,
                            Rcpp::IntegerVector subject,
                            Rcpp::IntegerVector method,
                            Rcpp::IntegerVector time,
                            int n_subjects,
                            int n_methods,
                            int n_times,
                            bool homogeneous,
                            bool constrain_vc,
                            double conf_level,
                            int n_threads = 1) {
  (void)n_threads;

  if (n_subjects < 2 || n_methods < 2 || n_times < 2) {
    Rcpp::stop(kCiaRmBalancedError);
  }

  const std::vector<double> cube =
    build_cube(y, subject, method, time, n_subjects, n_methods, n_times);
  const double zcrit =
    R::qnorm5(1.0 - (1.0 - conf_level) / 2.0, 0.0, 1.0, true, false);

  Rcpp::NumericVector se(
    static_cast<R_xlen_t>(n_methods) * n_methods * n_times,
    NA_REAL
  );
  Rcpp::NumericVector lwr(
    static_cast<R_xlen_t>(n_methods) * n_methods * n_times,
    NA_REAL
  );
  Rcpp::NumericVector upr(
    static_cast<R_xlen_t>(n_methods) * n_methods * n_times,
    NA_REAL
  );
  se.attr("dim") = Rcpp::IntegerVector::create(n_methods, n_methods, n_times);
  lwr.attr("dim") = Rcpp::IntegerVector::create(n_methods, n_methods, n_times);
  upr.attr("dim") = Rcpp::IntegerVector::create(n_methods, n_methods, n_times);

  Rcpp::NumericMatrix common_se(n_methods, n_methods);
  Rcpp::NumericMatrix common_lwr(n_methods, n_methods);
  Rcpp::NumericMatrix common_upr(n_methods, n_methods);
  std::fill(common_se.begin(), common_se.end(), NA_REAL);
  std::fill(common_lwr.begin(), common_lwr.end(), NA_REAL);
  std::fill(common_upr.begin(), common_upr.end(), NA_REAL);

  for (int j = 0; j < n_methods; ++j) {
    common_se(j, j) = 0.0;
    common_lwr(j, j) = 1.0;
    common_upr(j, j) = 1.0;
    for (int k = 0; k < n_times; ++k) {
      se[est_index(j, j, k, n_methods, n_times)] = 0.0;
      lwr[est_index(j, j, k, n_methods, n_times)] = 1.0;
      upr[est_index(j, j, k, n_methods, n_times)] = 1.0;
    }
  }

  for (int a = 0; a < n_methods - 1; ++a) {
    for (int b = a + 1; b < n_methods; ++b) {
      std::vector<double> D_rows(
        static_cast<std::size_t>(n_subjects * n_times), 0.0
      );
      std::vector<double> d(static_cast<std::size_t>(n_times), 0.0);
      std::vector<double> M(
        static_cast<std::size_t>(n_times * n_times), 0.0
      );

      for (int s = 0; s < n_subjects; ++s) {
        for (int k = 0; k < n_times; ++k) {
          const double value =
            cube[static_cast<std::size_t>(cube_index(s, a, k, n_methods, n_times))] -
            cube[static_cast<std::size_t>(cube_index(s, b, k, n_methods, n_times))];
          D_rows[static_cast<std::size_t>(s * n_times + k)] = value;
          d[static_cast<std::size_t>(k)] += value;
        }
      }

      for (int k = 0; k < n_times; ++k) {
        d[static_cast<std::size_t>(k)] /= static_cast<double>(n_subjects);
      }
      for (int s = 0; s < n_subjects; ++s) {
        const double* row_ptr =
          D_rows.data() + static_cast<std::size_t>(s * n_times);
        for (int r = 0; r < n_times; ++r) {
          for (int c = 0; c < n_times; ++c) {
            M[static_cast<std::size_t>(r + n_times * c)] +=
              row_ptr[r] * row_ptr[c];
          }
        }
      }
      for (double& value : M) {
        value /= static_cast<double>(n_subjects);
      }

      std::vector<double> m;
      m.reserve(static_cast<std::size_t>(n_times + n_times * n_times));
      m.insert(m.end(), d.begin(), d.end());
      m.insert(m.end(), M.begin(), M.end());

      const PairMomentPayload center =
        pair_payload_from_moments(
          m,
          n_times,
          n_subjects,
          constrain_vc,
          homogeneous
        );

      const std::vector<double> grad_est = numerical_gradient_matrix(
        [&](const std::vector<double>& mm) {
          return pair_payload_from_moments(
            mm,
            n_times,
            n_subjects,
            constrain_vc,
            homogeneous
          ).est;
        },
        m,
        n_times
      );

      for (int k = 0; k < n_times; ++k) {
        std::vector<double> grad_k(
          grad_est.begin() + static_cast<std::size_t>(k * m.size()),
          grad_est.begin() + static_cast<std::size_t>((k + 1) * m.size())
        );
        const double var_psi =
          projection_variance_quadratic(D_rows, n_subjects, n_times, grad_k);
        const double se_k =
          R_finite(var_psi) ? std::sqrt(std::max(var_psi, 0.0)) : NA_REAL;
        double lwr_k =
          center.est[static_cast<std::size_t>(k)] - zcrit * se_k;
        double upr_k =
          center.est[static_cast<std::size_t>(k)] + zcrit * se_k;
        apply_bounds_scalar(lwr_k, upr_k, constrain_vc);
        se[est_index(a, b, k, n_methods, n_times)] =
          se[est_index(b, a, k, n_methods, n_times)] = se_k;
        lwr[est_index(a, b, k, n_methods, n_times)] =
          lwr[est_index(b, a, k, n_methods, n_times)] = lwr_k;
        upr[est_index(a, b, k, n_methods, n_times)] =
          upr[est_index(b, a, k, n_methods, n_times)] = upr_k;
      }

      if (homogeneous && R_finite(center.common)) {
        const std::vector<double> grad_common = numerical_gradient_matrix(
          [&](const std::vector<double>& mm) {
            return std::vector<double>(
              1,
              pair_payload_from_moments(
                mm,
                n_times,
                n_subjects,
                constrain_vc,
                true
              ).common
            );
          },
          m,
          1
        );
        const double var_common =
          projection_variance_quadratic(D_rows, n_subjects, n_times, grad_common);
        const double se_common =
          R_finite(var_common) ? std::sqrt(std::max(var_common, 0.0)) : NA_REAL;
        double lwr_common = center.common - zcrit * se_common;
        double upr_common = center.common + zcrit * se_common;
        apply_bounds_scalar(lwr_common, upr_common, constrain_vc);
        common_se(a, b) = common_se(b, a) = se_common;
        common_lwr(a, b) = common_lwr(b, a) = lwr_common;
        common_upr(a, b) = common_upr(b, a) = upr_common;
      }
    }
  }

  Rcpp::NumericVector overall_se(n_times, NA_REAL);
  Rcpp::NumericVector overall_lwr(n_times, NA_REAL);
  Rcpp::NumericVector overall_upr(n_times, NA_REAL);
  double overall_common_se = NA_REAL;
  double overall_common_lwr = NA_REAL;
  double overall_common_upr = NA_REAL;

  if (n_methods >= 3) {
    const int P = n_methods * n_times;
    std::vector<double> X_rows(
      static_cast<std::size_t>(n_subjects * P), 0.0
    );
    std::vector<double> mu(static_cast<std::size_t>(P), 0.0);
    std::vector<double> M(static_cast<std::size_t>(P * P), 0.0);

    for (int s = 0; s < n_subjects; ++s) {
      for (int k = 0; k < n_times; ++k) {
        for (int j = 0; j < n_methods; ++j) {
          const int idx = j + n_methods * k;
          const double value =
            cube[static_cast<std::size_t>(cube_index(s, j, k, n_methods, n_times))];
          X_rows[static_cast<std::size_t>(s * P + idx)] = value;
          mu[static_cast<std::size_t>(idx)] += value;
        }
      }
    }
    for (double& value : mu) {
      value /= static_cast<double>(n_subjects);
    }
    for (int s = 0; s < n_subjects; ++s) {
      const double* row_ptr =
        X_rows.data() + static_cast<std::size_t>(s * P);
      for (int r = 0; r < P; ++r) {
        for (int c = 0; c < P; ++c) {
          M[static_cast<std::size_t>(r + P * c)] += row_ptr[r] * row_ptr[c];
        }
      }
    }
    for (double& value : M) {
      value /= static_cast<double>(n_subjects);
    }

    std::vector<double> m;
    m.reserve(static_cast<std::size_t>(P + P * P));
    m.insert(m.end(), mu.begin(), mu.end());
    m.insert(m.end(), M.begin(), M.end());

    const OverallMomentPayload center =
      overall_payload_from_moments(
        m,
        n_methods,
        n_times,
        n_subjects,
        homogeneous,
        constrain_vc
      );

    const std::vector<double> grad_overall = numerical_gradient_matrix(
      [&](const std::vector<double>& mm) {
        return overall_payload_from_moments(
          mm,
          n_methods,
          n_times,
          n_subjects,
          homogeneous,
          constrain_vc
        ).overall;
      },
      m,
      n_times
    );

    for (int k = 0; k < n_times; ++k) {
      std::vector<double> grad_k(
        grad_overall.begin() + static_cast<std::size_t>(k * m.size()),
        grad_overall.begin() + static_cast<std::size_t>((k + 1) * m.size())
      );
      const double var_psi =
        projection_variance_quadratic(X_rows, n_subjects, P, grad_k);
      const double se_k =
        R_finite(var_psi) ? std::sqrt(std::max(var_psi, 0.0)) : NA_REAL;
      double lwr_k =
        center.overall[static_cast<std::size_t>(k)] - zcrit * se_k;
      double upr_k =
        center.overall[static_cast<std::size_t>(k)] + zcrit * se_k;
      apply_bounds_scalar(lwr_k, upr_k, constrain_vc);
      overall_se[static_cast<std::size_t>(k)] = se_k;
      overall_lwr[static_cast<std::size_t>(k)] = lwr_k;
      overall_upr[static_cast<std::size_t>(k)] = upr_k;
    }

    if (homogeneous && R_finite(center.overall_common)) {
      const std::vector<double> grad_common = numerical_gradient_matrix(
        [&](const std::vector<double>& mm) {
          return std::vector<double>(
            1,
            overall_payload_from_moments(
              mm,
              n_methods,
              n_times,
              n_subjects,
              true,
              constrain_vc
            ).overall_common
          );
        },
        m,
        1
      );
      const double var_common =
        projection_variance_quadratic(X_rows, n_subjects, P, grad_common);
      overall_common_se =
        R_finite(var_common) ? std::sqrt(std::max(var_common, 0.0)) : NA_REAL;
      overall_common_lwr = center.overall_common - zcrit * overall_common_se;
      overall_common_upr = center.overall_common + zcrit * overall_common_se;
      apply_bounds_scalar(overall_common_lwr, overall_common_upr, constrain_vc);
    }
  }

  return Rcpp::List::create(
    Rcpp::_["se"] = se,
    Rcpp::_["lwr"] = lwr,
    Rcpp::_["upr"] = upr,
    Rcpp::_["common_se"] = common_se,
    Rcpp::_["common_lwr"] = common_lwr,
    Rcpp::_["common_upr"] = common_upr,
    Rcpp::_["overall_se"] = overall_se,
    Rcpp::_["overall_lwr"] = overall_lwr,
    Rcpp::_["overall_upr"] = overall_upr,
    Rcpp::_["overall_common_se"] = overall_common_se,
    Rcpp::_["overall_common_lwr"] = overall_common_lwr,
    Rcpp::_["overall_common_upr"] = overall_common_upr
  );
}
