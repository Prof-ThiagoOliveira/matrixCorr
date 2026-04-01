// Thiago de Paula Oliveira
// Repeated-measures Bland-Altman via stabilized EM/GLS with optional AR(1)
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>

// bring helpers
#include "matrixCorr_detail.h"

using namespace Rcpp;
using namespace arma;

// ---- use selected helpers from matrixCorr_detail ----
using matrixCorr_detail::clamp_policy::nan_preserve;
using matrixCorr_detail::linalg::inv_sympd_safe;
using matrixCorr_detail::linalg::logdet_spd_safe;
using matrixCorr_detail::linalg::solve_sympd_safe;
using matrixCorr_detail::indexing::reindex;

namespace {
constexpr const char* BA_RM_IDENTIFIABILITY_ERROR =
  "The repeated-measures Bland-Altman mixed model is not separately identifiable on the supplied paired data because there is insufficient within-subject replication after pairing to separate the residual and subject-level variance components.";
constexpr const char* BA_RM_SLOPE_SCALE_ERROR =
  "The proportional-bias slope is not estimable because the paired means are near-degenerate on the observed data scale.";
constexpr const char* BA_RM_CONVERGENCE_ERROR =
  "The repeated-measures Bland-Altman model is conceptually estimable on the supplied paired data, but the EM/GLS algorithm failed to converge to admissible finite variance-component estimates.";
constexpr double BA_RM_AR1_RIDGE = 1e-10;
} // namespace

static inline bool all_finite(const arma::vec& v) {
  for (arma::uword i = 0; i < v.n_elem; ++i) if (!std::isfinite(v[i])) return false;
  return true;
}

struct BaRMSlopeScaleInfo {
  double s_sd = NA_REAL;
  double s_iqr = NA_REAL;
  double s_mad = NA_REAL;
  double s_ref = NA_REAL;
  double s_m_star = NA_REAL;
  double tau = std::sqrt(std::numeric_limits<double>::epsilon());
  std::string source;
};

static inline bool ba_rm_is_finite_positive(double x) {
  return std::isfinite(x) && x > 0.0;
}

static std::vector<double> ba_rm_collect_finite_means(const std::vector<double>& mean_values) {
  std::vector<double> sorted;
  sorted.reserve(mean_values.size());
  for (double value : mean_values) {
    if (std::isfinite(value)) sorted.push_back(value);
  }
  return sorted;
}

static double ba_rm_iqr_scale(std::vector<double> sorted) {
  if (sorted.empty()) return NA_REAL;
  std::sort(sorted.begin(), sorted.end());
  const double q1 = matrixCorr_detail::standardise_bicor::quantile_sorted_std(sorted, 0.25);
  const double q3 = matrixCorr_detail::standardise_bicor::quantile_sorted_std(sorted, 0.75);
  return (q3 - q1) / 1.349;
}

static double ba_rm_mad_scale(const std::vector<double>& sorted) {
  if (sorted.empty()) return NA_REAL;
  std::vector<double> med_work = sorted;
  const double med = matrixCorr_detail::standardise_bicor::median_inplace(med_work);
  if (!std::isfinite(med)) return NA_REAL;

  std::vector<double> dev = sorted;
  for (double& value : dev) value = std::fabs(value - med);
  return 1.4826 * matrixCorr_detail::standardise_bicor::median_inplace(dev);
}

static bool ba_rm_is_near_zero(double scale, double s_ref, double tau) {
  return !ba_rm_is_finite_positive(scale) ||
    !ba_rm_is_finite_positive(s_ref) ||
    scale <= tau * s_ref;
}

static BaRMSlopeScaleInfo ba_rm_choose_slope_scale(const std::vector<double>& mean_values,
                                                   double tau = std::sqrt(std::numeric_limits<double>::epsilon())) {
  BaRMSlopeScaleInfo info;
  info.tau = tau;

  std::vector<double> sorted = ba_rm_collect_finite_means(mean_values);
  if (sorted.size() >= 2u) {
    arma::vec mean_vec(sorted.size(), arma::fill::zeros);
    for (arma::uword i = 0; i < mean_vec.n_elem; ++i) mean_vec[i] = sorted[i];
    info.s_sd = arma::stddev(mean_vec);
  }

  info.s_iqr = ba_rm_iqr_scale(sorted);
  info.s_mad = ba_rm_mad_scale(sorted);

  const double scales[] = {info.s_sd, info.s_iqr, info.s_mad};
  for (double scale : scales) {
    if (!ba_rm_is_finite_positive(scale)) continue;
    info.s_ref = ba_rm_is_finite_positive(info.s_ref) ? std::max(info.s_ref, scale) : scale;
  }

  if (!ba_rm_is_near_zero(info.s_sd, info.s_ref, info.tau)) {
    info.s_m_star = info.s_sd;
    info.source = "sd";
    return info;
  }

  if (!ba_rm_is_near_zero(info.s_iqr, info.s_ref, info.tau)) {
    info.s_m_star = info.s_iqr;
    info.source = "iqr";
    return info;
  }

  if (!ba_rm_is_near_zero(info.s_mad, info.s_ref, info.tau)) {
    info.s_m_star = info.s_mad;
    info.source = "mad";
    return info;
  }

  stop(BA_RM_SLOPE_SCALE_ERROR);
}

// [[Rcpp::export]]
Rcpp::List ba_rm_slope_scale_cpp(Rcpp::NumericVector mean_values) {
  std::vector<double> mean_vec(mean_values.begin(), mean_values.end());
  BaRMSlopeScaleInfo info = ba_rm_choose_slope_scale(mean_vec);
  return Rcpp::List::create(
    _["s_sd"] = info.s_sd,
    _["s_iqr"] = info.s_iqr,
    _["s_mad"] = info.s_mad,
    _["s_ref"] = info.s_ref,
    _["s_m_star"] = info.s_m_star,
    _["tau"] = info.tau,
    _["source"] = info.source
  );
}

static double ba_rm_slope_scale_denom(const std::vector<double>& mean_values,
                                      double tau = std::sqrt(std::numeric_limits<double>::epsilon())) {
  return ba_rm_choose_slope_scale(mean_values, tau).s_m_star;
}

struct PairData {
  std::vector<double> d;
  std::vector<double> mean;
  std::vector<int> subj;
  std::vector<int> time;
};

static inline uint64_t make_pair_key(int subj, int time) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(subj)) << 32) | static_cast<uint32_t>(time);
}

static PairData make_pairs(const NumericVector& y,
                           const IntegerVector& subject,
                           const IntegerVector& method,
                           const IntegerVector& time) {
  const int n = y.size();
  if (subject.size() != n || method.size() != n || time.size() != n)
    stop("lengths of y, subject, method, time must match.");

  struct V {
    bool has1 = false;
    bool has2 = false;
    double y1 = NA_REAL;
    double y2 = NA_REAL;
    int s = 0;
    int t = 0;
  };

  std::unordered_map<uint64_t, V> H;
  H.reserve(static_cast<size_t>(n));

  for (int i = 0; i < n; ++i) {
    if (IntegerVector::is_na(subject[i]) || IntegerVector::is_na(method[i]) || IntegerVector::is_na(time[i])) continue;
    if (NumericVector::is_na(y[i])) continue;
    const int s = subject[i];
    const int t = time[i];
    const int m = method[i];
    if (m != 1 && m != 2) continue;

    const uint64_t key = make_pair_key(s, t);
    V& v = H[key];
    v.s = s;
    v.t = t;
    if (m == 1) {
      v.has1 = true;
      v.y1 = y[i];
    } else {
      v.has2 = true;
      v.y2 = y[i];
    }
  }

  PairData P;
  P.d.reserve(H.size());
  P.mean.reserve(H.size());
  P.subj.reserve(H.size());
  P.time.reserve(H.size());

  for (auto& kv : H) {
    const V& v = kv.second;
    if (v.has1 && v.has2 && std::isfinite(v.y1) && std::isfinite(v.y2)) {
      P.d.push_back(v.y2 - v.y1);
      P.mean.push_back(0.5 * (v.y1 + v.y2));
      P.subj.push_back(v.s);
      P.time.push_back(v.t);
    }
  }

  if (P.d.empty()) stop("No complete subject-time pairs (both methods present).");
  return P;
}

// [[Rcpp::export]]
int ba_rm_complete_pairs_cpp(const NumericVector& y,
                             const IntegerVector& subject,
                             const IntegerVector& method,
                             const IntegerVector& time) {
  const int n = y.size();
  if (subject.size() != n || method.size() != n || time.size() != n) {
    stop("lengths of y, subject, method, time must match.");
  }

  struct PairSeen {
    bool has1 = false;
    bool has2 = false;
  };

  std::unordered_map<uint64_t, PairSeen> seen;
  seen.reserve(static_cast<size_t>(n));

  for (int i = 0; i < n; ++i) {
    if (NumericVector::is_na(y[i]) ||
        IntegerVector::is_na(subject[i]) ||
        IntegerVector::is_na(method[i]) ||
        IntegerVector::is_na(time[i])) {
      continue;
    }

    const int m = method[i];
    if (m != 1 && m != 2) continue;

    PairSeen& entry = seen[make_pair_key(subject[i], time[i])];
    if (m == 1) entry.has1 = true;
    else        entry.has2 = true;
  }

  int n_pairs = 0;
  for (const auto& kv : seen) {
    if (kv.second.has1 && kv.second.has2) ++n_pairs;
  }

  return n_pairs;
}

struct BlockSegment {
  int start = 0;
  int len = 0;
};

struct BlockData {
  int n_i = 0;
  std::vector<double> y;
  std::vector<double> x;
  std::vector<int> time;
  std::vector<BlockSegment> segments;
};

static std::vector<BlockData> build_subject_blocks(const PairData& P,
                                                   const std::vector<int>& subj_idx,
                                                   int m,
                                                   bool include_slope,
                                                   double x_mean,
                                                   double x_scale) {
  std::vector<std::vector<int> > rows(m);
  std::vector<std::vector<int> > tim(m);
  for (size_t i = 0; i < subj_idx.size(); ++i) {
    rows[subj_idx[i]].push_back(static_cast<int>(i));
    tim[subj_idx[i]].push_back(P.time[i]);
  }

  std::vector<BlockData> blocks(m);
  for (int i = 0; i < m; ++i) {
    auto& r = rows[i];
    auto& t = tim[i];
    if (r.empty()) continue;

    std::vector<int> ord(r.size());
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) {
      const int ta = t[a];
      const int tb = t[b];
      if (ta < 0 && tb < 0) return a < b;
      if (ta < 0) return false;
      if (tb < 0) return true;
      return ta < tb;
    });

    BlockData B;
    B.n_i = static_cast<int>(r.size());
    B.y.resize(B.n_i);
    B.time.resize(B.n_i);
    if (include_slope) B.x.resize(B.n_i);

    for (int k = 0; k < B.n_i; ++k) {
      const int idx = r[ord[k]];
      B.y[k] = P.d[idx];
      B.time[k] = P.time[idx];
      if (include_slope) B.x[k] = (P.mean[idx] - x_mean) / x_scale;
    }

    int s = 0;
    while (s < B.n_i) {
      if (B.time[s] < 0) {
        B.segments.push_back(BlockSegment{s, 1});
        ++s;
        continue;
      }
      int e = s;
      while (e + 1 < B.n_i && B.time[e + 1] >= 0 && B.time[e + 1] == B.time[e] + 1) ++e;
      B.segments.push_back(BlockSegment{s, e - s + 1});
      s = e + 1;
    }

    blocks[i] = std::move(B);
  }

  return blocks;
}

struct BlockStats {
  int n_i = 0;
  double UCU = 0.0;
  double UTCy = 0.0;
  double UTCx = 0.0;
  double yCy = 0.0;
  double xCy = 0.0;
  double xCx = 0.0;
  double logdetCinv = 0.0;
};

static inline bool ba_rm_at_su2_boundary(double su2, double tol = 1e-10) {
  return std::isfinite(su2) && su2 <= tol;
}

static bool ba_rm_accumulate_segment_ar1(const BlockData& B,
                                         const BlockSegment& seg,
                                         bool include_slope,
                                         double rho,
                                         BlockStats& out) {
  const int start = seg.start;
  const int len = seg.len;
  if (len <= 0) return true;

  const double r2 = rho * rho;
  const double denom = std::max(1.0 - r2, 1e-12);
  const double diag_end = (len == 1 ? 1.0 : 1.0 / denom) + BA_RM_AR1_RIDGE;
  const double diag_mid = (1.0 + r2) / denom + BA_RM_AR1_RIDGE;
  const double off = (len >= 2 ? (-rho / denom) : 0.0);

  double delta_prev = 0.0;
  for (int j = 0; j < len; ++j) {
    const int idx = start + j;
    const double yj = B.y[idx];
    const bool endpoint = (j == 0 || j == len - 1);
    const double dj = (len == 1 || endpoint) ? diag_end : diag_mid;

    out.UCU += dj;
    out.UTCy += dj * yj;
    out.yCy += dj * yj * yj;

    if (include_slope) {
      const double xj = B.x[idx];
      out.UTCx += dj * xj;
      out.xCy += dj * xj * yj;
      out.xCx += dj * xj * xj;
    }

    double delta = dj;
    if (j > 0) delta -= (off * off) / delta_prev;
    if (!std::isfinite(delta) || delta <= 0.0) return false;
    out.logdetCinv += std::log(delta);
    delta_prev = delta;
  }

  for (int j = 0; j < len - 1; ++j) {
    const int idx = start + j;
    const double yj = B.y[idx];
    const double yk = B.y[idx + 1];
    out.UCU += 2.0 * off;
    out.UTCy += off * (yj + yk);
    out.yCy += 2.0 * off * yj * yk;

    if (include_slope) {
      const double xj = B.x[idx];
      const double xk = B.x[idx + 1];
      out.UTCx += off * (xj + xk);
      out.xCy += off * (xj * yk + xk * yj);
      out.xCx += 2.0 * off * xj * xk;
    }
  }

  return std::isfinite(out.UCU) &&
    std::isfinite(out.UTCy) &&
    std::isfinite(out.yCy) &&
    std::isfinite(out.logdetCinv) &&
    (!include_slope || (std::isfinite(out.UTCx) && std::isfinite(out.xCy) && std::isfinite(out.xCx)));
}

static bool ba_rm_compute_block_stats(const BlockData& B,
                                      bool include_slope,
                                      bool use_ar1,
                                      double rho,
                                      BlockStats& out) {
  out = BlockStats();
  out.n_i = B.n_i;
  if (B.n_i == 0) return true;

  if (!use_ar1) {
    out.UCU = static_cast<double>(B.n_i);
    for (int j = 0; j < B.n_i; ++j) {
      const double yj = B.y[j];
      out.UTCy += yj;
      out.yCy += yj * yj;
      if (include_slope) {
        const double xj = B.x[j];
        out.UTCx += xj;
        out.xCy += xj * yj;
        out.xCx += xj * xj;
      }
    }
    return true;
  }

  for (const BlockSegment& seg : B.segments) {
    if (!ba_rm_accumulate_segment_ar1(B, seg, include_slope, rho, out)) return false;
  }
  return true;
}

static bool ba_rm_compute_all_block_stats(const std::vector<BlockData>& blocks,
                                          bool include_slope,
                                          bool use_ar1,
                                          double rho,
                                          std::vector<BlockStats>& out) {
  out.resize(blocks.size());
  for (size_t i = 0; i < blocks.size(); ++i) {
    if (!ba_rm_compute_block_stats(blocks[i], include_slope, use_ar1, rho, out[i])) return false;
  }
  return true;
}

struct RhoStatsCacheEntry {
  double rho = NA_REAL;
  std::vector<BlockStats> stats;
};

struct RhoStatsCache {
  const std::vector<BlockData>& blocks;
  bool include_slope = false;
  bool use_ar1 = false;
  std::vector<BlockStats> iid_stats;
  bool iid_ready = false;
  std::vector<RhoStatsCacheEntry> entries;

  explicit RhoStatsCache(const std::vector<BlockData>& blocks_,
                         bool include_slope_,
                         bool use_ar1_)
    : blocks(blocks_), include_slope(include_slope_), use_ar1(use_ar1_) {}

  bool get(double rho, const std::vector<BlockStats>*& stats_ptr) {
    if (!use_ar1) {
      if (!iid_ready) {
        if (!ba_rm_compute_all_block_stats(blocks, include_slope, false, 0.0, iid_stats)) return false;
        iid_ready = true;
      }
      stats_ptr = &iid_stats;
      return true;
    }

    for (RhoStatsCacheEntry& entry : entries) {
      if (entry.rho == rho) {
        stats_ptr = &entry.stats;
        return true;
      }
    }

    RhoStatsCacheEntry entry;
    entry.rho = rho;
    if (!ba_rm_compute_all_block_stats(blocks, include_slope, true, rho, entry.stats)) return false;
    entries.push_back(std::move(entry));
    stats_ptr = &entries.back().stats;
    return true;
  }
};

struct FitOut {
  arma::vec beta;
  arma::mat beta_vcov;
  double su2 = NA_REAL;
  double se2 = NA_REAL;
  int iter = 0;
  bool converged = false;
  bool boundary_su2_zero = false;
  std::string warn;
};

struct GLSAccumulator {
  int p = 1;
  double A00 = 0.0;
  double A01 = 0.0;
  double A11 = 0.0;
  double b0 = 0.0;
  double b1 = 0.0;
};

static inline void ba_rm_accumulate_gls_block(const BlockStats& B,
                                              int p,
                                              double su2,
                                              double se2,
                                              GLSAccumulator& acc) {
  const double eps = 1e-12;
  const double inv_se = 1.0 / std::max(se2, eps);
  const bool at_su0 = ba_rm_at_su2_boundary(su2);

  const double s0 = inv_se * B.UCU;
  const double sy = inv_se * B.UTCy;
  double Minv = 0.0;
  if (!at_su0) {
    const double inv_su = 1.0 / su2;
    const double M = std::max(inv_su + s0, 1e-12);
    Minv = 1.0 / M;
  }

  acc.A00 += s0 - s0 * Minv * s0;
  acc.b0 += sy - s0 * Minv * sy;

  if (p == 2) {
    const double s1 = inv_se * B.UTCx;
    acc.A01 += inv_se * B.UTCx - s0 * Minv * s1;
    acc.A11 += inv_se * B.xCx - s1 * Minv * s1;
    acc.b1 += inv_se * B.xCy - s1 * Minv * sy;
  }
}

static bool ba_rm_solve_gls_system(const GLSAccumulator& acc,
                                   arma::vec& beta,
                                   arma::mat& beta_vcov,
                                   double* logdet_xtvix = nullptr) {
  if (acc.p == 1) {
    if (!std::isfinite(acc.A00) || acc.A00 <= 0.0 || !std::isfinite(acc.b0)) return false;
    beta.set_size(1);
    beta[0] = acc.b0 / acc.A00;
    beta_vcov.set_size(1, 1);
    beta_vcov(0, 0) = 1.0 / acc.A00;
    if (logdet_xtvix) *logdet_xtvix = std::log(acc.A00);
    return std::isfinite(beta[0]) && std::isfinite(beta_vcov(0, 0));
  }

  arma::mat A(2, 2, arma::fill::zeros);
  A(0, 0) = acc.A00;
  A(0, 1) = acc.A01;
  A(1, 0) = acc.A01;
  A(1, 1) = acc.A11;
  arma::mat Ainv;
  if (!inv_sympd_safe(Ainv, A) || !Ainv.is_finite()) return false;

  arma::vec b(2, arma::fill::zeros);
  b[0] = acc.b0;
  b[1] = acc.b1;

  beta = Ainv * b;
  beta_vcov = Ainv;
  if (logdet_xtvix) *logdet_xtvix = logdet_spd_safe(A);
  return beta.is_finite() && beta_vcov.is_finite();
}

static inline double ba_rm_dot_xtviy_beta(const GLSAccumulator& acc,
                                          const arma::vec& beta) {
  if (beta.n_elem == 1) return acc.b0 * beta[0];
  return acc.b0 * beta[0] + acc.b1 * beta[1];
}

static inline void ba_rm_residual_quadratics(const BlockStats& B,
                                             const arma::vec& beta,
                                             double& UTCr,
                                             double& q) {
  const double b0 = beta[0];
  UTCr = B.UTCy - b0 * B.UCU;
  q = B.yCy - 2.0 * b0 * B.UTCy + b0 * b0 * B.UCU;

  if (beta.n_elem == 2) {
    const double b1 = beta[1];
    UTCr -= b1 * B.UTCx;
    q += -2.0 * b1 * B.xCy + 2.0 * b0 * b1 * B.UTCx + b1 * b1 * B.xCx;
  }
}

static FitOut fit_diff_boundary_su0(int p,
                                    const std::vector<BlockStats>& stats) {
  GLSAccumulator acc;
  acc.p = p;
  double yCy = 0.0;
  int n = 0;

  for (const BlockStats& B : stats) {
    if (B.n_i == 0) continue;
    acc.A00 += B.UCU;
    acc.b0 += B.UTCy;
    if (p == 2) {
      acc.A01 += B.UTCx;
      acc.A11 += B.xCx;
      acc.b1 += B.xCy;
    }
    yCy += B.yCy;
    n += B.n_i;
  }

  arma::vec beta;
  arma::mat XtCX_inv;
  if (!ba_rm_solve_gls_system(acc, beta, XtCX_inv)) {
    stop("Boundary fit with sigma2_subject = 0 failed because X'CX was not invertible.");
  }

  const double yPy = yCy - ba_rm_dot_xtviy_beta(acc, beta);
  const double se2 = std::max(1e-12, yPy / std::max(1, n - p));

  FitOut out;
  out.beta = beta;
  out.beta_vcov = se2 * XtCX_inv;
  out.su2 = 0.0;
  out.se2 = se2;
  out.iter = 0;
  out.converged = true;
  out.boundary_su2_zero = true;
  out.warn = "Constrained boundary fit used with sigma2_subject fixed at 0.";
  return out;
}

static double ba_rm_variance_all_y(const std::vector<BlockData>& blocks) {
  double sum = 0.0;
  double sumsq = 0.0;
  int n = 0;
  for (const BlockData& B : blocks) {
    for (double value : B.y) {
      sum += value;
      sumsq += value * value;
      ++n;
    }
  }
  if (n <= 1) return 0.0;
  const double mean = sum / static_cast<double>(n);
  return std::max(0.0, (sumsq - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1));
}

static FitOut fit_diff_em(int p,
                          const std::vector<BlockData>& blocks,
                          const std::vector<BlockStats>& stats,
                          int max_iter,
                          double tol) {
  const int m = static_cast<int>(blocks.size());
  int n = 0;
  for (const BlockData& B : blocks) n += B.n_i;

  const double vref = ba_rm_variance_all_y(blocks);
  const double EPS = std::max(1e-12, vref * 1e-12);
  const double MAXV = std::max(10.0 * vref, 1.0);

  std::vector<double> mu_s(m, 0.0);
  std::vector<int> cnt_s(m, 0);
  int used_mu = 0;
  int n_replicated_subjects = 0;
  double cnt_sum = 0.0;

  for (int i = 0; i < m; ++i) {
    cnt_s[i] = blocks[i].n_i;
    if (cnt_s[i] <= 0) continue;
    mu_s[i] = std::accumulate(blocks[i].y.begin(), blocks[i].y.end(), 0.0) / static_cast<double>(cnt_s[i]);
    ++used_mu;
    cnt_sum += static_cast<double>(cnt_s[i]);
    if (cnt_s[i] >= 2) ++n_replicated_subjects;
  }

  if (used_mu < 2 || n_replicated_subjects < 1) {
    stop(BA_RM_IDENTIFIABILITY_ERROR);
  }

  arma::vec mu_s_vec(used_mu, arma::fill::zeros);
  {
    int k = 0;
    for (int i = 0; i < m; ++i) if (cnt_s[i] > 0) mu_s_vec[k++] = mu_s[i];
  }

  const double var_mu = (used_mu >= 2 ? arma::var(mu_s_vec, /*unbiased*/ true) : 0.0);
  double num_w = 0.0;
  double den_w = 0.0;
  bool se2_init_fallback = false;

  for (const BlockData& B : blocks) {
    if (B.n_i < 2) continue;
    const double sum = std::accumulate(B.y.begin(), B.y.end(), 0.0);
    double sumsq = 0.0;
    for (double value : B.y) sumsq += value * value;
    const double mean = sum / static_cast<double>(B.n_i);
    const double vsi = (sumsq - static_cast<double>(B.n_i) * mean * mean) / static_cast<double>(B.n_i - 1);
    if (std::isfinite(vsi) && vsi > 0.0) {
      num_w += static_cast<double>(B.n_i - 1) * vsi;
      den_w += static_cast<double>(B.n_i - 1);
    }
  }

  double se2_init = NA_REAL;
  if (den_w > 0.5) {
    se2_init = num_w / den_w;
  } else {
    se2_init_fallback = true;
    se2_init = std::max(EPS, 0.5 * std::max(vref, 0.0));
  }

  double su2_init = std::max(
    0.0,
    var_mu - se2_init / std::max(1.0, cnt_sum / static_cast<double>(used_mu))
  );

  se2_init = nan_preserve(se2_init, EPS, MAXV);
  su2_init = nan_preserve(su2_init, 0.0, MAXV);

  auto damp_to_ratio = [](double oldv, double newv, double rmax) {
    if (!std::isfinite(newv)) return oldv;
    if (oldv <= 0.0) return nan_preserve(newv, 0.0, rmax);
    const double lo = oldv / rmax;
    const double hi = oldv * rmax;
    return nan_preserve(newv, std::min(lo, hi), std::max(lo, hi));
  };

  double su2 = su2_init;
  double se2 = se2_init;
  arma::vec beta(p, arma::fill::zeros);
  arma::mat beta_vcov(p, p, arma::fill::zeros);
  bool converged = false;
  int it = -1;
  std::string init_warn;
  if (se2_init_fallback) {
    init_warn = "Residual-variance initialization used 0.5 * v_ref only as a positive EM starting heuristic because no finite positive pooled within-subject variance could be formed after pairing.";
  }

  for (it = 0; it < max_iter; ++it) {
    GLSAccumulator acc;
    acc.p = p;
    for (const BlockStats& B : stats) {
      if (B.n_i == 0) continue;
      ba_rm_accumulate_gls_block(B, p, su2, se2, acc);
    }

    if (!ba_rm_solve_gls_system(acc, beta, beta_vcov) || !all_finite(beta) || !beta_vcov.is_finite()) break;

    double su_acc = 0.0;
    double se_num = 0.0;
    for (const BlockStats& B : stats) {
      if (B.n_i == 0) continue;
      const double inv_su = 1.0 / std::max(su2, EPS);
      const double M = std::max(inv_su + (1.0 / std::max(se2, EPS)) * B.UCU, 1e-12);
      const double Minv = 1.0 / M;

      double UTCr = 0.0;
      double q = 0.0;
      ba_rm_residual_quadratics(B, beta, UTCr, q);

      const double Utr = (1.0 / std::max(se2, EPS)) * UTCr;
      const double b_i = Minv * Utr;
      const double Eu2 = b_i * b_i + Minv;
      su_acc += Eu2;

      const double q_cond = q - 2.0 * b_i * UTCr + b_i * b_i * B.UCU;
      const double num_i = q_cond + B.UCU * Minv;
      se_num += num_i;
    }

    double su2_new = su_acc / std::max(1, m);
    double se2_new = se_num / std::max(1, n);
    su2_new = damp_to_ratio(su2, su2_new, 3.0);
    se2_new = damp_to_ratio(se2, se2_new, 3.0);
    su2_new = nan_preserve(su2_new, 0.0, MAXV);
    se2_new = nan_preserve(se2_new, EPS, MAXV);

    const bool admissible =
      std::isfinite(su2_new) &&
      std::isfinite(se2_new) &&
      su2_new >= 0.0 &&
      se2_new >= EPS;
    if (!admissible) break;

    const double diff = std::fabs(su2_new - su2) + std::fabs(se2_new - se2);
    su2 = su2_new;
    se2 = se2_new;
    if (diff < tol) {
      converged = true;
      break;
    }
  }

  if (!converged || !all_finite(beta) || !beta_vcov.is_finite() ||
      !std::isfinite(su2) || !std::isfinite(se2) ||
      su2 < 0.0 || se2 < EPS) {
    stop(BA_RM_CONVERGENCE_ERROR);
  }

  FitOut out;
  out.beta = beta;
  out.beta_vcov = beta_vcov;
  out.su2 = su2;
  out.se2 = se2;
  out.iter = it + 1;
  out.converged = true;
  out.warn = init_warn;
  return out;
}

static inline double loa_var_subject_equal_weight(const std::vector<BlockData>& blocks,
                                                  double mu0) {
  double acc = 0.0;
  int used = 0;
  for (const BlockData& B : blocks) {
    if (B.n_i == 0) continue;
    double ss_i = 0.0;
    for (double value : B.y) {
      const double dev = value - mu0;
      ss_i += dev * dev;
    }
    acc += ss_i / static_cast<double>(B.n_i);
    ++used;
  }
  if (used == 0) stop("No data per subject to compute LoA variance.");
  return acc / static_cast<double>(used);
}

static double estimate_rho_moments(const FitOut& fit,
                                   const std::vector<BlockData>& blocks) {
  const arma::vec beta = fit.beta;
  const double EPS = 1e-12;

  double z_sum = 0.0;
  double w_sum = 0.0;

  for (const BlockData& B : blocks) {
    for (const BlockSegment& seg : B.segments) {
      const int L = seg.len;
      if (L < 3) continue;

      arma::vec r(L, arma::fill::zeros);
      for (int k = 0; k < L; ++k) {
        const int idx = seg.start + k;
        double fitted = beta[0];
        if (beta.n_elem == 2) fitted += beta[1] * B.x[idx];
        r[k] = B.y[idx] - fitted;
      }

      arma::vec t(L, arma::fill::zeros);
      for (int k = 0; k < L; ++k) t[k] = static_cast<double>(k);
      arma::mat Z(L, 2, arma::fill::zeros);
      Z.col(0).ones();
      Z.col(1) = t - arma::mean(t);
      const arma::mat ZZ = Z.t() * Z;
      const arma::vec Zy = Z.t() * r;
      const arma::vec b = solve_sympd_safe(ZZ, Zy);
      const arma::vec u = r - Z * b;

      if (L <= 3) continue;

      const arma::vec u1 = u.subvec(0, L - 2);
      const arma::vec u2 = u.subvec(1, L - 1);
      const double den = arma::dot(u1, u1);
      if (den <= EPS) continue;
      double rho = arma::dot(u1, u2) / den;
      rho = nan_preserve(rho, -0.999, 0.999);

      const double adj = (1.0 - rho * rho) / std::max(3, L);
      const double rho_bc = nan_preserve(rho + adj, -0.999, 0.999);
      const double w = std::max(1.0, static_cast<double>(L) - 3.0);
      z_sum += w * std::atanh(rho_bc);
      w_sum += w;
    }
  }

  if (w_sum <= 0.0) return 0.0;
  return nan_preserve(std::tanh(z_sum / w_sum), -0.999, 0.999);
}

static bool ba_rm_beta_gls_given_theta(const std::vector<BlockStats>& stats,
                                       int p,
                                       double su2,
                                       double se2,
                                       arma::vec& beta,
                                       arma::mat& beta_vcov) {
  GLSAccumulator acc;
  acc.p = p;
  for (const BlockStats& B : stats) {
    if (B.n_i == 0) continue;
    ba_rm_accumulate_gls_block(B, p, su2, se2, acc);
  }
  return ba_rm_solve_gls_system(acc, beta, beta_vcov);
}

static double ba_rm_reml_loglik_fixed(const std::vector<BlockStats>& stats,
                                      const arma::vec& beta,
                                      double su2,
                                      double se2) {
  const int p = beta.n_elem;
  int n = 0;

  const double eps = 1e-12;
  const double se2_use = std::max(se2, eps);
  const double inv_se = 1.0 / se2_use;
  const bool at_su0 = ba_rm_at_su2_boundary(su2);

  GLSAccumulator acc;
  acc.p = p;
  double sum_logdetV = 0.0;
  double yViY = 0.0;

  for (const BlockStats& B : stats) {
    if (B.n_i == 0) continue;
    n += B.n_i;

    const double s0 = inv_se * B.UCU;
    const double sy = inv_se * B.UTCy;
    double Minv = 0.0;
    if (!at_su0) {
      const double su2_use = std::max(su2, eps);
      const double inv_su = 1.0 / su2_use;
      const double M = std::max(inv_su + s0, 1e-12);
      Minv = 1.0 / M;
    }

    acc.A00 += s0 - s0 * Minv * s0;
    acc.b0 += sy - s0 * Minv * sy;
    if (p == 2) {
      const double s1 = inv_se * B.UTCx;
      acc.A01 += inv_se * B.UTCx - s0 * Minv * s1;
      acc.A11 += inv_se * B.xCx - s1 * Minv * s1;
      acc.b1 += inv_se * B.xCy - s1 * Minv * sy;
    }

    double UTCr = 0.0;
    double q = 0.0;
    ba_rm_residual_quadratics(B, beta, UTCr, q);
    const double quad = inv_se * q - (inv_se * inv_se) * UTCr * UTCr * Minv;
    if (!std::isfinite(quad)) return -std::numeric_limits<double>::infinity();
    yViY += quad;

    double logdetV_i = static_cast<double>(B.n_i) * std::log(se2_use) - B.logdetCinv;
    if (!at_su0) logdetV_i += std::log1p((su2 / se2_use) * B.UCU);
    if (!std::isfinite(logdetV_i)) return -std::numeric_limits<double>::infinity();
    sum_logdetV += logdetV_i;
  }

  arma::vec beta_tmp;
  arma::mat beta_vcov_tmp;
  double logdetXtViX = NA_REAL;
  if (!ba_rm_solve_gls_system(acc, beta_tmp, beta_vcov_tmp, &logdetXtViX)) {
    return -std::numeric_limits<double>::infinity();
  }

  const double yPy = yViY - ba_rm_dot_xtviy_beta(acc, beta);
  if (!std::isfinite(logdetXtViX) || !std::isfinite(yPy)) {
    return -std::numeric_limits<double>::infinity();
  }

  const double two_pi = 2.0 * M_PI;
  return -0.5 * (
    (static_cast<double>(n) - static_cast<double>(p)) * std::log(two_pi) +
      sum_logdetV +
      logdetXtViX +
      yPy
  );
}

static double ba_rm_reml_loglik(const std::vector<BlockStats>& stats,
                                const FitOut& fit) {
  return ba_rm_reml_loglik_fixed(stats, fit.beta, fit.su2, fit.se2);
}

static double estimate_rho_profile(const std::vector<BlockData>& blocks,
                                   bool include_slope,
                                   int p,
                                   int max_iter,
                                   double tol,
                                   double rho_seed) {
  RhoStatsCache cache(blocks, include_slope, /*use_ar1*/ true);

  auto eval_rho = [&](double rho, double& loglik_out) -> bool {
    try {
      const double rho_use = nan_preserve(rho, -0.95, 0.95);
      const std::vector<BlockStats>* stats_ptr = nullptr;
      if (!cache.get(rho_use, stats_ptr) || stats_ptr == nullptr) return false;
      FitOut fit = fit_diff_em(p, blocks, *stats_ptr, max_iter, tol);
      const double loglik = ba_rm_reml_loglik(*stats_ptr, fit);
      if (!std::isfinite(loglik)) return false;
      loglik_out = loglik;
      return true;
    } catch (...) {
      return false;
    }
  };

  std::vector<double> candidates;
  for (int k = 0; k <= 8; ++k) candidates.push_back(-0.8 + 0.2 * static_cast<double>(k));
  if (std::isfinite(rho_seed)) candidates.push_back(nan_preserve(rho_seed, -0.95, 0.95));

  double best_rho = 0.0;
  double best_loglik = -std::numeric_limits<double>::infinity();

  auto scan_candidates = [&](const std::vector<double>& grid) {
    for (double rho : grid) {
      double loglik = NA_REAL;
      if (!eval_rho(rho, loglik)) continue;
      if (loglik > best_loglik) {
        best_loglik = loglik;
        best_rho = rho;
      }
    }
  };

  scan_candidates(candidates);

  double half_width = 0.2;
  for (int pass = 0; pass < 2; ++pass) {
    const double lo = std::max(-0.95, best_rho - half_width);
    const double hi = std::min(0.95, best_rho + half_width);
    std::vector<double> grid;
    grid.reserve(7);
    if (hi <= lo) {
      grid.push_back(lo);
    } else {
      for (int k = 0; k <= 6; ++k) {
        grid.push_back(lo + (hi - lo) * static_cast<double>(k) / 6.0);
      }
    }
    scan_candidates(grid);
    half_width /= 3.0;
  }

  return nan_preserve(best_rho, -0.95, 0.95);
}

struct ProfileThetaEval {
  arma::vec beta;
  arma::mat beta_vcov;
  double su2 = NA_REAL;
  double se2 = NA_REAL;
  double rho = NA_REAL;
  double mu0 = NA_REAL;
  double sd_loa = NA_REAL;
  double loa_lower = NA_REAL;
  double loa_upper = NA_REAL;
  double loglik = -std::numeric_limits<double>::infinity();
  bool ok = false;
};

static inline double ba_rm_pos_from_eta(double eta, double floor_val = 1e-12) {
  return std::max(std::exp(eta), floor_val);
}

static inline double ba_rm_rho_from_z(double z) {
  return 0.95 * std::tanh(z);
}

static inline double ba_rm_z_from_rho(double rho) {
  double x = rho / 0.95;
  x = std::max(-0.999999, std::min(0.999999, x));
  return std::atanh(x);
}

static bool ba_rm_eval_profile_theta(const arma::vec& theta_t,
                                     RhoStatsCache& cache,
                                     int p,
                                     bool use_ar1,
                                     bool rho_free,
                                     double rho_fixed,
                                     double loa_multiplier,
                                     ProfileThetaEval& out) {
  const int min_q = (rho_free ? 2 : 1);
  if (static_cast<int>(theta_t.n_elem) < min_q) return false;

  int pos = 0;
  double su2 = 0.0;
  double se2 = NA_REAL;
  double rho = 0.0;

  const bool su2_free = (theta_t.n_elem == (rho_free ? 3 : 2));

  if (su2_free) {
    su2 = ba_rm_pos_from_eta(theta_t[pos++]);
    if (su2 <= 1e-10) su2 = 0.0;
  } else {
    su2 = 0.0;
  }

  se2 = ba_rm_pos_from_eta(theta_t[pos++]);

  if (use_ar1) {
    rho = rho_free ? ba_rm_rho_from_z(theta_t[pos]) : rho_fixed;
  }

  const std::vector<BlockStats>* stats_ptr = nullptr;
  if (!cache.get(use_ar1 ? rho : 0.0, stats_ptr) || stats_ptr == nullptr) return false;

  arma::vec beta;
  arma::mat beta_vcov;
  if (!ba_rm_beta_gls_given_theta(*stats_ptr, p, su2, se2, beta, beta_vcov)) return false;

  const double loglik = ba_rm_reml_loglik_fixed(*stats_ptr, beta, su2, se2);
  if (!std::isfinite(loglik)) return false;

  const double V_loa = std::max(0.0, su2 + se2);
  const double sd_loa = std::sqrt(V_loa);
  const double mu0 = beta[0];

  out.beta = beta;
  out.beta_vcov = beta_vcov;
  out.su2 = su2;
  out.se2 = se2;
  out.rho = rho;
  out.mu0 = mu0;
  out.sd_loa = sd_loa;
  out.loa_lower = mu0 - loa_multiplier * sd_loa;
  out.loa_upper = mu0 + loa_multiplier * sd_loa;
  out.loglik = loglik;
  out.ok = true;
  return true;
}

template <typename Fn>
static bool ba_rm_central_gradient(const arma::vec& x,
                                   const arma::vec& h,
                                   Fn fn,
                                   arma::vec& grad) {
  const int q = x.n_elem;
  grad.set_size(q);

  for (int i = 0; i < q; ++i) {
    arma::vec xp = x, xm = x;
    xp[i] += h[i];
    xm[i] -= h[i];

    double fp = NA_REAL, fm = NA_REAL;
    if (!fn(xp, fp) || !fn(xm, fm)) return false;

    grad[i] = (fp - fm) / (2.0 * h[i]);
  }

  return grad.is_finite();
}

template <typename Fn>
static bool ba_rm_central_hessian(const arma::vec& x,
                                  const arma::vec& h,
                                  Fn fn,
                                  arma::mat& H) {
  const int q = x.n_elem;
  H.set_size(q, q);

  double f0 = NA_REAL;
  if (!fn(x, f0)) return false;

  for (int i = 0; i < q; ++i) {
    arma::vec xp = x, xm = x;
    xp[i] += h[i];
    xm[i] -= h[i];

    double fp = NA_REAL, fm = NA_REAL;
    if (!fn(xp, fp) || !fn(xm, fm)) return false;

    H(i, i) = (fp - 2.0 * f0 + fm) / (h[i] * h[i]);

    for (int j = i + 1; j < q; ++j) {
      arma::vec xpp = x, xpm = x, xmp = x, xmm = x;
      xpp[i] += h[i]; xpp[j] += h[j];
      xpm[i] += h[i]; xpm[j] -= h[j];
      xmp[i] -= h[i]; xmp[j] += h[j];
      xmm[i] -= h[i]; xmm[j] -= h[j];

      double fpp = NA_REAL, fpm = NA_REAL, fmp = NA_REAL, fmm = NA_REAL;
      if (!fn(xpp, fpp) || !fn(xpm, fpm) || !fn(xmp, fmp) || !fn(xmm, fmm)) return false;

      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * h[i] * h[j]);
      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return H.is_finite();
}

static bool ba_rm_invert_observed_info(const arma::mat& Iobs, arma::mat& Sigma) {
  arma::mat A = 0.5 * (Iobs + Iobs.t());
  if (inv_sympd_safe(Sigma, A) && Sigma.is_finite()) return true;

  const double tr = arma::trace(A);
  double ridge = std::max(1e-10, 1e-8 * std::max(1.0, std::fabs(tr) / std::max(1.0, static_cast<double>(A.n_rows))));

  for (int k = 0; k < 8; ++k) {
    arma::mat Ar = A + ridge * arma::eye<arma::mat>(A.n_rows, A.n_cols);
    if (inv_sympd_safe(Sigma, Ar) && Sigma.is_finite()) return true;
    ridge *= 10.0;
  }
  return false;
}

// [[Rcpp::export]]
Rcpp::List bland_altman_repeated_em_ext_cpp(
    Rcpp::NumericVector y,
    Rcpp::IntegerVector subject,
    Rcpp::IntegerVector method,
    Rcpp::IntegerVector time,
    bool include_slope = false,
    bool use_ar1 = false,
    double ar1_rho = NA_REAL,
    int max_iter = 200,
    double tol = 1e-6,
    double conf_level = 0.95,
    double loa_multiplier_arg = NA_REAL,
    bool use_cov_su_se = true
) {
  (void) use_cov_su_se;

  if (y.size() == 0) stop("Empty input.");
  if (use_ar1 && Rcpp::NumericVector::is_na(ar1_rho) == false && std::fabs(ar1_rho) >= 0.999)
    stop("ar1_rho must be in (-0.999, 0.999).");

  PairData P = make_pairs(y, subject, method, time);

  std::vector<int> subj_idx;
  int m = 0;
  reindex(P.subj, subj_idx, m);

  const int p = include_slope ? 2 : 1;
  double x2_mean = 0.0;
  double x2_scale = 1.0;
  if (include_slope) {
    arma::vec x2(P.mean.data(), P.mean.size(), /*copy_aux_mem*/ false, /*strict*/ true);
    x2_mean = arma::mean(x2);
    x2_scale = ba_rm_slope_scale_denom(P.mean);
  }

  std::vector<BlockData> blocks =
    build_subject_blocks(P, subj_idx, m, include_slope, x2_mean, x2_scale);

  bool ar1_estimated = false;
  double rho_used = NA_REAL;
  FitOut fit;

  auto select_best_fit = [&](const std::vector<BlockStats>& stats_local,
                             const FitOut* fit_interior_or_null) -> FitOut {
    FitOut fit_boundary = fit_diff_boundary_su0(p, stats_local);
    const double ll_boundary = ba_rm_reml_loglik(stats_local, fit_boundary);

    if (fit_interior_or_null == nullptr) return fit_boundary;

    const double ll_interior = ba_rm_reml_loglik(stats_local, *fit_interior_or_null);
    if (!std::isfinite(ll_interior)) return fit_boundary;
    if (!std::isfinite(ll_boundary)) return *fit_interior_or_null;
    if (ll_boundary >= ll_interior - 1e-8) return fit_boundary;
    return *fit_interior_or_null;
  };

  if (use_ar1 && Rcpp::NumericVector::is_na(ar1_rho)) {
    std::vector<BlockStats> stats_iid;
    if (!ba_rm_compute_all_block_stats(blocks, include_slope, false, 0.0, stats_iid)) {
      stop("Failed to build iid sufficient statistics for paired blocks.");
    }

    FitOut fit_iid_interior;
    FitOut* fit_iid_ptr = nullptr;
    try {
      fit_iid_interior = fit_diff_em(p, blocks, stats_iid, max_iter, tol);
      fit_iid_ptr = &fit_iid_interior;
    } catch (...) {
      fit_iid_ptr = nullptr;
    }
    FitOut fit_iid = select_best_fit(stats_iid, fit_iid_ptr);

    const double rho_mom = estimate_rho_moments(fit_iid, blocks);
    const double rho_hat = estimate_rho_profile(blocks, include_slope, p, max_iter, tol, rho_mom);
    rho_used = rho_hat;
    ar1_estimated = true;

    std::vector<BlockStats> stats_ar1;
    if (!ba_rm_compute_all_block_stats(blocks, include_slope, true, rho_hat, stats_ar1)) {
      stop("Failed to build AR(1) sufficient statistics for paired blocks.");
    }

    FitOut fit_ar1_interior;
    FitOut* fit_ar1_ptr = nullptr;
    try {
      fit_ar1_interior = fit_diff_em(p, blocks, stats_ar1, max_iter, tol);
      fit_ar1_ptr = &fit_ar1_interior;
    } catch (...) {
      fit_ar1_ptr = nullptr;
    }
    fit = select_best_fit(stats_ar1, fit_ar1_ptr);
  } else {
    if (use_ar1) rho_used = nan_preserve(ar1_rho, -0.999, 0.999);

    std::vector<BlockStats> stats;
    if (!ba_rm_compute_all_block_stats(blocks, include_slope, use_ar1, (use_ar1 ? rho_used : 0.0), stats)) {
      stop("Failed to build sufficient statistics for paired blocks.");
    }

    FitOut fit_interior;
    FitOut* fit_ptr = nullptr;
    try {
      fit_interior = fit_diff_em(p, blocks, stats, max_iter, tol);
      fit_ptr = &fit_interior;
    } catch (...) {
      fit_ptr = nullptr;
    }
    fit = select_best_fit(stats, fit_ptr);
  }

  arma::vec beta = fit.beta;
  arma::mat beta_vcov = fit.beta_vcov;

  const double beta_center = beta[0];
  double beta0_orig = beta[0];
  double beta1_orig = (p == 2 ? beta[1] : NA_REAL);

  if (include_slope) {
    beta1_orig = beta[1] / x2_scale;
    beta0_orig = beta[0] - beta1_orig * x2_mean;
  }

  const double su2 = fit.su2;
  const double se2 = fit.se2;
  const double mu0 = beta_center;

  const double alpha = 1.0 - std::min(std::max(conf_level, 0.0), 1.0);
  const double z = R::qnorm(1.0 - 0.5 * alpha, 0.0, 1.0, 1, 0);
  double loa_multiplier = loa_multiplier_arg;
  if (!std::isfinite(loa_multiplier) || loa_multiplier <= 0.0) loa_multiplier = z;

  const double V_loa_model = nan_preserve(su2 + se2, 0.0, 1e12);
  const double sd_loa = std::sqrt(V_loa_model);
  const double loa_lower = mu0 - loa_multiplier * sd_loa;
  const double loa_upper = mu0 + loa_multiplier * sd_loa;

  const double V_loa_empirical =
    nan_preserve(loa_var_subject_equal_weight(blocks, mu0), 0.0, 1e12);
  const double sd_loa_empirical = std::sqrt(V_loa_empirical);

  const double eps_theta = std::max(1e-10, 1e-10 * std::max(1.0, ba_rm_variance_all_y(blocks)));
  const bool rho_free = (use_ar1 && ar1_estimated);
  const bool su2_boundary = fit.boundary_su2_zero;

  arma::vec theta_hat((su2_boundary ? 1 : 2) + (rho_free ? 1 : 0), arma::fill::zeros);
  int th_pos = 0;
  if (!su2_boundary) theta_hat[th_pos++] = std::log(std::max(su2, eps_theta));
  theta_hat[th_pos++] = std::log(std::max(se2, eps_theta));
  if (rho_free) theta_hat[th_pos++] = ba_rm_z_from_rho(rho_used);

  arma::vec step(theta_hat.n_elem, arma::fill::zeros);
  for (arma::uword j = 0; j < theta_hat.n_elem; ++j) {
    step[j] = 1e-4 * std::max(1.0, std::fabs(theta_hat[j]));
  }

  RhoStatsCache profile_cache(blocks, include_slope, use_ar1);

  auto loglik_fn = [&](const arma::vec& th, double& out) -> bool {
    ProfileThetaEval E;
    if (!ba_rm_eval_profile_theta(th, profile_cache, p, use_ar1, rho_free, rho_used, loa_multiplier, E)) return false;
    out = E.loglik;
    return std::isfinite(out);
  };

  auto mu_fn = [&](const arma::vec& th, double& out) -> bool {
    ProfileThetaEval E;
    if (!ba_rm_eval_profile_theta(th, profile_cache, p, use_ar1, rho_free, rho_used, loa_multiplier, E)) return false;
    out = E.mu0;
    return std::isfinite(out);
  };

  auto lower_fn = [&](const arma::vec& th, double& out) -> bool {
    ProfileThetaEval E;
    if (!ba_rm_eval_profile_theta(th, profile_cache, p, use_ar1, rho_free, rho_used, loa_multiplier, E)) return false;
    out = E.loa_lower;
    return std::isfinite(out);
  };

  auto upper_fn = [&](const arma::vec& th, double& out) -> bool {
    ProfileThetaEval E;
    if (!ba_rm_eval_profile_theta(th, profile_cache, p, use_ar1, rho_free, rho_used, loa_multiplier, E)) return false;
    out = E.loa_upper;
    return std::isfinite(out);
  };

  arma::mat H_theta;
  if (!ba_rm_central_hessian(theta_hat, step, loglik_fn, H_theta)) {
    stop("Profile-REML CI calculation failed while evaluating the numerical Hessian.");
  }

  arma::mat I_theta = -0.5 * (H_theta + H_theta.t());
  arma::mat Sigma_theta;
  if (!ba_rm_invert_observed_info(I_theta, Sigma_theta)) {
    stop("Profile-REML CI calculation failed because the observed information matrix was not invertible.");
  }

  arma::vec g_mu, g_lower, g_upper;
  if (!ba_rm_central_gradient(theta_hat, step, mu_fn, g_mu)) {
    stop("Profile-REML CI calculation failed while evaluating the gradient for the bias.");
  }
  if (!ba_rm_central_gradient(theta_hat, step, lower_fn, g_lower)) {
    stop("Profile-REML CI calculation failed while evaluating the gradient for the lower LoA.");
  }
  if (!ba_rm_central_gradient(theta_hat, step, upper_fn, g_upper)) {
    stop("Profile-REML CI calculation failed while evaluating the gradient for the upper LoA.");
  }

  const double var_mu_cond = (beta_vcov.n_rows >= 1 ? std::max(0.0, beta_vcov(0, 0)) : 0.0);
  const double var_lower_cond = var_mu_cond;
  const double var_upper_cond = var_mu_cond;

  const double var_mu_prof = std::max(0.0, arma::as_scalar(g_mu.t() * Sigma_theta * g_mu));
  const double var_lower_prof = std::max(0.0, arma::as_scalar(g_lower.t() * Sigma_theta * g_lower));
  const double var_upper_prof = std::max(0.0, arma::as_scalar(g_upper.t() * Sigma_theta * g_upper));

  const double var_mu0 = var_mu_cond + var_mu_prof;
  const double var_loa_lower = var_lower_cond + var_lower_prof;
  const double var_loa_upper = var_upper_cond + var_upper_prof;

  const double se_bias = std::sqrt(std::max(0.0, var_mu0));
  const double se_loa_lower = std::sqrt(std::max(0.0, var_loa_lower));
  const double se_loa_upper = std::sqrt(std::max(0.0, var_loa_upper));

  const double bias_lwr = mu0 - z * se_bias;
  const double bias_upr = mu0 + z * se_bias;
  const double loa_lower_lwr = loa_lower - z * se_loa_lower;
  const double loa_lower_upr = loa_lower + z * se_loa_lower;
  const double loa_upper_lwr = loa_upper - z * se_loa_upper;
  const double loa_upper_upr = loa_upper + z * se_loa_upper;

  return List::create(
    _["n_pairs"] = static_cast<int>(P.d.size()),
    _["n_subjects"] = m,
    _["pairs_mean"] = Rcpp::NumericVector(P.mean.begin(), P.mean.end()),
    _["pairs_diff"] = Rcpp::NumericVector(P.d.begin(), P.d.end()),
    _["bias_mu0"] = mu0,
    _["bias_se"] = se_bias,
    _["bias_lwr"] = bias_lwr,
    _["bias_upr"] = bias_upr,
    _["beta_intercept"] = beta0_orig,
    _["beta_slope"] = (include_slope ? beta1_orig : NA_REAL),
    _["beta_center"] = beta_center,
    _["beta_center_reference_mean"] = (include_slope ? x2_mean : NA_REAL),
    _["sigma2_subject"] = su2,
    _["sigma2_resid"] = se2,
    _["sd_loa"] = sd_loa,
    _["loa_var_model"] = V_loa_model,
    _["loa_lower"] = loa_lower,
    _["loa_upper"] = loa_upper,
    _["loa_lower_lwr"] = loa_lower_lwr,
    _["loa_lower_upr"] = loa_lower_upr,
    _["loa_upper_lwr"] = loa_upper_lwr,
    _["loa_upper_upr"] = loa_upper_upr,
    _["loa_var_empirical"] = V_loa_empirical,
    _["sd_loa_empirical"] = sd_loa_empirical,
    _["use_ar1"] = use_ar1,
    _["ar1_rho"] = (use_ar1 ? rho_used : NA_REAL),
    _["ar1_estimated"] = (use_ar1 ? ar1_estimated : false),
    _["loa_multiplier"] = loa_multiplier,
    _["conf_level"] = conf_level,
    _["converged"] = fit.converged,
    _["iter"] = fit.iter,
    _["warn"] = fit.warn
  );
}
