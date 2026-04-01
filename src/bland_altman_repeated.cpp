// Thiago de Paula Oliveira
// Repeated-measures Bland–Altman via stabilized EM/GLS with optional AR(1)
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <sstream>
#include <cmath>
#include <numeric>
#include <algorithm>

// bring helpers
#include "matrixCorr_detail.h"

using namespace Rcpp;
using namespace arma;

// ---- use selected helpers from matrixCorr_detail ----
using matrixCorr_detail::clamp_policy::nan_preserve;
using matrixCorr_detail::linalg::inv_sympd_safe;
using matrixCorr_detail::linalg::logdet_spd_safe;
using matrixCorr_detail::linalg::solve_sympd_safe;
using matrixCorr_detail::moments::sample_var;
using matrixCorr_detail::moments::sample_cov;
using matrixCorr_detail::indexing::reindex;

namespace {
constexpr const char* BA_RM_IDENTIFIABILITY_ERROR =
  "The repeated-measures Bland-Altman mixed model is not separately identifiable on the supplied paired data because there is insufficient within-subject replication after pairing to separate the residual and subject-level variance components.";
constexpr const char* BA_RM_SLOPE_SCALE_ERROR =
  "The proportional-bias slope is not estimable because the paired means are near-degenerate on the observed data scale.";
constexpr const char* BA_RM_CONVERGENCE_ERROR =
  "The repeated-measures Bland-Altman model is conceptually estimable on the supplied paired data, but the EM/GLS algorithm failed to converge to admissible finite variance-component estimates.";
} // namespace

// ---------- local helpers ----------
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
    arma::vec mean_vec(sorted.data(), sorted.size(), /*copy_aux_mem*/ true);
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

// ---------- Section 1: Efficient Pair Construction (integer-based, no string keys) ----------
struct PairData {
  std::vector<double> d;     // y2 - y1
  std::vector<double> mean;  // (y1 + y2)/2
  std::vector<int>    subj;
  std::vector<int>    time;
};

// Composite key for subject-time, using 64-bit integer for uniqueness and speed
static inline uint64_t make_pair_key(int subj, int time) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(subj)) << 32) | static_cast<uint32_t>(time);
}

static PairData make_pairs(const NumericVector& y,
                           const IntegerVector& subject,
                           const IntegerVector& method,
                           const IntegerVector& time) {
  const int n = y.size();
  if (subject.size()!=n || method.size()!=n || time.size()!=n)
    stop("lengths of y, subject, method, time must match.");
  struct V { bool has1=false, has2=false; double y1=NA_REAL, y2=NA_REAL; int s, t; };
  std::unordered_map<uint64_t, V> H; H.reserve((size_t)n);
  for (int i=0;i<n;++i) {
    if (IntegerVector::is_na(subject[i]) || IntegerVector::is_na(method[i]) || IntegerVector::is_na(time[i])) continue;
    if (NumericVector::is_na(y[i])) continue;
    int s = subject[i], t = time[i], m = method[i];
    if (m!=1 && m!=2) continue;     // only the two target methods for this pair
    uint64_t key = make_pair_key(s, t);
    auto &v = H[key]; v.s = s; v.t = t;
    if (m==1) { v.has1=true; v.y1=y[i]; }
    else      { v.has2=true; v.y2=y[i]; }
  }
  PairData P;
  P.d.reserve(H.size()); P.mean.reserve(H.size());
  P.subj.reserve(H.size()); P.time.reserve(H.size());
  for (auto &kv : H) {
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

// ---------- Section 2: Fast Subject/Block Indexing (preallocated, sorted) ----------
struct BySubjBA {
  std::vector< std::vector<int> > rows;
  std::vector< std::vector<int> > tim;
};
static BySubjBA index_by_subject_ba(const std::vector<int>& subj_idx,
                                    const std::vector<int>& time) {
  const int m = *std::max_element(subj_idx.begin(), subj_idx.end()) + 1;
  BySubjBA S; S.rows.assign(m,{}); S.tim.assign(m,{});
  for (size_t i=0;i<subj_idx.size();++i) {
    int j = subj_idx[i];
    S.rows[j].push_back((int)i);
    S.tim[j].push_back(time[i]);
  }
  // Sort by time, but avoid unnecessary copying
  for (int i=0;i<m;++i) {
    auto &r = S.rows[i]; auto &t = S.tim[i];
    if (r.size() <= 1) continue;
    std::vector<int> ord(r.size()); std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){
      int ta=t[a], tb=t[b];
      if (ta<0 && tb<0) return a<b;
      if (ta<0) return false;
      if (tb<0) return true;
      return ta<tb;
    });
    std::vector<int> r2(r.size()), t2(t.size());
    for (size_t k=0; k<ord.size(); ++k) { r2[k]=r[ord[k]]; t2[k]=t[ord[k]]; }
    r.swap(r2); t.swap(t2);
  }
  return S;
}

// ---------- Section 3: Per-Subject Precomputation (preallocated, AR(1) tridiagonal) ----------
struct Precomp {
  int n_i=0;
  mat X_i;   // n_i x p
  vec y_i;   // n_i
  mat Cinv;  // n_i x n_i  (I or AR1 precision)

  // sufficient stats with Cinv (scale-free):
  mat XTCX;  // p x p
  vec XTCy;  // p
  vec UTCX;  // p (1^T C X)
  double UTCy = 0.0;   // 1^T C y
  double UCU  = 0.0;   // 1^T C 1
};

// Efficient AR(1) precision: tridiagonal, in-place
static inline void ar1_precision_tridiag(const std::vector<int>& tim, double rho, mat& Cinv) {
  const int n = (int)tim.size();
  Cinv.zeros(n, n);
  if (n == 0) return;
  const double r2 = rho * rho;
  const double denom = std::max(1.0 - r2, 1e-12);
  int s = 0;
  while (s < n) {
    if (tim[s] < 0) { Cinv(s,s) += 1.0; ++s; continue; }  // singleton/NA -> iid
    int e = s;
    while (e + 1 < n && tim[e+1] >= 0 && tim[e+1] == tim[e] + 1) ++e; // contiguous only
    const int L = e - s + 1;
    if (L == 1) {
      Cinv(s,s) += 1.0;
    } else {
      Cinv(s,s) += 1.0 / denom;
      Cinv(e,e) += 1.0 / denom;
      for (int t = s+1; t <= e-1; ++t) Cinv(t,t) += (1.0 + r2) / denom;
      for (int t = s; t <= e-1; ++t) {
        Cinv(t, t+1) += -rho / denom;
        Cinv(t+1, t) += -rho / denom;
      }
    }
    s = e + 1;
  }
  Cinv.diag() += 1e-10; // tiny ridge
}

static std::vector<Precomp> precompute_blocks(const mat& X, const vec& y,
                                              const BySubjBA& S,
                                              bool use_ar1, double rho) {
  const int p = X.n_cols;
  const int m = (int)S.rows.size();
  std::vector<Precomp> out(m);
  for (int i=0;i<m;++i) {
    const auto& rows = S.rows[i];
    const auto& tim  = S.tim[i];
    const int n_i = (int)rows.size();
    if (n_i==0) continue;
    Precomp P;
    P.n_i = n_i;
    P.X_i.set_size(n_i,p); P.y_i.set_size(n_i);
    std::vector<int> tim_ord(n_i,-1);
    for (int k=0;k<n_i;++k) { int g = rows[k]; P.X_i.row(k)=X.row(g); P.y_i[k]=y[g]; tim_ord[k]=tim[k]; }
    if (use_ar1) ar1_precision_tridiag(tim_ord, rho, P.Cinv);
    else P.Cinv.eye(n_i, n_i);

    vec ones(n_i, fill::ones);

    // C-side sufficient stats
    mat CX = P.Cinv * P.X_i;
    P.XTCX = P.X_i.t() * CX;
    P.XTCy = P.X_i.t() * (P.Cinv * P.y_i);
    vec CU = P.Cinv * ones;
    P.UTCX = (ones.t() * CX).t();
    P.UTCy = arma::as_scalar(ones.t() * (P.Cinv * P.y_i));
    P.UCU  = arma::as_scalar(ones.t() * CU);

    out[i] = std::move(P);
  }
  return out;
}

// ---------- Section 4: EM/GLS for Differences (blockwise, memory-efficient) ----------
struct FitOut {
  arma::vec beta;
  arma::mat beta_vcov;   // GLS covariance of beta conditional on fitted variance params
  double su2 = NA_REAL, se2 = NA_REAL;
  arma::vec su_term;     // E[u_i^2|y] per subject (for delta)
  arma::vec se_term;     // per-subject contribution for se2 variance (normalized)
  int iter = 0;
  bool converged = false;
  bool boundary_su2_zero = false;
  std::string warn;
};

static inline bool ba_rm_at_su2_boundary(double su2, double tol = 1e-10) {
  return std::isfinite(su2) && su2 <= tol;
}

// Boundary fit: sigma2_subject = 0
static FitOut fit_diff_boundary_su0(const mat& X,
                                    const vec& y,
                                    const std::vector<Precomp>& PC) {
  const int p = X.n_cols;
  int n = 0;
  arma::mat XtCX(p, p, arma::fill::zeros);
  arma::vec XtCy(p, arma::fill::zeros);
  double yCy = 0.0;
  for (const Precomp& P : PC) {
    if (P.n_i == 0) continue;
    XtCX += P.XTCX;
    XtCy += P.XTCy;
    yCy  += arma::as_scalar(P.y_i.t() * (P.Cinv * P.y_i));
    n    += P.n_i;
  }
  arma::mat XtCX_inv;
  if (!inv_sympd_safe(XtCX_inv, XtCX) || !XtCX_inv.is_finite()) {
    stop("Boundary fit with sigma2_subject = 0 failed because X'CX was not invertible.");
  }
  arma::vec beta = XtCX_inv * XtCy;
  const double yPy = yCy - arma::as_scalar(XtCy.t() * XtCX_inv * XtCy);
  const double se2 = std::max(1e-12, yPy / std::max(1, n - p));
  FitOut out;
  out.beta = beta;
  out.beta_vcov = se2 * XtCX_inv;
  out.su2 = 0.0;
  out.se2 = se2;
  out.su_term = arma::vec(PC.size(), arma::fill::zeros);
  out.se_term = arma::vec(PC.size(), arma::fill::zeros);
  out.iter = 0;
  out.converged = true;
  out.boundary_su2_zero = true;
  out.warn = "Constrained boundary fit used with sigma2_subject fixed at 0.";
  return out;
}

// Interior fit: EM for random subject intercept + (optional) AR(1)
static FitOut fit_diff_em(const mat& X, const vec& y,
                          const BySubjBA& S,
                          const std::vector<Precomp>& PC,
                          int max_iter, double tol) {
  const int p = X.n_cols;
  const int m = (int)S.rows.size();
  const int n = y.n_rows;
  const double vref = arma::var(y, /*unbiased*/1);
  const double EPS  = std::max(1e-12, vref * 1e-12);
  const double MAXV = std::max(10.0 * vref, 1.0);
  std::vector<double> mu_s(m, 0.0);
  std::vector<int> cnt_s(m, 0);
  for (int i = 0; i < m; ++i) {
    for (int r : S.rows[i]) {
      mu_s[i] += y[r];
      ++cnt_s[i];
    }
  }
  int used_mu = 0;
  for (int i = 0; i < m; ++i) {
    if (cnt_s[i] > 0) {
      mu_s[i] /= static_cast<double>(cnt_s[i]);
      ++used_mu;
    }
  }
  int n_replicated_subjects = 0;
  for (int i = 0; i < m; ++i) if (cnt_s[i] >= 2) ++n_replicated_subjects;
  if (used_mu < 2 || n_replicated_subjects < 1) {
    stop(BA_RM_IDENTIFIABILITY_ERROR);
  }
  arma::vec mu_s_vec(used_mu, arma::fill::zeros);
  {
    int k = 0;
    for (int i = 0; i < m; ++i) if (cnt_s[i] > 0) mu_s_vec[k++] = mu_s[i];
  }
  double var_mu = (used_mu >= 2 ? arma::var(mu_s_vec, /*unbiased*/true) : 0.0);
  double num_w = 0.0, den_w = 0.0;
  bool se2_init_fallback = false;
  for (int i = 0; i < m; ++i) if ((int)S.rows[i].size() >= 2) {
    arma::vec yi(S.rows[i].size());
    for (size_t k = 0; k < S.rows[i].size(); ++k) yi[k] = y[S.rows[i][k]];
    double vsi = arma::var(yi, /*unbiased*/true);
    if (std::isfinite(vsi) && vsi > 0.0) {
      num_w += ((int)yi.n_elem - 1) * vsi;
      den_w += ((int)yi.n_elem - 1);
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
    var_mu - se2_init / std::max(1.0, arma::mean(arma::conv_to<arma::vec>::from(cnt_s)))
  );
  se2_init = nan_preserve(se2_init, EPS, MAXV);
  su2_init = nan_preserve(su2_init, 0.0, MAXV);
  auto damp_to_ratio = [](double oldv, double newv, double rmax) {
    if (!std::isfinite(newv)) return oldv;
    if (oldv <= 0.0) return nan_preserve(newv, 0.0, rmax);
    double lo = oldv / rmax, hi = oldv * rmax;
    return nan_preserve(newv, std::min(lo, hi), std::max(lo, hi));
  };
  double su2 = su2_init, se2 = se2_init;
  vec beta(p, fill::zeros);
  mat beta_vcov(p, p, fill::zeros);
  vec su_term(m, fill::zeros);
  vec se_term(m, fill::zeros);
  bool converged = false;
  int it = -1;
  std::string init_warn;
  if (se2_init_fallback) {
    init_warn = "Residual-variance initialization used 0.5 * v_ref only as a positive EM starting heuristic because no finite positive pooled within-subject variance could be formed after pairing.";
  }
  for (it = 0; it < max_iter; ++it) {
    mat XtViX(p, p, fill::zeros);
    vec XtViy(p, fill::zeros);
    const double inv_se = 1.0 / std::max(se2, EPS);
    for (int i = 0; i < m; ++i) {
      const Precomp& P = PC[i];
      if (P.n_i == 0) continue;
      double M = (1.0 / std::max(su2, EPS)) + inv_se * P.UCU;
      M = std::max(M, 1e-12);
      const double Minv = 1.0 / M;
      double S_uy = inv_se * P.UTCy;
      double Z_y  = Minv * S_uy;
      vec S_ux = inv_se * P.UTCX;
      vec Z_X  = Minv * S_ux;
      mat XTRinvX = inv_se * P.XTCX;
      vec XTRinvY = inv_se * P.XTCy;
      XtViy += XTRinvY - S_ux * Z_y;
      XtViX += XTRinvX - S_ux * Z_X.t();
    }
    mat XtViX_inv;
    if (!inv_sympd_safe(XtViX_inv, XtViX) || !XtViX_inv.is_finite()) break;
    beta = XtViX_inv * XtViy;
    beta_vcov = XtViX_inv;
    if (!all_finite(beta)) break;
    double su_acc = 0.0;
    double se_num = 0.0;
    vec r_global = y - X * beta;
    for (int i = 0; i < m; ++i) {
      const Precomp& P = PC[i];
      if (P.n_i == 0) continue;
      vec r_i(P.n_i);
      for (int k = 0; k < P.n_i; ++k) r_i[k] = r_global[S.rows[i][k]];
      const double inv_su = 1.0 / std::max(su2, EPS);
      double M = inv_su + (1.0 / std::max(se2, EPS)) * P.UCU;
      M = std::max(M, 1e-12);
      const double Minv = 1.0 / M;
      double Utr = (1.0 / std::max(se2, EPS)) *
        arma::as_scalar(arma::ones<vec>(P.n_i).t() * (P.Cinv * r_i));
      double b_i = Minv * Utr;
      double Eu2 = b_i * b_i + Minv;
      su_acc += Eu2;
      su_term[i] = Eu2;
      vec e = r_i - arma::ones<vec>(P.n_i) * b_i;
      double q = arma::as_scalar(e.t() * (P.Cinv * e));
      double num_i = q + P.UCU * Minv;
      se_num += num_i;
      se_term[i] = num_i / std::max(1, P.n_i);
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
    double diff = std::fabs(su2_new - su2) + std::fabs(se2_new - se2);
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
  out.su_term = su_term;
  out.se_term = se_term;
  out.iter = it + 1;
  out.converged = true;
  out.warn = init_warn;
  return out;
}

// subject-equal bias (equal weight to each subject's mean difference)
static inline double bias_subject_equal_weight(const std::vector<double>& d,
                                               const std::vector<int>& subj_idx,
                                               int m,
                                               arma::vec& subj_means_out) {
  std::vector<double> sum(m,0.0); std::vector<int> cnt(m,0);
  for (size_t i=0;i<d.size();++i){ sum[subj_idx[i]] += d[i]; ++cnt[subj_idx[i]]; }
  subj_means_out.set_size(m);
  int used=0; double acc=0.0;
  for (int i=0;i<m;++i) {
    if (cnt[i]>0) { double mu_i = sum[i]/(double)cnt[i]; subj_means_out[i]=mu_i; acc += mu_i; ++used; }
    else subj_means_out[i]=0.0;
  }
  if (used==0) stop("No data per subject to compute bias.");
  return acc / (double)used;
}

static inline double loa_var_subject_equal_weight(const arma::vec& d,
                                                  const BySubjBA& S,
                                                  double mu0) {
  double acc = 0.0;
  int used = 0;
  for (const auto& rows_i : S.rows) {
    if (rows_i.empty()) continue;
    double ss_i = 0.0;
    for (int row : rows_i) {
      const double dev = d[row] - mu0;
      ss_i += dev * dev;
    }
    acc += ss_i / static_cast<double>(rows_i.size());
    ++used;
  }
  if (used == 0) stop("No data per subject to compute LoA variance.");
  return acc / static_cast<double>(used);
}

// ----- AR(1) rho from residuals (robust for short L)
static double estimate_rho_moments(const FitOut& fit,
                                   const arma::mat& X,
                                   const arma::vec& y,
                                   const BySubjBA& S,
                                   const std::vector<Precomp>& /*unused*/) {
  const arma::vec beta = fit.beta;
  const double EPS = 1e-12;

  double z_sum = 0.0, w_sum = 0.0;   // Fisher-z

  for (size_t i = 0; i < S.rows.size(); ++i) {
    const auto& rows = S.rows[i];
    const auto& tim  = S.tim[i];
    const int n_i = (int)rows.size();
    if (n_i <= 2) continue;

    int s = 0;
    while (s < n_i) {
      if (tim[s] < 0) { ++s; continue; }
      int e = s;
      while (e + 1 < n_i && tim[e+1] == tim[e] + 1) ++e;
      const int L = e - s + 1;
      if (L >= 3) {
        arma::vec r(L);
        for (int k = 0; k < L; ++k)
          r[k] = y[ rows[s+k] ] - arma::as_scalar( X.row( rows[s+k] ) * beta );

        // Detrend by (intercept + linear time) within the block
        arma::vec t(L); for (int k = 0; k < L; ++k) t[k] = (double)k;
        arma::mat Z(L, 2); Z.col(0).ones(); Z.col(1) = t - arma::mean(t);
        arma::mat ZZ = Z.t() * Z;
        arma::vec Zy = Z.t() * r;
        arma::vec b  = solve_sympd_safe(ZZ, Zy);
        arma::vec u  = r - Z * b;

        if (L <= 3) { s = e + 1; continue; } // need at least 4 after detrend to be stable

        // Lag-1 on detrended residuals
        arma::vec u1 = u.subvec(0, L - 2);
        arma::vec u2 = u.subvec(1, L - 1);
        double den = arma::dot(u1, u1);
        if (den <= EPS) { s = e + 1; continue; }
        double rho = arma::dot(u1, u2) / den;
        rho = nan_preserve(rho, -0.999, 0.999);

        // small-sample adjustment
        double adj = (1.0 - rho * rho) / std::max(3, L);
        double rho_bc = nan_preserve(rho + adj, -0.999, 0.999);

        // Fisher-z pool; effective weight approx. L - 3
        double w = std::max(1.0, (double)L - 3.0);
        z_sum += w * std::atanh(rho_bc);
        w_sum += w;
      }
      s = e + 1;
    }
  }
  if (w_sum <= 0.0) return 0.0;
  return nan_preserve(std::tanh(z_sum / w_sum), -0.999, 0.999);
}

static double ba_rm_reml_loglik(const arma::mat& X,
                                const arma::vec& y,
                                const std::vector<Precomp>& PC,
                                const FitOut& fit) {
  const int p = X.n_cols;
  const int n = y.n_rows;
  const double eps = 1e-12;
  const double su2 = std::max(fit.su2, 0.0);
  const double se2 = std::max(fit.se2, eps);
  const double inv_se = 1.0 / se2;
  const bool at_su0 = ba_rm_at_su2_boundary(su2);

  arma::mat XtViX(p, p, arma::fill::zeros);
  arma::vec XtViy(p, arma::fill::zeros);
  double sum_logdetV = 0.0;
  double yViY = 0.0;

  for (const Precomp& P : PC) {
    if (P.n_i == 0) continue;

    double Minv = 0.0;
    if (!at_su0) {
      const double inv_su = 1.0 / su2;
      const double M = std::max(inv_su + inv_se * P.UCU, 1e-12);
      Minv = 1.0 / M;
    }

    const arma::vec S_ux = inv_se * P.UTCX;
    const arma::vec XTRinvY = inv_se * P.XTCy;
    const arma::mat XTRinvX = inv_se * P.XTCX;
    const double S_uy = inv_se * P.UTCy;
    const double Z_y = Minv * S_uy;
    const arma::vec Z_X = Minv * S_ux;

    XtViy += XTRinvY - S_ux * Z_y;
    XtViX += XTRinvX - S_ux * Z_X.t();

    const arma::vec r_i = P.y_i - P.X_i * fit.beta;
    const double UTCr = P.UTCy - arma::dot(P.UTCX, fit.beta);
    const double rTCr = arma::as_scalar(r_i.t() * (P.Cinv * r_i));
    const double quad = inv_se * rTCr - (inv_se * inv_se) * UTCr * UTCr * Minv;
    if (!std::isfinite(quad)) return -std::numeric_limits<double>::infinity();
    yViY += quad;

    double logdetV_i = P.n_i * std::log(se2) - logdet_spd_safe(P.Cinv);
    if (!at_su0) {
      logdetV_i += std::log1p((su2 / se2) * P.UCU);
    }
    if (!std::isfinite(logdetV_i)) return -std::numeric_limits<double>::infinity();
    sum_logdetV += logdetV_i;
  }

  arma::mat XtViX_inv;
  if (!inv_sympd_safe(XtViX_inv, XtViX) || !XtViX_inv.is_finite()) {
    return -std::numeric_limits<double>::infinity();
  }

  const double logdetXtViX = logdet_spd_safe(XtViX);
  const double yPy = yViY - arma::as_scalar(XtViy.t() * (XtViX_inv * XtViy));
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

static double estimate_rho_profile(const arma::mat& X,
                                   const arma::vec& y,
                                   const BySubjBA& S,
                                   int max_iter,
                                   double tol,
                                   double rho_seed) {
  auto eval_rho = [&](double rho, double& loglik_out) -> bool {
    try {
      const double rho_use = nan_preserve(rho, -0.95, 0.95);
      std::vector<Precomp> PC = precompute_blocks(X, y, S, /*use_ar1*/ true, rho_use);
      FitOut fit = fit_diff_em(X, y, S, PC, max_iter, tol);
      const double loglik = ba_rm_reml_loglik(X, y, PC, fit);
      if (!std::isfinite(loglik)) return false;
      loglik_out = loglik;
      return true;
    } catch (...) {
      return false;
    }
  };

  std::vector<double> candidates;
  for (int k = 0; k <= 8; ++k) {
    candidates.push_back(-0.8 + 0.2 * static_cast<double>(k));
  }
  if (std::isfinite(rho_seed)) {
    candidates.push_back(nan_preserve(rho_seed, -0.95, 0.95));
  }

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
    const double hi = std::min( 0.95, best_rho + half_width);
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


// ---------- Section 5: Profile-REML CI helpers (efficient, blockwise) ----------
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
  // Match the profile-search support used elsewhere in this file.
  return 0.95 * std::tanh(z);
}

static inline double ba_rm_z_from_rho(double rho) {
  double x = rho / 0.95;
  x = std::max(-0.999999, std::min(0.999999, x));
  return std::atanh(x);
}

// GLS for beta at fixed variance parameters (blockwise, memory-efficient)
static bool ba_rm_beta_gls_given_theta(const std::vector<Precomp>& PC,
                                       int p,
                                       double su2,
                                       double se2,
                                       arma::vec& beta,
                                       arma::mat& beta_vcov) {
  const double eps = 1e-12;
  const double inv_se = 1.0 / std::max(se2, eps);
  const bool at_su0 = ba_rm_at_su2_boundary(su2);
  arma::mat XtViX(p, p, arma::fill::zeros);
  arma::vec XtViy(p, arma::fill::zeros);
  for (const Precomp& P : PC) {
    if (P.n_i == 0) continue;
    double Minv = 0.0;
    if (!at_su0) {
      const double inv_su = 1.0 / su2;
      const double M = std::max(inv_su + inv_se * P.UCU, 1e-12);
      Minv = 1.0 / M;
    }
    const arma::vec S_ux = inv_se * P.UTCX;
    const double    S_uy = inv_se * P.UTCy;
    const arma::mat XTRinvX = inv_se * P.XTCX;
    const arma::vec XTRinvY = inv_se * P.XTCy;
    const arma::vec Z_X = Minv * S_ux;
    const double    Z_y = Minv * S_uy;
    XtViX += XTRinvX - S_ux * Z_X.t();
    XtViy += XTRinvY - S_ux * Z_y;
  }
  arma::mat XtViX_inv;
  if (!inv_sympd_safe(XtViX_inv, XtViX) || !XtViX_inv.is_finite()) return false;
  beta = XtViX_inv * XtViy;
  beta_vcov = XtViX_inv;
  return beta.is_finite() && beta_vcov.is_finite();
}

static double ba_rm_reml_loglik_fixed(const std::vector<Precomp>& PC,
                                      const arma::vec& beta,
                                      double su2,
                                      double se2) {
  const int p = beta.n_elem;
  int n = 0;

  const double eps = 1e-12;
  const double se2_use = std::max(se2, eps);
  const double inv_se = 1.0 / se2_use;
  const bool at_su0 = ba_rm_at_su2_boundary(su2);

  arma::mat XtViX(p, p, arma::fill::zeros);
  arma::vec XtViy(p, arma::fill::zeros);
  double sum_logdetV = 0.0;
  double yViY = 0.0;

  for (const Precomp& P : PC) {
    if (P.n_i == 0) continue;
    n += P.n_i;

    double Minv = 0.0;
    if (!at_su0) {
      const double su2_use = std::max(su2, eps);
      const double inv_su = 1.0 / su2_use;
      const double M = std::max(inv_su + inv_se * P.UCU, 1e-12);
      Minv = 1.0 / M;
    }

    const arma::vec S_ux = inv_se * P.UTCX;
    const double    S_uy = inv_se * P.UTCy;

    const arma::mat XTRinvX = inv_se * P.XTCX;
    const arma::vec XTRinvY = inv_se * P.XTCy;

    const arma::vec Z_X = Minv * S_ux;
    const double    Z_y = Minv * S_uy;

    XtViX += XTRinvX - S_ux * Z_X.t();
    XtViy += XTRinvY - S_ux * Z_y;

    const arma::vec r_i = P.y_i - P.X_i * beta;
    const double UTCr = P.UTCy - arma::dot(P.UTCX, beta);
    const double rTCr = arma::as_scalar(r_i.t() * (P.Cinv * r_i));

    const double quad =
      inv_se * rTCr -
      (inv_se * inv_se) * UTCr * UTCr * Minv;

    if (!std::isfinite(quad)) return -std::numeric_limits<double>::infinity();
    yViY += quad;

    double logdetV_i =
      P.n_i * std::log(se2_use) -
      logdet_spd_safe(P.Cinv);

    if (!at_su0) {
      logdetV_i += std::log1p((su2 / se2_use) * P.UCU);
    }

    if (!std::isfinite(logdetV_i)) return -std::numeric_limits<double>::infinity();
    sum_logdetV += logdetV_i;
  }

  arma::mat XtViX_inv;
  if (!inv_sympd_safe(XtViX_inv, XtViX) || !XtViX_inv.is_finite()) {
    return -std::numeric_limits<double>::infinity();
  }

  const double logdetXtViX = logdet_spd_safe(XtViX);
  const double yPy = yViY - arma::as_scalar(XtViy.t() * (XtViX_inv * XtViy));

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

static bool ba_rm_eval_profile_theta(const arma::vec& theta_t,
                                     const arma::mat& X,
                                     const arma::vec& y,
                                     const BySubjBA& S,
                                     bool use_ar1,
                                     bool rho_free,
                                     double rho_fixed,
                                     double loa_multiplier,
                                     ProfileThetaEval& out) {
  const int min_q = (rho_free ? 2 : 1);
  if ((int)theta_t.n_elem < min_q) return false;

  int pos = 0;
  double su2 = 0.0;
  double se2 = NA_REAL;
  double rho = 0.0;

  // If theta has 2 variance parameters before rho, use interior form.
  // If theta has 1 variance parameter before rho, treat su2 as fixed at 0.
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

  std::vector<Precomp> PC = precompute_blocks(X, y, S, use_ar1, rho);

  arma::vec beta;
  arma::mat beta_vcov;
  if (!ba_rm_beta_gls_given_theta(PC, X.n_cols, su2, se2, beta, beta_vcov)) {
    return false;
  }

  const double loglik = ba_rm_reml_loglik_fixed(PC, beta, su2, se2);
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
  double ridge = std::max(1e-10, 1e-8 * std::max(1.0, std::fabs(tr) / std::max(1.0, (double)A.n_rows)));

  for (int k = 0; k < 8; ++k) {
    arma::mat Ar = A + ridge * arma::eye<arma::mat>(A.n_rows, A.n_cols);
    if (inv_sympd_safe(Sigma, Ar) && Sigma.is_finite()) return true;
    ridge *= 10.0;
  }
  return false;
}

// [[Rcpp::export]]
Rcpp::List bland_altman_repeated_em_ext_cpp(
    Rcpp::NumericVector y,          // stacked responses
    Rcpp::IntegerVector subject,    // subject id (any integers)
    Rcpp::IntegerVector method,     // MUST be 1 or 2 for the selected pair
    Rcpp::IntegerVector time,       // pairing key (replicate/time)
    bool include_slope = false,     // proportional bias: slope vs pair mean
    bool use_ar1 = false,           // AR(1) within subject on residuals
    double ar1_rho = NA_REAL,       // if NA and use_ar1=TRUE, estimate it
    int    max_iter = 200,
    double tol = 1e-6,
    double conf_level = 0.95,
    double loa_multiplier_arg = NA_REAL,
    bool   use_cov_su_se = true
) {
  if (y.size() == 0) stop("Empty input.");
  if (use_ar1 && Rcpp::NumericVector::is_na(ar1_rho) == false && std::fabs(ar1_rho) >= 0.999)
    stop("ar1_rho must be in (-0.999, 0.999).");

  // 1) Build paired differences
  PairData P = make_pairs(y, subject, method, time);

  // 2) Subject indexing over PAIRS (not raw rows)
  std::vector<int> subj_idx;
  int m = 0;
  reindex(P.subj, subj_idx, m);
  BySubjBA S = index_by_subject_ba(subj_idx, P.time);

  // 3) Fixed-effects design: [Intercept, (optional) scaled pair mean]
  const int p = include_slope ? 2 : 1;
  mat X(P.d.size(), p, fill::zeros);
  for (size_t i = 0; i < P.d.size(); ++i) X(i, 0) = 1.0;

  vec x2;
  double x2_mean = 0.0, x2_scale = 1.0;
  if (include_slope) {
    x2 = vec(P.mean.data(), P.mean.size(), /*copy*/ false);
    x2_mean = arma::mean(x2);
    x2_scale = ba_rm_slope_scale_denom(P.mean);
    vec x2c = (x2 - x2_mean) / x2_scale;
    X.col(1) = x2c;
  }
  vec ydiff(P.d.data(), P.d.size(), /*copy*/ false);

  // Strategy:
  //   - Build the requested covariance structure.
  //   - Try the interior EM fit.
  //   - Also compute the exact boundary fit with sigma2_subject = 0.
  //   - Select the better fit by profiled REML.
  //   - If AR1 requested & rho is NA, estimate rho from iid first, then repeat.
  bool ar1_estimated = false;
  double rho_used = NA_REAL;
  FitOut fit;

  auto select_best_fit = [&](const std::vector<Precomp>& PC_local,
                             const FitOut* fit_interior_or_null) -> FitOut {
                               FitOut fit_boundary = fit_diff_boundary_su0(X, ydiff, PC_local);
                               const double ll_boundary = ba_rm_reml_loglik(X, ydiff, PC_local, fit_boundary);

                               if (fit_interior_or_null == nullptr) {
                                 return fit_boundary;
                               }

                               const double ll_interior =
                                 ba_rm_reml_loglik(X, ydiff, PC_local, *fit_interior_or_null);

                               if (!std::isfinite(ll_interior)) {
                                 return fit_boundary;
                               }
                               if (!std::isfinite(ll_boundary)) {
                                 return *fit_interior_or_null;
                               }

                               // prefer boundary if it is at least as good up to a tiny tolerance
                               if (ll_boundary >= ll_interior - 1e-8) {
                                 return fit_boundary;
                               }
                               return *fit_interior_or_null;
                             };

  if (use_ar1 && Rcpp::NumericVector::is_na(ar1_rho)) {
    // First pass: iid
    std::vector<Precomp> PC_iid =
      precompute_blocks(X, ydiff, S, /*use_ar1*/ false, 0.0);

    FitOut fit_iid_interior;
    FitOut* fit_iid_ptr = nullptr;
    try {
      fit_iid_interior = fit_diff_em(X, ydiff, S, PC_iid, max_iter, tol);
      fit_iid_ptr = &fit_iid_interior;
    } catch (...) {
      fit_iid_ptr = nullptr;
    }
    FitOut fit_iid = select_best_fit(PC_iid, fit_iid_ptr);

    // Estimate rho from the selected iid fit
    double rho_mom = estimate_rho_moments(fit_iid, X, ydiff, S, PC_iid);
    double rho_hat = estimate_rho_profile(X, ydiff, S, max_iter, tol, rho_mom);
    rho_used = rho_hat;
    ar1_estimated = true;

    // Second pass: AR1 with rho_hat
    std::vector<Precomp> PC_ar1 =
      precompute_blocks(X, ydiff, S, /*use_ar1*/ true, rho_hat);

    FitOut fit_ar1_interior;
    FitOut* fit_ar1_ptr = nullptr;
    try {
      fit_ar1_interior = fit_diff_em(X, ydiff, S, PC_ar1, max_iter, tol);
      fit_ar1_ptr = &fit_ar1_interior;
    } catch (...) {
      fit_ar1_ptr = nullptr;
    }
    fit = select_best_fit(PC_ar1, fit_ar1_ptr);

  } else {
    if (use_ar1) {
      rho_used = nan_preserve(ar1_rho, -0.999, 0.999);
    }

    std::vector<Precomp> PC =
      precompute_blocks(X, ydiff, S, use_ar1, (use_ar1 ? rho_used : 0.0));

    FitOut fit_interior;
    FitOut* fit_ptr = nullptr;
    try {
      fit_interior = fit_diff_em(X, ydiff, S, PC, max_iter, tol);
      fit_ptr = &fit_interior;
    } catch (...) {
      fit_ptr = nullptr;
    }
    fit = select_best_fit(PC, fit_ptr);
  }

  vec beta = fit.beta;
  mat beta_vcov = fit.beta_vcov;

  // Keep both parameterisations:
  // - beta_center: fitted mean at the centred pair mean (x2c = 0)
  // - beta0_orig / beta1_orig: intercept/slope on the original pair-mean scale
  const double beta_center = beta[0];
  double beta0_orig = beta[0];
  double beta1_orig = (p == 2 ? beta[1] : NA_REAL);

  if (include_slope) {
    beta1_orig = beta[1] / x2_scale;
    beta0_orig = beta[0] - beta1_orig * x2_mean;
  }

  const double su2 = fit.su2;
  const double se2 = fit.se2;

  // --- Model-based BA centre ---
  // For include_slope = FALSE, this is the model intercept.
  // For include_slope = TRUE, this is the fitted mean difference at the average paired mean,
  // because x2 was centred before fitting.
  const double mu0 = beta_center;

  // 7) LoA setup
  const double alpha = 1.0 - std::min(std::max(conf_level, 0.0), 1.0);
  const double z = R::qnorm(1.0 - 0.5 * alpha, 0.0, 1.0, 1, 0);
  double loa_multiplier = loa_multiplier_arg;
  if (!std::isfinite(loa_multiplier) || loa_multiplier <= 0.0) {
    loa_multiplier = z;
  }

  // --- Main LoA point estimate: model-based predictive variance ---
  const double V_loa_model = nan_preserve(su2 + se2, 0.0, 1e12);
  const double sd_loa = std::sqrt(V_loa_model);
  const double loa_lower = mu0 - loa_multiplier * sd_loa;
  const double loa_upper = mu0 + loa_multiplier * sd_loa;

  // --- Empirical paired-difference LoA variance retained as diagnostic only ---
  const double V_loa_empirical =
    nan_preserve(loa_var_subject_equal_weight(ydiff, S, mu0), 0.0, 1e12);
  const double sd_loa_empirical = std::sqrt(V_loa_empirical);

  // 8) Profile-REML / delta-method CIs
  //
  // Decomposition:
  //   Var{T(hat)} = E[ Var{T(hat) | theta_hat} ] + Var{ E[T(hat) | theta_hat] }
  //
  // Here:
  //   - conditional beta uncertainty comes from beta_vcov at fixed variance parameters;
  //   - profile uncertainty in theta = (log su2, log se2, [z_rho]) comes from the
  //     observed Hessian of the profiled REML log-likelihood;
  //   - the second term is propagated numerically by central finite differences on the
  //     derived quantities mu0, LoA_lower, LoA_upper.

  const double eps_theta = std::max(1e-10, 1e-10 * std::max(1.0, arma::var(ydiff, /*unbiased*/1)));
  const bool rho_free = (use_ar1 && ar1_estimated);

  const bool su2_boundary = fit.boundary_su2_zero;

  arma::vec theta_hat(
      (su2_boundary ? 1 : 2) + (rho_free ? 1 : 0),
      arma::fill::zeros
  );

  int th_pos = 0;
  if (!su2_boundary) {
    theta_hat[th_pos++] = std::log(std::max(su2, eps_theta));
  }
  theta_hat[th_pos++] = std::log(std::max(se2, eps_theta));
  if (rho_free) {
    theta_hat[th_pos++] = ba_rm_z_from_rho(rho_used);
  }

  arma::vec step(theta_hat.n_elem, arma::fill::zeros);
  for (arma::uword j = 0; j < theta_hat.n_elem; ++j) {
    step[j] = 1e-4 * std::max(1.0, std::fabs(theta_hat[j]));
  }

  auto loglik_fn = [&](const arma::vec& th, double& out) -> bool {
    ProfileThetaEval E;
    if (!ba_rm_eval_profile_theta(th, X, ydiff, S, use_ar1, rho_free, rho_used, loa_multiplier, E)) return false;
    out = E.loglik;
    return std::isfinite(out);
  };

  auto mu_fn = [&](const arma::vec& th, double& out) -> bool {
    ProfileThetaEval E;
    if (!ba_rm_eval_profile_theta(th, X, ydiff, S, use_ar1, rho_free, rho_used, loa_multiplier, E)) return false;
    out = E.mu0;
    return std::isfinite(out);
  };

  auto lower_fn = [&](const arma::vec& th, double& out) -> bool {
    ProfileThetaEval E;
    if (!ba_rm_eval_profile_theta(th, X, ydiff, S, use_ar1, rho_free, rho_used, loa_multiplier, E)) return false;
    out = E.loa_lower;
    return std::isfinite(out);
  };

  auto upper_fn = [&](const arma::vec& th, double& out) -> bool {
    ProfileThetaEval E;
    if (!ba_rm_eval_profile_theta(th, X, ydiff, S, use_ar1, rho_free, rho_used, loa_multiplier, E)) return false;
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

  // Conditional beta uncertainty at fixed theta:
  // mu0 = beta_center = beta[0], and each LoA endpoint depends on beta only through beta[0].
  const double var_mu_cond = (beta_vcov.n_rows >= 1 ? std::max(0.0, beta_vcov(0, 0)) : 0.0);
  const double var_lower_cond = var_mu_cond;
  const double var_upper_cond = var_mu_cond;

  // Profile uncertainty in theta, propagated through the derived quantities.
  const double var_mu_prof =
    std::max(0.0, arma::as_scalar(g_mu.t() * Sigma_theta * g_mu));
  const double var_lower_prof =
    std::max(0.0, arma::as_scalar(g_lower.t() * Sigma_theta * g_lower));
  const double var_upper_prof =
    std::max(0.0, arma::as_scalar(g_upper.t() * Sigma_theta * g_upper));

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

    // model-based BA centre and CI
    _["bias_mu0"] = mu0,
    _["bias_se"] = se_bias,
    _["bias_lwr"] = bias_lwr,
    _["bias_upr"] = bias_upr,

    // fixed effects on original pair-mean scale
    _["beta_intercept"] = beta0_orig,
    _["beta_slope"] = (include_slope ? beta1_orig : NA_REAL),
    _["beta_center"] = beta_center,
    _["beta_center_reference_mean"] = (include_slope ? x2_mean : NA_REAL),

    // variance components
    _["sigma2_subject"] = su2,
    _["sigma2_resid"] = se2,

    // LoA (model-based main result)
    _["sd_loa"] = sd_loa,
    _["loa_var_model"] = V_loa_model,
    _["loa_lower"] = loa_lower,
    _["loa_upper"] = loa_upper,
    _["loa_lower_lwr"] = loa_lower_lwr,
    _["loa_lower_upr"] = loa_lower_upr,
    _["loa_upper_lwr"] = loa_upper_lwr,
    _["loa_upper_upr"] = loa_upper_upr,

    // empirical paired-difference diagnostics
    _["loa_var_empirical"] = V_loa_empirical,
    _["sd_loa_empirical"] = sd_loa_empirical,

    // AR(1)
    _["use_ar1"] = use_ar1,
    _["ar1_rho"] = (use_ar1 ? rho_used : NA_REAL),
    _["ar1_estimated"] = (use_ar1 ? ar1_estimated : false),

    // misc
    _["loa_multiplier"] = loa_multiplier,
    _["conf_level"] = conf_level,
    _["converged"] = fit.converged,
    _["iter"] = fit.iter,
    _["warn"] = fit.warn
  );
}
