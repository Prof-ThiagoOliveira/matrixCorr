// Thiago de Paula Oliveira
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R_ext/Applic.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <limits>
#include <vector>

#include "matrixCorr_detail.h"

using matrixCorr_detail::clamp_policy::nan_preserve;
using matrixCorr_detail::norm1::Phi;
using matrixCorr_detail::norm1::phi;
using matrixCorr_detail::norm1::qnorm01;
using matrixCorr_detail::brent::optimize;
using matrixCorr_detail::bvn_fast::Phi2;
using matrixCorr_detail::bvn_fast::rect_prob;

using namespace Rcpp;

namespace {

struct PolyserialPrepared {
  std::vector<double> x;
  std::vector<int> y;
};

struct PolyserialFitDetails {
  std::vector<double> z;
  std::vector<int> y;
  std::vector<double> par;
  double estimate{NA_REAL};
  double fmin{NA_REAL};
  int ncat{0};
  int n_obs{0};
};

struct PolyserialOptimData {
  const std::vector<double>* z;
  const std::vector<int>* y;
  double maxcor;
};

inline double polyserial_penalty() {
  return std::numeric_limits<double>::infinity();
}

bool polyserial_prepare_xy(const NumericVector& x,
                           const IntegerVector& y,
                           bool pairwise_complete,
                           PolyserialPrepared& out) {
  const int n = x.size();
  if (n != y.size()) return false;

  out.x.clear();
  out.y.clear();
  out.x.reserve(n);
  out.y.reserve(n);

  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    const int yi = y[i];
    const bool bad_x = NumericVector::is_na(xi) || !std::isfinite(xi);
    const bool bad_y = IntegerVector::is_na(yi);
    if (bad_x || bad_y) {
      if (pairwise_complete) continue;
      return false;
    }
    out.x.push_back(xi);
    out.y.push_back(yi);
  }

  const int n_complete = static_cast<int>(out.x.size());
  if (n_complete < 2) return false;

  std::vector<int> levels = out.y;
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  if (levels.size() < 2) return false;

  for (int& yi : out.y) {
    yi = static_cast<int>(std::lower_bound(levels.begin(), levels.end(), yi) - levels.begin()) + 1;
  }

  return true;
}

double polyserial_negloglik_core(const std::vector<double>& z,
                                 const std::vector<int>& y,
                                 const double* pars,
                                 const int ncat,
                                 const double maxcor) {
  if (ncat < 2) return NA_REAL;

  double rho = pars[0];
  if (std::abs(rho) > maxcor) rho = (rho > 0.0 ? maxcor : -maxcor);
  if (std::abs(rho) >= 1.0) return polyserial_penalty();

  const double denom_sq = 1.0 - rho * rho;
  if (!(denom_sq > 0.0) || !std::isfinite(denom_sq)) return polyserial_penalty();
  const double denom = std::sqrt(denom_sq);

  double out = 0.0;
  for (int i = 0, n = static_cast<int>(z.size()); i < n; ++i) {
    const int yi = y[i];
    if (yi < 1 || yi > ncat) return polyserial_penalty();

    double lower_cut = (yi == 1) ? -INFINITY : pars[yi - 1];
    double upper_cut = (yi == ncat) ? INFINITY : pars[yi];
    if (upper_cut < lower_cut) return polyserial_penalty();

    const double upper = (upper_cut - rho * z[i]) / denom;
    const double lower = (lower_cut - rho * z[i]) / denom;
    const double p = Phi(upper) - Phi(lower);
    if (!(p > 0.0) || !std::isfinite(p)) return polyserial_penalty();

    const double dens = phi(z[i]);
    if (!(dens > 0.0) || !std::isfinite(dens)) return polyserial_penalty();

    out -= std::log(dens * p);
  }

  return out;
}

double polyserial_optim_fn(int n, double* par, void* ex) {
  PolyserialOptimData* data = static_cast<PolyserialOptimData*>(ex);
  return polyserial_negloglik_core(*data->z, *data->y, par, n, data->maxcor);
}

double polyserial_fit_full_mle(const std::vector<double>& x,
                               const std::vector<int>& y,
                               const double maxcor = 0.9999);

bool polyserial_fit_full_details(const std::vector<double>& x,
                                 const std::vector<int>& y,
                                 PolyserialFitDetails& fit,
                                 const double maxcor = 0.9999) {
  const int n = static_cast<int>(x.size());
  fit = PolyserialFitDetails{};
  if (n != static_cast<int>(y.size()) || n < 2) return false;

  double mx = 0.0;
  double my = 0.0;
  for (int i = 0; i < n; ++i) {
    mx += x[i];
    my += y[i];
  }
  mx /= n;
  my /= n;

  double sxx = 0.0;
  double syy = 0.0;
  double sxy = 0.0;
  std::vector<double> z(n);
  for (int i = 0; i < n; ++i) {
    const double dx = x[i] - mx;
    const double dy = y[i] - my;
    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
    z[i] = dx;
  }

  if (!(sxx > 0.0) || !(syy > 0.0)) return false;

  const double sx = std::sqrt(sxx / (n - 1));
  const double sy = std::sqrt(syy / (n - 1));
  if (!(sx > 0.0) || !(sy > 0.0)) return false;

  for (double& zi : z) zi /= sx;

  const int ncat = *std::max_element(y.begin(), y.end());
  if (ncat < 2) return false;

  std::vector<double> freq(ncat, 0.0);
  for (int yi : y) freq[yi - 1] += 1.0;
  if (std::count_if(freq.begin(), freq.end(), [](double f) { return f > 0.0; }) < 2) {
    return false;
  }

  std::vector<double> par(ncat);
  double acc = 0.0;
  for (int k = 0; k < ncat - 1; ++k) {
    acc += freq[k] / n;
    par[k + 1] = qnorm01(nan_preserve(acc, 1e-12, 1.0 - 1e-12));
  }

  double denom = 0.0;
  for (int k = 1; k < ncat; ++k) denom += phi(par[k]);
  if (!(denom > 0.0) || !std::isfinite(denom)) return false;

  double rho = std::sqrt(static_cast<double>(n - 1) / n) * sy * (sxy / std::sqrt(sxx * syy)) / denom;
  if (std::abs(rho) > maxcor) rho = (rho > 0.0 ? maxcor : -maxcor);
  par[0] = rho;

  std::vector<double> opt_par = par;
  double fmin = polyserial_penalty();
  int fail = 0;
  int fncount = 0;
  PolyserialOptimData opt_data{&z, &y, maxcor};
  nmmin(
    ncat,
    par.data(),
    opt_par.data(),
    &fmin,
    polyserial_optim_fn,
    &fail,
    R_NegInf,
    std::sqrt(std::numeric_limits<double>::epsilon()),
    &opt_data,
    1.0,
    0.5,
    2.0,
    0,
    &fncount,
    500
  );

  if (fail != 0 || !std::isfinite(opt_par[0])) return false;

  fit.z = std::move(z);
  fit.y = y;
  fit.par = std::move(opt_par);
  fit.estimate = nan_preserve(fit.par[0], -maxcor, maxcor);
  fit.fmin = fmin;
  fit.ncat = ncat;
  fit.n_obs = n;
  fit.par[0] = fit.estimate;
  return true;
}

double polyserial_fit_full_mle(const std::vector<double>& x,
                               const std::vector<int>& y,
                               const double maxcor) {
  PolyserialFitDetails fit;
  if (!polyserial_fit_full_details(x, y, fit, maxcor)) return NA_REAL;
  return fit.estimate;
}

} // namespace

static inline NumericMatrix drop_zero_margins(const NumericMatrix& tab) {
  const int nr = tab.nrow();
  const int nc = tab.ncol();
  std::vector<int> keep_r;
  std::vector<int> keep_c;
  keep_r.reserve(nr);
  keep_c.reserve(nc);

  for (int i = 0; i < nr; ++i) {
    double rs = 0.0;
    for (int j = 0; j < nc; ++j) rs += tab(i, j);
    if (rs > 0.0) keep_r.push_back(i);
  }
  for (int j = 0; j < nc; ++j) {
    double cs = 0.0;
    for (int i = 0; i < nr; ++i) cs += tab(i, j);
    if (cs > 0.0) keep_c.push_back(j);
  }

  NumericMatrix out(keep_r.size(), keep_c.size());
  for (int i = 0; i < static_cast<int>(keep_r.size()); ++i) {
    for (int j = 0; j < static_cast<int>(keep_c.size()); ++j) {
      out(i, j) = tab(keep_r[i], keep_c[j]);
    }
  }
  return out;
}

static inline void cutpoints_from_table(const NumericMatrix& tab,
                                        std::vector<double>& rc,
                                        std::vector<double>& cc) {
  const int nr = tab.nrow();
  const int nc = tab.ncol();
  rc.clear();
  cc.clear();
  if (nr < 2 || nc < 2) return;

  NumericVector rsum(nr), csum(nc);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      rsum[i] += tab(i, j);
      csum[j] += tab(i, j);
    }
  }

  double acc = 0.0;
  rc.reserve(nr - 1);
  for (int i = 0; i < nr - 1; ++i) {
    acc += rsum[i];
    rc.push_back(qnorm01(nan_preserve(acc, 1e-12, 1.0 - 1e-12)));
  }

  acc = 0.0;
  cc.reserve(nc - 1);
  for (int j = 0; j < nc - 1; ++j) {
    acc += csum[j];
    cc.push_back(qnorm01(nan_preserve(acc, 1e-12, 1.0 - 1e-12)));
  }
}

static inline void cutpoints_from_rowmajor(const std::vector<double>& tab,
                                           int nr, int nc,
                                           std::vector<double>& rc,
                                           std::vector<double>& cc) {
  rc.clear();
  cc.clear();
  if (nr < 2 || nc < 2) return;

  std::vector<double> rsum(nr, 0.0), csum(nc, 0.0);
  for (int i = 0; i < nr; ++i) {
    const int off = i * nc;
    for (int j = 0; j < nc; ++j) {
      const double v = tab[off + j];
      rsum[i] += v;
      csum[j] += v;
    }
  }

  double acc = 0.0;
  rc.reserve(nr - 1);
  for (int i = 0; i < nr - 1; ++i) {
    acc += rsum[i];
    rc.push_back(qnorm01(nan_preserve(acc, 1e-12, 1.0 - 1e-12)));
  }

  acc = 0.0;
  cc.reserve(nc - 1);
  for (int j = 0; j < nc - 1; ++j) {
    acc += csum[j];
    cc.push_back(qnorm01(nan_preserve(acc, 1e-12, 1.0 - 1e-12)));
  }
}

struct PolyWorkspace {
  int nr{0};
  int nc{0};
  std::vector<double> rc;
  std::vector<double> cc;
  std::vector<double> rc_phi;
  std::vector<double> cc_phi;
  std::vector<double> prob;

  void reset(const std::vector<double>& rc_in, const std::vector<double>& cc_in) {
    rc = rc_in;
    cc = cc_in;
    nr = static_cast<int>(rc.size()) + 1;
    nc = static_cast<int>(cc.size()) + 1;
    rc_phi.resize(rc.size());
    cc_phi.resize(cc.size());
    for (int i = 0; i < static_cast<int>(rc.size()); ++i) rc_phi[i] = Phi(rc[i]);
    for (int j = 0; j < static_cast<int>(cc.size()); ++j) cc_phi[j] = Phi(cc[j]);
    prob.assign(static_cast<std::size_t>(nr * nc), 0.0);
  }
};

static inline double tetra_nll_fixed_thresholds(double rho,
                                                double a, double b, double c, double d,
                                                double rc, double cc) {
  const double p11 = Phi2(rc, cc, rho);
  const double p21 = Phi(rc) - p11;
  const double p12 = Phi(cc) - p11;
  const double p22 = 1.0 - Phi(rc) - p12;

  if ((a > 0.0 && (!(p11 > 0.0) || !std::isfinite(p11))) ||
      (b > 0.0 && (!(p21 > 0.0) || !std::isfinite(p21))) ||
      (c > 0.0 && (!(p12 > 0.0) || !std::isfinite(p12))) ||
      (d > 0.0 && (!(p22 > 0.0) || !std::isfinite(p22)))) {
    return R_PosInf;
  }

  return -(a * std::log(p11) + b * std::log(p21) +
    c * std::log(p12) + d * std::log(p22));
}

static inline void poly_prob_from_thresholds(double rho, PolyWorkspace& ws) {
  const int nr = ws.nr;
  const int nc = ws.nc;
  std::fill(ws.prob.begin(), ws.prob.end(), 0.0);

  if (nr == 2 && nc == 2) {
    const double p11 = Phi2(ws.rc[0], ws.cc[0], rho);
    ws.prob[0] = p11;
    ws.prob[1] = ws.rc_phi[0] - p11;
    ws.prob[2] = ws.cc_phi[0] - p11;
    ws.prob[3] = 1.0 - ws.rc_phi[0] - ws.prob[2];
    return;
  }

  for (int i = 0; i < nr - 1; ++i) {
    const double al = (i == 0) ? -INFINITY : ws.rc[i - 1];
    const double au = ws.rc[i];
    const int off = i * nc;
    for (int j = 0; j < nc - 1; ++j) {
      const double bl = (j == 0) ? -INFINITY : ws.cc[j - 1];
      const double bu = ws.cc[j];
      ws.prob[off + j] = rect_prob(al, au, bl, bu, rho);
    }
  }

  double sum_row0 = 0.0;
  for (int j = 0; j < nc - 1; ++j) sum_row0 += ws.prob[j];
  ws.prob[nc - 1] = ws.rc_phi[0] - sum_row0;

  double sum_col0 = 0.0;
  for (int i = 0; i < nr - 1; ++i) sum_col0 += ws.prob[i * nc];
  ws.prob[(nr - 1) * nc] = ws.cc_phi[0] - sum_col0;

  for (int i = 1; i < nr - 1; ++i) {
    double row_sum = 0.0;
    const int off = i * nc;
    for (int j = 0; j < nc - 1; ++j) row_sum += ws.prob[off + j];
    ws.prob[off + nc - 1] = ws.rc_phi[i] - ws.rc_phi[i - 1] - row_sum;
  }

  for (int j = 1; j < nc - 1; ++j) {
    double col_sum = 0.0;
    for (int i = 0; i < nr - 1; ++i) col_sum += ws.prob[i * nc + j];
    ws.prob[(nr - 1) * nc + j] = ws.cc_phi[j] - ws.cc_phi[j - 1] - col_sum;
  }

  double last_row_sum = 0.0;
  const int last_row = (nr - 1) * nc;
  for (int j = 0; j < nc - 1; ++j) last_row_sum += ws.prob[last_row + j];
  ws.prob[last_row + nc - 1] = 1.0 - ws.rc_phi[nr - 2] - last_row_sum;
}

static inline double poly_nll_fixed_thresholds(double rho,
                                               const double* tab,
                                               PolyWorkspace& ws) {
  poly_prob_from_thresholds(rho, ws);

  double out = 0.0;
  const int n = ws.nr * ws.nc;
  for (int idx = 0; idx < n; ++idx) {
    const double nij = tab[idx];
    if (nij <= 0.0) continue;
    const double pij = ws.prob[idx];
    if (!(pij > 0.0) || !std::isfinite(pij)) continue;
    out -= nij * std::log(pij);
  }
  return out;
}

// [[Rcpp::export]]
double matrixCorr_tetrachoric_mle_cpp(Rcpp::NumericMatrix tab, double correct = 0.5) {
  if (tab.nrow() != 2 || tab.ncol() != 2) return NA_REAL;

  double a = tab(0, 0), b = tab(1, 0), c = tab(0, 1), d = tab(1, 1);
  if (a + b + c + d <= 0.0) return NA_REAL;

  const bool any_zero = (a <= 0.0) || (b <= 0.0) || (c <= 0.0) || (d <= 0.0);
  if (correct > 0.0 && any_zero) {
    if (a <= 0.0) a += correct;
    if (b <= 0.0) b += correct;
    if (c <= 0.0) c += correct;
    if (d <= 0.0) d += correct;
  }

  const double tot = a + b + c + d;
  a /= tot; b /= tot; c /= tot; d /= tot;
  const double rc = qnorm01(nan_preserve(a + c, 1e-12, 1.0 - 1e-12));
  const double cc = qnorm01(nan_preserve(a + b, 1e-12, 1.0 - 1e-12));

  auto nll = [&](double rho) { return tetra_nll_fixed_thresholds(rho, a, b, c, d, rc, cc); };
  const double est = optimize(nll, -1.0, 1.0, 1e-8, 250);
  return nan_preserve(est, -1.0, 1.0);
}

// [[Rcpp::export]]
double matrixCorr_tetrachoric_fixed_cpp(Rcpp::NumericMatrix tab, double rc, double cc, double correct = 0.5) {
  if (tab.nrow() != 2 || tab.ncol() != 2) return NA_REAL;

  double a = tab(0, 0), b = tab(1, 0), c = tab(0, 1), d = tab(1, 1);
  if (a + b + c + d <= 0.0) return NA_REAL;

  const bool any_zero = (a <= 0.0) || (b <= 0.0) || (c <= 0.0) || (d <= 0.0);
  if (correct > 0.0 && any_zero) {
    if (a <= 0.0) a += correct;
    if (b <= 0.0) b += correct;
    if (c <= 0.0) c += correct;
    if (d <= 0.0) d += correct;
  }

  auto nll = [&](double rho) { return tetra_nll_fixed_thresholds(rho, a, b, c, d, rc, cc); };
  const double est = optimize(nll, -1.0, 1.0, 1e-8, 250);
  return nan_preserve(est, -1.0, 1.0);
}

// [[Rcpp::export]]
double matrixCorr_polychoric_mle_cpp(NumericMatrix tab, double correct = 0.5) {
  if (tab.nrow() < 2 || tab.ncol() < 2) return NA_REAL;

  NumericMatrix N = clone(tab);
  double tot = 0.0;
  for (int i = 0; i < N.nrow(); ++i) {
    for (int j = 0; j < N.ncol(); ++j) {
      tot += N(i, j);
    }
  }
  if (tot <= 0.0) return NA_REAL;

  for (int i = 0; i < N.nrow(); ++i) {
    for (int j = 0; j < N.ncol(); ++j) {
      N(i, j) /= tot;
    }
  }
  N = drop_zero_margins(N);
  if (N.nrow() < 2 || N.ncol() < 2) return NA_REAL;

  if (correct > 0.0) {
    const double add = correct / tot;
    for (int i = 0; i < N.nrow(); ++i) {
      for (int j = 0; j < N.ncol(); ++j) {
        if (N(i, j) <= 0.0) N(i, j) = add;
      }
    }
  }

  std::vector<double> rc, cc;
  cutpoints_from_table(N, rc, cc);
  if (rc.empty() || cc.empty()) return NA_REAL;

  PolyWorkspace ws;
  ws.reset(rc, cc);
  std::vector<double> tab_vec(static_cast<std::size_t>(N.nrow() * N.ncol()));
  for (int i = 0; i < N.nrow(); ++i) {
    for (int j = 0; j < N.ncol(); ++j) {
      tab_vec[static_cast<std::size_t>(i * N.ncol() + j)] = N(i, j);
    }
  }

  auto nll = [&](double rho) { return poly_nll_fixed_thresholds(rho, tab_vec.data(), ws); };
  const double est = optimize(nll, -1.0, 1.0, 1e-8, 250);
  return nan_preserve(est, -1.0, 1.0);
}

// [[Rcpp::export]]
double matrixCorr_polychoric_fixed_cpp(NumericMatrix tab, NumericVector rc_in, NumericVector cc_in) {
  if (tab.nrow() < 2 || tab.ncol() < 2) return NA_REAL;

  NumericMatrix N = clone(tab);
  double tot = 0.0;
  for (int i = 0; i < N.nrow(); ++i) {
    for (int j = 0; j < N.ncol(); ++j) {
      tot += N(i, j);
    }
  }
  if (tot <= 0.0) return NA_REAL;

  for (int i = 0; i < N.nrow(); ++i) {
    for (int j = 0; j < N.ncol(); ++j) {
      N(i, j) /= tot;
    }
  }

  if (rc_in.size() != N.nrow() - 1 || cc_in.size() != N.ncol() - 1) {
    return NA_REAL;
  }
  std::vector<double> rc(rc_in.size()), cc(cc_in.size());
  std::copy(rc_in.begin(), rc_in.end(), rc.begin());
  std::copy(cc_in.begin(), cc_in.end(), cc.begin());

  PolyWorkspace ws;
  ws.reset(rc, cc);
  std::vector<double> tab_vec(static_cast<std::size_t>(N.nrow() * N.ncol()));
  for (int i = 0; i < N.nrow(); ++i) {
    for (int j = 0; j < N.ncol(); ++j) {
      tab_vec[static_cast<std::size_t>(i * N.ncol() + j)] = N(i, j);
    }
  }

  auto nll = [&](double rho) { return poly_nll_fixed_thresholds(rho, tab_vec.data(), ws); };
  const double est = optimize(nll, -1.0, 1.0, 1e-8, 250);
  return nan_preserve(est, -1.0, 1.0);
}

static inline double tetra_from_counts_fixed(const double* tab,
                                             double rc, double cc,
                                             double correct) {
  double a = tab[0], b = tab[2], c = tab[1], d = tab[3];
  if (a + b + c + d <= 0.0) return NA_REAL;

  const bool any_zero = (a <= 0.0) || (b <= 0.0) || (c <= 0.0) || (d <= 0.0);
  if (correct > 0.0 && any_zero) {
    if (a <= 0.0) a += correct;
    if (b <= 0.0) b += correct;
    if (c <= 0.0) c += correct;
    if (d <= 0.0) d += correct;
  }

  auto nll = [&](double rho) { return tetra_nll_fixed_thresholds(rho, a, b, c, d, rc, cc); };
  const double est = optimize(nll, -1.0, 1.0, 1e-8, 250);
  return nan_preserve(est, -1.0, 1.0);
}

struct PolyPairBuffers {
  std::vector<double> row_sum;
  std::vector<double> col_sum;
  std::vector<int> keep_r;
  std::vector<int> keep_c;
  std::vector<double> tab;
  std::vector<double> rc;
  std::vector<double> cc;
  PolyWorkspace ws;
};

static inline double poly_from_counts_local(const double* raw, int nr, int nc,
                                            double correct, PolyPairBuffers& buf) {
  double tot = 0.0;
  buf.row_sum.assign(nr, 0.0);
  buf.col_sum.assign(nc, 0.0);
  for (int i = 0; i < nr; ++i) {
    const int off = i * nc;
    for (int j = 0; j < nc; ++j) {
      const double v = raw[off + j];
      tot += v;
      buf.row_sum[i] += v;
      buf.col_sum[j] += v;
    }
  }
  if (tot <= 0.0) return NA_REAL;

  buf.keep_r.clear();
  buf.keep_c.clear();
  for (int i = 0; i < nr; ++i) if (buf.row_sum[i] > 0.0) buf.keep_r.push_back(i);
  for (int j = 0; j < nc; ++j) if (buf.col_sum[j] > 0.0) buf.keep_c.push_back(j);
  const int nr2 = static_cast<int>(buf.keep_r.size());
  const int nc2 = static_cast<int>(buf.keep_c.size());
  if (nr2 < 2 || nc2 < 2) return NA_REAL;

  buf.tab.assign(static_cast<std::size_t>(nr2 * nc2), 0.0);
  for (int i = 0; i < nr2; ++i) {
    const int src_off = buf.keep_r[i] * nc;
    const int dst_off = i * nc2;
    for (int j = 0; j < nc2; ++j) {
      buf.tab[dst_off + j] = raw[src_off + buf.keep_c[j]] / tot;
    }
  }

  if (correct > 0.0) {
    const double add = correct / tot;
    for (double& v : buf.tab) {
      if (v <= 0.0) v = add;
    }
  }

  cutpoints_from_rowmajor(buf.tab, nr2, nc2, buf.rc, buf.cc);
  if (buf.rc.empty() || buf.cc.empty()) return NA_REAL;

  buf.ws.reset(buf.rc, buf.cc);
  auto nll = [&](double rho) { return poly_nll_fixed_thresholds(rho, buf.tab.data(), buf.ws); };
  const double est = optimize(nll, -1.0, 1.0, 1e-8, 250);
  return nan_preserve(est, -1.0, 1.0);
}

static inline double poly_from_counts_fixed(const double* raw, int nr, int nc,
                                            const double* rc, const double* cc) {
  double tot = 0.0;
  for (int idx = 0; idx < nr * nc; ++idx) tot += raw[idx];
  if (tot <= 0.0) return NA_REAL;

  std::vector<double> tab(static_cast<std::size_t>(nr * nc));
  for (int idx = 0; idx < nr * nc; ++idx) tab[idx] = raw[idx] / tot;
  std::vector<double> rc_vec(rc, rc + nr - 1);
  std::vector<double> cc_vec(cc, cc + nc - 1);
  PolyWorkspace ws;
  ws.reset(rc_vec, cc_vec);
  auto nll = [&](double rho) { return poly_nll_fixed_thresholds(rho, tab.data(), ws); };
  const double est = optimize(nll, -1.0, 1.0, 1e-8, 250);
  return nan_preserve(est, -1.0, 1.0);
}

// [[Rcpp::export]]
NumericMatrix matrixCorr_tetrachoric_matrix_cpp(IntegerMatrix x, NumericVector tau,
                                                double correct = 0.5,
                                                bool pairwise_complete = true) {
  const int n = x.nrow();
  const int p = x.ncol();
  NumericMatrix out(p, p);
  std::fill(out.begin(), out.end(), NA_REAL);

  for (int j = 0; j < p; ++j) {
    bool seen0 = false, seen1 = false;
    for (int i = 0; i < n; ++i) {
      const int v = x(i, j);
      if (v == NA_INTEGER) continue;
      if (v == 1) seen0 = true;
      else if (v == 2) seen1 = true;
      if (seen0 && seen1) break;
    }
    out(j, j) = (seen0 && seen1) ? 1.0 : NA_REAL;
  }

  double tab[4];
  for (int j = 0; j < p - 1; ++j) {
    for (int k = j + 1; k < p; ++k) {
      std::fill(tab, tab + 4, 0.0);
      bool seen_x0 = false, seen_x1 = false, seen_y0 = false, seen_y1 = false;
      for (int i = 0; i < n; ++i) {
        const int xv = x(i, j);
        const int yv = x(i, k);
        if (xv == NA_INTEGER || yv == NA_INTEGER) {
          if (pairwise_complete) continue;
          continue;
        }
        if (xv < 1 || xv > 2 || yv < 1 || yv > 2) continue;
        const int xi = xv - 1;
        const int yi = yv - 1;
        tab[xi * 2 + yi] += 1.0;
        seen_x0 = seen_x0 || (xi == 0);
        seen_x1 = seen_x1 || (xi == 1);
        seen_y0 = seen_y0 || (yi == 0);
        seen_y1 = seen_y1 || (yi == 1);
      }
      if (!(seen_x0 && seen_x1 && seen_y0 && seen_y1)) continue;
      const double est = tetra_from_counts_fixed(tab, tau[k], tau[j], correct);
      out(j, k) = est;
      out(k, j) = est;
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix matrixCorr_polychoric_matrix_cpp(IntegerMatrix x, IntegerVector n_levels,
                                               NumericMatrix tau_mat, bool global_all = false,
                                               double correct = 0.5,
                                               bool pairwise_complete = true) {
  const int n = x.nrow();
  const int p = x.ncol();
  NumericMatrix out(p, p);
  std::fill(out.begin(), out.end(), NA_REAL);

  int max_levels = 0;
  for (int j = 0; j < p; ++j) max_levels = std::max(max_levels, n_levels[j]);
  std::vector<double> raw(static_cast<std::size_t>(max_levels * max_levels), 0.0);
  PolyPairBuffers buf;

  for (int j = 0; j < p; ++j) {
    std::vector<int> seen(static_cast<std::size_t>(n_levels[j]), 0);
    int n_seen = 0;
    for (int i = 0; i < n; ++i) {
      const int v = x(i, j);
      if (v == NA_INTEGER || v < 1 || v > n_levels[j]) continue;
      if (!seen[static_cast<std::size_t>(v - 1)]) {
        seen[static_cast<std::size_t>(v - 1)] = 1;
        n_seen++;
      }
      if (n_seen >= 2) break;
    }
    out(j, j) = (n_seen >= 2) ? 1.0 : NA_REAL;
  }

  for (int j = 0; j < p - 1; ++j) {
    for (int k = j + 1; k < p; ++k) {
      const int nr = n_levels[j];
      const int nc = n_levels[k];
      const int m = nr * nc;
      std::fill(raw.begin(), raw.begin() + m, 0.0);
      std::vector<int> seen_r(static_cast<std::size_t>(nr), 0);
      std::vector<int> seen_c(static_cast<std::size_t>(nc), 0);
      int nr_seen = 0, nc_seen = 0;

      for (int i = 0; i < n; ++i) {
        const int xv = x(i, j);
        const int yv = x(i, k);
        if (xv == NA_INTEGER || yv == NA_INTEGER) {
          if (pairwise_complete) continue;
          continue;
        }
        if (xv < 1 || xv > nr || yv < 1 || yv > nc) continue;
        raw[(xv - 1) * nc + (yv - 1)] += 1.0;
        if (!seen_r[static_cast<std::size_t>(xv - 1)]) {
          seen_r[static_cast<std::size_t>(xv - 1)] = 1;
          nr_seen++;
        }
        if (!seen_c[static_cast<std::size_t>(yv - 1)]) {
          seen_c[static_cast<std::size_t>(yv - 1)] = 1;
          nc_seen++;
        }
      }

      if (nr_seen < 2 || nc_seen < 2) continue;

      double est = NA_REAL;
      if (global_all) {
        est = poly_from_counts_fixed(
          raw.data(), nr, nc,
          &tau_mat(0, j), &tau_mat(0, k)
        );
      } else {
        est = poly_from_counts_local(raw.data(), nr, nc, correct, buf);
      }
      out(j, k) = est;
      out(k, j) = est;
    }
  }

  return out;
}

// [[Rcpp::export]]
double matrixCorr_biserial_latent_cpp(NumericVector x, LogicalVector y) {
  const int n = x.size();
  if (n != y.size() || n < 2) return NA_REAL;

  double mx = 0.0;
  for (int i = 0; i < n; ++i) mx += x[i];
  mx /= n;

  double v = 0.0;
  for (int i = 0; i < n; ++i) {
    const double d = x[i] - mx;
    v += d * d;
  }
  if (v <= 0.0) return NA_REAL;
  const double sx = std::sqrt(v / (n - 1));

  int n1 = 0, n0 = 0;
  double s1 = 0.0, s0 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (y[i]) {
      s1 += x[i];
      n1++;
    } else {
      s0 += x[i];
      n0++;
    }
  }
  if (n1 == 0 || n0 == 0) return NA_REAL;

  const double mean1 = s1 / n1;
  const double mean0 = s0 / n0;
  const double p = nan_preserve(static_cast<double>(n1) / n, 1e-12, 1.0 - 1e-12);
  const double q = 1.0 - p;
  const double zp = phi(qnorm01(p));
  const double r = ((mean1 - mean0) / sx) * (p * q) / zp;
  return nan_preserve(r, -1.0, 1.0);
}

// [[Rcpp::export]]
double matrixCorr_polyserial_mle_cpp(NumericVector x, IntegerVector y) {
  PolyserialPrepared prep;
  if (!polyserial_prepare_xy(x, y, false, prep)) return NA_REAL;
  return polyserial_fit_full_mle(prep.x, prep.y);
}

// [[Rcpp::export]]
double matrixCorr_polyserial_negloglik_cpp(NumericVector z, IntegerVector y,
                                           NumericVector pars, double maxcor = 0.9999) {
  const int n = z.size();
  const int ncat = pars.size();
  if (n != y.size() || ncat < 2) return NA_REAL;

  std::vector<double> z_std(n);
  std::vector<int> y_int(n);
  for (int i = 0; i < n; ++i) {
    z_std[i] = z[i];
    y_int[i] = y[i];
  }

  return polyserial_negloglik_core(z_std, y_int, REAL(pars), ncat, maxcor);
}

namespace {

constexpr double kLatentMaxcor = 0.9999;

inline double rho_from_theta(double theta) {
  return kLatentMaxcor * std::tanh(theta);
}

inline double drho_dtheta(double theta) {
  const double th = std::tanh(theta);
  return kLatentMaxcor * (1.0 - th * th);
}

double latent_safe_step(const std::vector<double>& theta,
                        int idx,
                        const std::function<double(const std::vector<double>&)>& fn) {
  double h = std::pow(std::numeric_limits<double>::epsilon(), 0.25) *
    std::max(1.0, std::abs(theta[static_cast<std::size_t>(idx)]));
  for (int it = 0; it < 12; ++it) {
    std::vector<double> tp = theta;
    std::vector<double> tm = theta;
    tp[static_cast<std::size_t>(idx)] += h;
    tm[static_cast<std::size_t>(idx)] -= h;
    const double fp = fn(tp);
    const double fm = fn(tm);
    if (std::isfinite(fp) && std::isfinite(fm)) return h;
    h *= 0.5;
  }
  return NA_REAL;
}

arma::mat latent_numeric_hessian(const std::vector<double>& theta,
                                 const std::function<double(const std::vector<double>&)>& fn) {
  const int p = static_cast<int>(theta.size());
  arma::mat H(p, p, arma::fill::zeros);
  const double f0 = fn(theta);
  if (!std::isfinite(f0)) {
    H.fill(NA_REAL);
    return H;
  }

  std::vector<double> step(static_cast<std::size_t>(p), NA_REAL);
  for (int i = 0; i < p; ++i) {
    step[static_cast<std::size_t>(i)] = latent_safe_step(theta, i, fn);
    if (!std::isfinite(step[static_cast<std::size_t>(i)])) {
      H.fill(NA_REAL);
      return H;
    }
  }

  for (int i = 0; i < p; ++i) {
    const double hi = step[static_cast<std::size_t>(i)];
    std::vector<double> tp = theta;
    std::vector<double> tm = theta;
    tp[static_cast<std::size_t>(i)] += hi;
    tm[static_cast<std::size_t>(i)] -= hi;
    const double fp = fn(tp);
    const double fm = fn(tm);
    if (!std::isfinite(fp) || !std::isfinite(fm)) {
      H.fill(NA_REAL);
      return H;
    }
    H(i, i) = (fp - 2.0 * f0 + fm) / (hi * hi);

    for (int j = i + 1; j < p; ++j) {
      const double hj = step[static_cast<std::size_t>(j)];
      std::vector<double> tpp = theta;
      std::vector<double> tpm = theta;
      std::vector<double> tmp = theta;
      std::vector<double> tmm = theta;
      tpp[static_cast<std::size_t>(i)] += hi;
      tpp[static_cast<std::size_t>(j)] += hj;
      tpm[static_cast<std::size_t>(i)] += hi;
      tpm[static_cast<std::size_t>(j)] -= hj;
      tmp[static_cast<std::size_t>(i)] -= hi;
      tmp[static_cast<std::size_t>(j)] += hj;
      tmm[static_cast<std::size_t>(i)] -= hi;
      tmm[static_cast<std::size_t>(j)] -= hj;
      const double fpp = fn(tpp);
      const double fpm = fn(tpm);
      const double fmp = fn(tmp);
      const double fmm = fn(tmm);
      if (!std::isfinite(fpp) || !std::isfinite(fpm) ||
          !std::isfinite(fmp) || !std::isfinite(fmm)) {
        H.fill(NA_REAL);
        return H;
      }
      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * hi * hj);
      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return 0.5 * (H + H.t());
}

List latent_wald_result(double estimate,
                        double jacobian_rho,
                        const arma::mat& hessian,
                        int n_obs,
                        double conf_level,
                        const std::string& ci_method,
                        const std::string& inference_method,
                        double nll) {
  double se = NA_REAL;
  double statistic = NA_REAL;
  double p_value = NA_REAL;
  double lwr = NA_REAL;
  double upr = NA_REAL;

  if (hessian.n_rows >= 1 && hessian.n_cols >= 1 &&
      hessian.is_finite() &&
      std::isfinite(hessian(0, 0)) &&
      hessian(0, 0) > 0.0) {
    arma::mat vcov_theta;
    bool ok = arma::inv_sympd(vcov_theta, hessian);
    if (!ok) ok = arma::inv(vcov_theta, hessian);
    if (ok && vcov_theta.n_rows >= 1 && vcov_theta.n_cols >= 1 &&
        std::isfinite(vcov_theta(0, 0)) && vcov_theta(0, 0) >= 0.0) {
      const double var_rho = jacobian_rho * jacobian_rho * vcov_theta(0, 0);
      if (std::isfinite(var_rho) && var_rho >= 0.0) {
        se = std::sqrt(var_rho);
      }
    }
  }

  if (std::isfinite(estimate) && std::isfinite(se) && se > 0.0) {
    statistic = estimate / se;
    p_value = 2.0 * R::pnorm5(std::abs(statistic), 0.0, 1.0, false, false);
    const double zcrit = R::qnorm5(0.5 * (1.0 + conf_level), 0.0, 1.0, true, false);
    lwr = std::max(-1.0, estimate - zcrit * se);
    upr = std::min(1.0, estimate + zcrit * se);
  }

  return List::create(
    _["estimate"] = estimate,
    _["se"] = se,
    _["statistic"] = statistic,
    _["parameter"] = NA_REAL,
    _["p_value"] = p_value,
    _["lwr"] = lwr,
    _["upr"] = upr,
    _["n_obs"] = n_obs,
    _["conf.level"] = conf_level,
    _["ci.method"] = ci_method,
    _["inference.method"] = inference_method,
    _["nll"] = nll
  );
}

List tetrachoric_inference_impl(const NumericMatrix& tab,
                                double correct,
                                double conf_level) {
  if (tab.nrow() != 2 || tab.ncol() != 2) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              0, conf_level,
                              "wald_information_tetrachoric",
                              "wald_z_tetrachoric",
                              NA_REAL);
  }

  double a = tab(0, 0), b = tab(1, 0), c = tab(0, 1), d = tab(1, 1);
  const double n_obs = a + b + c + d;
  if (!(n_obs > 0.0)) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              0, conf_level,
                              "wald_information_tetrachoric",
                              "wald_z_tetrachoric",
                              NA_REAL);
  }

  const bool any_zero = (a <= 0.0) || (b <= 0.0) || (c <= 0.0) || (d <= 0.0);
  if (correct > 0.0 && any_zero) {
    if (a <= 0.0) a += correct;
    if (b <= 0.0) b += correct;
    if (c <= 0.0) c += correct;
    if (d <= 0.0) d += correct;
  }
  const double tot = a + b + c + d;
  if (!(tot > 0.0)) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              static_cast<int>(n_obs), conf_level,
                              "wald_information_tetrachoric",
                              "wald_z_tetrachoric",
                              NA_REAL);
  }

  const double pa = a / tot;
  const double pb = b / tot;
  const double pc = c / tot;
  const double pd = d / tot;
  const double rc0 = qnorm01(nan_preserve(pa + pb, 1e-12, 1.0 - 1e-12));
  const double cc0 = qnorm01(nan_preserve(pa + pc, 1e-12, 1.0 - 1e-12));
  auto nll_fixed = [&](double rho) { return tetra_nll_fixed_thresholds(rho, pa, pb, pc, pd, rc0, cc0); };
  const double est = nan_preserve(optimize(nll_fixed, -1.0, 1.0, 1e-8, 250), -1.0, 1.0);
  const double est_clamped = nan_preserve(est, -kLatentMaxcor, kLatentMaxcor);
  if (!std::isfinite(est_clamped)) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              static_cast<int>(n_obs), conf_level,
                              "wald_information_tetrachoric",
                              "wald_z_tetrachoric",
                              NA_REAL);
  }

  std::vector<double> theta0(1);
  theta0[0] = std::atanh(est_clamped / kLatentMaxcor);
  auto fn = [&](const std::vector<double>& theta) {
    return tetra_nll_fixed_thresholds(rho_from_theta(theta[0]), a, b, c, d, rc0, cc0);
  };
  const arma::mat H = latent_numeric_hessian(theta0, fn);
  return latent_wald_result(est_clamped, drho_dtheta(theta0[0]), H, static_cast<int>(n_obs), conf_level,
                            "wald_information_tetrachoric",
                            "wald_z_tetrachoric",
                            fn(theta0));
}

List polychoric_inference_impl(const NumericMatrix& tab,
                               double correct,
                               double conf_level) {
  if (tab.nrow() < 2 || tab.ncol() < 2) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              0, conf_level,
                              "wald_information_polychoric",
                              "wald_z_polychoric",
                              NA_REAL);
  }

  NumericMatrix N = clone(tab);
  double n_obs = 0.0;
  for (int i = 0; i < N.nrow(); ++i) {
    for (int j = 0; j < N.ncol(); ++j) n_obs += N(i, j);
  }
  if (!(n_obs > 0.0)) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              0, conf_level,
                              "wald_information_polychoric",
                              "wald_z_polychoric",
                              NA_REAL);
  }

  for (int i = 0; i < N.nrow(); ++i) {
    for (int j = 0; j < N.ncol(); ++j) N(i, j) /= n_obs;
  }

  N = drop_zero_margins(N);
  if (N.nrow() < 2 || N.ncol() < 2) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              static_cast<int>(n_obs), conf_level,
                              "wald_information_polychoric",
                              "wald_z_polychoric",
                              NA_REAL);
  }

  if (correct > 0.0) {
    for (int i = 0; i < N.nrow(); ++i) {
      for (int j = 0; j < N.ncol(); ++j) {
        if (N(i, j) <= 0.0) N(i, j) = correct / n_obs;
      }
    }
  }

  const int nr = N.nrow();
  const int nc = N.ncol();
  std::vector<double> probs(static_cast<std::size_t>(nr * nc), 0.0);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      probs[static_cast<std::size_t>(i * nc + j)] = N(i, j);
    }
  }

  std::vector<double> rc0, cc0;
  cutpoints_from_table(N, rc0, cc0);
  if (rc0.empty() || cc0.empty()) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              static_cast<int>(n_obs), conf_level,
                              "wald_information_polychoric",
                              "wald_z_polychoric",
                              NA_REAL);
  }

  PolyWorkspace ws0;
  ws0.reset(rc0, cc0);
  auto nll_fixed = [&](double rho) { return poly_nll_fixed_thresholds(rho, probs.data(), ws0); };
  const double est = nan_preserve(optimize(nll_fixed, -1.0, 1.0, 1e-8, 250), -1.0, 1.0);
  const double est_clamped = nan_preserve(est, -kLatentMaxcor, kLatentMaxcor);
  if (!std::isfinite(est_clamped)) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              static_cast<int>(n_obs), conf_level,
                              "wald_information_polychoric",
                              "wald_z_polychoric",
                              NA_REAL);
  }

  for (double& v : probs) v *= n_obs;
  std::vector<double> theta0(1);
  theta0[0] = std::atanh(est_clamped / kLatentMaxcor);
  PolyWorkspace ws;
  ws.reset(rc0, cc0);
  auto fn = [&](const std::vector<double>& theta) {
    return poly_nll_fixed_thresholds(rho_from_theta(theta[0]), probs.data(), ws);
  };
  const arma::mat H = latent_numeric_hessian(theta0, fn);
  return latent_wald_result(est_clamped, drho_dtheta(theta0[0]), H, static_cast<int>(n_obs), conf_level,
                            "wald_information_polychoric",
                            "wald_z_polychoric",
                            fn(theta0));
}

List polyserial_inference_impl(const NumericVector& x,
                               const IntegerVector& y,
                               double conf_level) {
  PolyserialPrepared prep;
  if (!polyserial_prepare_xy(x, y, false, prep)) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              0, conf_level,
                              "wald_information_polyserial",
                              "wald_z_polyserial",
                              NA_REAL);
  }

  PolyserialFitDetails fit;
  if (!polyserial_fit_full_details(prep.x, prep.y, fit, kLatentMaxcor)) {
    return latent_wald_result(NA_REAL, NA_REAL, arma::mat(1, 1, arma::fill::value(NA_REAL)),
                              static_cast<int>(prep.x.size()), conf_level,
                              "wald_information_polyserial",
                              "wald_z_polyserial",
                              NA_REAL);
  }

  auto fn = [&](const std::vector<double>& par) {
    return polyserial_negloglik_core(fit.z, fit.y, par.data(), fit.ncat, kLatentMaxcor);
  };
  const arma::mat H = latent_numeric_hessian(fit.par, fn);
  return latent_wald_result(fit.estimate, 1.0, H, fit.n_obs, conf_level,
                            "wald_information_polyserial",
                            "wald_z_polyserial",
                            fit.fmin);
}

} // namespace

// [[Rcpp::export]]
List matrixCorr_tetrachoric_inference_cpp(NumericMatrix tab,
                                          double correct = 0.5,
                                          double conf_level = 0.95) {
  return tetrachoric_inference_impl(tab, correct, conf_level);
}

// [[Rcpp::export]]
List matrixCorr_polychoric_inference_cpp(NumericMatrix tab,
                                         double correct = 0.5,
                                         double conf_level = 0.95) {
  return polychoric_inference_impl(tab, correct, conf_level);
}

// [[Rcpp::export]]
List matrixCorr_polyserial_inference_cpp(NumericVector x,
                                         IntegerVector y,
                                         double conf_level = 0.95) {
  return polyserial_inference_impl(x, y, conf_level);
}

// [[Rcpp::export]]
double matrixCorr_polydi_mle_cpp(NumericMatrix tab, double correct = 0.5) {
  const int R = tab.nrow(), C = tab.ncol();
  if (R < 2 || C != 2) return NA_REAL;
  return matrixCorr_polychoric_mle_cpp(tab, correct);
}
