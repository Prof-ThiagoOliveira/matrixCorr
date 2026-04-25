// kendall_corr.cpp
// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <string>
#include <vector>
#include "matrixCorr_detail.h"
#include "matrixCorr_omp.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using matrixCorr_detail::order_stats::getMs_ll;
using matrixCorr_detail::order_stats::insertion_sort_range;
using matrixCorr_detail::order_stats::inv_count_inplace;
using matrixCorr_detail::order_stats::tau_two_vectors_fast;

namespace {

struct KendallColumnOrder {
  std::vector<int> ord;
  std::vector<std::pair<int, int>> tie_runs;
  long long m1 = 0;
};

inline KendallColumnOrder build_kendall_column_order(const long long* x, const int n) {
  KendallColumnOrder state;
  state.ord.resize(n);
  std::iota(state.ord.begin(), state.ord.end(), 0);
  std::sort(state.ord.begin(), state.ord.end(),
            [&](int a, int b) { return x[a] < x[b]; });

  state.tie_runs.reserve(64);
  for (int s = 0; s < n; ) {
    int e = s + 1;
    const long long x0 = x[state.ord[s]];
    while (e < n && x[state.ord[e]] == x0) ++e;
    const int len = e - s;
    if (len > 1) {
      state.tie_runs.emplace_back(s, e);
      state.m1 += 1LL * len * (len - 1) / 2;
    }
    s = e;
  }
  return state;
}

inline double kendall_tau_from_order(const long long* y,
                                     const KendallColumnOrder& state) {
  const int n = static_cast<int>(state.ord.size());
  thread_local std::vector<long long> ybuf, mrg;
  if (static_cast<int>(ybuf.capacity()) < n) {
    ybuf.reserve(n);
    mrg.reserve(n);
  }
  ybuf.resize(n);
  mrg.resize(n);

  for (int k = 0; k < n; ++k) ybuf[k] = y[state.ord[k]];

  long long s_acc = 0;
  for (const auto& run : state.tie_runs) {
    const int s = run.first;
    const int e = run.second;
    const int len = e - s;
    if (len <= 32) insertion_sort_range(ybuf.data(), s, e);
    else           std::sort(ybuf.begin() + s, ybuf.begin() + e);
    s_acc += getMs_ll(ybuf.data() + s, len);
  }

  const long long inv = inv_count_inplace<long long>(ybuf.data(), mrg.data(), n);
  const long long m2  = getMs_ll(ybuf.data(), n);
  const long long n0  = 1LL * n * (n - 1) / 2LL;
  const double s = double(n0) - double(state.m1) - double(m2)
    - 2.0 * double(inv) + double(s_acc);
  const double den1 = double(n0 - state.m1);
  const double den2 = double(n0 - m2);
  if (den1 <= 0.0 || den2 <= 0.0) return NA_REAL;
  return s / std::sqrt(den1 * den2);
}

inline double clamp_tau(double x) {
  if (!std::isfinite(x)) return NA_REAL;
  if (x > 1.0) return 1.0;
  if (x < -1.0) return -1.0;
  return x;
}

inline double bounded_atanh(double x) {
  constexpr double eps = 1e-15;
  const double xc = std::min(1.0 - eps, std::max(-1.0 + eps, x));
  return std::atanh(xc);
}

template <class Func>
double bisect_root(Func&& g, double lo, double hi,
                   double f_lo, double f_hi,
                   const double tol = 1e-10,
                   const int max_iter = 200) {
  if (std::abs(f_lo) <= tol) return lo;
  if (std::abs(f_hi) <= tol) return hi;

  double left = lo;
  double right = hi;
  double f_left = f_lo;

  for (int iter = 0; iter < max_iter; ++iter) {
    const double mid = 0.5 * (left + right);
    const double f_mid = g(mid);
    if (!std::isfinite(f_mid) || std::abs(f_mid) <= tol ||
        0.5 * (right - left) <= tol) {
      return mid;
    }
    if ((f_left <= 0.0 && f_mid <= 0.0) || (f_left >= 0.0 && f_mid >= 0.0)) {
      left = mid;
      f_left = f_mid;
    } else {
      right = mid;
    }
  }

  return 0.5 * (left + right);
}

inline bool kendall_fieller_ci_core(const double tau,
                                    const int n,
                                    const double conf_level,
                                    double& lwr,
                                    double& upr) {
  lwr = NA_REAL;
  upr = NA_REAL;
  if (!std::isfinite(tau) || n <= 4) return false;

  const double z = bounded_atanh(tau);
  const double se = std::sqrt(0.437 / (static_cast<double>(n) - 4.0));
  const double alpha = 1.0 - conf_level;
  const double crit = R::qnorm(1.0 - 0.5 * alpha, 0.0, 1.0, 1, 0);

  lwr = clamp_tau(std::tanh(z - crit * se));
  upr = clamp_tau(std::tanh(z + crit * se));
  return std::isfinite(lwr) && std::isfinite(upr);
}

struct KendallBBCell {
  int xr = 0;
  int yr = 0;
  int count = 0;
  double pi_c = 0.0;
  double pi_d = 0.0;
};

class FenwickCounts {
 public:
  explicit FenwickCounts(const int n = 0) : tree_(static_cast<std::size_t>(n + 1), 0.0) {}

  void reset() {
    std::fill(tree_.begin(), tree_.end(), 0.0);
  }

  void add(int idx, double value) {
    ++idx;
    while (idx < static_cast<int>(tree_.size())) {
      tree_[static_cast<std::size_t>(idx)] += value;
      idx += idx & -idx;
    }
  }

  double sum_prefix(int idx) const {
    if (idx < 0) return 0.0;
    double out = 0.0;
    ++idx;
    while (idx > 0) {
      out += tree_[static_cast<std::size_t>(idx)];
      idx -= idx & -idx;
    }
    return out;
  }

 private:
  std::vector<double> tree_;
};

inline bool kendall_brown_benedetti_ci_core(const arma::vec& x,
                                            const arma::vec& y,
                                            const double conf_level,
                                            double& tau,
                                            double& lwr,
                                            double& upr) {
  lwr = NA_REAL;
  upr = NA_REAL;
  tau = NA_REAL;

  const arma::uword n = x.n_elem;
  if (n != y.n_elem || n < 2u) return false;

  std::vector<double> ux(x.begin(), x.end());
  std::vector<double> uy(y.begin(), y.end());
  std::sort(ux.begin(), ux.end());
  ux.erase(std::unique(ux.begin(), ux.end()), ux.end());
  std::sort(uy.begin(), uy.end());
  uy.erase(std::unique(uy.begin(), uy.end()), uy.end());

  const int nr = static_cast<int>(ux.size());
  const int nc = static_cast<int>(uy.size());
  if (nr < 2 || nc < 2) return false;

  std::vector<std::pair<int, int>> coords;
  coords.reserve(n);
  for (arma::uword i = 0; i < n; ++i) {
    const int xr = static_cast<int>(std::lower_bound(ux.begin(), ux.end(), x[i]) - ux.begin());
    const int yr = static_cast<int>(std::lower_bound(uy.begin(), uy.end(), y[i]) - uy.begin());
    coords.emplace_back(xr, yr);
  }
  std::sort(coords.begin(), coords.end());

  std::vector<KendallBBCell> cells;
  cells.reserve(coords.size());
  std::vector<double> row_counts(static_cast<std::size_t>(nr), 0.0);
  std::vector<double> col_counts(static_cast<std::size_t>(nc), 0.0);

  for (std::size_t i = 0; i < coords.size(); ) {
    std::size_t j = i + 1u;
    while (j < coords.size() && coords[j] == coords[i]) ++j;

    KendallBBCell cell;
    cell.xr = coords[i].first;
    cell.yr = coords[i].second;
    cell.count = static_cast<int>(j - i);
    cells.push_back(cell);
    row_counts[static_cast<std::size_t>(cell.xr)] += static_cast<double>(cell.count);
    col_counts[static_cast<std::size_t>(cell.yr)] += static_cast<double>(cell.count);
    i = j;
  }

  FenwickCounts bit(nc);
  double total_before = 0.0;
  for (std::size_t s = 0; s < cells.size(); ) {
    std::size_t e = s + 1u;
    while (e < cells.size() && cells[e].xr == cells[s].xr) ++e;

    for (std::size_t idx = s; idx < e; ++idx) {
      const int yr = cells[idx].yr;
      const double sw = bit.sum_prefix(yr - 1);
      const double nw = total_before - bit.sum_prefix(yr);
      cells[idx].pi_c = sw;
      cells[idx].pi_d = nw;
    }

    for (std::size_t idx = s; idx < e; ++idx) {
      bit.add(cells[idx].yr, static_cast<double>(cells[idx].count));
      total_before += static_cast<double>(cells[idx].count);
    }
    s = e;
  }

  bit.reset();
  double total_after = 0.0;
  for (std::size_t e = cells.size(); e > 0u; ) {
    std::size_t s = e - 1u;
    while (s > 0u && cells[s - 1u].xr == cells[e - 1u].xr) --s;

    for (std::size_t idx = s; idx < e; ++idx) {
      const int yr = cells[idx].yr;
      const double se = bit.sum_prefix(yr - 1);
      const double ne = total_after - bit.sum_prefix(yr);
      cells[idx].pi_c += ne;
      cells[idx].pi_d += se;
    }

    for (std::size_t idx = s; idx < e; ++idx) {
      bit.add(cells[idx].yr, static_cast<double>(cells[idx].count));
      total_after += static_cast<double>(cells[idx].count);
    }
    e = s;
  }

  double C = 0.0;
  double D = 0.0;
  for (const auto& cell : cells) {
    C += cell.pi_c * static_cast<double>(cell.count);
    D += cell.pi_d * static_cast<double>(cell.count);
  }
  C *= 0.5;
  D *= 0.5;

  const double nd = static_cast<double>(n);
  const double n0 = nd * (nd - 1.0) / 2.0;
  double n1 = 0.0;
  double n2 = 0.0;
  double rowsum_sq = 0.0;
  double colsum_sq = 0.0;
  for (double rc : row_counts) {
    n1 += rc * (rc - 1.0) / 2.0;
    const double p = rc / nd;
    rowsum_sq += p * p;
  }
  for (double cc : col_counts) {
    n2 += cc * (cc - 1.0) / 2.0;
    const double p = cc / nd;
    colsum_sq += p * p;
  }

  const double den_tau = std::sqrt((n0 - n1) * (n0 - n2));
  if (!(den_tau > 0.0)) return false;
  tau = clamp_tau((C - D) / den_tau);
  if (!std::isfinite(tau)) return false;

  const double delta1 = std::sqrt(std::max(0.0, 1.0 - rowsum_sq));
  const double delta2 = std::sqrt(std::max(0.0, 1.0 - colsum_sq));
  if (!(delta1 > 0.0) || !(delta2 > 0.0)) return false;

  const double Pdiff = 2.0 * (C - D) / (nd * nd);
  double sum_pi_tauphi = 0.0;
  double sum_pi_tauphi2 = 0.0;
  for (const auto& cell : cells) {
    const double pi = static_cast<double>(cell.count) / nd;
    const double rowprob = row_counts[static_cast<std::size_t>(cell.xr)] / nd;
    const double colprob = col_counts[static_cast<std::size_t>(cell.yr)] / nd;
    const double pdiff = (cell.pi_c - cell.pi_d) / nd;
    const double tauphi =
      (2.0 * pdiff + Pdiff * colprob) * delta2 * delta1 +
      (Pdiff * rowprob * delta2) / delta1;
    sum_pi_tauphi += pi * tauphi;
    sum_pi_tauphi2 += pi * tauphi * tauphi;
  }

  double sigma2 =
    (sum_pi_tauphi2 - sum_pi_tauphi * sum_pi_tauphi) /
    std::pow(delta1 * delta2, 4.0) / nd;
  if (sigma2 < std::numeric_limits<double>::epsilon() * 10.0) sigma2 = 0.0;
  if (!(sigma2 >= 0.0) || !std::isfinite(sigma2)) return false;

  const double pr2 = 1.0 - (1.0 - conf_level) / 2.0;
  const double crit = R::qnorm(pr2, 0.0, 1.0, 1, 0);
  const double half_width = crit * std::sqrt(sigma2);
  lwr = clamp_tau(std::max(tau - half_width, -1.0));
  upr = clamp_tau(std::min(tau + half_width, 1.0));
  return std::isfinite(lwr) && std::isfinite(upr);
}

inline bool kendall_ifel_pseudo_obs(const arma::vec& x,
                                    const arma::vec& y,
                                    const double tau_input,
                                    arma::vec& w,
                                    double& tau_hat) {
  const arma::uword n = x.n_elem;
  if (n != y.n_elem || n < 3u || !std::isfinite(tau_input)) return false;

  arma::vec s_sum(n, arma::fill::zeros);
  arma::vec a_sum(n, arma::fill::zeros);
  arma::vec b_sum(n, arma::fill::zeros);

  double s_total = 0.0;
  double a_total = 0.0;
  double b_total = 0.0;

  for (arma::uword i = 0; i + 1u < n; ++i) {
    for (arma::uword j = i + 1u; j < n; ++j) {
      const double dx = x[i] - x[j];
      const double dy = y[i] - y[j];
      const int sx = (dx > 0.0) - (dx < 0.0);
      const int sy = (dy > 0.0) - (dy < 0.0);
      const double s_ij = static_cast<double>(sx * sy);
      const double a_ij = (sx != 0) ? 1.0 : 0.0;
      const double b_ij = (sy != 0) ? 1.0 : 0.0;

      s_sum[i] += s_ij;
      s_sum[j] += s_ij;
      a_sum[i] += a_ij;
      a_sum[j] += a_ij;
      b_sum[i] += b_ij;
      b_sum[j] += b_ij;

      s_total += s_ij;
      a_total += a_ij;
      b_total += b_ij;
    }
  }

  if (!(a_total > 0.0) || !(b_total > 0.0)) return false;

  tau_hat = s_total / std::sqrt(a_total * b_total);
  tau_hat = clamp_tau(tau_hat);
  if (!std::isfinite(tau_hat)) return false;

  const double n0 = 0.5 * static_cast<double>(n) * static_cast<double>(n - 1u);
  const double s_bar = s_total / n0;
  const double a_bar = a_total / n0;
  const double b_bar = b_total / n0;
  if (!(a_bar > 0.0) || !(b_bar > 0.0)) return false;

  const double den = std::sqrt(a_bar * b_bar);
  const double inv_nm1 = 1.0 / static_cast<double>(n - 1u);

  w.set_size(n);
  for (arma::uword i = 0; i < n; ++i) {
    const double s_i = s_sum[i] * inv_nm1;
    const double a_i = a_sum[i] * inv_nm1;
    const double b_i = b_sum[i] * inv_nm1;
    const double psi_i = 2.0 * (
      (s_i - s_bar) / den -
        0.5 * tau_hat * ((a_i - a_bar) / a_bar + (b_i - b_bar) / b_bar)
    );
    w[i] = tau_hat + psi_i;
  }

  return true;
}

inline bool kendall_el_ratio_stat(const arma::vec& w,
                                  const double theta,
                                  double& stat) {
  constexpr double tol = 1e-12;
  constexpr double domain_eps = 1e-12;

  stat = std::numeric_limits<double>::infinity();
  if (w.n_elem == 0u || !std::isfinite(theta)) return false;

  const arma::vec v = w - theta;
  const double vmax = v.max();
  const double vmin = v.min();

  if (std::max(std::abs(vmax), std::abs(vmin)) <= tol) {
    stat = 0.0;
    return true;
  }

  if (!(vmin < 0.0 && vmax > 0.0)) return false;

  double lo = -std::numeric_limits<double>::infinity();
  double hi = std::numeric_limits<double>::infinity();
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    const double vi = v[i];
    if (vi > 0.0) {
      lo = std::max(lo, -1.0 / vi);
    } else if (vi < 0.0) {
      hi = std::min(hi, -1.0 / vi);
    }
  }
  lo += domain_eps;
  hi -= domain_eps;
  if (!(lo < hi) || !std::isfinite(lo) || !std::isfinite(hi)) return false;

  const auto g = [&](double lambda) -> double {
    double acc = 0.0;
    for (arma::uword i = 0; i < v.n_elem; ++i) {
      const double den = 1.0 + lambda * v[i];
      if (!(den > 0.0) || !std::isfinite(den)) {
        return (lambda > 0.0)
          ? std::numeric_limits<double>::infinity()
          : -std::numeric_limits<double>::infinity();
      }
      acc += v[i] / den;
    }
    return acc;
  };

  const double g0 = arma::accu(v);
  double lambda = 0.0;
  if (std::abs(g0) > tol) {
    const double g_lo = g(lo);
    const double g_hi = g(hi);
    if (!std::isfinite(g_lo) || !std::isfinite(g_hi) || g_lo <= 0.0 || g_hi >= 0.0) {
      return false;
    }
    lambda = bisect_root(g, lo, hi, g_lo, g_hi);
  }

  double llr = 0.0;
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    const double den = 1.0 + lambda * v[i];
    if (!(den > 0.0) || !std::isfinite(den)) return false;
    llr += 2.0 * std::log(den);
  }
  stat = llr;
  return std::isfinite(stat);
}

inline bool kendall_ifel_ci_core(const arma::vec& x,
                                 const arma::vec& y,
                                 const double tau_input,
                                 const double conf_level,
                                 double& lwr,
                                 double& upr) {
  constexpr double theta_eps = 1e-10;

  lwr = NA_REAL;
  upr = NA_REAL;

  arma::vec w;
  double tau_hat = tau_input;
  if (!kendall_ifel_pseudo_obs(x, y, tau_input, w, tau_hat)) return false;

  const double spread = arma::mean(arma::square(w - tau_hat));
  if (!std::isfinite(spread)) return false;
  if (spread <= 1e-14) {
    lwr = tau_hat;
    upr = tau_hat;
    return true;
  }

  const double crit = R::qchisq(conf_level, 1.0, 1, 0);
  const auto f = [&](double theta) -> double {
    double stat = 0.0;
    if (!kendall_el_ratio_stat(w, theta, stat)) {
      return std::numeric_limits<double>::infinity();
    }
    return stat - crit;
  };

  const double w_min = w.min();
  const double w_max = w.max();
  const double left_edge = std::max(-1.0, w_min);
  const double right_edge = std::min(1.0, w_max);

  if (tau_hat <= left_edge + theta_eps) {
    lwr = left_edge;
  } else {
    const double left_start = std::min(tau_hat, left_edge + theta_eps);
    const double f_left = f(left_start);
    const double f_tau = f(tau_hat);
    if (!std::isfinite(f_left) || f_left > 0.0) {
      lwr = bisect_root(f, left_start, tau_hat, f_left, f_tau);
    } else {
      lwr = left_edge;
    }
  }

  if (tau_hat >= right_edge - theta_eps) {
    upr = right_edge;
  } else {
    const double right_start = std::max(tau_hat, right_edge - theta_eps);
    const double f_tau = f(tau_hat);
    const double f_right = f(right_start);
    if (!std::isfinite(f_right) || f_right > 0.0) {
      upr = bisect_root(f, tau_hat, right_start, f_tau, f_right);
    } else {
      upr = right_edge;
    }
  }

  lwr = clamp_tau(lwr);
  upr = clamp_tau(upr);
  return std::isfinite(lwr) && std::isfinite(upr);
}

} // namespace

// ---------------------------- Matrix version ---------------------------------
// [[Rcpp::export]]
Rcpp::NumericMatrix kendall_matrix_cpp(Rcpp::NumericMatrix mat){
  const int n = mat.nrow();
  const int p = mat.ncol();

  // --- Scalar fast path (p == 2): raw doubles, no discretisation
  if (p == 2) {
    const double* x = &mat(0, 0);
    const double* y = &mat(0, 1);
    const double tau = tau_two_vectors_fast(x, y, n);
    Rcpp::NumericMatrix out(2, 2);
    out(0,0) = out(1,1) = 1.0;
    out(0,1) = out(1,0) = tau;
    return out;
  }

  // --- Discretise once per column for p >= 3 (matrix path)
  std::vector<long long> cols(static_cast<std::size_t>(p) * static_cast<std::size_t>(n));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if(p > 1 && omp_get_max_threads() > 1)
#endif
  for (int j = 0; j < p; ++j) {
    const double* cj = &mat(0, j);
    long long* col_j = cols.data() + static_cast<std::size_t>(j) * static_cast<std::size_t>(n);
    double max_abs = 0.0;
    for (int i = 0; i < n; ++i)
      max_abs = std::max(max_abs, std::abs(cj[i]));

    double scale = 1e8;
    if (max_abs > 0.0) {
      const double max_ll = static_cast<double>(std::numeric_limits<long long>::max()) / 4.0;
      double safe_scale = std::floor(max_ll / max_abs);
      if (!std::isfinite(safe_scale) || safe_scale < 1.0) safe_scale = 1.0;
      scale = std::min(1e8, safe_scale);
    }

    for (int i = 0; i < n; ++i)
      col_j[i] = static_cast<long long>(std::floor(cj[i] * scale));
  }

  Rcpp::NumericMatrix out(p, p);
  for (int j = 0; j < p; ++j) out(j,j) = 1.0;

  auto fill_row = [&](int i) {
    const long long* col_i = cols.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(n);
    const KendallColumnOrder state = build_kendall_column_order(col_i, n);
    for (int j = i + 1; j < p; ++j) {
      const long long* col_j = cols.data() + static_cast<std::size_t>(j) * static_cast<std::size_t>(n);
      const double tau = kendall_tau_from_order(col_j, state);
      out(i, j) = out(j, i) = tau;
    }
  };

  // --- Avoid OpenMP overhead when p is tiny
#if defined(_OPENMP)
  if (p >= 3 && omp_get_max_threads() > 1) {
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < p - 1; ++i) {
      fill_row(i);
    }
  } else
#endif
{
  // single-threaded path
  for (int i = 0; i < p - 1; ++i) {
    fill_row(i);
  }
}

  return out;
}

// ---------------------------- Two-vector wrapper ------------------------------
// [[Rcpp::export]]
double kendall_tau2_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  const R_xlen_t n = x.size();
  if (n != y.size() || n < 2) return NA_REAL;
  const double* px = x.begin();
  const double* py = y.begin();
  return matrixCorr_detail::order_stats::tau_two_vectors_fast(px, py, static_cast<int>(n));
}


// [[Rcpp::export]]
double kendall_tau2_from_mat_cpp(Rcpp::NumericMatrix mat) {
  if (mat.ncol() != 2 || mat.nrow() < 2) return NA_REAL;
  const int n = mat.nrow();
  const double* x = &mat(0, 0);
  const double* y = &mat(0, 1);
  return matrixCorr_detail::order_stats::tau_two_vectors_fast(x, y, n);
}

// [[Rcpp::export]]
Rcpp::List kendall_matrix_pairwise_cpp(SEXP X_,
                                       const bool return_ci = false,
                                       const double conf_level = 0.95,
                                       const std::string ci_method = "fieller") {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");
  if (return_ci && !(conf_level > 0.0 && conf_level < 1.0))
    Rcpp::stop("conf_level must be in (0,1).");
  if (return_ci &&
      ci_method != "fieller" &&
      ci_method != "if_el" &&
      ci_method != "brown_benedetti")
    Rcpp::stop("ci_method must be one of 'fieller', 'if_el', or 'brown_benedetti'.");

  const arma::uword n = static_cast<arma::uword>(Rf_nrows(X_));
  const arma::uword p = static_cast<arma::uword>(Rf_ncols(X_));
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  arma::mat est(p, p, arma::fill::none);
  est.fill(arma::datum::nan);
  arma::Mat<int> n_complete(p, p, arma::fill::zeros);

  arma::mat lwr;
  arma::mat upr;
  if (return_ci) {
    lwr.set_size(p, p);
    upr.set_size(p, p);
    lwr.fill(arma::datum::nan);
    upr.fill(arma::datum::nan);
  }

  std::vector<arma::uvec> finite_idx(p);
  for (arma::uword j = 0; j < p; ++j) {
    finite_idx[j] = arma::find_finite(X.col(j));
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(p >= 4u && omp_get_max_threads() > 1)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const arma::uvec& idx_j = finite_idx[uj];
    const arma::uword n_idx_j = idx_j.n_elem;
    const arma::uword* idx_j_ptr = idx_j.memptr();
    const double* colj_ptr = X.colptr(uj);

    static thread_local arma::vec xbuf;
    static thread_local arma::vec ybuf;

    for (arma::uword k = uj + 1u; k < p; ++k) {
      const arma::uvec& idx_k = finite_idx[k];
      const arma::uword n_idx_k = idx_k.n_elem;
      const std::size_t possible = std::min(n_idx_j, n_idx_k);
      if (possible < 2u) continue;

      if (xbuf.n_elem < possible) {
        xbuf.set_size(static_cast<arma::uword>(possible));
        ybuf.set_size(static_cast<arma::uword>(possible));
      }
      const arma::uword* idx_k_ptr = idx_k.memptr();
      const double* colk_ptr = X.colptr(k);
      arma::uword ia = 0u, ib = 0u;
      arma::uword overlap_n = 0u;
      while (ia < n_idx_j && ib < n_idx_k) {
        const arma::uword a = idx_j_ptr[ia];
        const arma::uword b = idx_k_ptr[ib];
        if (a == b) {
          xbuf[overlap_n] = colj_ptr[a];
          ybuf[overlap_n] = colk_ptr[b];
          ++overlap_n;
          ++ia;
          ++ib;
        } else if (a < b) {
          ++ia;
        } else {
          ++ib;
        }
      }

      const int overlap_n_i = static_cast<int>(overlap_n);
      n_complete(uj, k) = overlap_n_i;
      n_complete(k, uj) = overlap_n_i;
      if (overlap_n < 2u) continue;

      arma::vec x_view(xbuf.memptr(), overlap_n, false, true);
      arma::vec y_view(ybuf.memptr(), overlap_n, false, true);
      const double tau = tau_two_vectors_fast(x_view.memptr(), y_view.memptr(),
                                              static_cast<int>(overlap_n));
      est(uj, k) = tau;
      est(k, uj) = tau;

      if (return_ci && std::isfinite(tau)) {
        double lo = arma::datum::nan;
        double hi = arma::datum::nan;
        bool ok = false;
        if (ci_method == "fieller") {
          ok = kendall_fieller_ci_core(tau, overlap_n_i, conf_level, lo, hi);
        } else if (ci_method == "brown_benedetti") {
          double tau_bb = arma::datum::nan;
          ok = kendall_brown_benedetti_ci_core(x_view, y_view, conf_level, tau_bb, lo, hi);
          if (ok && std::isfinite(tau_bb)) {
            est(uj, k) = tau_bb;
            est(k, uj) = tau_bb;
          }
        } else {
          ok = kendall_ifel_ci_core(x_view, y_view, tau, conf_level, lo, hi);
        }
        if (ok) {
          lwr(uj, k) = lo;
          lwr(k, uj) = lo;
          upr(uj, k) = hi;
          upr(k, uj) = hi;
        }
      }
    }
  }

  for (arma::uword j = 0; j < p; ++j) {
    n_complete(j, j) = static_cast<int>(finite_idx[j].n_elem);
    est(j, j) = (finite_idx[j].n_elem >= 2u) ? 1.0 : arma::datum::nan;
    if (return_ci) {
      lwr(j, j) = arma::datum::nan;
      upr(j, j) = arma::datum::nan;
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["est"] = est,
    Rcpp::_["n_complete"] = n_complete
  );
  if (return_ci) {
    out["lwr"] = lwr;
    out["upr"] = upr;
    out["conf_level"] = conf_level;
    out["ci_method"] = ci_method;
  }
  return out;
}
