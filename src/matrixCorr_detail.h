#pragma once
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <functional>

namespace matrixCorr_detail {

//------------------------------------------------------------------------------
// ---------- clamp variants ----------
//------------------------------------------------------------------------------
namespace clamp_policy {
// NaN-preserving (comparisons with NaN are false -> returns x i.e., NaN)
inline double nan_preserve(double x, double lo, double hi){
  return (x < lo) ? lo : (x > hi) ? hi : x;
}
// Saturate NaN to midpoint of [lo, hi] (alternative behavior)
inline double saturate_nan(double x, double lo, double hi){
  if (std::isnan(x)) return 0.5 * (lo + hi);
  return (x < lo) ? lo : (x > hi) ? hi : x;
}
}
//------------------------------------------------------------------------------
// ---------- standard normal helpers ----------
//------------------------------------------------------------------------------
namespace norm1 {
inline double Phi(double x)     { return R::pnorm(x, 0.0, 1.0, 1, 0); }
inline double phi(double x)     { return R::dnorm(x, 0.0, 1.0, 0); }
inline double qnorm01(double p) { return R::qnorm(p, 0.0, 1.0, 1, 0); }
}

//------------------------------------------------------------------------------
// ---------- Brent minimizer (scalar 1D) ----------
//------------------------------------------------------------------------------
namespace brent {
inline double optimize(std::function<double(double)> f,
                       double a, double b,
                       double tol = 1e-6, int max_iter = 100){
  const double golden = 0.3819660;
  double x = a + golden * (b - a), w = x, v = x;
  double fx = f(x), fw = fx, fv = fx;
  double d = 0.0, e = 0.0;

  for (int iter = 0; iter < max_iter; ++iter){
    double m = 0.5 * (a + b);
    double tol1 = tol * std::abs(x) + 1e-10;
    double tol2 = 2.0 * tol1;

    if (std::abs(x - m) <= tol2 - 0.5 * (b - a)) break;

    double p = 0, q = 0, r = 0;
    if (std::abs(e) > tol1){
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0) p = -p;
      q = std::abs(q);
      double temp = e;
      e = d;

      if (std::abs(p) >= std::abs(0.5 * q * temp) || p <= q * (a - x) || p >= q * (b - x)){
        e = (x < m) ? b - x : a - x;
        d = golden * e;
      } else {
        d = p / q;
        double u = x + d;
        if (u - a < tol2 || b - u < tol2)
          d = tol1 * ((m - x >= 0) ? 1 : -1);
      }
    } else {
      e = (x < m) ? b - x : a - x;
      d = golden * e;
    }

    double u = x + ((std::abs(d) >= tol1) ? d : tol1 * ((d > 0) ? 1 : -1));
    double fu = f(u);

    if (fu <= fx){
      if (u < x) b = x; else a = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x){
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w){
        v = u; fv = fu;
      }
    }
  }
  return x;
}
}
//------------------------------------------------------------------------------
// ---------- BVN CDF (adaptive 1D integration) ----------
//------------------------------------------------------------------------------
namespace bvn_adaptive {
using matrixCorr_detail::norm1::Phi;
using matrixCorr_detail::norm1::qnorm01;

struct BVNParams { double b, rho, sigma; };

inline double bvn_integrand(double t, const BVNParams& P){
  if (t <= 0.0) return 0.0;
  if (t >= 1.0) t = std::nextafter(1.0, 0.0);
  double z = qnorm01(t);
  double arg = (P.b - P.rho * z) / P.sigma;
  return Phi(arg);
}

inline double simpson_rec(std::function<double(double)> f, double a, double b,
                          double fa, double fb, double fm, double S, int depth){
  double m = 0.5 * (a + b);
  double lm = 0.5 * (a + m);
  double rm = 0.5 * (m + b);
  double flm = f(lm);
  double frm = f(rm);
  double Sleft  = (m - a) * (fa + 4.0 * flm + fm) / 6.0;
  double Sright = (b - m) * (fm + 4.0 * frm + fb) / 6.0;
  double S2 = Sleft + Sright;
  if (depth <= 0 || std::abs(S2 - S) < 1e-9 * (1.0 + std::abs(S2))) return S2;
  return simpson_rec(f, a, m, fa, fm, flm, Sleft,  depth - 1) +
    simpson_rec(f, m, b, fm, fb, frm, Sright, depth - 1);
}

inline double integrate_adaptive(std::function<double(double)> f, double a, double b){
  double fa = f(a), fb = f(b), fm = f(0.5 * (a + b));
  double S  = (b - a) * (fa + 4.0 * fm + fb) / 6.0;
  return simpson_rec(f, a, b, fa, fb, fm, S, 20);
}

inline double Phi2(double a, double b, double rho){
  using matrixCorr_detail::clamp_policy::nan_preserve; // we’ll clamp rho directly here
  if (std::isinf(a) && a < 0) return 0.0;
  if (std::isinf(b) && b < 0) return 0.0;
  if (std::isinf(a) && a > 0 && std::isinf(b) && b > 0) return 1.0;
  if (std::isinf(a) && a > 0) return Phi(b);
  if (std::isinf(b) && b > 0) return Phi(a);

  rho = nan_preserve(rho, -0.999999, 0.999999);
  double pa = Phi(a);
  if (pa <= 0.0) return 0.0;
  double sigma = std::sqrt(1.0 - rho * rho);

  if (rho > 0.99999)   return Phi(std::min(a, b));
  if (rho < -0.99999)  return std::max(0.0, Phi(a) + Phi(b) - 1.0);

  BVNParams P{b, rho, sigma};
  auto f = [&](double t){ return bvn_integrand(t, P); };
  // integrate t in [0, Φ(a)]
  double pa_clamped = std::min(1.0, std::max(0.0, pa));
  return integrate_adaptive(f, 0.0, pa_clamped);
}

inline double rect_prob(double al, double au, double bl, double bu, double rho){
  double A = Phi2(au, bu, rho);
  double B = Phi2(al, bu, rho);
  double C = Phi2(au, bl, rho);
  double D = Phi2(al, bl, rho);
  double p = A - B - C + D;
  return (p <= 0.0) ? 0.0 : p;
}
}
//------------------------------------------------------------------------------
// ---------- ranking utilities ----------
//------------------------------------------------------------------------------
namespace ranking {

// Average rank with tie handling
inline arma::vec rank_vector(const arma::vec& x) {
  const arma::uword n = x.n_elem;
  arma::uvec idx = arma::sort_index(x);
  arma::vec ranks(n);

  arma::uword i = 0;
  while (i < n) {
    arma::uword j = i + 1;
    const double xi = x(idx[i]);
    while (j < n && x(idx[j]) == xi) ++j;
    const double avg_rank = (static_cast<double>(i + j - 1) * 0.5) + 1.0;
    for (arma::uword k = i; k < j; ++k) ranks(idx[k]) = avg_rank;
    i = j;
  }
  return ranks;
}

// Invert standard deviations safely (avoid div by zero)
inline arma::vec safe_inv_stddev(const arma::vec& s) {
  arma::vec inv_s(s.n_elem, arma::fill::zeros);
  arma::uvec nz = arma::find(s > 0.0);
  inv_s.elem(nz) = 1.0 / s.elem(nz);
  return inv_s;
}
}
//------------------------------------------------------------------------------
// ---------- order-statistics + pairwise ----------
//------------------------------------------------------------------------------
namespace order_stats {

// stable integer discretisation (no NA handling here; caller decides policy)
inline std::vector<long long> discretise(const arma::vec& x, double scale){
  const arma::uword n = x.n_elem;
  std::vector<long long> out(n);
  for (arma::uword i = 0; i < n; ++i) out[i] = (long long) std::floor(x[i] * scale);
  return out;
}

inline bool has_ties(const std::vector<long long>& v){
  std::vector<long long> s = v;
  std::sort(s.begin(), s.end());
  return std::adjacent_find(s.begin(), s.end()) != s.end();
}

// ----- merge-sort inversion counter (for tau-a fast path) -----
inline long long merge_count(std::vector<int>& a, std::vector<int>& tmp, int L, int M, int R){
  int i = L, j = M, k = L;
  long long inv = 0;
  while (i < M && j <= R){
    if (a[i] <= a[j]) tmp[k++] = a[i++];
    else { tmp[k++] = a[j++]; inv += (M - i); }
  }
  while (i < M)  tmp[k++] = a[i++];
  while (j <= R) tmp[k++] = a[j++];
  for (int t = L; t <= R; ++t) a[t] = tmp[t];
  return inv;
}

inline long long sort_count(std::vector<int>& a, std::vector<int>& tmp, int L, int R){
  if (R - L < 1) return 0LL;
  int M = L + (R - L) / 2;
  long long inv = 0;
  inv += sort_count(a, tmp, L, M);
  inv += sort_count(a, tmp, M + 1, R);
  inv += merge_count(a, tmp, L, M + 1, R);
  return inv;
}

// ----- Fenwick (BIT) for tau-b numerator -----
struct Fenwick {
  std::vector<long long> t;
  int n;
  explicit Fenwick(int n): t(n + 1, 0), n(n) {}
  void add(int i, long long v){ for (; i <= n; i += i & -i) t[i] += v; }
  long long sum(int i) const { long long s = 0; for (; i > 0; i -= i & -i) s += t[i]; return s; }
};

// coordinate compression of arbitrary long long vector (returns 1-based ranks)
inline std::vector<int> compress_1based(const std::vector<long long>& v){
  std::vector<long long> u = v;
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  std::vector<int> r(v.size());
  for (size_t i = 0; i < v.size(); ++i){
    r[i] = int(std::lower_bound(u.begin(), u.end(), v[i]) - u.begin()) + 1;
  }
  return r;
}

} // namespace order_stats

namespace pairwise {

// Generic pairwise matrix builder for any column-wise metric:
//   metric(col_i, col_j) -> double
template <class Metric>
inline arma::mat apply(const arma::mat& X, Metric metric){
  const arma::uword p = X.n_cols;
  arma::mat M(p, p, arma::fill::ones);
  for (arma::uword j = 0; j < p; ++j){
    for (arma::uword k = j + 1; k < p; ++k){
      const double v = metric(X.col(j), X.col(k));
      M(j,k) = v; M(k,j) = v;
    }
  }
  return M;
}

} // namespace pairwise

} // namespace matrixCorr_detail
