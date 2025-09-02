#pragma once
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <limits>

namespace matrixCorr_detail {

//------------------------------------------------------------------------------
// ---------- clamp variants ----------
//------------------------------------------------------------------------------
namespace clamp_policy {
// NaN-preserving (comparisons with NaN are false -> returns x i.e., NaN)
inline double nan_preserve(double x, double lo, double hi) noexcept {
  return (x < lo) ? lo : (x > hi) ? hi : x;
}
// Saturate NaN to midpoint of [lo, hi] (alternative behavior)
inline double saturate_nan(double x, double lo, double hi) noexcept {
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

template<class F>
inline double optimize(F&& f, double a, double b,
                       double tol = 1e-8, int max_iter = 100)
  noexcept(noexcept(std::declval<F&>()(std::declval<double>())))
  {
    auto&& func = f; // non-const reference to the callable
    using std::abs;
    constexpr double golden = 0.3819660112501051; // more precise
    const double eps = std::numeric_limits<double>::epsilon();

    double x = a + golden * (b - a), w = x, v = x;
    double fx = func(x), fw = fx, fv = fx;
    double d = 0.0, e = 0.0;

    for (int iter = 0; iter < max_iter; ++iter){
      double m = 0.5 * (a + b);
      double tol1 = tol * abs(x) + eps;
      double tol2 = 2.0 * tol1;

      if (abs(x - m) <= tol2 - 0.5 * (b - a)) break;

      double p = 0, q = 0, r = 0;
      if (abs(e) > tol1){
        r = (x - w) * (fx - fv);
        q = (x - v) * (fx - fw);
        p = (x - v) * q - (x - w) * r;
        q = 2.0 * (q - r);
        if (q > 0) p = -p;
        q = abs(q);
        const double etemp = e;
        e = d;

        if (abs(p) >= abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)){
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

      double u = x + ((abs(d) >= tol1) ? d : tol1 * ((d > 0) ? 1 : -1));
      double fu = func(u);

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
} // namespace brent


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

template<class F>
inline double simpson_rec(const F& f, double a, double b,
                          double fa, double fb, double fm, double S, int depth)
  noexcept(noexcept(f(a)))
  {
    double m  = 0.5 * (a + b);
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

template<class F>
inline double integrate_adaptive(F&& f, double a, double b)
  noexcept(noexcept(std::declval<F&>()(a)))
  {
    const F& func = f;                           // bind once
    double fa = func(a), fb = func(b), fm = func(0.5 * (a + b));
    double S  = (b - a) * (fa + 4.0 * fm + fb) / 6.0;
    return simpson_rec(func, a, b, fa, fb, fm, S, 20);
  }

inline double Phi2(double a, double b, double rho){
  using matrixCorr_detail::clamp_policy::nan_preserve; // clamp rho here
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
} // namespace bvn_adaptive

//------------------------------------------------------------------------------
// ---------- ranking utilities ----------
//------------------------------------------------------------------------------
namespace ranking {

// Close-enough comparator: treat a and b as a tie if |a-b| <= abs_eps + rel_eps*max(|a|,|b|)
inline bool is_close(double a, double b, double abs_eps, double rel_eps) noexcept {
  const double diff = std::abs(a - b);
  return diff <= (abs_eps + rel_eps * std::max(std::abs(a), std::abs(b)));
}

// tolerance-capable ranker (average ranks)
inline arma::vec rank_vector_eps(const arma::vec& x,
                                 double abs_eps = 0.0,
                                 double rel_eps = 0.0,
                                 bool   stable  = true)
{
  const arma::uword n = x.n_elem;
  // stable sort gives deterministic grouping for near-ties
  arma::uvec idx = stable ? arma::stable_sort_index(x) : arma::sort_index(x);
  // precompute sorted values to avoid x(idx[...]) in the loop
  arma::vec xs   = x.elem(idx);
  arma::vec ranks(n);

  arma::uword i = 0;
  while (i < n) {
    arma::uword j = i + 1;
    const double xi = xs[i];
    while (j < n && is_close(xs[j], xi, abs_eps, rel_eps)) ++j;
    // 1-based average
    const double avg_rank = 0.5 * static_cast<double>(i + j - 1) + 1.0;
    for (arma::uword k = i; k < j; ++k) ranks(idx[k]) = avg_rank;
    i = j;
  }
  return ranks;
}

// Backwards-compatible wrapper: exact ties only (your original behavior)
inline arma::vec rank_vector(const arma::vec& x) {
  return rank_vector_eps(x, /*abs_eps=*/0.0, /*rel_eps=*/0.0, /*stable=*/true);
}

// Invert standard deviations safely (unchanged)
inline arma::vec safe_inv_stddev(const arma::vec& s) {
  arma::vec inv_s(s.n_elem, arma::fill::zeros);
  arma::uvec nz = arma::find(s > 0.0);
  inv_s.elem(nz) = 1.0 / s.elem(nz);
  return inv_s;
}

} // namespace ranking
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

//------------------------------------------------------------------------------
// ---------- pairwise ----------
//------------------------------------------------------------------------------

namespace pairwise {

// Generic pairwise matrix builder for any column-wise metric:
//   metric(col_i, col_j) -> double
template <class Metric>
inline arma::mat apply_ptr(const arma::mat& X, Metric&& metric){
  const arma::uword p = X.n_cols, n = X.n_rows;
  arma::mat M(p, p, arma::fill::ones);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j){
    const double* cj = X.colptr(static_cast<arma::uword>(j));
    for (arma::uword k = static_cast<arma::uword>(j) + 1; k < p; ++k){
      const double* ck = X.colptr(k);
      const double v = metric(cj, ck, n);
      M(j,k) = v; M(k,j) = v;
    }
  }
  return M;
}

// double(const arma::vec&, const arma::vec&)
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

//------------------------------------------------------------------------------
// ---------- linear algebra helpers ----------
//------------------------------------------------------------------------------
namespace linalg {

// X'X (upper triangle via BLAS SYRK) then symmetrize.
// Returns p x p matrix XtX = X'X.
inline arma::mat crossprod_no_copy(const arma::mat& X) {
  const arma::uword n = X.n_rows, p = X.n_cols;
  arma::mat XtX(p, p);
#if defined(ARMA_USE_BLAS)
{
  XtX.zeros();
  const arma::blas_int N = static_cast<arma::blas_int>(p);
  const arma::blas_int K = static_cast<arma::blas_int>(n);
  const double alpha = 1.0, beta = 0.0;
  const char uplo  = 'U';
  const char trans = 'T';
  arma::blas::syrk<double>(&uplo, &trans, &N, &K,
                           &alpha, X.memptr(), &K,
                           &beta,  XtX.memptr(), &N);
                           XtX = arma::symmatu(XtX);
}
#else
XtX = X.t() * X;
#endif
return XtX;
}

// M := M - n * mu * mu' (operate on upper triangle) then symmetrize.
// mu is 1 x p.
inline void subtract_n_outer_mu(arma::mat& M, const arma::rowvec& mu, double n) {
  const arma::uword p = M.n_cols;
  const double scale = -n;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double muj = mu[uj];
    const double fj  = scale * muj; // = -n * mu_j
    for (arma::uword i = 0; i <= uj; ++i) {
      M(i, uj) += fj * mu[i];       // -= n * mu_i * mu_j
    }
  }
  M = arma::symmatu(M);
}

// Ensure positive definiteness by geometric diagonal jitter until chol succeeds.
inline void make_pd_inplace(arma::mat& S, double& jitter, const double max_jitter = 1e-2) {
  if (jitter < 0) jitter = 0.0;
  arma::mat C;
  for (;;) {
    if (arma::chol(C, S, "upper")) return;
    if (jitter == 0.0) jitter = 1e-8; else jitter *= 10.0;
    if (jitter > max_jitter)
      Rcpp::stop("Covariance not positive definite; jitter exceeded limit.");
    S.diag() += jitter;
  }
}

} // namespace linalg

//------------------------------------------------------------------------------
// ---------- covariance shrinkage (OAS) ----------
//------------------------------------------------------------------------------
namespace cov_shrinkage {

// Chen–Wiesel–Hero (2010) OAS shrinkage to mu*I.
// Input: cov_mle = (1/n) * (X - mu)'(X - mu), with n samples.
// Returns Sigma = (1 - rho)*S + rho*mu*I, rho in [0,1]; writes rho to rho_out.
inline arma::mat oas_shrink(const arma::mat& cov_mle, double n, double& rho_out) {
  const arma::uword p = cov_mle.n_cols;
  const double trS  = arma::trace(cov_mle);
  const double trS2 = arma::accu(cov_mle % cov_mle); // Frobenius^2 = tr(S^2)
  const double mu   = trS / static_cast<double>(p);

  const double pd  = static_cast<double>(p);
  const double num = (1.0 - 2.0 / pd) * trS2 + trS * trS;
  const double den = (n + 1.0 - 2.0 / pd) * (trS2 - (trS * trS) / pd);

  double rho = (den > 0.0) ? (num / den) : 1.0;
  rho = std::max(0.0, std::min(1.0, rho));
  rho_out = rho;

  arma::mat Sigma = (1.0 - rho) * cov_mle;
  Sigma.diag() += rho * mu; // add rho*mu*I
  return Sigma;
}

} // namespace cov_shrinkage

//------------------------------------------------------------------------------
// ---------- moments ----------
//------------------------------------------------------------------------------
namespace moments {

// Means and population variances (uses arma::var * ((n-1)/n) exactly as in your code)
inline void col_means_vars_pop(const arma::mat& X,
                               arma::vec& means,
                               arma::vec& vars_pop) {
  const arma::uword p = X.n_cols;
  const double n = static_cast<double>(X.n_rows);
  means.set_size(p);
  vars_pop.set_size(p);
  for (arma::uword j = 0; j < p; ++j) {
    const auto col = X.col(j);
    means[j]    = arma::mean(col);
    vars_pop[j] = arma::var(col) * ((n - 1.0) / n);
  }
}


// Population covariance via manual loop (matches ccc_cpp evaluation route)
inline double cov_xy_pop_manual(const arma::vec& x, const arma::vec& y,
                                double mean_x, double mean_y) {
  const arma::uword n = x.n_elem;
  double s = 0.0;
  for (arma::uword k = 0; k < n; ++k) s += (x[k] - mean_x) * (y[k] - mean_y);
  return s / static_cast<double>(n);
}

// Population covariance from arma::cov (unbiased) scaled to population
inline double cov_xy_pop_arma(const arma::vec& x, const arma::vec& y) {
  const double n = static_cast<double>(x.n_elem);
  return arma::as_scalar(arma::cov(x, y)) * ((n - 1.0) / n);
}

// Correlation from cov and population variances
inline double corr_from_covvar(double cov_xy, double var_x, double var_y) {
  const double denom = std::sqrt(var_x * var_y);
  return (denom > 0.0) ? (cov_xy / denom) : std::numeric_limits<double>::quiet_NaN();
}

} // namespace moments

//------------------------------------------------------------------------------
// ---------- ccc helpers ----------
//------------------------------------------------------------------------------

namespace ccc_bits {

// EXACT evaluation order as your current code (r -> sxy -> p).
inline double ccc_from_stats_via_r(double mean_x, double mean_y,
                                   double var_x,  double var_y,
                                   double cov_xy) {
  const double r   = cov_xy / std::sqrt(var_x * var_y);
  const double sxy = r * std::sqrt(var_x * var_y);  // equals cov_xy, but keeps order
  const double dmu = mean_x - mean_y;
  const double den = var_x + var_y + dmu * dmu;
  return (den > 0.0) ? (2.0 * sxy / den) : std::numeric_limits<double>::quiet_NaN();
}

// Algebraically equivalent, slightly fewer ops (use only if tiny FP diffs are acceptable)
inline double ccc_from_stats_via_cov(double mean_x, double mean_y,
                                     double var_x,  double var_y,
                                     double cov_xy) {
  const double dmu = mean_x - mean_y;
  const double den = var_x + var_y + dmu * dmu;
  return (den > 0.0) ? (2.0 * cov_xy / den) : std::numeric_limits<double>::quiet_NaN();
}

} // namespace ccc_bits

//------------------------------------------------------------------------------
// ---------- fisherz ----------
//------------------------------------------------------------------------------

namespace fisherz {
using matrixCorr_detail::clamp_policy::nan_preserve;

inline double z(double r) noexcept {
  const double one_minus = std::nextafter(1.0, 0.0);       // < 1.0
  const double neg_one_plus = std::nextafter(-1.0, 0.0);   // > -1.0
  r = nan_preserve(r, neg_one_plus, one_minus);
  return std::atanh(r);
}
inline double inv(double z) noexcept {
  double r = std::tanh(z);
  const double one_minus = std::nextafter(1.0, 0.0);
  if (r >  one_minus) r = one_minus;
  if (r < -one_minus) r = -one_minus;
  return r;
}

inline void ci_from_z(double p, double se_t, double zcrit,
                      double& lci, double& uci) noexcept {
  const double t   = z(p);
  const double lzt = t - zcrit * se_t;
  const double uzt = t + zcrit * se_t;
  lci = inv(lzt);
  uci = inv(uzt);
}
} // namespace fisherz

//------------------------------------------------------------------------------
// ---------- ccc_se ----------
//------------------------------------------------------------------------------

namespace ccc_se {
inline double se_delta(double r, double p, double u, int n) {
  if (n <= 2) return std::numeric_limits<double>::quiet_NaN();
  const double eps = 1e-12;
  if (!std::isfinite(r) || std::abs(r) < eps)
    return std::numeric_limits<double>::infinity();

  const double r2 = r * r;
  const double p2 = p * p;
  const double term =
    ((1.0 - r2) * p2 * (1.0 - p2) / r2)
    + (2.0 * p * p * p * (1.0 - p) * u * u / r)
    - (0.5 * p2 * p2 * u * u * u * u / r2);
    return std::sqrt(term / (static_cast<double>(n) - 2.0));
}
} // namespace ccc_se

//------------------------------------------------------------------------------
// ---------- symm ----------
//------------------------------------------------------------------------------
namespace symm {

// Assign both [i,j] and [j,i]
inline void put(arma::mat& M, arma::uword i, arma::uword j, double v) {
  M(i, j) = v; M(j, i) = v;
}

// Iterate (i<j) pairs; body(i,j) is a callable.
// Parallelise outer loop safely if desired.
template <class Body>
inline void for_upper_pairs(arma::uword p, Body body) {
  for (arma::uword i = 0; i < p; ++i)
    for (arma::uword j = i + 1; j < p; ++j)
      body(i, j);
}

} // namespace symm

//------------------------------------------------------------------------------
// ---------- covmat ----------
//------------------------------------------------------------------------------
namespace covmat {

// Returns population covariance: S = (X'X - n * mu * mu') / n
// population covariance: S = (X'X - n * mu * mu') / n
// Default = fast BLAS-3 two-pass path (your current implementation).
// Optionally: stable_centered=true -> center into a temporary then do SYRK.
// Optionally: method = "welford" -> one-pass, very stable, BLAS-2; good when p is small/moderate.
inline arma::mat cov_pop(const arma::mat& X,
                         arma::rowvec& mu_out,
                         bool stable_centered = false)
{
  const arma::uword n_rows = X.n_rows;
  const arma::uword p      = X.n_cols;

  if (n_rows == 0) {
    mu_out.set_size(p); mu_out.fill(arma::datum::nan);
    arma::mat out(p, p);
    out.fill(arma::datum::nan);
    return out;
  }

  const double n = static_cast<double>(n_rows);
  mu_out = arma::mean(X, 0);

  if (stable_centered) {
    // More numerically stable when |mu| is large relative to the spread:
    // compute S = (Xc'Xc)/n with Xc = X - 1*mu. Costs one p-by-n copy.
    arma::mat Xc = X.each_row() - mu_out;
    arma::mat XtXc = linalg::crossprod_no_copy(Xc);   // (X - mu)'(X - mu)
    XtXc /= n;
    return XtXc;
  }

  // Fast two-pass (your original): BLAS SYRK on X, then subtract n * mu * mu'
  arma::mat XtX = linalg::crossprod_no_copy(X);       // X'X
  linalg::subtract_n_outer_mu(XtX, mu_out, n);        // XtX -= n * mu * mu'
  XtX /= n;                                           // population scaling
  return XtX;
}

// One-pass Welford/Chan-Golub-LeVeque style population covariance.
// Very stable; updates mean and M2 with each row. Complexity O(n p^2).
// Prefer when p is small/medium.
inline arma::mat cov_pop_welford(const arma::mat& X, arma::rowvec& mu_out)
{
  const arma::uword n_rows = X.n_rows;
  const arma::uword p      = X.n_cols;

  if (n_rows == 0) {
    mu_out.set_size(p); mu_out.fill(arma::datum::nan);
    arma::mat out(p, p);
    out.fill(arma::datum::nan);
    return out;
  }

  arma::rowvec mean(p, arma::fill::zeros);
  arma::mat    M2(p, p, arma::fill::zeros);

  double k = 0.0;
  for (arma::uword i = 0; i < n_rows; ++i) {
    const arma::rowvec xi = X.row(i);
    const arma::rowvec delta  = xi - mean;
    k += 1.0;
    const arma::rowvec delta2 = delta / k;
    mean += delta2;
    // rank-1 update: M2 += (xi - mean_old)' * (xi - mean_new)
    // which equals (delta)' * (xi - mean) after the mean update.
    M2  += arma::trans(delta) * (xi - mean);
  }

  mu_out = mean;
  // population covariance
  M2 /= k;
  // enforce exact symmetry (accumulated FP noise)
  M2 = arma::symmatu(M2);
  return M2;
}

} // namespace covmat

//------------------------------------------------------------------------------
// ---------- quantile_utils ----------
//------------------------------------------------------------------------------
namespace quantile_utils {
// Weighted quantile for sorted values 's'
// If not sorted, sort both together before calling this function.
inline double weighted_quantile_sorted(const arma::vec &s, const arma::vec &ws, double p) {
  const std::size_t n = s.n_elem;
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  if (p <= 0.0) return s.front();
  if (p >= 1.0) return s.back();
  const double W = arma::accu(ws);
  if (!(W > 0.0)) return std::numeric_limits<double>::quiet_NaN();

  // cumulative weight target
  const double target = p * W;
  double csum = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    csum += ws[i];
    if (csum >= target) {
      // Optional linear interpolation to the previous point
      if (i == 0) return s[0];
      const double csum_prev = csum - ws[i];
      const double wgap = ws[i];
      if (wgap <= 0.0) return s[i];
      const double frac = (target - csum_prev) / wgap; // in [0,1]
      return s[i-1] + frac * (s[i] - s[i-1]);
    }
  }
  return s.back();
}

inline double weighted_median(const arma::vec &x, const arma::vec &w) {
  arma::uvec ord = arma::stable_sort_index(x);       // stable to keep determinism
  arma::vec xs = x.elem(ord), ws = w.elem(ord);
  return weighted_quantile_sorted(xs, ws, 0.5);
}

inline double weighted_mad(const arma::vec &x, const arma::vec &w, double med) {
  arma::vec dev = arma::abs(x - med);
  return weighted_median(dev, w);
}

// Linear interpolation for quantiles on a sorted vector (p in [0,1]).
inline double quantile_sorted(const arma::vec &s, double p) {
  const std::size_t n = s.n_elem;
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  if (p <= 0.0) return s.front();
  if (p >= 1.0) return s.back();
  double pos = p * (n - 1);
  std::size_t lo = static_cast<std::size_t>(std::floor(pos));
  std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
  double frac = pos - lo;
  return (1.0 - frac) * s[lo] + frac * s[hi];
}
} // namespace quantile_utils

//------------------------------------------------------------------------------
// ---------- standardise_bicor ----------
//------------------------------------------------------------------------------
namespace standardise_bicor {

// Standardise one column with *weights* (no NA, already subset if needed).
inline void standardise_bicor_column_weighted(const arma::vec& x,
                                              const arma::vec& wobs,
                                              arma::vec& z,
                                              int   pearson_fallback_mode,
                                              double c_const,
                                              double maxPOutliers,
                                              bool& col_is_valid)
{
  const std::size_t n = x.n_elem;
  if (z.n_elem != n) z.set_size(n);
  z.zeros();
  col_is_valid = false;

  if (n == 0) { z.fill(arma::datum::nan); return; }

  // Degenerate if all observation weights are zero
  const double Wtot = arma::accu(wobs);
  if (!(Wtot > 0.0)) { z.fill(arma::datum::nan); return; }

  // Force Pearson
  if (pearson_fallback_mode == 2) {
    const double mu = arma::dot(x, wobs) / Wtot;
    arma::vec centered = x - mu;
    const double denom2 = arma::dot(centered % wobs, centered % wobs);
    if (denom2 > 0.0) {
      z = (centered % wobs) / std::sqrt(denom2);
      col_is_valid = true;
    } else {
      z.fill(arma::datum::nan);
    }
    return;
  }

  // Weighted median and MAD
  const double med = matrixCorr_detail::quantile_utils::weighted_median(x, wobs);
  const double mad = matrixCorr_detail::quantile_utils::weighted_mad(x, wobs, med);

  if (!(mad > 0.0)) {
    if (pearson_fallback_mode == 0) { z.fill(arma::datum::nan); return; }
    // Pearson fallback with observation weights
    const double mu = arma::dot(x, wobs) / Wtot;
    arma::vec centered = x - mu;
    const double denom2 = arma::dot(centered % wobs, centered % wobs);
    if (!(denom2 > 0.0)) { z.fill(arma::datum::nan); return; }
    z = (centered % wobs) / std::sqrt(denom2);
    col_is_valid = true;
    return;
  }

  // Side-cap quantiles for rescaling
  double scale_neg = 1.0, scale_pos = 1.0;
  if (maxPOutliers < 1.0) {
    arma::uvec ord = arma::stable_sort_index(x);
    arma::vec xs = x.elem(ord), ws = wobs.elem(ord);
    const double qL = matrixCorr_detail::quantile_utils::weighted_quantile_sorted(xs, ws, maxPOutliers);
    const double qU = matrixCorr_detail::quantile_utils::weighted_quantile_sorted(xs, ws, 1.0 - maxPOutliers);
    const double uL = (qL - med) / (c_const * mad);
    const double uU = (qU - med) / (c_const * mad);
    if (std::abs(uL) > 1.0) scale_neg = std::abs(uL);
    if (std::abs(uU) > 1.0) scale_pos = std::abs(uU);
  }

  arma::vec xm = x - med;
  arma::vec u  = xm / (c_const * mad);
  if (maxPOutliers < 1.0 && (scale_neg > 1.0 || scale_pos > 1.0)) {
    for (std::size_t i = 0; i < n; ++i) {
      if (xm[i] < 0.0)      u[i] /= scale_neg;
      else if (xm[i] > 0.0) u[i] /= scale_pos;
    }
  }

  arma::vec wt(n, arma::fill::zeros);
  for (std::size_t i = 0; i < n; ++i) {
    const double a = u[i];
    if (std::abs(a) < 1.0) {
      const double t = (1.0 - a*a);
      wt[i] = t * t;
    }
  }
  wt %= wobs; // combine Tukey and observation weights

  arma::vec r = xm % wt;
  const double denom2 = arma::dot(r, r);
  if (!(denom2 > 0.0)) {
    if (pearson_fallback_mode >= 1) {
      const double mu = arma::dot(x, wobs) / Wtot;
      arma::vec centered = x - mu;
      const double d2 = arma::dot(centered % wobs, centered % wobs);
      if (d2 > 0.0) {
        z = (centered % wobs) / std::sqrt(d2);
        col_is_valid = true;
        return;
      }
    }
    z.fill(arma::datum::nan);
    return;
  }

  z = r / std::sqrt(denom2);
  col_is_valid = true;
}

// Core unweighted standardiser (exactly same maths/flow as your .cpp)
inline void standardise_bicor_column(const arma::vec& x,
                                     arma::vec& z,
                                     int   pearson_fallback_mode,
                                     double c_const,
                                     double maxPOutliers,
                                     bool& col_is_valid)
{
  const std::size_t n = x.n_elem;
  if (z.n_elem != n) z.set_size(n);
  z.zeros();
  col_is_valid = false;

  // Force Pearson for this column
  if (pearson_fallback_mode == 2) {
    const double mu = arma::mean(x);
    arma::vec centered = x - mu;
    const double denom2 = arma::dot(centered, centered);
    if (denom2 > 0.0) {
      z = centered / std::sqrt(denom2);
      col_is_valid = true;
    } else {
      z.fill(arma::datum::nan);
    }
    return;
  }

  // Median and MAD (use sorted copies for determinism)
  arma::vec xc = arma::sort(x);
  const double med = arma::median(xc);
  arma::vec absdev = arma::abs(x - med);
  const double mad = arma::median(arma::sort(absdev));

  // If MAD == 0, either NA or Pearson fallback
  if (!(mad > 0.0)) {
    if (pearson_fallback_mode == 0) { z.fill(arma::datum::nan); return; }
    const double mu = arma::mean(x);
    arma::vec centered = x - mu;
    const double denom2 = arma::dot(centered, centered);
    if (!(denom2 > 0.0)) { z.fill(arma::datum::nan); return; }
    z = centered / std::sqrt(denom2);
    col_is_valid = true;
    return;
  }

  // Side-cap quantiles so chosen tails map to |u| = 1 (Langfelder & Horvath)
  double scale_neg = 1.0, scale_pos = 1.0;
  if (maxPOutliers < 1.0) {
    const double qL = matrixCorr_detail::quantile_utils::quantile_sorted(xc, maxPOutliers);
    const double qU = matrixCorr_detail::quantile_utils::quantile_sorted(xc, 1.0 - maxPOutliers);
    const double uL = (qL - med) / (c_const * mad);
    const double uU = (qU - med) / (c_const * mad);
    if (std::abs(uL) > 1.0) scale_neg = std::abs(uL);
    if (std::abs(uU) > 1.0) scale_pos = std::abs(uU);
  }

  arma::vec xm = x - med;
  arma::vec u  = xm / (c_const * mad);
  if (maxPOutliers < 1.0 && (scale_neg > 1.0 || scale_pos > 1.0)) {
    for (std::size_t i = 0; i < n; ++i) {
      if (xm[i] < 0.0)      u[i] /= scale_neg;
      else if (xm[i] > 0.0) u[i] /= scale_pos;
    }
  }

  arma::vec w(n, arma::fill::zeros);
  for (std::size_t i = 0; i < n; ++i) {
    const double a = u[i];
    if (std::abs(a) < 1.0) {
      const double t = (1.0 - a*a);
      w[i] = t * t;
    }
  }
  arma::vec r = xm % w;

  const double denom2 = arma::dot(r, r);
  if (!(denom2 > 0.0)) {
    if (pearson_fallback_mode >= 1) {
      const double mu = arma::mean(x);
      arma::vec centered = x - mu;
      const double d2 = arma::dot(centered, centered);
      if (d2 > 0.0) {
        z = centered / std::sqrt(d2);
        col_is_valid = true;
        return;
      }
    }
    z.fill(arma::datum::nan);
    return;
  }

  z = r / std::sqrt(denom2);
  col_is_valid = true;
}

// Compatibility overload keeping your old trailing argument.
// Safe no-op forwarder; avoids touching existing call-sites.
inline void standardise_bicor_column(const arma::vec& x,
                                     arma::vec& z,
                                     int   pearson_fallback_mode,
                                     double c_const,
                                     double maxPOutliers,
                                     bool& col_is_valid,
                                     int /*n_threads_unused*/)
{
  standardise_bicor_column(x, z, pearson_fallback_mode, c_const, maxPOutliers,
                           col_is_valid);
}

} // namespace standardise_bicor

} // namespace matrixCorr_detail
