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

// Tiny insertion sort for short ranges [s, e) (used within x-tie runs).
template <typename T>
inline void insertion_sort_range(T* a, int s, int e) noexcept {
  for (int i = s + 1; i < e; ++i) {
    T v = a[i]; int j = i - 1;
    while (j >= s && a[j] > v) { a[j + 1] = a[j]; --j; }
    a[j + 1] = v;
  }
}

// Count ties (“Ms”) on a sorted buffer (generic).
template <typename T>
inline long long getMs_T(const T* dat, int len) noexcept {
  long long Ms = 0, tie = 0;
  for (int i = 1; i < len; ++i) {
    if (dat[i] == dat[i - 1]) ++tie;
    else if (tie) { Ms += tie * (tie + 1) / 2; tie = 0; }
  }
  if (tie) Ms += tie * (tie + 1) / 2;
  return Ms;
}

// Convenience aliases.
inline long long getMs_ll(const long long* dat, int len) noexcept {
  return getMs_T<long long>(dat, len);
}
inline long long getMs_double(const double* dat, int len) noexcept {
  return getMs_T<double>(dat, len);
}

// Non-allocating inversion counter (bottom-up mergesort with insertion warm-up).
// Sorts a[0..n) in place (nondecreasing) and returns the inversion count.
// Requires external buffer buf[0..n).
template <typename T>
inline long long inv_count_inplace(T* __restrict__ a,
                                   T* __restrict__ buf,
                                   const int n) noexcept {
  if (n < 2) return 0LL;
  long long inv = 0;
  constexpr int TH = 96; // tune 64–128 if desired

  auto insertion_block = [&](int L, int R) noexcept {
    long long v_inv = 0;
    for (int i = L + 1; i <= R; ++i) {
      T v = a[i]; int j = i - 1;
      while (j >= L && a[j] > v) { a[j + 1] = a[j]; --j; ++v_inv; }
      a[j + 1] = v;
    }
    return v_inv;
  };

  // Warm-up small blocks with insertion sort (cache-friendly).
  for (int L = 0; L < n; L += TH) {
    const int R = std::min(L + TH - 1, n - 1);
    inv += insertion_block(L, R);
  }

  // Bottom-up merge; count cross inversions.
  for (int width = TH; width < n; width <<= 1) {
    for (int L = 0; L < n; L += 2 * width) {
      const int M = std::min(L + width, n);
      const int R = std::min(L + 2 * width, n);
      int i = L, j = M, k = L;
      while (i < M && j < R) {
        if (a[i] <= a[j]) buf[k++] = a[i++];
        else { buf[k++] = a[j++]; inv += (M - i); }
      }
      while (i < M) buf[k++] = a[i++];
      while (j < R) buf[k++] = a[j++];
      for (int t = L; t < R; ++t) a[t] = buf[t];
    }
  }
  return inv;
}

// Scalar fast path using raw doubles (no discretisation).
// Knight O(n log n): order by x, sort y within x-runs (Ms in runs),
// global inversion count on y|x + Ms on globally sorted y.
inline double tau_two_vectors_raw(const double* x,
                                  const double* y,
                                  int n) {
  if (n < 2) return NA_REAL;

  // Per-thread reusable buffers
  thread_local std::vector<int>     ord;
  thread_local std::vector<double>  ybuf, mrg;

  ord.resize(n);
  ybuf.resize(n);
  mrg.resize(n);

  // ord := argsort(x)
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) noexcept { return x[a] < x[b]; });

  // ybuf := y in x-sorted order
  for (int k = 0; k < n; ++k) ybuf[k] = y[ord[k]];

  // tie-runs in x (on x[ord[*]]); sort y within runs; accumulate m1 and within-run Ms
  long long m1 = 0, s_acc = 0;
  for (int s = 0; s < n; ) {
    int e = s + 1;
    const double xi = x[ord[s]];
    while (e < n && x[ord[e]] == xi) ++e;
    const int L = e - s;
    if (L > 1) {
      m1 += 1LL * L * (L - 1) / 2;
      // small runs: insertion sort; larger: std::sort
      if (L <= 32) insertion_sort_range<double>(ybuf.data(), s, e);
      else         std::sort(ybuf.begin() + s, ybuf.begin() + e);
      s_acc += getMs_double(ybuf.data() + s, L);
    }
    s = e;
  }

  // Global inversion count on y|x, and global Ms on sorted y
  const long long inv = inv_count_inplace(ybuf.data(), mrg.data(), n);
  const long long m2  = getMs_double(ybuf.data(), n);

  // Kendall tau-b in double (counts are 64-bit)
  const long long n0  = 1LL * n * (n - 1) / 2LL;
  const double    S   = double(n0) - double(m1) - double(m2)
    - 2.0 * double(inv) + double(s_acc);
  const double    den1 = double(n0 - m1);
  const double    den2 = double(n0 - m2);
  if (den1 <= 0.0 || den2 <= 0.0) return NA_REAL;
  return S / std::sqrt(den1 * den2);
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

// Robust inverse of SPD/SPSD with geometric diagonal jitter; SVD fallback.
// Returns true if a numerically finite inverse is obtained.
inline bool inv_sympd_safe(arma::mat& out,
                           const arma::mat& A,
                           double init_jitter = 1e-8,
                           double max_jitter  = 1e-2)
{
  if (arma::inv_sympd(out, A)) return true;

  arma::mat Aj = A;
  double base = 1.0;
  if (A.n_rows > 0) {
    const double tr = arma::trace(A);
    if (std::isfinite(tr) && tr > 0.0) base = std::max(1.0, tr / A.n_rows);
  }
  double lam = std::max(1e-12, init_jitter * base);

  for (; lam <= max_jitter; lam *= 10.0) {
    Aj = A; Aj.diag() += lam;
    if (arma::inv_sympd(out, Aj)) return true;
  }
  out = arma::pinv(A);
  return out.is_finite();
}

// Solve A * X = B for SPD/SPSD A with the same jitter/SVD strategy.
inline arma::mat solve_sympd_safe(const arma::mat& A,
                                  const arma::mat& B,
                                  double init_jitter = 1e-8,
                                  double max_jitter  = 1e-2)
{
  arma::mat X;
  // Fast path: try Cholesky/solve first
  {
    arma::mat L;
    if (arma::chol(L, A, "lower")) {
      X = arma::solve(arma::trimatl(L), B, arma::solve_opts::fast);
      X = arma::solve(arma::trimatu(L.t()), X, arma::solve_opts::fast);
      return X;
    }
  }
  // Fallback via inv_sympd_safe
  arma::mat Ai;
  if (!inv_sympd_safe(Ai, A, init_jitter, max_jitter)) {
    // As a last resort, pseudo-inverse
    Ai = arma::pinv(A);
  }
  return Ai * B;
}

// Log|A| for SPD/SPSD with jitter; SVD fallback (sums log positive singulars).
inline double logdet_spd_safe(const arma::mat& A,
                              double init_jitter = 1e-8,
                              double max_jitter  = 1e-2)
{
  arma::mat L;
  if (arma::chol(L, A, "lower"))
    return 2.0 * arma::sum(arma::log(L.diag()));

  arma::mat Aj = A;
  double base = 1.0;
  if (A.n_rows > 0) {
    const double tr = arma::trace(A);
    if (std::isfinite(tr) && tr > 0.0) base = std::max(1.0, tr / A.n_rows);
  }
  double lam = std::max(1e-12, init_jitter * base);

  for (; lam <= max_jitter; lam *= 10.0) {
    Aj = A; Aj.diag() += lam;
    if (arma::chol(L, Aj, "lower"))
      return 2.0 * arma::sum(arma::log(L.diag()));
  }
  arma::vec s;
  arma::mat U, V;
  arma::svd(U, s, V, A);
  double acc = 0.0;
  for (arma::uword i = 0; i < s.n_elem; ++i) if (s[i] > 0.0) acc += std::log(s[i]);
  return acc;
}

// Copy a subset of rows (by integer indices) into 'out'
inline void rows_take_to(const arma::mat& M,
                         const std::vector<int>& rows,
                         arma::mat& out)
{
  const int n_i = static_cast<int>(rows.size());
  const int q   = static_cast<int>(M.n_cols);
  out.set_size(n_i, q);
  for (int k = 0; k < n_i; ++k) out.row(k) = M.row(rows[k]);
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

// Unbiased sample variance (n-1 denominator); returns 0 for n<2 (as a safe default).
inline double sample_var(const arma::vec& v) {
  const arma::uword n = v.n_elem;
  if (n < 2u) return 0.0;
  const double mu = arma::mean(v);
  double acc = 0.0;
  for (arma::uword i = 0; i < n; ++i) { const double d = v[i] - mu; acc += d*d; }
  return acc / static_cast<double>(n - 1u);
}

inline double sample_cov(const arma::vec& a, const arma::vec& b) {
  const arma::uword n = std::min(a.n_elem, b.n_elem);
  if (n < 2u) return 0.0;
  const double ma = arma::mean(a);
  const double mb = arma::mean(b);
  double acc = 0.0;
  for (arma::uword i = 0; i < n; ++i) acc += (a[i] - ma) * (b[i] - mb);
  return acc / (double)(n - 1);
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

namespace timeseries {
namespace ar1 {

// κ_T = mean of all entries of the TxT AR(1) correlation matrix (used for time-averaged factors)
inline double kappa_T(double rho, int T) {
  if (T <= 1) return 1.0;
  const double r = rho;
  double acc = static_cast<double>(T);
  double rpow = r;
  for (int k = 1; k <= T - 1; ++k) { acc += 2.0 * static_cast<double>(T - k) * rpow; rpow *= r; }
  const double TT = static_cast<double>(T) * static_cast<double>(T);
  return std::max(acc / TT, 1e-12);
}

// Log|R| when residuals are block-diagonal AR(1) sequences separated by negatives in tim_i.
// se = σ_e, so log|Σ| = log|R| + n_i*log(se^2) = log|R| + 2*n_i*log(se) (we return the latter form).
inline double logdet_R_blocks(const std::vector<int>& tim_i,
                              double rho,
                              double se,
                              double eps = 1e-12)
{
  const int n_i = static_cast<int>(tim_i.size());
  const double lg_se = std::log(std::max(se, eps));
  const double lg_one_minus_r2 = std::log(std::max(1.0 - rho*rho, 1e-12));
  double out = 0.0;
  int s = 0;
  while (s < n_i) {
    if (tim_i[s] < 0) { out += lg_se; ++s; continue; }
    int e = s; while (e + 1 < n_i && tim_i[e + 1] >= 0) ++e;
    const int L = e - s + 1;
    out += L * lg_se;
    if (L >= 2) out += (L - 1) * lg_one_minus_r2;
    s = e + 1;
  }
  return out;
}

// Build the block-tri-diagonal R^{-1} for AR(1) sequences separated by negatives in tim_i.
// The matrix returned is symmetric and slightly jittered on the diagonal for safety.
inline void make_Cinv(const std::vector<int>& tim_i,
                      double rho,
                      arma::mat& Cinv)
{
  const int n_i = static_cast<int>(tim_i.size());
  Cinv.zeros(n_i, n_i);
  if (n_i == 0) return;

  const double r2 = rho * rho;
  const double denom = std::max(1.0 - r2, 1e-12);

  int s = 0;
  while (s < n_i) {
    if (tim_i[s] < 0) { Cinv(s, s) += 1.0; ++s; continue; }
    int e = s;
    while (e + 1 < n_i && tim_i[e + 1] >= 0) ++e;
    const int L = e - s + 1;

    if (L == 1) {
      Cinv(s, s) += 1.0;
    } else {
      Cinv(s, s) += 1.0 / denom;
      Cinv(e, e) += 1.0 / denom;
      for (int t = s + 1; t <= e - 1; ++t) Cinv(t, t) += (1.0 + r2) / denom;
      for (int t = s; t <= e - 1; ++t) {
        Cinv(t, t + 1) += -rho / denom;
        Cinv(t + 1, t) += -rho / denom;
      }
    }
    s = e + 1;
  }
  Cinv.diag() += 1e-10;
}

} // namespace matrixCorr_detail::timeseries::ar1
} // namespace matrixCorr_detail::timeseries

namespace indexing {

struct BySubject {
  std::vector<std::vector<int>> rows; // row indices per subject
  std::vector<std::vector<int>> met;  // method codes per row (−1 if missing)
  std::vector<std::vector<int>> tim;  // time   codes per row (−1 if missing)
};

// Reindex an integer label vector to 0..(m-1). NA handling optional.
// For Rcpp::IntegerVector pass na_value = NA_INTEGER; otherwise choose a sentinel (e.g., INT_MIN).
template <class IntVec>
inline void reindex(const IntVec& ids,
                    std::vector<int>& out_idx,
                    int& m_out,
                    int na_value = std::numeric_limits<int>::min())
{
  const int n = static_cast<int>(ids.size());
  out_idx.resize(n);
  std::unordered_map<int,int> map;
  map.reserve(static_cast<size_t>(n));
  int next = 0;

  for (int i = 0; i < n; ++i) {
    const int s = ids[i];
    if (s == na_value) Rcpp::stop("reindex(): input contains NA/sentinel");
    auto it = map.find(s);
    if (it == map.end()) { map.emplace(s, next); out_idx[i] = next; ++next; }
    else out_idx[i] = it->second;
  }
  m_out = next;
}

// Group rows by subject (already reindexed 0..m-1). Method/time may be empty vectors.
template <class MethodVec, class TimeVec>
inline BySubject group_by_subject(const std::vector<int>& subj_idx,
                                  const MethodVec& method,
                                  const TimeVec&  time,
                                  int m)
{
  BySubject S;
  S.rows.assign(m, {});
  S.met.assign(m, {});
  S.tim.assign(m, {});
  const int n = static_cast<int>(subj_idx.size());

  for (int i = 0; i < n; ++i) {
    const int j = subj_idx[i];
    S.rows[j].push_back(i);
    if (static_cast<int>(method.size()) > 0) {
      const int v = method[i];
      S.met[j].push_back(v == NA_INTEGER ? -1 : (v - 1));
    }
    if (static_cast<int>(time.size()) > 0) {
      const int v = time[i];
      S.tim[j].push_back(v == NA_INTEGER ? -1 : (v - 1));
    }
  }
  return S;
}

} // namespace matrixCorr_detail::indexing

namespace design {

// Build the "base" subject-level design U = [1 | method dummies | time dummies]
inline void build_U_base(const std::vector<int>& met_i,   // size n_i, −1 allowed
                         const std::vector<int>& tim_i,   // size n_i, −1 allowed
                         int nm,                          // #method levels used in U
                         int nt,                          // #time   levels used in U
                         arma::mat& U)                    // n_i x r_base
{
  const int n_i   = static_cast<int>(tim_i.size());
  const int rbase = 1 + (nm > 0 ? nm : 0) + (nt > 0 ? nt : 0);
  U.zeros(n_i, rbase);

  for (int t = 0; t < n_i; ++t) U(t, 0) = 1.0;

  int col = 1;
  if (nm > 0) {
    for (int t = 0; t < n_i; ++t) { const int l = met_i[t]; if (l >= 0) U(t, col + l) = 1.0; }
    col += nm;
  }
  if (nt > 0) {
    for (int t = 0; t < n_i; ++t) { const int tt = tim_i[t]; if (tt >= 0) U(t, col + tt) = 1.0; }
  }
}

// Gram matrix UtU for U = [1 | method dummies | time dummies]
inline void gram_UtU(const std::vector<int>& met_i,
                     const std::vector<int>& tim_i,
                     int n_i, int nm, int nt,
                     arma::mat& UtU)                      // r x r
{
  const int r = 1 + (nm > 0 ? nm : 0) + (nt > 0 ? nt : 0);
  UtU.zeros(r, r);
  UtU(0, 0) = n_i;

  if (nm > 0) {
    arma::vec cm(nm, arma::fill::zeros);
    for (int v : met_i) if (v >= 0) cm[v] += 1.0;
    for (int l = 0; l < nm; ++l) {
      UtU(0, 1 + l) = UtU(1 + l, 0) = cm[l];
      UtU(1 + l, 1 + l) = cm[l];
    }
    if (nt > 0) {
      arma::vec ct(nt, arma::fill::zeros);
      for (int v : tim_i) if (v >= 0) ct[v] += 1.0;
      for (int t = 0; t < nt; ++t) {
        const int jt = 1 + nm + t;
        UtU(0, jt) = UtU(jt, 0) = ct[t];
        UtU(jt, jt) = ct[t];
      }
      arma::mat cmt(nm, nt, arma::fill::zeros);
      for (int k = 0; k < n_i; ++k) {
        const int l = met_i[k], t = tim_i[k];
        if (l >= 0 && t >= 0) cmt(l, t) += 1.0;
      }
      for (int l = 0; l < nm; ++l)
        for (int t = 0; t < nt; ++t) {
          const int il = 1 + l;
          const int jt = 1 + nm + t;
          UtU(il, jt) = UtU(jt, il) = cmt(l, t);
        }
    }
  } else if (nt > 0) {
    arma::vec ct(nt, arma::fill::zeros);
    for (int v : tim_i) if (v >= 0) ct[v] += 1.0;
    for (int t = 0; t < nt; ++t) {
      const int jt = 1 + t;
      UtU(0, jt) = UtU(jt, 0) = ct[t];
      UtU(jt, jt) = ct[t];
    }
  }
}

// Accumulate U^T * v where v is accessed via a callable v(row_index).
template <class Accessor>
inline void accumulate_Ut_vec(const std::vector<int>& rows_i,
                              const std::vector<int>& met_i,
                              const std::vector<int>& tim_i,
                              int nm, int nt,
                              Accessor v, arma::vec& Utv)         // length r
{
  const int n_i = static_cast<int>(rows_i.size());
  const int r_expected = 1 + (nm > 0 ? nm : 0) + (nt > 0 ? nt : 0);
  Utv.zeros(r_expected);

  double s0 = 0.0;
  for (int ridx : rows_i) s0 += v(ridx);
  Utv[0] = s0;

  if (nm > 0) {
    for (int l = 0; l < nm; ++l) {
      double sm = 0.0;
      for (int k = 0; k < n_i; ++k) if (met_i[k] == l) sm += v(rows_i[k]);
      Utv[1 + l] = sm;
    }
  }
  if (nt > 0) {
    const int off = 1 + (nm > 0 ? nm : 0);
    for (int t = 0; t < nt; ++t) {
      double st = 0.0;
      for (int k = 0; k < n_i; ++k) if (tim_i[k] == t) st += v(rows_i[k]);
      Utv[off + t] = st;
    }
  }
}

// out += U * a (U as base design with optional method/time columns)
inline void add_U_times(const std::vector<int>& rows_i,
                        const std::vector<int>& met_i,
                        const std::vector<int>& tim_i,
                        int nm, int nt,
                        const arma::vec& a, arma::vec& out)        // out length n_i
{
  const int n_i = static_cast<int>(rows_i.size());
  const double a0 = a[0];

  std::vector<double> am(nm > 0 ? nm : 0), at(nt > 0 ? nt : 0);
  if (nm > 0) for (int l = 0; l < nm; ++l) am[l] = a[1 + l];
  if (nt > 0) for (int t = 0; t < nt; ++t) at[t] = a[1 + (nm > 0 ? nm : 0) + t];

  for (int k = 0; k < n_i; ++k) {
    double val = a0;
    if (nm > 0 && met_i[k] >= 0) val += am[ met_i[k] ];
    if (nt > 0 && tim_i[k] >= 0) val += at[ tim_i[k] ];
    out[k] += val;
  }
}

// out += Uextra * a_extra (no-op if empty)
inline void add_extra_times(const arma::mat& Uextra,
                            const arma::vec& a_extra,
                            arma::vec& out)
{
  if (Uextra.n_rows == 0 || Uextra.n_cols == 0) return;
  out += Uextra * a_extra;
}

} // namespace matrixCorr_detail::design
} // namespace matrixCorr_detail
