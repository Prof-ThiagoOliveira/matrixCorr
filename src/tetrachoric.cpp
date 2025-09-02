// Thiago de Paula Oliveira
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <functional>

#include "matrixCorr_detail.h"   // <-- new header

using namespace Rcpp;

namespace { // explicit policy choices for THIS translation unit

// choose clamp variant
inline double clamp(double x, double lo, double hi){
  return matrixCorr_detail::clamp_policy::nan_preserve(x, lo, hi);
}
// choose normal helpers
using matrixCorr_detail::norm1::Phi;
using matrixCorr_detail::norm1::phi;
using matrixCorr_detail::norm1::qnorm01;

// choose Brent implementation
inline double optimize_brent(std::function<double(double)> f,
                             double a, double b, double tol = 1e-6, int max_iter = 100){
  return matrixCorr_detail::brent::optimize(std::move(f), a, b, tol, max_iter);
}

// choose BVN policy
inline double Phi2(double a, double b, double rho){
  return matrixCorr_detail::bvn_adaptive::Phi2(a, b, rho);
}
inline double rect_prob(double al, double au, double bl, double bu, double rho){
  return matrixCorr_detail::bvn_adaptive::rect_prob(al, au, bl, bu, rho);
}

// --- log-likelihood helpers (use chosen policies) ---
double tetra_loglik(double rho, double a, double b, double c, double d){
  const double N = a + b + c + d;
  const double eps = 1e-12;
  double p_row1 = clamp((a + b) / N, eps, 1.0 - eps);
  double p_col1 = clamp((a + c) / N, eps, 1.0 - eps);
  double q1 = qnorm01(p_row1);
  double q2 = qnorm01(p_col1);

  double p11 = std::max(Phi2(q1, q2, rho), 1e-16);
  double p1_ = Phi(q1), _1p = Phi(q2);
  double p10 = std::max(p1_ - p11, 1e-16);
  double p01 = std::max(_1p - p11, 1e-16);
  double p00 = std::max(1.0 - p1_ - _1p + p11, 1e-16);

  return a*std::log(p11) + b*std::log(p10) + c*std::log(p01) + d*std::log(p00);
}

void build_cutpoints(const NumericMatrix& N, std::vector<double>& alpha, std::vector<double>& beta){
  int R = N.nrow(), C = N.ncol();
  NumericVector rowp(R), colp(C);
  double tot = 0.0;
  for (int i=0;i<R;++i) for (int j=0;j<C;++j){ rowp[i] += N(i,j); colp[j] += N(i,j); tot += N(i,j); }
  if (tot <= 0) stop("Empty table");

  rowp = rowp / tot; colp = colp / tot;
  const double eps = 1e-12;

  NumericVector rowcum(R-1), colcum(C-1);
  double acc = 0.0;
  for (int i=0;i<R-1;++i){ acc += rowp[i]; rowcum[i] = clamp(acc, eps, 1.0 - eps); }
  acc = 0.0;
  for (int j=0;j<C-1;++j){ acc += colp[j]; colcum[j] = clamp(acc, eps, 1.0 - eps); }

  alpha.assign(R+1, 0.0);
  beta.assign(C+1, 0.0);
  alpha[0] = -INFINITY; alpha[R] = INFINITY;
  beta[0]  = -INFINITY; beta[C]  = INFINITY;
  for (int i=1;i<=R-1;++i) alpha[i] = qnorm01(rowcum[i-1]);
  for (int j=1;j<=C-1;++j) beta[j]  = qnorm01(colcum[j-1]);
}

double polychoric_loglik(double rho, const NumericMatrix& N,
                         const std::vector<double>& alpha, const std::vector<double>& beta){
  int R = N.nrow(), C = N.ncol();
  double ll = 0.0;
  for (int i=1;i<=R; ++i){
    for (int j=1;j<=C; ++j){
      double pij = std::max(rect_prob(alpha[i-1], alpha[i], beta[j-1], beta[j], rho), 1e-16);
      ll += N(i-1, j-1) * std::log(pij);
    }
  }
  return ll;
}

}

// ---------- exports----------

// [[Rcpp::export]]
double matrixCorr_tetrachoric_mle_cpp(NumericMatrix tab, double correct = 0.5){
  if (tab.nrow() != 2 || tab.ncol() != 2) return NA_REAL;

  double a = tab(0,0), b = tab(0,1), c = tab(1,0), d = tab(1,1);
  if (a <= 0) a += correct;
  if (b <= 0) b += correct;
  if (c <= 0) c += correct;
  if (d <= 0) d += correct;

  auto nll = [&](double rho){ return -tetra_loglik(rho, a, b, c, d); };
  double est = optimize_brent(nll, -0.99, 0.99);
  return clamp(est, -1.0, 1.0);
}

// [[Rcpp::export]]
double matrixCorr_polychoric_mle_cpp(NumericMatrix tab, double correct = 0.5){
  int R = tab.nrow(), C = tab.ncol();
  if (R < 2 || C < 2) return NA_REAL;

  NumericMatrix N = clone(tab);
  for (int i=0;i<R;++i) for (int j=0;j<C;++j) if (N(i,j) <= 0.0) N(i,j) += correct;

  std::vector<double> alpha, beta;
  build_cutpoints(N, alpha, beta);

  auto nll = [&](double rho){ return -polychoric_loglik(rho, N, alpha, beta); };
  double est = optimize_brent(nll, -0.99, 0.99);
  return clamp(est, -1.0, 1.0);
}

// [[Rcpp::export]]
double matrixCorr_biserial_latent_cpp(NumericVector x, LogicalVector y){
  int n = x.size();
  if (n != y.size() || n < 2) return NA_REAL;

  double mx = 0.0; for (int i=0;i<n;++i) mx += x[i]; mx /= n;
  double v  = 0.0; for (int i=0;i<n;++i){ double d = x[i]-mx; v += d*d; }
  if (v <= 0.0) return NA_REAL;
  double sx = std::sqrt(v / (n - 1));

  int n1 = 0, n0 = 0; double s1 = 0.0, s0 = 0.0;
  for (int i=0;i<n;++i){
    if (y[i]){ s1 += x[i]; n1++; }
    else     { s0 += x[i]; n0++; }
  }
  if (n1 == 0 || n0 == 0) return NA_REAL;

  double mean1 = s1 / n1, mean0 = s0 / n0;
  double p = clamp(double(n1)/n, 1e-12, 1.0 - 1e-12), q = 1.0 - p;
  double zp = phi(qnorm01(p));
  double r = ((mean1 - mean0) / sx) * std::sqrt(p*q) / zp;
  return clamp(r, -1.0, 1.0);
}

// [[Rcpp::export]]
double matrixCorr_polyserial_mle_cpp(NumericVector x, IntegerVector y){
  int n = x.size();
  if (n != y.size() || n < 2) return NA_REAL;

  // standardize x
  double mx = 0.0; for (int i=0;i<n;++i) mx += x[i]; mx /= n;
  double v  = 0.0; for (int i=0;i<n;++i){ double d = x[i]-mx; v += d*d; }
  if (v <= 0.0) return NA_REAL;
  double sx = std::sqrt(v / (n - 1));
  std::vector<double> xs(n);
  for (int i=0;i<n;++i) xs[i] = (x[i] - mx) / sx;

  // map y to 1..K and build cutpoints from marginals
  IntegerVector yu = clone(y);
  int ymin = yu[0], ymax = yu[0];
  for (int i=1;i<n;++i){ if (yu[i]<ymin) ymin=yu[i]; if (yu[i]>ymax) ymax=yu[i]; }
  for (int i=0;i<n;++i) yu[i] = yu[i] - ymin + 1;
  int K = ymax - ymin + 1;
  if (K < 2) return NA_REAL;

  std::vector<double> beta(K+1);
  beta[0] = -INFINITY; beta[K] = INFINITY;
  std::vector<double> freq(K, 0.0);
  for (int i=0;i<n;++i) freq[ yu[i]-1 ] += 1.0;
  const double eps = 1e-12;
  double acc = 0.0;
  for (int k=1; k<=K-1; ++k){
    acc += freq[k-1] / n;
    acc = clamp(acc, eps, 1.0 - eps);
    beta[k] = qnorm01(acc);
  }

  auto nll = [&](double rho){
    rho = clamp(rho, -0.999999, 0.999999);
    double sig = std::sqrt(1.0 - rho*rho);
    double negLL = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:negLL)
#endif
    for (int i=0; i<n; ++i){
      int k = yu[i];
      double xz = xs[i];
      double u = (beta[k]   - rho * xz) / sig;
      double l = (beta[k-1] - rho * xz) / sig;
      double pk = clamp(Phi(u) - Phi(l), 1e-16, 1.0);
      negLL += -std::log(pk);
    }
    return negLL;
  };

  double est = optimize_brent(nll, -0.99, 0.99);
  return clamp(est, -1.0, 1.0);
}

// [[Rcpp::export]]
double matrixCorr_polydi_mle_cpp(NumericMatrix tab, double correct = 0.5){
  int R = tab.nrow(), C = tab.ncol();
  if (R < 2 || C != 2) return NA_REAL;
  return matrixCorr_polychoric_mle_cpp(tab, correct);
}
