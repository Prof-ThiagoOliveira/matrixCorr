// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

namespace {

inline double clamp01(double x) {
  if (!std::isfinite(x)) return NA_REAL;
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

inline double pnorm_std(double x) {
  return R::pnorm(x, 0.0, 1.0, 1, 0);
}

inline double dnorm_std(double x) {
  return R::dnorm(x, 0.0, 1.0, 0);
}

inline double inv_link(double eta, int link_code) {
  if (link_code == 0) return eta;
  if (link_code == 1) {
    if (eta >= 0.0) {
      const double z = std::exp(-eta);
      return 1.0 / (1.0 + z);
    }
    const double z = std::exp(eta);
    return z / (1.0 + z);
  }
  return pnorm_std(eta);
}

inline double deriv_link(double eta, int link_code) {
  if (link_code == 0) return 1.0;
  if (link_code == 1) {
    const double p = inv_link(eta, link_code);
    return p * (1.0 - p);
  }
  return dnorm_std(eta);
}

inline double theta_from_diff(double diff, double se, double lower, double upper) {
  if (!std::isfinite(diff) || !std::isfinite(se) ||
      !std::isfinite(lower) || !std::isfinite(upper)) {
    return NA_REAL;
  }
  if (se == 0.0) {
    return (diff >= lower && diff <= upper) ? 1.0 : 0.0;
  }
  const double z_l = (lower - diff) / se;
  const double z_u = (upper - diff) / se;
  return clamp01(pnorm_std(z_u) - pnorm_std(z_l));
}

inline double se_from_var(double var) {
  if (!std::isfinite(var)) return NA_REAL;
  if (var < 0.0) {
    if (var > -1e-12) var = 0.0;
    else return NA_REAL;
  }
  return std::sqrt(std::max(0.0, var));
}

inline double variance_quad(const arma::vec& g1,
                            const arma::vec& g2,
                            const arma::mat& vcov1,
                            const arma::mat& vcov2,
                            const arma::mat& cov12) {
  return arma::as_scalar(g1.t() * vcov1 * g1) +
    arma::as_scalar(g2.t() * vcov2 * g2) -
    2.0 * arma::as_scalar(g1.t() * cov12 * g2);
}

inline double theta_from_betas(const arma::rowvec& xi,
                               const arma::vec& beta1,
                               const arma::vec& beta2,
                               const arma::mat& vcov1,
                               const arma::mat& vcov2,
                               const arma::mat& cov12,
                               int link_code,
                               double lower,
                               double upper) {
  const double eta1 = arma::as_scalar(xi * beta1);
  const double eta2 = arma::as_scalar(xi * beta2);
  const double est1 = inv_link(eta1, link_code);
  const double est2 = inv_link(eta2, link_code);
  arma::vec g1 = (deriv_link(eta1, link_code) * xi).t();
  arma::vec g2 = (deriv_link(eta2, link_code) * xi).t();
  const double se = se_from_var(variance_quad(g1, g2, vcov1, vcov2, cov12));
  return theta_from_diff(est1 - est2, se, lower, upper);
}

struct GlmFit {
  arma::vec beta;
  arma::mat vcov;
  int n_obs;
  bool converged;
  int iterations;
};

inline GlmFit fit_binomial_curve(const arma::vec& y,
                                 const arma::vec& x,
                                 int link_code,
                                 int max_iter,
                                 double tol) {
  const int n = y.n_elem;
  arma::mat X(n, 2);
  X.col(0).ones();
  X.col(1) = x;
  arma::vec beta(2, arma::fill::zeros);
  const double ybar = std::min(1.0 - 1e-6, std::max(1e-6, arma::mean(y)));
  beta[0] = link_code == 1 ? std::log(ybar / (1.0 - ybar)) : R::qnorm(ybar, 0.0, 1.0, 1, 0);
  bool converged = false;
  int iter = 0;
  arma::mat xtwx(2, 2, arma::fill::zeros);

  for (iter = 0; iter < max_iter; ++iter) {
    arma::vec eta = X * beta;
    arma::vec z(n), w(n);
    for (int i = 0; i < n; ++i) {
      double mu = inv_link(eta[i], link_code);
      mu = std::min(1.0 - 1e-8, std::max(1e-8, mu));
      double deta = deriv_link(eta[i], link_code);
      deta = std::max(1e-10, deta);
      w[i] = (deta * deta) / (mu * (1.0 - mu));
      z[i] = eta[i] + (y[i] - mu) / deta;
    }
    arma::mat Xw = X.each_col() % w;
    xtwx = X.t() * Xw;
    arma::vec xtwz = X.t() * (w % z);
    arma::vec beta_new;
    bool ok = arma::solve(beta_new, xtwx, xtwz, arma::solve_opts::likely_sympd);
    if (!ok || !beta_new.is_finite()) {
      ok = arma::solve(beta_new, xtwx, xtwz);
    }
    if (!ok || !beta_new.is_finite()) {
      stop("Binomial reliability model fit failed because the information matrix is singular.");
    }
    if (arma::max(arma::abs(beta_new - beta)) < tol) {
      beta = beta_new;
      converged = true;
      ++iter;
      break;
    }
    beta = beta_new;
  }

  arma::vec eta = X * beta;
  arma::vec w(n);
  for (int i = 0; i < n; ++i) {
    double mu = inv_link(eta[i], link_code);
    mu = std::min(1.0 - 1e-8, std::max(1e-8, mu));
    double deta = deriv_link(eta[i], link_code);
    deta = std::max(1e-10, deta);
    w[i] = (deta * deta) / (mu * (1.0 - mu));
  }
  arma::mat Xw = X.each_col() % w;
  xtwx = X.t() * Xw;
  arma::mat vcov;
  bool ok = arma::inv_sympd(vcov, xtwx);
  if (!ok || !vcov.is_finite()) {
    vcov = arma::pinv(xtwx);
  }

  return GlmFit{beta, vcov, n, converged, iter};
}

} // namespace

// [[Rcpp::export]]
Rcpp::List prob_agree_fit_cpp(const Rcpp::NumericVector& response,
                              const Rcpp::NumericVector& predictor,
                              const Rcpp::IntegerVector& group,
                              const Rcpp::NumericVector& eval_predictor,
                              int link_code,
                              const Rcpp::NumericVector& lower,
                              const Rcpp::NumericVector& upper,
                              bool ci,
                              int ci_method,
                              double conf_level,
                              int max_iter,
                              double tol) {
  const int n = response.size();
  if (predictor.size() != n || group.size() != n) {
    stop("response, predictor, and group must have the same length.");
  }
  std::vector<double> y1, y2, x1, x2;
  y1.reserve(n); y2.reserve(n); x1.reserve(n); x2.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(response[i]) || !std::isfinite(predictor[i])) continue;
    if (response[i] != 0.0 && response[i] != 1.0) {
      stop("response must contain only 0/1 values after removing missing values.");
    }
    if (group[i] == 1) {
      y1.push_back(response[i]); x1.push_back(predictor[i]);
    } else if (group[i] == 2) {
      y2.push_back(response[i]); x2.push_back(predictor[i]);
    }
  }
  if (y1.size() < 3 || y2.size() < 3) {
    stop("Each group must contain at least three complete observations.");
  }
  arma::vec ay1(y1), ay2(y2), ax1(x1), ax2(x2);
  if (arma::min(ay1) == arma::max(ay1) || arma::min(ay2) == arma::max(ay2)) {
    stop("Each group must contain both response outcomes (0 and 1) to fit a reliability curve.");
  }
  GlmFit fit1 = fit_binomial_curve(ay1, ax1, link_code, max_iter, tol);
  GlmFit fit2 = fit_binomial_curve(ay2, ax2, link_code, max_iter, tol);

  const int m = eval_predictor.size();
  if (lower.size() != m || upper.size() != m) {
    stop("Tolerance vectors must have length equal to eval_predictor.");
  }
  NumericVector est1(m), est2(m), diff(m), se_diff(m), se_theta(m, NA_REAL),
    theta(m), lwr(m, NA_REAL), upr(m, NA_REAL);
  const double alpha = 1.0 - conf_level;
  const double crit = R::qnorm(1.0 - alpha / 2.0, 0.0, 1.0, 1, 0);
  arma::mat cov12(2, 2, arma::fill::zeros);

  for (int i = 0; i < m; ++i) {
    arma::rowvec xi(2);
    xi[0] = 1.0;
    xi[1] = eval_predictor[i];
    const double eta1 = arma::as_scalar(xi * fit1.beta);
    const double eta2 = arma::as_scalar(xi * fit2.beta);
    est1[i] = inv_link(eta1, link_code);
    est2[i] = inv_link(eta2, link_code);
    diff[i] = est1[i] - est2[i];
    arma::vec g1 = (deriv_link(eta1, link_code) * xi).t();
    arma::vec g2 = (deriv_link(eta2, link_code) * xi).t();
    se_diff[i] = se_from_var(variance_quad(g1, g2, fit1.vcov, fit2.vcov, cov12));
    theta[i] = theta_from_diff(diff[i], se_diff[i], lower[i], upper[i]);

    if (ci_method == 1 && std::isfinite(theta[i])) {
      arma::vec grad(4, arma::fill::zeros);
      for (int j = 0; j < 2; ++j) {
        const double step1 = std::max(1e-6, std::abs(fit1.beta[j]) * 1e-5);
        arma::vec bp = fit1.beta;
        arma::vec bm = fit1.beta;
        bp[j] += step1;
        bm[j] -= step1;
        const double tp = theta_from_betas(xi, bp, fit2.beta, fit1.vcov, fit2.vcov, cov12, link_code, lower[i], upper[i]);
        const double tm = theta_from_betas(xi, bm, fit2.beta, fit1.vcov, fit2.vcov, cov12, link_code, lower[i], upper[i]);
        grad[j] = (tp - tm) / (2.0 * step1);

        const double step2 = std::max(1e-6, std::abs(fit2.beta[j]) * 1e-5);
        bp = fit2.beta;
        bm = fit2.beta;
        bp[j] += step2;
        bm[j] -= step2;
        const double tp2 = theta_from_betas(xi, fit1.beta, bp, fit1.vcov, fit2.vcov, cov12, link_code, lower[i], upper[i]);
        const double tm2 = theta_from_betas(xi, fit1.beta, bm, fit1.vcov, fit2.vcov, cov12, link_code, lower[i], upper[i]);
        grad[2 + j] = (tp2 - tm2) / (2.0 * step2);
      }
      arma::mat sigma(4, 4, arma::fill::zeros);
      sigma.submat(0, 0, 1, 1) = fit1.vcov;
      sigma.submat(2, 2, 3, 3) = fit2.vcov;
      se_theta[i] = se_from_var(arma::as_scalar(grad.t() * sigma * grad));
      if (ci && std::isfinite(se_theta[i])) {
        lwr[i] = clamp01(theta[i] - crit * se_theta[i]);
        upr[i] = clamp01(theta[i] + crit * se_theta[i]);
      }
    }
  }

  List out = List::create(
    _["predictor"] = clone(eval_predictor),
    _["prob_agree"] = theta
  );
  if (ci) {
    out["lwr.ci"] = lwr;
    out["upr.ci"] = upr;
  }
  out["beta1"] = wrap(fit1.beta);
  out["beta2"] = wrap(fit2.beta);
  out["vcov1"] = wrap(fit1.vcov);
  out["vcov2"] = wrap(fit2.vcov);
  out["n1"] = fit1.n_obs;
  out["n2"] = fit2.n_obs;
  out["converged"] = LogicalVector::create(fit1.converged, fit2.converged);
  out["iterations"] = IntegerVector::create(fit1.iterations, fit2.iterations);
  return out;
}
