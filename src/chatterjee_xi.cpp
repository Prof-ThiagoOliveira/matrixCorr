// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
#include "matrixCorr_omp.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

namespace {

inline double na_nan() {
  return NA_REAL;
}

inline bool is_finite_d(const double x) {
  return std::isfinite(x);
}

void shuffle_slice(std::vector<int>& ord, const int first, const int last) {
  for (int i = last - 1; i > first; --i) {
    const double u = R::runif(0.0, 1.0);
    const int j = first + static_cast<int>(std::floor(u * static_cast<double>(i - first + 1)));
    std::swap(ord[static_cast<std::size_t>(i)], ord[static_cast<std::size_t>(j)]);
  }
}

void sort_order_by_x(const std::vector<double>& x,
                     std::vector<int>& ord,
                     const std::string& tie_method) {
  const int n = static_cast<int>(x.size());
  ord.resize(static_cast<std::size_t>(n));
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](const int a, const int b) {
    return x[static_cast<std::size_t>(a)] < x[static_cast<std::size_t>(b)];
  });

  if (tie_method != "random") return;

  int start = 0;
  while (start < n) {
    int end = start + 1;
    const double v = x[static_cast<std::size_t>(ord[static_cast<std::size_t>(start)])];
    while (end < n && x[static_cast<std::size_t>(ord[static_cast<std::size_t>(end)])] == v) {
      ++end;
    }
    if (end - start > 1) {
      shuffle_slice(ord, start, end);
    }
    start = end;
  }
}

double xi_raw_complete(const std::vector<double>& x,
                       const std::vector<double>& y,
                       const std::string& tie_method) {
  const int n = static_cast<int>(x.size());
  if (n < 2 || y.size() != x.size()) return na_nan();

  double ymin = y[0];
  double ymax = y[0];
  for (int i = 1; i < n; ++i) {
    const double v = y[static_cast<std::size_t>(i)];
    if (v < ymin) ymin = v;
    if (v > ymax) ymax = v;
  }
  if (!(ymin < ymax)) return na_nan();

  std::vector<int> ord;
  sort_order_by_x(x, ord, tie_method);

  std::vector<double> y_sorted = y;
  std::sort(y_sorted.begin(), y_sorted.end());

  std::vector<double> r(static_cast<std::size_t>(n));
  double denom_sum = 0.0;
  for (int pos = 0; pos < n; ++pos) {
    const int idx = ord[static_cast<std::size_t>(pos)];
    const double yi = y[static_cast<std::size_t>(idx)];
    const double le = static_cast<double>(
      std::upper_bound(y_sorted.begin(), y_sorted.end(), yi) - y_sorted.begin()
    );
    const double lt = static_cast<double>(
      std::lower_bound(y_sorted.begin(), y_sorted.end(), yi) - y_sorted.begin()
    );
    const double li = static_cast<double>(n) - lt;
    r[static_cast<std::size_t>(pos)] = le;
    denom_sum += li * (static_cast<double>(n) - li);
  }

  if (!(denom_sum > 0.0) || !std::isfinite(denom_sum)) return na_nan();

  double diff_sum = 0.0;
  for (int i = 0; i < n - 1; ++i) {
    diff_sum += std::abs(r[static_cast<std::size_t>(i + 1)] - r[static_cast<std::size_t>(i)]);
  }

  const double numerator = static_cast<double>(n) * diff_sum;
  const double denominator = 2.0 * denom_sum;
  return 1.0 - numerator / denominator;
}

double xi_complete(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::string& tie_method,
                   const std::string& bias_correction) {
  const double raw = xi_raw_complete(x, y, tie_method);
  if (!std::isfinite(raw)) return raw;
  if (bias_correction != "upper_bound") return raw;

  const double upper = xi_raw_complete(y, y, tie_method);
  if (!(upper > 0.0) || !std::isfinite(upper)) return na_nan();
  const double corrected = raw / upper;
  return corrected < -1.0 ? -1.0 : corrected;
}

double xi_from_ptrs(const double* x,
                    const double* y,
                    const int n,
                    const std::string& tie_method,
                    const std::string& bias_correction) {
  std::vector<double> xv;
  std::vector<double> yv;
  xv.reserve(static_cast<std::size_t>(n));
  yv.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    const double yi = y[i];
    if (is_finite_d(xi) && is_finite_d(yi)) {
      xv.push_back(xi);
      yv.push_back(yi);
    }
  }
  return xi_complete(xv, yv, tie_method, bias_correction);
}

bool matrix_all_finite(const arma::mat& X) {
  const double* ptr = X.memptr();
  const arma::uword n_elem = X.n_elem;
  for (arma::uword i = 0; i < n_elem; ++i) {
    if (!is_finite_d(ptr[i])) return false;
  }
  return true;
}

bool columns_have_no_ties(const arma::mat& X) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  std::vector<double> values(static_cast<std::size_t>(n));
  for (arma::uword j = 0; j < p; ++j) {
    const double* col = X.colptr(j);
    for (arma::uword i = 0; i < n; ++i) {
      values[static_cast<std::size_t>(i)] = col[i];
    }
    std::sort(values.begin(), values.end());
    for (arma::uword i = 1; i < n; ++i) {
      if (values[static_cast<std::size_t>(i)] == values[static_cast<std::size_t>(i - 1)]) {
        return false;
      }
    }
  }
  return true;
}

struct XiPreparedColumn {
  std::vector<int> order;
  std::vector<double> rank_le;
  double denom = std::numeric_limits<double>::quiet_NaN();
  double self_raw = std::numeric_limits<double>::quiet_NaN();
  bool usable = false;
};

void prepare_order_first(const arma::mat& X,
                         const arma::uword column,
                         std::vector<int>& order) {
  const int n = static_cast<int>(X.n_rows);
  const double* x = X.colptr(column);
  order.resize(static_cast<std::size_t>(n));
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](const int a, const int b) {
    return x[a] < x[b];
  });
}

void prepare_response_ranks(const arma::mat& X,
                            const arma::uword column,
                            XiPreparedColumn& prep) {
  const int n = static_cast<int>(X.n_rows);
  const double* y = X.colptr(column);
  std::vector<double> y_sorted(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    y_sorted[static_cast<std::size_t>(i)] = y[i];
  }
  std::sort(y_sorted.begin(), y_sorted.end());

  prep.rank_le.resize(static_cast<std::size_t>(n));
  prep.denom = 0.0;
  for (int i = 0; i < n; ++i) {
    const double yi = y[i];
    const double le = static_cast<double>(
      std::upper_bound(y_sorted.begin(), y_sorted.end(), yi) - y_sorted.begin()
    );
    const double lt = static_cast<double>(
      std::lower_bound(y_sorted.begin(), y_sorted.end(), yi) - y_sorted.begin()
    );
    const double li = static_cast<double>(n) - lt;
    prep.rank_le[static_cast<std::size_t>(i)] = le;
    prep.denom += li * (static_cast<double>(n) - li);
  }
  prep.usable = prep.denom > 0.0 && std::isfinite(prep.denom);
}

double xi_raw_from_prepared(const std::vector<int>& order,
                            const XiPreparedColumn& response,
                            const int n) {
  if (!response.usable) return na_nan();
  double diff_sum = 0.0;
  const std::vector<double>& ranks = response.rank_le;
  for (int i = 0; i < n - 1; ++i) {
    const int idx0 = order[static_cast<std::size_t>(i)];
    const int idx1 = order[static_cast<std::size_t>(i + 1)];
    diff_sum += std::abs(
      ranks[static_cast<std::size_t>(idx1)] - ranks[static_cast<std::size_t>(idx0)]
    );
  }
  const double numerator = static_cast<double>(n) * diff_sum;
  const double denominator = 2.0 * response.denom;
  return 1.0 - numerator / denominator;
}

arma::mat chatterjee_xi_matrix_prepared_cpp(const arma::mat& X,
                                            const std::string& bias_correction,
                                            const int n_threads) {
  const int n = static_cast<int>(X.n_rows);
  const arma::uword p = X.n_cols;
  std::vector<XiPreparedColumn> prep(static_cast<std::size_t>(p));

  for (arma::uword j = 0; j < p; ++j) {
    prepare_order_first(X, j, prep[static_cast<std::size_t>(j)].order);
    prepare_response_ranks(X, j, prep[static_cast<std::size_t>(j)]);
  }

  for (arma::uword j = 0; j < p; ++j) {
    XiPreparedColumn& col = prep[static_cast<std::size_t>(j)];
    col.self_raw = xi_raw_from_prepared(col.order, col, n);
  }

  arma::mat out(p, p, arma::fill::none);
  out.fill(arma::datum::nan);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads) if(n_threads > 1)
#endif
  for (arma::sword ii = 0; ii < static_cast<arma::sword>(p); ++ii) {
    const arma::uword i = static_cast<arma::uword>(ii);
    const std::vector<int>& order = prep[static_cast<std::size_t>(i)].order;
    for (arma::uword j = 0; j < p; ++j) {
      const XiPreparedColumn& response = prep[static_cast<std::size_t>(j)];
      const double raw = xi_raw_from_prepared(order, response, n);
      if (!std::isfinite(raw) || bias_correction != "upper_bound") {
        out(i, j) = raw;
        continue;
      }
      const double upper = response.self_raw;
      if (!(upper > 0.0) || !std::isfinite(upper)) {
        out(i, j) = na_nan();
        continue;
      }
      const double corrected = raw / upper;
      out(i, j) = corrected < -1.0 ? -1.0 : corrected;
    }
  }

  return out;
}

void sample_without_replacement(const int n, const int m, std::vector<int>& out) {
  out.resize(static_cast<std::size_t>(n));
  std::iota(out.begin(), out.end(), 0);
  for (int i = 0; i < m; ++i) {
    const double u = R::runif(0.0, 1.0);
    const int j = i + static_cast<int>(std::floor(u * static_cast<double>(n - i)));
    std::swap(out[static_cast<std::size_t>(i)], out[static_cast<std::size_t>(j)]);
  }
  out.resize(static_cast<std::size_t>(m));
}

double sample_sd(const std::vector<double>& x) {
  const int n = static_cast<int>(x.size());
  if (n < 2) return na_nan();
  double mean = 0.0;
  for (double v : x) mean += v;
  mean /= static_cast<double>(n);
  double ss = 0.0;
  for (double v : x) {
    const double d = v - mean;
    ss += d * d;
  }
  return std::sqrt(ss / static_cast<double>(n - 1));
}

double quantile_type7(std::vector<double> x, const double prob) {
  x.erase(std::remove_if(x.begin(), x.end(), [](double v) { return !std::isfinite(v); }), x.end());
  const int n = static_cast<int>(x.size());
  if (n == 0) return na_nan();
  if (n == 1) return x[0];
  std::sort(x.begin(), x.end());
  const double h = 1.0 + (static_cast<double>(n) - 1.0) * prob;
  const int lo = static_cast<int>(std::floor(h));
  const double frac = h - static_cast<double>(lo);
  if (lo <= 1) return x[0] + frac * (x[1] - x[0]);
  if (lo >= n) return x[static_cast<std::size_t>(n - 1)];
  const double x0 = x[static_cast<std::size_t>(lo - 1)];
  const double x1 = x[static_cast<std::size_t>(lo)];
  return x0 + frac * (x1 - x0);
}

int choose_m(const int n, const int user_m, const std::string& method) {
  if (n < 3) return NA_INTEGER;
  int m = user_m;
  if (m == NA_INTEGER || m <= 0) {
    if (method == "n_choose_m") {
      m = static_cast<int>(std::round(2.0 * std::sqrt(static_cast<double>(n))));
    } else {
      m = static_cast<int>(std::floor(std::sqrt(static_cast<double>(n))));
    }
  }
  if (m < 2) m = 2;
  if (m > n - 1) m = n - 1;
  return m;
}

void bootstrap_ci_pair(const std::vector<double>& x,
                       const std::vector<double>& y,
                       const double xi_n,
                       const double conf_level,
                       const std::string& method,
                       const int bootstrap_reps,
                       const int m,
                       const std::string& tie_method,
                       const std::string& bias_correction,
                       double& lwr,
                       double& upr,
                       double& se) {
  lwr = na_nan();
  upr = na_nan();
  se = na_nan();
  const int n = static_cast<int>(x.size());
  if (!std::isfinite(xi_n) || n < 3 || m < 2 || m >= n || bootstrap_reps < 2) return;

  std::vector<int> sampled;
  std::vector<double> xb(static_cast<std::size_t>(m));
  std::vector<double> yb(static_cast<std::size_t>(m));
  std::vector<double> boot;
  boot.reserve(static_cast<std::size_t>(bootstrap_reps));

  for (int b = 0; b < bootstrap_reps; ++b) {
    sample_without_replacement(n, m, sampled);
    for (int k = 0; k < m; ++k) {
      const int idx = sampled[static_cast<std::size_t>(k)];
      xb[static_cast<std::size_t>(k)] = x[static_cast<std::size_t>(idx)];
      yb[static_cast<std::size_t>(k)] = y[static_cast<std::size_t>(idx)];
    }
    const double val = xi_complete(xb, yb, tie_method, bias_correction);
    if (std::isfinite(val)) boot.push_back(val);
  }

  if (boot.size() < 2u) return;
  const double alpha = 1.0 - conf_level;
  const double z = R::qnorm(1.0 - alpha / 2.0, 0.0, 1.0, 1, 0);

  if (method == "dette_kroll") {
    const double sd_boot = sample_sd(boot);
    if (!std::isfinite(sd_boot)) return;
    const double sigma_hat = std::sqrt(static_cast<double>(m)) * sd_boot;
    se = sigma_hat / std::sqrt(static_cast<double>(n));
    lwr = xi_n - z * se;
    upr = xi_n + z * se;
    return;
  }

  std::vector<double> tstar;
  tstar.reserve(boot.size());
  const double scale_m = std::sqrt(static_cast<double>(m));
  for (double val : boot) {
    tstar.push_back(scale_m * (val - xi_n));
  }
  const double q_hi = quantile_type7(tstar, 1.0 - alpha / 2.0);
  const double q_lo = quantile_type7(tstar, alpha / 2.0);
  if (!std::isfinite(q_hi) || !std::isfinite(q_lo)) return;

  // Finite m-out-of-n bootstrap samples can put both quantiles on the same
  // side of zero. Anchor the inversion at zero so the reported interval keeps
  // the observed estimate inside without truncating to the parameter range.
  const double q_hi_anchor = std::max(q_hi, 0.0);
  const double q_lo_anchor = std::min(q_lo, 0.0);
  const double scale_n = std::sqrt(static_cast<double>(n));
  lwr = xi_n - q_hi_anchor / scale_n;
  upr = xi_n - q_lo_anchor / scale_n;
  se = sample_sd(tstar) / scale_n;
}

void validate_options(const std::string& tie_method,
                      const std::string& bias_correction) {
  if (tie_method != "random" && tie_method != "first") {
    Rcpp::stop("tie_method must be 'random' or 'first'.");
  }
  if (bias_correction != "none" && bias_correction != "upper_bound") {
    Rcpp::stop("bias_correction must be 'none' or 'upper_bound'.");
  }
}

} // namespace

// [[Rcpp::export]]
double chatterjee_xi_vec_cpp(Rcpp::NumericVector x,
                             Rcpp::NumericVector y,
                             std::string tie_method = "random",
                             std::string bias_correction = "none") {
  validate_options(tie_method, bias_correction);
  if (x.size() != y.size()) Rcpp::stop("x and y must have the same length.");
  return xi_from_ptrs(REAL(x), REAL(y), x.size(), tie_method, bias_correction);
}

// [[Rcpp::export]]
arma::mat chatterjee_xi_matrix_cpp(const arma::mat& X,
                                   std::string tie_method = "random",
                                   std::string bias_correction = "none",
                                   int n_threads = 1) {
  validate_options(tie_method, bias_correction);
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (n < 2u || p < 2u) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  if (matrix_all_finite(X) && (tie_method == "first" || columns_have_no_ties(X))) {
    return chatterjee_xi_matrix_prepared_cpp(X, bias_correction, n_threads);
  }

  arma::mat out(p, p, arma::fill::none);
  out.fill(arma::datum::nan);

  const bool can_parallel = tie_method == "first" && n_threads > 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads) if(can_parallel)
#endif
  for (arma::sword i = 0; i < static_cast<arma::sword>(p); ++i) {
    const arma::uword ui = static_cast<arma::uword>(i);
    const double* xi = X.colptr(ui);
    for (arma::uword j = 0; j < p; ++j) {
      const double* yj = X.colptr(j);
      out(ui, j) = xi_from_ptrs(
        xi,
        yj,
        static_cast<int>(n),
        tie_method,
        bias_correction
      );
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List chatterjee_xi_matrix_pairwise_cpp(SEXP X_,
                                             bool return_ci = false,
                                             double conf_level = 0.95,
                                             std::string ci_method = "auto",
                                             int bootstrap_reps = 999,
                                             Rcpp::Nullable<Rcpp::IntegerVector> m = R_NilValue,
                                             int large_sample_cutoff = 1000,
                                             std::string tie_method = "random",
                                             std::string bias_correction = "none",
                                             int n_threads = 1) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_)) Rcpp::stop("Numeric double matrix required.");
  validate_options(tie_method, bias_correction);
  if (ci_method != "auto" && ci_method != "dette_kroll" && ci_method != "n_choose_m") {
    Rcpp::stop("ci_method must be 'auto', 'dette_kroll', or 'n_choose_m'.");
  }
  if (return_ci && !(conf_level > 0.0 && conf_level < 1.0)) {
    Rcpp::stop("conf_level must be in (0,1).");
  }
  if (return_ci && bootstrap_reps < 2) {
    Rcpp::stop("bootstrap_reps must be >= 2 when return_ci is TRUE.");
  }

  const arma::uword n = static_cast<arma::uword>(Rf_nrows(X_));
  const arma::uword p = static_cast<arma::uword>(Rf_ncols(X_));
  if (n < 2u || p < 2u) Rcpp::stop("Need >= 2 rows and >= 2 columns.");
  arma::mat X(REAL(X_), n, p, false, true);

  int user_m = NA_INTEGER;
  if (m.isNotNull()) {
    Rcpp::IntegerVector mv(m);
    if (mv.size() > 0) user_m = mv[0];
  }

  arma::mat est(p, p, arma::fill::none);
  est.fill(arma::datum::nan);
  arma::Mat<int> n_complete(p, p, arma::fill::zeros);
  arma::mat lwr(p, p, arma::fill::none);
  arma::mat upr(p, p, arma::fill::none);
  arma::mat se(p, p, arma::fill::none);
  arma::Mat<int> m_mat(p, p, arma::fill::zeros);
  Rcpp::CharacterMatrix method_mat(p, p);
  lwr.fill(arma::datum::nan);
  upr.fill(arma::datum::nan);
  se.fill(arma::datum::nan);

  std::vector<arma::uvec> finite_idx(static_cast<std::size_t>(p));
  for (arma::uword j = 0; j < p; ++j) {
    finite_idx[static_cast<std::size_t>(j)] = arma::find_finite(X.col(j));
  }

  const bool can_parallel = !return_ci && tie_method == "first" && n_threads > 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads) if(can_parallel)
#endif
  for (arma::sword ii = 0; ii < static_cast<arma::sword>(p); ++ii) {
    const arma::uword i = static_cast<arma::uword>(ii);
    for (arma::uword j = 0; j < p; ++j) {
      const arma::uvec& idx_i = finite_idx[static_cast<std::size_t>(i)];
      const arma::uvec& idx_j = finite_idx[static_cast<std::size_t>(j)];
      const arma::uword* pi = idx_i.memptr();
      const arma::uword* pj = idx_j.memptr();
      const double* xi = X.colptr(i);
      const double* yj = X.colptr(j);

      std::vector<double> xv;
      std::vector<double> yv;
      xv.reserve(static_cast<std::size_t>(std::min(idx_i.n_elem, idx_j.n_elem)));
      yv.reserve(static_cast<std::size_t>(std::min(idx_i.n_elem, idx_j.n_elem)));

      arma::uword ai = 0u;
      arma::uword aj = 0u;
      while (ai < idx_i.n_elem && aj < idx_j.n_elem) {
        const arma::uword ri = pi[ai];
        const arma::uword rj = pj[aj];
        if (ri == rj) {
          xv.push_back(xi[ri]);
          yv.push_back(yj[rj]);
          ++ai;
          ++aj;
        } else if (ri < rj) {
          ++ai;
        } else {
          ++aj;
        }
      }

      const int nc = static_cast<int>(xv.size());
      n_complete(i, j) = nc;
      if (nc < 2) continue;

      const double val = xi_complete(xv, yv, tie_method, bias_correction);
      est(i, j) = val;

      if (return_ci) {
        std::string method_use = ci_method;
        if (method_use == "auto") {
          method_use = (nc > large_sample_cutoff) ? "n_choose_m" : "dette_kroll";
        }
        const int mi = choose_m(nc, user_m, method_use);
        m_mat(i, j) = mi;
        method_mat(i, j) = method_use;
        double lo = na_nan();
        double hi = na_nan();
        double sei = na_nan();
        bootstrap_ci_pair(
          xv, yv, val, conf_level, method_use, bootstrap_reps, mi,
          tie_method, bias_correction, lo, hi, sei
        );
        lwr(i, j) = lo;
        upr(i, j) = hi;
        se(i, j) = sei;
      }
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
    out["ci_method"] = method_mat;
    out["se"] = se;
    out["m"] = m_mat;
    out["bootstrap_reps"] = bootstrap_reps;
  }
  return out;
}
