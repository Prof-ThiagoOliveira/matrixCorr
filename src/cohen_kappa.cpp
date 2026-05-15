// Thiago de Paula Oliveira
// src-cohen_kappa.cpp
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>
#include "matrixCorr_omp.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

namespace {

struct KappaStats {
  double est = NA_REAL;
  int n_complete = 0;
  double po = NA_REAL;
  double pe = NA_REAL;
  double se = NA_REAL;
  double lwr = NA_REAL;
  double upr = NA_REAL;
  double statistic = NA_REAL;
  double p_value = NA_REAL;
};

inline bool valid_level(const int code, const int n_levels) {
  return code >= 1 && code <= n_levels;
}

KappaStats compute_kappa_from_counts(const std::vector<int>& counts,
                                     const int n_levels,
                                     const bool return_inference,
                                     const double conf_level) {
  KappaStats out;
  long long n_complete = 0;
  std::vector<double> q(static_cast<std::size_t>(std::max(0, n_levels)), 0.0);
  std::vector<double> r(static_cast<std::size_t>(std::max(0, n_levels)), 0.0);
  double agree = 0.0;

  for (int a = 0; a < n_levels; ++a) {
    for (int b = 0; b < n_levels; ++b) {
      const double cij = static_cast<double>(counts[static_cast<std::size_t>(a * n_levels + b)]);
      n_complete += static_cast<long long>(counts[static_cast<std::size_t>(a * n_levels + b)]);
      q[static_cast<std::size_t>(a)] += cij;
      r[static_cast<std::size_t>(b)] += cij;
      if (a == b) {
        agree += cij;
      }
    }
  }

  out.n_complete = static_cast<int>(n_complete);
  if (n_complete < 2 || n_levels <= 0) {
    return out;
  }

  const double n = static_cast<double>(n_complete);
  for (int a = 0; a < n_levels; ++a) {
    q[static_cast<std::size_t>(a)] /= n;
    r[static_cast<std::size_t>(a)] /= n;
  }

  const double A = agree / n;
  double E = 0.0;
  for (int a = 0; a < n_levels; ++a) {
    E += q[static_cast<std::size_t>(a)] * r[static_cast<std::size_t>(a)];
  }

  const double D = 1.0 - E;
  out.po = A;
  out.pe = E;
  if (!std::isfinite(D) || D <= 1e-15) {
    return out;
  }

  const double kappa = (A - E) / D;
  if (!std::isfinite(kappa)) {
    return out;
  }
  out.est = kappa;

  if (!return_inference) {
    return out;
  }

  double mean_g = 0.0;
  double mean_g2 = 0.0;
  const double D2 = D * D;
  for (int a = 0; a < n_levels; ++a) {
    for (int b = 0; b < n_levels; ++b) {
      const double pab = static_cast<double>(counts[static_cast<std::size_t>(a * n_levels + b)]) / n;
      if (pab <= 0.0) {
        continue;
      }
      const double qab = r[static_cast<std::size_t>(a)] + q[static_cast<std::size_t>(b)];
      const double gab = ((a == b) ? D : 0.0) + (A - 1.0) * qab;
      const double g = gab / D2;
      mean_g += pab * g;
      mean_g2 += pab * g * g;
    }
  }

  const double var = (mean_g2 - mean_g * mean_g) / n;
  if (!std::isfinite(var) || var <= 0.0) {
    return out;
  }

  const double se = std::sqrt(std::max(var, 0.0));
  if (!std::isfinite(se) || se <= 0.0) {
    return out;
  }

  out.se = se;
  out.statistic = kappa / se;
  out.p_value = 2.0 * R::pnorm5(std::fabs(out.statistic), 0.0, 1.0, 0, 0);
  const double z = R::qnorm5(0.5 * (1.0 + conf_level), 0.0, 1.0, 1, 0);
  out.lwr = std::max(-1.0, kappa - z * se);
  out.upr = std::min(1.0, kappa + z * se);
  return out;
}

KappaStats compute_kappa_pair(const Rcpp::IntegerVector& x,
                              const Rcpp::IntegerVector& y,
                              const int n_levels,
                              const bool return_inference,
                              const double conf_level) {
  std::vector<int> counts(static_cast<std::size_t>(std::max(0, n_levels * n_levels)), 0);
  const R_xlen_t n = x.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    const int xi = x[i];
    const int yi = y[i];
    if (xi == NA_INTEGER || yi == NA_INTEGER) {
      continue;
    }
    if (!valid_level(xi, n_levels) || !valid_level(yi, n_levels)) {
      continue;
    }
    counts[static_cast<std::size_t>((xi - 1) * n_levels + (yi - 1))] += 1;
  }
  return compute_kappa_from_counts(counts, n_levels, return_inference, conf_level);
}

KappaStats compute_kappa_pair_matrix(const Rcpp::IntegerMatrix& X,
                                     const int col_j,
                                     const int col_k,
                                     const int n_levels,
                                     const bool return_inference,
                                     const double conf_level) {
  std::vector<int> counts(static_cast<std::size_t>(std::max(0, n_levels * n_levels)), 0);
  const int n = X.nrow();
  for (int i = 0; i < n; ++i) {
    const int xij = X(i, col_j);
    const int xik = X(i, col_k);
    if (xij == NA_INTEGER || xik == NA_INTEGER) {
      continue;
    }
    if (!valid_level(xij, n_levels) || !valid_level(xik, n_levels)) {
      continue;
    }
    counts[static_cast<std::size_t>((xij - 1) * n_levels + (xik - 1))] += 1;
  }
  return compute_kappa_from_counts(counts, n_levels, return_inference, conf_level);
}

KappaStats compute_kappa_diag(const Rcpp::IntegerMatrix& X,
                              const int col_j,
                              const int n_levels,
                              const bool return_inference,
                              const double conf_level) {
  KappaStats out;
  std::vector<int> freq(static_cast<std::size_t>(std::max(0, n_levels)), 0);
  const int n = X.nrow();
  for (int i = 0; i < n; ++i) {
    const int xij = X(i, col_j);
    if (xij == NA_INTEGER || !valid_level(xij, n_levels)) {
      continue;
    }
    freq[static_cast<std::size_t>(xij - 1)] += 1;
    out.n_complete += 1;
  }

  if (out.n_complete < 2) {
    return out;
  }

  out.est = 1.0;
  out.po = 1.0;
  const double n_complete = static_cast<double>(out.n_complete);
  double pe = 0.0;
  for (int a = 0; a < n_levels; ++a) {
    const double pa = static_cast<double>(freq[static_cast<std::size_t>(a)]) / n_complete;
    pe += pa * pa;
  }
  out.pe = pe;

  if (return_inference) {
    out.lwr = 1.0;
    out.upr = 1.0;
    out.se = NA_REAL;
    out.statistic = NA_REAL;
    out.p_value = NA_REAL;
    if (!std::isfinite(conf_level)) {
      out.lwr = NA_REAL;
      out.upr = NA_REAL;
    }
  }
  return out;
}

inline void store_sym_double(Rcpp::NumericMatrix& M,
                             const int j,
                             const int k,
                             const double value) {
  M(j, k) = value;
  M(k, j) = value;
}

inline void store_sym_int(Rcpp::IntegerMatrix& M,
                          const int j,
                          const int k,
                          const int value) {
  M(j, k) = value;
  M(k, j) = value;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List cohen_kappa_pair_cpp(Rcpp::IntegerVector x,
                                Rcpp::IntegerVector y,
                                int n_levels,
                                bool return_inference = false,
                                double conf_level = 0.95) {
  KappaStats stats = compute_kappa_pair(x, y, std::max(0, n_levels), return_inference, conf_level);
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["est"] = stats.est,
    Rcpp::_["n_complete"] = stats.n_complete,
    Rcpp::_["po"] = stats.po,
    Rcpp::_["pe"] = stats.pe
  );
  if (return_inference) {
    out["se"] = stats.se;
    out["lwr"] = stats.lwr;
    out["upr"] = stats.upr;
    out["statistic"] = stats.statistic;
    out["p_value"] = stats.p_value;
    out["conf_level"] = conf_level;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List cohen_kappa_matrix_cpp(Rcpp::IntegerMatrix X,
                                  int n_levels,
                                  bool pairwise_complete = false,
                                  bool return_inference = false,
                                  double conf_level = 0.95,
                                  int n_threads = 1) {
  const int p = X.ncol();
  Rcpp::NumericMatrix est(p, p);
  Rcpp::IntegerMatrix n_complete(p, p);
  Rcpp::NumericMatrix po(p, p);
  Rcpp::NumericMatrix pe(p, p);
  std::fill(est.begin(), est.end(), NA_REAL);
  std::fill(n_complete.begin(), n_complete.end(), NA_INTEGER);
  std::fill(po.begin(), po.end(), NA_REAL);
  std::fill(pe.begin(), pe.end(), NA_REAL);

  Rcpp::NumericMatrix se;
  Rcpp::NumericMatrix lwr;
  Rcpp::NumericMatrix upr;
  Rcpp::NumericMatrix statistic;
  Rcpp::NumericMatrix p_value;
  if (return_inference) {
    se = Rcpp::NumericMatrix(p, p);
    lwr = Rcpp::NumericMatrix(p, p);
    upr = Rcpp::NumericMatrix(p, p);
    statistic = Rcpp::NumericMatrix(p, p);
    p_value = Rcpp::NumericMatrix(p, p);
    std::fill(se.begin(), se.end(), NA_REAL);
    std::fill(lwr.begin(), lwr.end(), NA_REAL);
    std::fill(upr.begin(), upr.end(), NA_REAL);
    std::fill(statistic.begin(), statistic.end(), NA_REAL);
    std::fill(p_value.begin(), p_value.end(), NA_REAL);
  }

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    KappaStats diag_stats = compute_kappa_diag(X, j, std::max(0, n_levels), return_inference, conf_level);
    est(j, j) = diag_stats.est;
    n_complete(j, j) = diag_stats.n_complete;
    po(j, j) = diag_stats.po;
    pe(j, j) = diag_stats.pe;
    if (return_inference) {
      se(j, j) = diag_stats.se;
      lwr(j, j) = diag_stats.lwr;
      upr(j, j) = diag_stats.upr;
      statistic(j, j) = diag_stats.statistic;
      p_value(j, j) = diag_stats.p_value;
    }

    for (int k = j + 1; k < p; ++k) {
      KappaStats stats = compute_kappa_pair_matrix(
        X,
        j,
        k,
        std::max(0, n_levels),
        return_inference,
        conf_level
      );
      store_sym_double(est, j, k, stats.est);
      store_sym_int(n_complete, j, k, stats.n_complete);
      store_sym_double(po, j, k, stats.po);
      store_sym_double(pe, j, k, stats.pe);
      if (return_inference) {
        store_sym_double(se, j, k, stats.se);
        store_sym_double(lwr, j, k, stats.lwr);
        store_sym_double(upr, j, k, stats.upr);
        store_sym_double(statistic, j, k, stats.statistic);
        store_sym_double(p_value, j, k, stats.p_value);
      }
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["est"] = est,
    Rcpp::_["n_complete"] = n_complete,
    Rcpp::_["po"] = po,
    Rcpp::_["pe"] = pe
  );
  if (return_inference) {
    out["se"] = se;
    out["lwr"] = lwr;
    out["upr"] = upr;
    out["statistic"] = statistic;
    out["p_value"] = p_value;
    out["conf_level"] = conf_level;
  }
  out["pairwise_complete"] = pairwise_complete;
  return out;
}

// [[Rcpp::export]]
Rcpp::List cohen_kappa_threshold_triplets_cpp(Rcpp::IntegerMatrix X,
                                              int n_levels,
                                              double threshold = 0.0,
                                              bool diag = true,
                                              int block_size = 256,
                                              int n_threads = 1) {
  if (!(threshold >= 0.0) || !std::isfinite(threshold)) {
    Rcpp::stop("threshold must be finite and >= 0.");
  }
  if (block_size < 1) {
    Rcpp::stop("block_size must be >= 1.");
  }

  const int p = X.ncol();
  const int n = X.nrow();
  if (p < 1 || n < 1) {
    return Rcpp::List::create(
      Rcpp::_["i"] = Rcpp::IntegerVector(),
      Rcpp::_["j"] = Rcpp::IntegerVector(),
      Rcpp::_["x"] = Rcpp::NumericVector()
    );
  }

  const int max_threads =
#ifdef _OPENMP
    std::max(1, n_threads);
#else
    1;
#endif
  std::vector< std::vector<int> > i_buf(static_cast<std::size_t>(max_threads));
  std::vector< std::vector<int> > j_buf(static_cast<std::size_t>(max_threads));
  std::vector< std::vector<double> > x_buf(static_cast<std::size_t>(max_threads));

#ifdef _OPENMP
  omp_set_num_threads(max_threads);
#pragma omp parallel
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    std::vector<int>& i_local = i_buf[static_cast<std::size_t>(tid)];
    std::vector<int>& j_local = j_buf[static_cast<std::size_t>(tid)];
    std::vector<double>& x_local = x_buf[static_cast<std::size_t>(tid)];

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int j = 0; j < p; ++j) {
      if (diag) {
        KappaStats d = compute_kappa_diag(X, j, std::max(0, n_levels), false, 0.95);
        if (std::isfinite(d.est) && std::fabs(d.est) >= threshold) {
          i_local.push_back(j + 1);
          j_local.push_back(j + 1);
          x_local.push_back(d.est);
        }
      }
      for (int k = j + 1; k < p; ++k) {
        KappaStats stats = compute_kappa_pair_matrix(X, j, k, std::max(0, n_levels), false, 0.95);
        if (std::isfinite(stats.est) && std::fabs(stats.est) >= threshold) {
          i_local.push_back(j + 1);
          j_local.push_back(k + 1);
          x_local.push_back(stats.est);
        }
      }
    }
  }

  std::size_t total = 0u;
  for (int t = 0; t < max_threads; ++t) {
    total += i_buf[static_cast<std::size_t>(t)].size();
  }

  std::vector<int> i_out;
  std::vector<int> j_out;
  std::vector<double> x_out;
  i_out.reserve(total);
  j_out.reserve(total);
  x_out.reserve(total);
  for (int t = 0; t < max_threads; ++t) {
    i_out.insert(
      i_out.end(),
      i_buf[static_cast<std::size_t>(t)].begin(),
      i_buf[static_cast<std::size_t>(t)].end()
    );
    j_out.insert(
      j_out.end(),
      j_buf[static_cast<std::size_t>(t)].begin(),
      j_buf[static_cast<std::size_t>(t)].end()
    );
    x_out.insert(
      x_out.end(),
      x_buf[static_cast<std::size_t>(t)].begin(),
      x_buf[static_cast<std::size_t>(t)].end()
    );
  }

  return Rcpp::List::create(
    Rcpp::_["i"] = Rcpp::wrap(i_out),
    Rcpp::_["j"] = Rcpp::wrap(j_out),
    Rcpp::_["x"] = Rcpp::wrap(x_out)
  );
}
