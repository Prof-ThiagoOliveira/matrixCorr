// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <vector>
#include "matrixCorr_detail.h"
#include "threshold_triplets.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using matrixCorr_detail::ranking::safe_inv_stddev;
using matrixCorr_detail::ranking::rank_vector;

namespace {

inline double clamp_corr(double x) {
  if (!std::isfinite(x)) return NA_REAL;
  if (x > 1.0) return 1.0;
  if (x < -1.0) return -1.0;
  return x;
}

inline double pearson_from_ranks(const arma::vec& rx, const arma::vec& ry) {
  const arma::uword n = rx.n_elem;
  if (n != ry.n_elem || n < 2u) return arma::datum::nan;

  const double* rx_ptr = rx.memptr();
  const double* ry_ptr = ry.memptr();

  double sum_xx = 0.0;
  double sum_yy = 0.0;
  double sum_xy = 0.0;
  for (arma::uword i = 0u; i < n; ++i) {
    const double vx = rx_ptr[i];
    const double vy = ry_ptr[i];
    sum_xx += vx * vx;
    sum_yy += vy * vy;
    sum_xy += vx * vy;
  }

  const double n_d = static_cast<double>(n);
  const double mean_rank = 0.5 * (n_d + 1.0);
  const double adj = n_d * mean_rank * mean_rank;
  const double sxx = sum_xx - adj;
  const double syy = sum_yy - adj;
  const double sxy = sum_xy - adj;
  if (!(sxx > 0.0) || !(syy > 0.0)) return arma::datum::nan;
  return clamp_corr(sxy / std::sqrt(sxx * syy));
}

inline double spearman_pair_core(const arma::vec& x,
                                 const arma::vec& y,
                                 arma::vec& rx,
                                 arma::vec& ry) {
  if (x.n_elem != y.n_elem || x.n_elem < 2u) return arma::datum::nan;
  rx.set_size(x.n_elem);
  ry.set_size(y.n_elem);
  rank_vector(x, rx);
  rank_vector(y, ry);
  return pearson_from_ranks(rx, ry);
}

inline double spearman_pair_core(const arma::vec& x, const arma::vec& y) {
  arma::vec rx, ry;
  return spearman_pair_core(x, y, rx, ry);
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

inline bool spearman_jel_ci_core(const arma::vec& x,
                                 const arma::vec& y,
                                 const double u,
                                 const double conf_level,
                                 double& lwr,
                                 double& upr) {
  const arma::uword n = x.n_elem;
  lwr = arma::datum::nan;
  upr = arma::datum::nan;
  if (!std::isfinite(u) || n < 3u) return false;

  arma::vec pseudo(n, arma::fill::none);
  arma::vec x_drop(n - 1u), y_drop(n - 1u), rx, ry;

  for (arma::uword i = 0; i < n; ++i) {
    arma::uword out = 0u;
    for (arma::uword j = 0; j < n; ++j) {
      if (j == i) continue;
      x_drop[out] = x[j];
      y_drop[out] = y[j];
      ++out;
    }
    const double u_drop = spearman_pair_core(x_drop, y_drop, rx, ry);
    if (!std::isfinite(u_drop)) return false;
    pseudo[i] = static_cast<double>(n) * u -
      static_cast<double>(n - 1u) * u_drop;
  }

  const double pseudo_var = arma::mean(arma::square(pseudo - u));
  if (!std::isfinite(pseudo_var)) return false;
  if (pseudo_var <= 1e-14) {
    lwr = u;
    upr = u;
    return true;
  }

  const double crit = R::qchisq(conf_level, 1.0, /*lower_tail*/ 1, /*log_p*/ 0);
  const auto g = [&](double theta) -> double {
    const arma::vec diff = pseudo - theta;
    const double s = arma::mean(arma::square(diff));
    const double num = static_cast<double>(n) * (u - theta) * (u - theta);
    if (!(s > 0.0)) {
      return (num <= 1e-18) ? -crit : std::numeric_limits<double>::infinity();
    }
    return num / s - crit;
  };

  if (u <= -1.0) {
    lwr = u;
  } else {
    const double g_lo = g(-1.0);
    const double g_u = g(u);
    if (!std::isfinite(g_lo) || g_lo <= 0.0) {
      lwr = -1.0;
    } else if (std::isfinite(g_u) && g_u <= 0.0) {
      lwr = bisect_root(g, -1.0, u, g_lo, g_u);
    } else {
      lwr = u;
    }
  }

  if (u >= 1.0) {
    upr = u;
  } else {
    const double g_hi = g(1.0);
    const double g_u = g(u);
    if (!std::isfinite(g_hi) || g_hi <= 0.0) {
      upr = 1.0;
    } else if (std::isfinite(g_u) && g_u <= 0.0) {
      upr = bisect_root(g, u, 1.0, g_u, g_hi);
    } else {
      upr = u;
    }
  }

  lwr = clamp_corr(lwr);
  upr = clamp_corr(upr);
  return std::isfinite(lwr) && std::isfinite(upr);
}

} // namespace

// [[Rcpp::export]]
arma::mat spearman_matrix_cpp(SEXP X_) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");

  const arma::uword n = static_cast<arma::uword>(Rf_nrows(X_));
  const arma::uword p = static_cast<arma::uword>(Rf_ncols(X_));
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  arma::mat R(n, p);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (p >= 4u && omp_get_max_threads() > 1)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    arma::vec out_col(R.colptr(uj), n, false, true);
    rank_vector(X.col(uj), out_col);
  }

  arma::mat RtR(p, p);
#if defined(ARMA_USE_BLAS)
  {
    RtR.zeros();
    const arma::blas_int N = static_cast<arma::blas_int>(p);
    const arma::blas_int K = static_cast<arma::blas_int>(n);
    const double alpha = 1.0, beta = 0.0;
    const char uplo  = 'U';
    const char trans = 'T';
    arma::blas::syrk<double>(&uplo, &trans, &N, &K,
                             &alpha, R.memptr(), &K,
                             &beta,  RtR.memptr(), &N);
    matrixCorr_detail::linalg::copy_upper_to_lower_inplace(RtR);
  }
#else
  RtR = R.t() * R;
#endif

  const double mu = 0.5 * (static_cast<double>(n) + 1.0);
  RtR -= static_cast<double>(n) * mu * mu;

  arma::vec s2 = RtR.diag() / static_cast<double>(n - 1);
  s2 = arma::clamp(s2, 0.0, std::numeric_limits<double>::infinity());
  arma::vec s  = arma::sqrt(s2);

  RtR /= static_cast<double>(n - 1);
  arma::vec inv_s = safe_inv_stddev(s);
  RtR.each_row() %= inv_s.t();
  RtR.each_col() %= inv_s;

  RtR.diag().ones();
  arma::uvec zero = arma::find(s == 0.0);
  for (arma::uword k = 0; k < zero.n_elem; ++k) {
    const arma::uword j = zero[k];
    RtR.row(j).fill(arma::datum::nan);
    RtR.col(j).fill(arma::datum::nan);
    RtR(j, j) = arma::datum::nan;
  }

  return RtR;
}

// [[Rcpp::export]]
Rcpp::List spearman_matrix_pairwise_cpp(SEXP X_,
                                        const bool return_ci = false,
                                        const double conf_level = 0.95) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");
  if (return_ci && !(conf_level > 0.0 && conf_level < 1.0))
    Rcpp::stop("conf_level must be in (0,1).");

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
#pragma omp parallel for schedule(dynamic) if (p >= 4u && omp_get_max_threads() > 1)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const arma::uvec& idx_j = finite_idx[uj];
    const arma::uword n_idx_j = idx_j.n_elem;
    const arma::uword* idx_j_ptr = idx_j.memptr();
    const double* colj_ptr = X.colptr(uj);

    static thread_local arma::vec xbuf;
    static thread_local arma::vec ybuf;
    static thread_local arma::vec rx;
    static thread_local arma::vec ry;

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
      const double rho = spearman_pair_core(x_view, y_view, rx, ry);
      est(uj, k) = rho;
      est(k, uj) = rho;

      if (return_ci && std::isfinite(rho)) {
        double lo = arma::datum::nan;
        double hi = arma::datum::nan;
        if (spearman_jel_ci_core(x_view, y_view, rho, conf_level, lo, hi)) {
          lwr(uj, k) = lo;
          lwr(k, uj) = lo;
          upr(uj, k) = hi;
          upr(k, uj) = hi;
        }
      }
    }
  }

  for (arma::uword j = 0; j < p; ++j) {
    const arma::uvec& idx = finite_idx[j];
    n_complete(j, j) = static_cast<int>(idx.n_elem);
    if (idx.n_elem < 2u) {
      est.row(j).fill(arma::datum::nan);
      est.col(j).fill(arma::datum::nan);
      if (return_ci) {
        lwr.row(j).fill(arma::datum::nan);
        lwr.col(j).fill(arma::datum::nan);
        upr.row(j).fill(arma::datum::nan);
        upr.col(j).fill(arma::datum::nan);
      }
      continue;
    }

    const double* col_ptr = X.colptr(j);
    double min_v = std::numeric_limits<double>::infinity();
    double max_v = -std::numeric_limits<double>::infinity();
    for (arma::uword t = 0; t < idx.n_elem; ++t) {
      const double v = col_ptr[idx[t]];
      if (v < min_v) min_v = v;
      if (v > max_v) max_v = v;
    }
    est(j, j) = (min_v < max_v) ? 1.0 : arma::datum::nan;
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["est"] = est,
    Rcpp::_["n_complete"] = n_complete
  );
  if (return_ci) {
    out["lwr"] = lwr;
    out["upr"] = upr;
    out["conf_level"] = conf_level;
  }
  return out;
}

// Complete-data Spearman upper-triplet kernel for thresholded outputs.
// [[Rcpp::export]]
Rcpp::List spearman_threshold_triplets_cpp(SEXP X_,
                                           const double threshold = 0.0,
                                           const bool diag = true,
                                           const int block_size = 256) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");
  if (!(threshold >= 0.0) || !std::isfinite(threshold))
    Rcpp::stop("threshold must be finite and >= 0.");
  if (block_size < 1)
    Rcpp::stop("block_size must be >= 1.");

  const arma::uword n = static_cast<arma::uword>(Rf_nrows(X_));
  const arma::uword p = static_cast<arma::uword>(Rf_ncols(X_));
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);
  const double mean_rank = 0.5 * (static_cast<double>(n) + 1.0);

  arma::mat Z(n, p, arma::fill::zeros);
  std::vector<unsigned char> valid(static_cast<std::size_t>(p), 0u);

  for (arma::uword j = 0; j < p; ++j) {
    arma::vec rj = rank_vector(X.col(j));
    rj -= mean_rank;
    const double ss = arma::dot(rj, rj);
    if (ss > 0.0 && std::isfinite(ss)) {
      Z.col(j) = rj / std::sqrt(ss);
      valid[static_cast<std::size_t>(j)] = 1u;
    }
  }

  const auto trip = matrixCorr_detail::threshold_triplets::collect_upper_triplets(
    static_cast<std::size_t>(p),
    static_cast<std::size_t>(block_size),
    diag,
    threshold,
    [&](std::size_t j0, std::size_t j1, std::size_t k0, std::size_t k1) -> arma::mat {
      arma::mat blk = Z.cols(static_cast<arma::uword>(j0), static_cast<arma::uword>(j1 - 1u)).t() *
        Z.cols(static_cast<arma::uword>(k0), static_cast<arma::uword>(k1 - 1u));

      for (std::size_t r = 0u; r < (j1 - j0); ++r) {
        const std::size_t gj = j0 + r;
        for (std::size_t c = 0u; c < (k1 - k0); ++c) {
          const std::size_t gk = k0 + c;
          double val = NA_REAL;
          if (valid[gj] && valid[gk]) {
            if (gj == gk) {
              val = 1.0;
            } else {
              val = clamp_corr(blk(
                static_cast<arma::uword>(r),
                static_cast<arma::uword>(c)
              ));
            }
          }
          blk(static_cast<arma::uword>(r), static_cast<arma::uword>(c)) = val;
        }
      }
      return blk;
    }
  );

  return matrixCorr_detail::threshold_triplets::as_list(trip);
}
