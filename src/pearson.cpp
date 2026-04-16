// Thiago de Paula Oliveira
// pearson.cpp -- upper-triangular Pearson kernel with a single mirror pass
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <vector>
#include "threshold_triplets.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

namespace {

inline double clamp_corr(double x) {
  if (!std::isfinite(x)) return NA_REAL;
  if (x > 1.0) return 1.0;
  if (x < -1.0) return -1.0;
  return x;
}

inline double fisher_z(double r) {
  const double one_minus = std::nextafter(1.0, 0.0);
  if (r > one_minus) r = one_minus;
  if (r < -one_minus) r = -one_minus;
  return std::atanh(r);
}

inline double fisher_z_inv(double z) {
  double r = std::tanh(z);
  const double one_minus = std::nextafter(1.0, 0.0);
  if (r > one_minus) r = one_minus;
  if (r < -one_minus) r = -one_minus;
  return r;
}

} // namespace

struct PearsonScaleInfo {
  arma::rowvec mu;
  arma::vec inv_s;
  std::vector<unsigned char> valid;
};

inline PearsonScaleInfo compute_pearson_scale_info(const arma::mat& X) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const double n_d = static_cast<double>(n);

  PearsonScaleInfo out{
    arma::rowvec(p, arma::fill::zeros),
    arma::vec(p, arma::fill::zeros),
    std::vector<unsigned char>(static_cast<std::size_t>(p), 0u)
  };

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double* xj = X.colptr(uj);

    double sum = 0.0;
    double sumsq = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      const double v = xj[i];
      sum += v;
      sumsq += v * v;
    }

    const double mu = sum / n_d;
    const double d = sumsq - n_d * mu * mu;

    out.mu[uj] = mu;
    if (d > 0.0 && std::isfinite(d)) {
      out.inv_s[uj] = 1.0 / std::sqrt(d);
      out.valid[static_cast<std::size_t>(uj)] = 1u;
    }
  }

  return out;
}

// X is n x p (double, no NA). Returns p x p Pearson correlation matrix.
// [[Rcpp::export]]
arma::mat pearson_matrix_cpp(SEXP X_) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");

  const arma::uword n = Rf_nrows(X_);
  const arma::uword p = Rf_ncols(X_);
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);
  const double n_d = static_cast<double>(n);

  const PearsonScaleInfo scale = compute_pearson_scale_info(X);
  const double* mu_ptr = scale.mu.memptr();
  const double* inv_ptr = scale.inv_s.memptr();
  const std::vector<unsigned char>& valid = scale.valid;

  arma::mat R(p, p, arma::fill::zeros);

#if defined(ARMA_USE_BLAS)
{
  const arma::blas_int N = static_cast<arma::blas_int>(p);
  const arma::blas_int K = static_cast<arma::blas_int>(n);
  const double alpha = 1.0;
  const double beta = 0.0;
  const char uplo = 'U';
  const char trans = 'T';
  arma::blas::syrk<double>(
    &uplo, &trans, &N, &K,
    &alpha, X.memptr(), &K,
    &beta, R.memptr(), &N
  );
}
#else
R = X.t() * X;
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
  const arma::uword uj = static_cast<arma::uword>(j);
  const bool valid_j = valid[static_cast<std::size_t>(uj)] != 0u;
  const double mu_j = mu_ptr[uj];
  const double inv_j = inv_ptr[uj];

  for (arma::uword i = 0; i < uj; ++i) {
    double val = NA_REAL;

    if (valid_j && valid[static_cast<std::size_t>(i)] != 0u) {
      const double num = R(i, uj) - n_d * mu_ptr[i] * mu_j;
      val = clamp_corr(num * inv_ptr[i] * inv_j);
    }

    R(i, uj) = val;
    R(uj, i) = val;
  }

  R(uj, uj) = valid_j ? 1.0 : NA_REAL;
}

return R;
}

// Complete-data Pearson upper-triplet kernel for thresholded outputs.
// [[Rcpp::export]]
Rcpp::List pearson_threshold_triplets_cpp(SEXP X_,
                                          const double threshold = 0.0,
                                          const bool diag = true,
                                          const int block_size = 256) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");
  if (!(threshold >= 0.0) || !std::isfinite(threshold))
    Rcpp::stop("threshold must be finite and >= 0.");
  if (block_size < 1)
    Rcpp::stop("block_size must be >= 1.");

  const arma::uword n = Rf_nrows(X_);
  const arma::uword p = Rf_ncols(X_);
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);
  const double n_d = static_cast<double>(n);

  const PearsonScaleInfo scale = compute_pearson_scale_info(X);
  const double* mu_ptr = scale.mu.memptr();
  const double* inv_ptr = scale.inv_s.memptr();
  const std::vector<unsigned char>& valid = scale.valid;

  matrixCorr_detail::threshold_triplets::TripletBuffer out;

  const std::size_t pz = static_cast<std::size_t>(p);
  const std::size_t total_upper = diag
  ? (pz * (pz + 1u)) / 2u
  : (pz * (pz - 1u)) / 2u;

  if (threshold <= 0.0) {
    out.i.reserve(total_upper);
    out.j.reserve(total_upper);
    out.x.reserve(total_upper);
  } else {
    const std::size_t initial =
      std::min<std::size_t>(total_upper, std::max<std::size_t>(1024u, 8u * pz));
    out.i.reserve(initial);
    out.j.reserve(initial);
    out.x.reserve(initial);
  }

  const std::size_t bs = std::max<std::size_t>(
    1u, static_cast<std::size_t>(block_size)
  );

  for (std::size_t j0 = 0u; j0 < pz; j0 += bs) {
    const std::size_t j1 = std::min<std::size_t>(pz, j0 + bs);
    const std::size_t bj = j1 - j0;

    for (std::size_t k0 = j0; k0 < pz; k0 += bs) {
      const std::size_t k1 = std::min<std::size_t>(pz, k0 + bs);
      const std::size_t bk = k1 - k0;

      arma::mat blk =
        X.cols(static_cast<arma::uword>(j0), static_cast<arma::uword>(j1 - 1u)).t() *
        X.cols(static_cast<arma::uword>(k0), static_cast<arma::uword>(k1 - 1u));

      if (j0 == k0) {
        for (std::size_t r = 0u; r < bj; ++r) {
          const std::size_t gj = j0 + r;
          if (valid[gj] == 0u) continue;

          const double mu_j = mu_ptr[gj];
          const double inv_j = inv_ptr[gj];
          const std::size_t c_start = diag ? r : (r + 1u);

          for (std::size_t c = c_start; c < bk; ++c) {
            const std::size_t gk = k0 + c;
            if (valid[gk] == 0u) continue;

            const double val = (gj == gk)
              ? 1.0
            : clamp_corr(
                (blk(static_cast<arma::uword>(r), static_cast<arma::uword>(c)) -
                  n_d * mu_j * mu_ptr[gk]) *
                  inv_j * inv_ptr[gk]
            );

            if (!matrixCorr_detail::threshold_triplets::retain_value(val, threshold))
              continue;

            out.i.push_back(static_cast<int>(gj + 1u));
            out.j.push_back(static_cast<int>(gk + 1u));
            out.x.push_back(val);
          }
        }
      } else {
        for (std::size_t r = 0u; r < bj; ++r) {
          const std::size_t gj = j0 + r;
          if (valid[gj] == 0u) continue;

          const double mu_j = mu_ptr[gj];
          const double inv_j = inv_ptr[gj];

          for (std::size_t c = 0u; c < bk; ++c) {
            const std::size_t gk = k0 + c;
            if (valid[gk] == 0u) continue;

            const double val = clamp_corr(
              (blk(static_cast<arma::uword>(r), static_cast<arma::uword>(c)) -
                n_d * mu_j * mu_ptr[gk]) *
                inv_j * inv_ptr[gk]
            );

            if (!matrixCorr_detail::threshold_triplets::retain_value(val, threshold))
              continue;

            out.i.push_back(static_cast<int>(gj + 1u));
            out.j.push_back(static_cast<int>(gk + 1u));
            out.x.push_back(val);
          }
        }
      }
    }
  }

  return matrixCorr_detail::threshold_triplets::as_list(out);
}

// Pairwise-complete Pearson matrix with optional Fisher-z confidence intervals.
// [[Rcpp::export]]
Rcpp::List pearson_matrix_pairwise_cpp(SEXP X_,
                                       const bool return_ci = false,
                                       const double conf_level = 0.95) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");
  if (return_ci && !(conf_level > 0.0 && conf_level < 1.0))
    Rcpp::stop("conf_level must be in (0,1).");

  const arma::uword n = Rf_nrows(X_);
  const arma::uword p = Rf_ncols(X_);
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  arma::mat R(p, p, arma::fill::none);
  R.fill(arma::datum::nan);

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

  const double zcrit = return_ci
    ? R::qnorm5(0.5 * (1.0 + conf_level), 0.0, 1.0, /*lower_tail*/ 1, /*log_p*/ 0)
    : NA_REAL;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const arma::uvec& idx_j = finite_idx[uj];
    const arma::uword n_idx_j = idx_j.n_elem;
    const arma::uword* idx_j_ptr = idx_j.memptr();
    const double* colj_ptr = X.colptr(uj);

    for (arma::uword k = uj + 1u; k < p; ++k) {
      const arma::uvec& idx_k = finite_idx[k];
      const arma::uword n_idx_k = idx_k.n_elem;
      const std::size_t possible = std::min(n_idx_j, n_idx_k);
      if (possible < 2u) continue;

      const arma::uword* idx_k_ptr = idx_k.memptr();
      const double* colk_ptr = X.colptr(k);
      arma::uword ia = 0u, ib = 0u;
      int overlap_n = 0;
      double mean_x = 0.0;
      double mean_y = 0.0;
      double sxx = 0.0;
      double syy = 0.0;
      double sxy = 0.0;
      while (ia < n_idx_j && ib < n_idx_k) {
        const arma::uword a = idx_j_ptr[ia];
        const arma::uword b = idx_k_ptr[ib];
        if (a == b) {
          const double x = colj_ptr[a];
          const double y = colk_ptr[b];
          ++overlap_n;
          const double inv_n = 1.0 / static_cast<double>(overlap_n);
          const double dx = x - mean_x;
          const double dy = y - mean_y;
          mean_x += dx * inv_n;
          mean_y += dy * inv_n;
          sxx += dx * (x - mean_x);
          syy += dy * (y - mean_y);
          sxy += dx * (y - mean_y);
          ++ia;
          ++ib;
        } else if (a < b) {
          ++ia;
        } else {
          ++ib;
        }
      }

      n_complete(uj, k) = overlap_n;
      n_complete(k, uj) = overlap_n;
      if (overlap_n < 2) continue;

      double r = arma::datum::nan;
      if (sxx > 0.0 && syy > 0.0) {
        r = clamp_corr(sxy / std::sqrt(sxx * syy));
      }
      R(uj, k) = r;
      R(k, uj) = r;

      if (return_ci && std::isfinite(r) && overlap_n > 3) {
        const double zr = fisher_z(r);
        const double se = 1.0 / std::sqrt(static_cast<double>(overlap_n) - 3.0);
        lwr(uj, k) = fisher_z_inv(zr - zcrit * se);
        upr(uj, k) = fisher_z_inv(zr + zcrit * se);
        lwr(k, uj) = lwr(uj, k);
        upr(k, uj) = upr(uj, k);
      }
    }
  }

  for (arma::uword j = 0; j < p; ++j) {
    const arma::uvec& idx = finite_idx[j];
    n_complete(j, j) = static_cast<int>(idx.n_elem);

    if (idx.n_elem < 2u) {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
      if (return_ci) {
        lwr.row(j).fill(arma::datum::nan);
        lwr.col(j).fill(arma::datum::nan);
        upr.row(j).fill(arma::datum::nan);
        upr.col(j).fill(arma::datum::nan);
      }
      continue;
    }

    const double* col_ptr = X.colptr(j);
    double sum = 0.0;
    for (arma::uword t = 0; t < idx.n_elem; ++t) {
      sum += col_ptr[idx[t]];
    }
    const double mean = sum / static_cast<double>(idx.n_elem);

    double ss = 0.0;
    for (arma::uword t = 0; t < idx.n_elem; ++t) {
      const double d = col_ptr[idx[t]] - mean;
      ss += d * d;
    }

    if (ss > 0.0) {
      R(j, j) = 1.0;
    } else {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["est"] = R,
    Rcpp::_["n_complete"] = n_complete
  );
  if (return_ci) {
    out["lwr"] = lwr;
    out["upr"] = upr;
    out["conf_level"] = conf_level;
  }
  return out;
}
