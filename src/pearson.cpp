// Thiago de Paula Oliveira
// pearson.cpp -- upper-triangular Pearson kernel with a single mirror pass
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <vector>
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

// X is n x p (double, no NA). Returns p x p Pearson correlation matrix.
// [[Rcpp::export]]
arma::mat pearson_matrix_cpp(SEXP X_) {
  if (!Rf_isReal(X_) || !Rf_isMatrix(X_))
    Rcpp::stop("Numeric double matrix required.");

  const arma::uword n = Rf_nrows(X_);
  const arma::uword p = Rf_ncols(X_);
  if (n < 2 || p < 2) Rcpp::stop("Need >= 2 rows and >= 2 columns.");

  // No-copy view over the input matrix.
  arma::mat X(REAL(X_), n, p, /*copy_aux_mem*/ false, /*strict*/ true);

  // Column means (1 x p).
  arma::rowvec mu = arma::sum(X, 0) / static_cast<double>(n);

  // XtX := X'X via SYRK (upper triangle); fallback to GEMM.
  arma::mat XtX(p, p, arma::fill::zeros);
#if defined(ARMA_USE_BLAS)
  {
    const arma::blas_int N = static_cast<arma::blas_int>(p);
    const arma::blas_int K = static_cast<arma::blas_int>(n);
    const double alpha = 1.0;
    const double beta = 0.0;
    const char uplo = 'U';
    const char trans = 'T';
    arma::blas::syrk<double>(&uplo, &trans, &N, &K,
                             &alpha, X.memptr(), &K,
                             &beta, XtX.memptr(), &N);
  }
#else
  XtX = X.t() * X;
#endif

  // XtX := XtX - n * mu * mu' on the stored upper triangle only.
  const double neg_n = -static_cast<double>(n);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double fj = neg_n * mu[uj];
    for (arma::uword i = 0; i <= uj; ++i) {
      XtX(i, uj) += fj * mu[i];
    }
  }

  // For correlation, the global covariance factor 1 / (n - 1) cancels.
  // Read the centred cross-product diagonal directly and normalise once.
  arma::vec centered_diag(p);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double d = XtX(uj, uj);
    centered_diag[uj] = d > 0.0 ? d : 0.0;
  }

  arma::vec inv_s = 1.0 / arma::sqrt(centered_diag);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
    const arma::uword uj = static_cast<arma::uword>(j);
    const double sj = inv_s[uj];
    if (!std::isfinite(sj) || sj == 0.0) continue;
    for (arma::uword i = 0; i < uj; ++i) {
      const double si = inv_s[i];
      if (!std::isfinite(si) || si == 0.0) continue;
      XtX(i, uj) *= (si * sj);
    }
    XtX(uj, uj) = 1.0;
  }

  XtX = arma::symmatu(XtX);

  // Zero-variance rows/cols are undefined.
  arma::uvec zero = arma::find(centered_diag == 0.0);
  for (arma::uword k = 0; k < zero.n_elem; ++k) {
    const arma::uword j = zero[k];
    XtX.row(j).fill(arma::datum::nan);
    XtX.col(j).fill(arma::datum::nan);
    XtX(j, j) = arma::datum::nan;
  }

  return XtX;
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
    const double* colj_ptr = X.colptr(uj);

    static thread_local std::vector<arma::uword> overlap_idx;

    for (arma::uword k = uj + 1u; k < p; ++k) {
      const arma::uvec& idx_k = finite_idx[k];
      const std::size_t possible = std::min(idx_j.n_elem, idx_k.n_elem);
      if (possible < 2u) continue;

      overlap_idx.clear();
      overlap_idx.reserve(possible);

      arma::uword ia = 0u, ib = 0u;
      while (ia < idx_j.n_elem && ib < idx_k.n_elem) {
        const arma::uword a = idx_j[ia];
        const arma::uword b = idx_k[ib];
        if (a == b) {
          overlap_idx.push_back(a);
          ++ia;
          ++ib;
        } else if (a < b) {
          ++ia;
        } else {
          ++ib;
        }
      }

      const int overlap_n = static_cast<int>(overlap_idx.size());
      n_complete(uj, k) = overlap_n;
      n_complete(k, uj) = overlap_n;
      if (overlap_n < 2) continue;

      const double* colk_ptr = X.colptr(k);
      double sum_x = 0.0;
      double sum_y = 0.0;
      for (std::size_t t = 0; t < overlap_idx.size(); ++t) {
        const arma::uword row = overlap_idx[t];
        sum_x += colj_ptr[row];
        sum_y += colk_ptr[row];
      }

      const double mean_x = sum_x / static_cast<double>(overlap_n);
      const double mean_y = sum_y / static_cast<double>(overlap_n);

      double sxx = 0.0;
      double syy = 0.0;
      double sxy = 0.0;
      for (std::size_t t = 0; t < overlap_idx.size(); ++t) {
        const arma::uword row = overlap_idx[t];
        const double dx = colj_ptr[row] - mean_x;
        const double dy = colk_ptr[row] - mean_y;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
      }

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
