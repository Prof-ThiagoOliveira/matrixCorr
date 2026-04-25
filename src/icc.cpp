// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

constexpr int ICC_FORM_1 = 0;
constexpr int ICC_FORM_2 = 1;
constexpr int ICC_FORM_3 = 2;

inline bool is_positive_finite(const double x) {
  return std::isfinite(x) && x > 0.0;
}

struct PairMoments {
  int n = 0;
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_x2 = 0.0;
  double sum_y2 = 0.0;
  double sum_xy = 0.0;
};

inline bool compute_overlap_moments(const arma::uvec& idx_j,
                                    const arma::uvec& idx_k,
                                    const double* colj_ptr,
                                    const double* colk_ptr,
                                    const int min_n,
                                    PairMoments& out) {
  const std::size_t possible = std::min(idx_j.n_elem, idx_k.n_elem);
  if (possible < static_cast<std::size_t>(min_n)) return false;

  out = PairMoments{};

  arma::uword ia = 0u, ib = 0u;
  while (ia < idx_j.n_elem && ib < idx_k.n_elem) {
    const arma::uword a = idx_j[ia];
    const arma::uword b = idx_k[ib];
    if (a == b) {
      const double x = colj_ptr[a];
      const double y = colk_ptr[b];
      out.n += 1;
      out.sum_x += x;
      out.sum_y += y;
      out.sum_x2 += x * x;
      out.sum_y2 += y * y;
      out.sum_xy += x * y;
      ++ia;
      ++ib;
    } else if (a < b) {
      ++ia;
    } else {
      ++ib;
    }
  }

  return out.n >= min_n;
}

inline bool icc_pair_from_moments(const PairMoments& moments,
                                  const int form_code,
                                  const bool average_unit,
                                  const bool return_ci,
                                  const double conf_level,
                                  double& estimate,
                                  double& lwr,
                                  double& upr) {
  estimate = NA_REAL;
  lwr = NA_REAL;
  upr = NA_REAL;

  const int n = moments.n;
  if (n < 2) return false;

  const double n_d = static_cast<double>(n);
  const double mean_x = moments.sum_x / n_d;
  const double mean_y = moments.sum_y / n_d;

  double ss_x = moments.sum_x2 - (moments.sum_x * moments.sum_x) / n_d;
  double ss_y = moments.sum_y2 - (moments.sum_y * moments.sum_y) / n_d;
  const double cov_num = moments.sum_xy - (moments.sum_x * moments.sum_y) / n_d;

  const double tol_x = 1e-12 * std::max(1.0, moments.sum_x2);
  const double tol_y = 1e-12 * std::max(1.0, moments.sum_y2);
  if (ss_x < 0.0 && std::abs(ss_x) <= tol_x) ss_x = 0.0;
  if (ss_y < 0.0 && std::abs(ss_y) <= tol_y) ss_y = 0.0;

  const double ssr = 0.5 * (ss_x + ss_y + 2.0 * cov_num);
  const double sse = 0.5 * (ss_x + ss_y - 2.0 * cov_num);
  const double mean_d = mean_x - mean_y;
  const double ssc = 0.5 * n_d * mean_d * mean_d;

  if (!is_positive_finite(ss_x) || !is_positive_finite(ss_y)) return false;

  const double df_subject = static_cast<double>(n - 1);
  const double df_error = static_cast<double>(n - 1);
  const double df_within = static_cast<double>(n);
  if (!(df_subject > 0.0) || !(df_error > 0.0) || !(df_within > 0.0)) return false;

  const double msr = ssr / df_subject;
  const double msc = ssc;
  const double mse = sse / df_error;
  const double msw = (ssc + sse) / df_within;
  const double k = 2.0;

  double single_est = NA_REAL;
  if (form_code == ICC_FORM_1) {
    const double denom = msr + (k - 1.0) * msw;
    if (!is_positive_finite(denom)) return false;
    single_est = (msr - msw) / denom;
  } else if (form_code == ICC_FORM_2) {
    const double denom = msr + (k - 1.0) * mse + k * (msc - mse) / n_d;
    if (!is_positive_finite(denom)) return false;
    single_est = (msr - mse) / denom;
  } else {
    const double denom = msr + (k - 1.0) * mse;
    if (!is_positive_finite(denom)) return false;
    single_est = (msr - mse) / denom;
  }

  estimate = average_unit
    ? (k * single_est) / (1.0 + (k - 1.0) * single_est)
    : single_est;

  if (!return_ci) return true;

  const double alpha = 1.0 - conf_level;
  if (!(alpha > 0.0 && alpha < 1.0)) return true;

  if (form_code == ICC_FORM_1) {
    if (!is_positive_finite(msw)) return true;
    const double F = msr / msw;
    if (!std::isfinite(F) || F < 0.0) return true;
    const double F_lo = F / R::qf(1.0 - alpha / 2.0, df_subject, static_cast<double>(n), 1, 0);
    const double F_hi = F * R::qf(1.0 - alpha / 2.0, static_cast<double>(n), df_subject, 1, 0);
    if (!std::isfinite(F_lo) || !std::isfinite(F_hi)) return true;
    lwr = (F_lo - 1.0) / (F_lo + (k - 1.0));
    upr = (F_hi - 1.0) / (F_hi + (k - 1.0));
  } else if (form_code == ICC_FORM_3) {
    if (!is_positive_finite(mse)) return true;
    const double F = msr / mse;
    if (!std::isfinite(F) || F < 0.0) return true;
    const double F_lo = F / R::qf(1.0 - alpha / 2.0, df_subject, df_error, 1, 0);
    const double F_hi = F * R::qf(1.0 - alpha / 2.0, df_error, df_subject, 1, 0);
    if (!std::isfinite(F_lo) || !std::isfinite(F_hi)) return true;
    lwr = (F_lo - 1.0) / (F_lo + (k - 1.0));
    upr = (F_hi - 1.0) / (F_hi + (k - 1.0));
  } else {
    if (!is_positive_finite(mse)) return true;
    const double Fj = msc / mse;
    if (!std::isfinite(Fj) || Fj < 0.0) return true;

    const double vn = (k - 1.0) * df_subject *
      std::pow(k * single_est * Fj + n_d * (1.0 + (k - 1.0) * single_est) - k * single_est, 2.0);
    const double vd =
      df_subject * k * k * single_est * single_est * Fj * Fj +
      std::pow(n_d * (1.0 + (k - 1.0) * single_est) - k * single_est, 2.0);
    if (!is_positive_finite(vn) || !is_positive_finite(vd)) return true;

    const double v = vn / vd;
    if (!is_positive_finite(v)) return true;

    const double F_hi = R::qf(1.0 - alpha / 2.0, df_subject, v, 1, 0);
    const double F_lo = R::qf(1.0 - alpha / 2.0, v, df_subject, 1, 0);
    if (!std::isfinite(F_lo) || !std::isfinite(F_hi)) return true;

    const double denom_lo = F_hi * (k * msc + (k * n_d - k - n_d) * mse) + n_d * msr;
    const double denom_hi = k * msc + (k * n_d - k - n_d) * mse + n_d * F_lo * msr;
    if (!is_positive_finite(denom_lo) || !is_positive_finite(denom_hi)) return true;

    lwr = n_d * (msr - F_hi * mse) / denom_lo;
    upr = n_d * (F_lo * msr - mse) / denom_hi;
  }

  if (average_unit && std::isfinite(lwr) && std::isfinite(upr)) {
    lwr = (k * lwr) / (1.0 + (k - 1.0) * lwr);
    upr = (k * upr) / (1.0 + (k - 1.0) * upr);
  }

  if (std::isfinite(lwr) && std::isfinite(upr) && lwr > upr) std::swap(lwr, upr);
  return true;
}

inline bool icc_pair_complete_core(const double* x,
                                   const double* y,
                                   const int n,
                                   const int form_code,
                                   const bool average_unit,
                                   const bool return_ci,
                                   const double conf_level,
                                   double& estimate,
                                   double& lwr,
                                   double& upr) {
  PairMoments moments;
  moments.n = n;
  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    const double yi = y[i];
    moments.sum_x += xi;
    moments.sum_y += yi;
    moments.sum_x2 += xi * xi;
    moments.sum_y2 += yi * yi;
    moments.sum_xy += xi * yi;
  }
  return icc_pair_from_moments(
    moments,
    form_code,
    average_unit,
    return_ci,
    conf_level,
    estimate,
    lwr,
    upr
  );
}

} // namespace

// [[Rcpp::export]]
Rcpp::List icc_matrix_cpp(const arma::mat& X,
                          const int form_code = 0,
                          const bool average_unit = false,
                          const bool pairwise_complete = false,
                          const bool return_ci = false,
                          const double conf_level = 0.95,
                          const int n_threads = 1) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (n < 2u || p < 2u) Rcpp::stop("Need >= 2 rows and >= 2 columns.");
  if (return_ci && !(conf_level > 0.0 && conf_level < 1.0)) {
    Rcpp::stop("conf_level must be in (0,1).");
  }

  arma::mat est(p, p, arma::fill::eye);
  arma::Mat<int> n_complete(p, p, arma::fill::zeros);

  arma::mat lwr;
  arma::mat upr;
  if (return_ci) {
    lwr.set_size(p, p);
    upr.set_size(p, p);
    lwr.fill(arma::datum::nan);
    upr.fill(arma::datum::nan);
  }

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#endif

  if (!pairwise_complete) {
    if (!X.is_finite()) {
      Rcpp::stop("X contains NA/NaN/Inf; please handle missingness upstream.");
    }

    const double n_d = static_cast<double>(n);
    const arma::vec col_means = arma::mean(X, 0).t();
    const arma::mat cross = X.t() * X;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
      const arma::uword uj = static_cast<arma::uword>(j);
      for (arma::uword k = uj + 1u; k < p; ++k) {
        PairMoments moments;
        moments.n = static_cast<int>(n);
        moments.sum_x = col_means[uj] * n_d;
        moments.sum_y = col_means[k] * n_d;
        moments.sum_x2 = cross(uj, uj);
        moments.sum_y2 = cross(k, k);
        moments.sum_xy = cross(uj, k);

        double estimate = NA_REAL, ci_l = NA_REAL, ci_u = NA_REAL;
        const bool ok = icc_pair_from_moments(
          moments,
          form_code,
          average_unit,
          return_ci,
          conf_level,
          estimate,
          ci_l,
          ci_u
        );
        n_complete(uj, k) = static_cast<int>(n);
        n_complete(k, uj) = static_cast<int>(n);
        est(uj, k) = ok ? estimate : NA_REAL;
        est(k, uj) = est(uj, k);
        if (return_ci) {
          lwr(uj, k) = ok ? ci_l : NA_REAL;
          lwr(k, uj) = lwr(uj, k);
          upr(uj, k) = ok ? ci_u : NA_REAL;
          upr(k, uj) = upr(uj, k);
        }
      }
    }

    for (arma::uword j = 0u; j < p; ++j) n_complete(j, j) = static_cast<int>(n);
  } else {
    std::vector<arma::uvec> finite_idx(p);
    for (arma::uword j = 0u; j < p; ++j) {
      finite_idx[j] = arma::find_finite(X.col(j));
      n_complete(j, j) = static_cast<int>(finite_idx[j].n_elem);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (arma::sword j = 0; j < static_cast<arma::sword>(p); ++j) {
      const arma::uword uj = static_cast<arma::uword>(j);
      const arma::uvec& idx_j = finite_idx[uj];
      const double* colj_ptr = X.colptr(uj);

      for (arma::uword k = uj + 1u; k < p; ++k) {
        const arma::uvec& idx_k = finite_idx[k];
        PairMoments moments;
        if (!compute_overlap_moments(idx_j, idx_k, colj_ptr, X.colptr(k), 2, moments)) {
          n_complete(uj, k) = 0;
          n_complete(k, uj) = 0;
          continue;
        }

        double estimate = NA_REAL, ci_l = NA_REAL, ci_u = NA_REAL;
        const bool ok = icc_pair_from_moments(
          moments,
          form_code,
          average_unit,
          return_ci,
          conf_level,
          estimate,
          ci_l,
          ci_u
        );
        n_complete(uj, k) = moments.n;
        n_complete(k, uj) = moments.n;
        est(uj, k) = ok ? estimate : NA_REAL;
        est(k, uj) = est(uj, k);
        if (return_ci) {
          lwr(uj, k) = ok ? ci_l : NA_REAL;
          lwr(k, uj) = lwr(uj, k);
          upr(uj, k) = ok ? ci_u : NA_REAL;
          upr(k, uj) = upr(uj, k);
        }
      }
    }
  }

  est.diag().fill(1.0);
  if (return_ci) {
    lwr.diag().fill(1.0);
    upr.diag().fill(1.0);
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["est"] = est,
    Rcpp::_["n_complete"] = n_complete
  );
  if (return_ci) {
    out["lwr.ci"] = lwr;
    out["upr.ci"] = upr;
    out["conf_level"] = conf_level;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List icc_overall_cpp(const arma::mat& X,
                           const bool return_ci = false,
                           const double conf_level = 0.95) {
  const arma::uword n = X.n_rows;
  const arma::uword k = X.n_cols;
  if (n < 2u || k < 2u) {
    Rcpp::stop("Need >= 2 rows and >= 2 columns.");
  }
  if (!X.is_finite()) {
    Rcpp::stop("X contains NA/NaN/Inf; please handle missingness upstream.");
  }
  if (return_ci && !(conf_level > 0.0 && conf_level < 1.0)) {
    Rcpp::stop("conf_level must be in (0,1).");
  }

  const double n_d = static_cast<double>(n);
  const double k_d = static_cast<double>(k);
  const double grand_mean = arma::mean(arma::vectorise(X));
  const arma::rowvec col_means = arma::mean(X, 0);
  const arma::colvec row_means = arma::mean(X, 1);

  const double ss_subject = k_d * arma::accu(arma::square(row_means - grand_mean));
  const double ss_rater = n_d * arma::accu(arma::square(col_means.t() - grand_mean));
  const double ss_total = arma::accu(arma::square(X - grand_mean));
  const double ss_error = ss_total - ss_subject - ss_rater;
  const double ss_within = ss_total - ss_subject;

  const double df_subject = n_d - 1.0;
  const double df_rater = k_d - 1.0;
  const double df_error = (n_d - 1.0) * (k_d - 1.0);
  const double df_within = n_d * (k_d - 1.0);

  const double ms_subject = (df_subject > 0.0) ? ss_subject / df_subject : NA_REAL;
  const double ms_rater = (df_rater > 0.0) ? ss_rater / df_rater : NA_REAL;
  const double ms_error = (df_error > 0.0) ? ss_error / df_error : NA_REAL;
  const double ms_within = (df_within > 0.0) ? ss_within / df_within : NA_REAL;

  const double icc1 = (is_positive_finite(ms_subject + (k_d - 1.0) * ms_within))
    ? (ms_subject - ms_within) / (ms_subject + (k_d - 1.0) * ms_within)
    : NA_REAL;
  const double icc2 = (is_positive_finite(ms_subject + (k_d - 1.0) * ms_error + k_d * (ms_rater - ms_error) / n_d))
    ? (ms_subject - ms_error) / (ms_subject + (k_d - 1.0) * ms_error + k_d * (ms_rater - ms_error) / n_d)
    : NA_REAL;
  const double icc3 = (is_positive_finite(ms_subject + (k_d - 1.0) * ms_error))
    ? (ms_subject - ms_error) / (ms_subject + (k_d - 1.0) * ms_error)
    : NA_REAL;
  const double icc1k = is_positive_finite(ms_subject) ? (ms_subject - ms_within) / ms_subject : NA_REAL;
  const double icc2k = is_positive_finite(ms_subject + (ms_rater - ms_error) / n_d)
    ? (ms_subject - ms_error) / (ms_subject + (ms_rater - ms_error) / n_d)
    : NA_REAL;
  const double icc3k = is_positive_finite(ms_subject) ? (ms_subject - ms_error) / ms_subject : NA_REAL;

  const double f11 = is_positive_finite(ms_within) ? ms_subject / ms_within : NA_REAL;
  const double f21 = is_positive_finite(ms_error) ? ms_subject / ms_error : NA_REAL;
  const double p11 = std::isfinite(f11)
    ? -std::expm1(R::pf(f11, df_subject, df_within, 1, 1))
    : NA_REAL;
  const double p21 = std::isfinite(f21)
    ? -std::expm1(R::pf(f21, df_subject, df_error, 1, 1))
    : NA_REAL;

  Rcpp::NumericVector lwr(6, NA_REAL);
  Rcpp::NumericVector upr(6, NA_REAL);

  if (return_ci) {
    const double alpha = 1.0 - conf_level;

    if (std::isfinite(f11) && df_subject > 0.0 && df_within > 0.0) {
      const double f1l = f11 / R::qf(1.0 - alpha / 2.0, df_subject, df_within, 1, 0);
      const double f1u = f11 * R::qf(1.0 - alpha / 2.0, df_within, df_subject, 1, 0);
      if (std::isfinite(f1l) && std::isfinite(f1u)) {
        lwr[0] = (f1l - 1.0) / (f1l + (k_d - 1.0));
        upr[0] = (f1u - 1.0) / (f1u + (k_d - 1.0));
        lwr[3] = 1.0 - 1.0 / f1l;
        upr[3] = 1.0 - 1.0 / f1u;
      }
    }

    if (std::isfinite(f21) && df_subject > 0.0 && df_error > 0.0) {
      const double f3l = f21 / R::qf(1.0 - alpha / 2.0, df_subject, df_error, 1, 0);
      const double f3u = f21 * R::qf(1.0 - alpha / 2.0, df_error, df_subject, 1, 0);
      if (std::isfinite(f3l) && std::isfinite(f3u)) {
        lwr[2] = (f3l - 1.0) / (f3l + (k_d - 1.0));
        upr[2] = (f3u - 1.0) / (f3u + (k_d - 1.0));
        lwr[5] = 1.0 - 1.0 / f3l;
        upr[5] = 1.0 - 1.0 / f3u;
      }
    }

    if (std::isfinite(icc2) && is_positive_finite(ms_error)) {
      const double fj = ms_rater / ms_error;
      const double vn = (k_d - 1.0) * (n_d - 1.0) *
        std::pow(k_d * icc2 * fj + n_d * (1.0 + (k_d - 1.0) * icc2) - k_d * icc2, 2.0);
      const double vd = (n_d - 1.0) * k_d * k_d * icc2 * icc2 * fj * fj +
        std::pow(n_d * (1.0 + (k_d - 1.0) * icc2) - k_d * icc2, 2.0);
      if (is_positive_finite(vn) && is_positive_finite(vd)) {
        const double v = vn / vd;
        if (is_positive_finite(v)) {
          const double f2u = R::qf(1.0 - alpha / 2.0, n_d - 1.0, v, 1, 0);
          const double f2l = R::qf(1.0 - alpha / 2.0, v, n_d - 1.0, 1, 0);
          if (std::isfinite(f2u) && std::isfinite(f2l)) {
            const double num_l = n_d * (ms_subject - f2u * ms_error);
            const double den_l = f2u * (k_d * ms_rater + (k_d * n_d - k_d - n_d) * ms_error) + n_d * ms_subject;
            const double num_u = n_d * (f2l * ms_subject - ms_error);
            const double den_u = k_d * ms_rater + (k_d * n_d - k_d - n_d) * ms_error + n_d * f2l * ms_subject;
            if (is_positive_finite(den_l) && is_positive_finite(den_u)) {
              lwr[1] = num_l / den_l;
              upr[1] = num_u / den_u;
              lwr[4] = lwr[1] * k_d / (1.0 + lwr[1] * (k_d - 1.0));
              upr[4] = upr[1] * k_d / (1.0 + upr[1] * (k_d - 1.0));
            }
          }
        }
      }
    }
  }

  Rcpp::CharacterVector coefficient =
    Rcpp::CharacterVector::create("ICC1", "ICC2", "ICC3", "ICC1k", "ICC2k", "ICC3k");
  Rcpp::CharacterVector label =
    Rcpp::CharacterVector::create("Single absolute", "Single random", "Single fixed",
                                  "Average absolute", "Average random", "Average fixed");
  Rcpp::NumericVector estimate =
    Rcpp::NumericVector::create(icc1, icc2, icc3, icc1k, icc2k, icc3k);
  Rcpp::NumericVector statistic =
    Rcpp::NumericVector::create(f11, f21, f21, f11, f21, f21);
  Rcpp::NumericVector df1 =
    Rcpp::NumericVector::create(df_subject, df_subject, df_subject, df_subject, df_subject, df_subject);
  Rcpp::NumericVector df2 =
    Rcpp::NumericVector::create(df_within, df_error, df_error, df_within, df_error, df_error);
  Rcpp::NumericVector p_value =
    Rcpp::NumericVector::create(p11, p21, p21, p11, p21, p21);

  Rcpp::List coefficient_list = Rcpp::List::create(
    Rcpp::_["coefficient"] = coefficient,
    Rcpp::_["label"] = label,
    Rcpp::_["estimate"] = estimate,
    Rcpp::_["statistic"] = statistic,
    Rcpp::_["df1"] = df1,
    Rcpp::_["df2"] = df2,
    Rcpp::_["p_value"] = p_value
  );
  if (return_ci) {
    coefficient_list["lwr"] = lwr;
    coefficient_list["upr"] = upr;
  }
  coefficient_list.attr("class") = "data.frame";
  coefficient_list.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -6);

  Rcpp::List anova = Rcpp::List::create(
    Rcpp::_["component"] = Rcpp::CharacterVector::create("subjects", "raters", "residual"),
    Rcpp::_["df"] = Rcpp::NumericVector::create(df_subject, df_rater, df_error),
    Rcpp::_["sum_sq"] = Rcpp::NumericVector::create(ss_subject, ss_rater, ss_error),
    Rcpp::_["mean_sq"] = Rcpp::NumericVector::create(ms_subject, ms_rater, ms_error),
    Rcpp::_["statistic"] = Rcpp::NumericVector::create(f11, is_positive_finite(ms_error) ? ms_rater / ms_error : NA_REAL, NA_REAL),
    Rcpp::_["p_value"] = Rcpp::NumericVector::create(p21, is_positive_finite(ms_error) ? -std::expm1(R::pf(ms_rater / ms_error, df_rater, df_error, 1, 1)) : NA_REAL, NA_REAL)
  );
  anova.attr("class") = "data.frame";
  anova.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -3);

  Rcpp::NumericVector mean_squares =
    Rcpp::NumericVector::create(ms_subject, ms_rater, ms_error, ms_within);
  mean_squares.names() =
    Rcpp::CharacterVector::create("ms_subject", "ms_rater", "ms_error", "ms_within");

  return Rcpp::List::create(
    Rcpp::_["coefficients"] = coefficient_list,
    Rcpp::_["anova"] = anova,
    Rcpp::_["mean_squares"] = mean_squares,
    Rcpp::_["n_subjects"] = static_cast<int>(n),
    Rcpp::_["n_raters"] = static_cast<int>(k)
  );
}
