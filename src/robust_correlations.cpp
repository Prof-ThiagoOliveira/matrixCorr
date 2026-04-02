// Thiago de Paula Oliveira
// robust_correlations.cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <random>
#include "matrixCorr_omp.h"

#include "matrixCorr_detail.h"

using namespace Rcpp;
using namespace arma;
using matrixCorr_detail::ranking::rank_vector;

namespace {

inline double clamp_corr(double x) {
  if (!std::isfinite(x)) return arma::datum::nan;
  if (x > 1.0) return 1.0;
  if (x < -1.0) return -1.0;
  return x;
}

inline double median_sorted_vec(const arma::vec& xs) {
  const arma::uword n = xs.n_elem;
  if (n == 0) return arma::datum::nan;
  if ((n % 2u) == 1u) return xs[n / 2u];
  return 0.5 * (xs[n / 2u - 1u] + xs[n / 2u]);
}

inline double median_sorted_std(const std::vector<double>& xs) {
  const std::size_t n = xs.size();
  if (n == 0u) return arma::datum::nan;
  if ((n % 2u) == 1u) return xs[n / 2u];
  return 0.5 * (xs[n / 2u - 1u] + xs[n / 2u]);
}

inline double order_stat_std(std::vector<double>& x, const std::size_t k) {
  std::nth_element(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(k), x.end());
  return x[k];
}

inline double order_stat_pair_std(std::vector<double>& x,
                                  const std::size_t k_hi,
                                  double& x_lo,
                                  double& x_hi) {
  std::nth_element(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(k_hi), x.end());
  x_hi = x[k_hi];
  x_lo = *std::max_element(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(k_hi));
  return x_hi;
}

inline double median_select_std(std::vector<double>& x) {
  const std::size_t n = x.size();
  if (n == 0u) return arma::datum::nan;
  if ((n % 2u) == 1u) {
    return order_stat_std(x, n / 2u);
  }
  double lo = 0.0, hi = 0.0;
  order_stat_pair_std(x, n / 2u, lo, hi);
  return 0.5 * (lo + hi);
}

inline double raw_mad(const arma::vec& x, const double med) {
  arma::vec dev = arma::abs(x - med);
  dev = arma::sort(dev);
  return median_sorted_vec(dev);
}

inline double idealf_iqr(const arma::vec& x) {
  const arma::uword n = x.n_elem;
  if (n < 2u) return 0.0;
  arma::vec y = arma::sort(x);
  const double j_d = std::floor(static_cast<double>(n) / 4.0 + 5.0 / 12.0);
  arma::uword j = static_cast<arma::uword>(std::max(1.0, j_d));
  if (j >= n) j = n - 1u;
  const double g = static_cast<double>(n) / 4.0 - static_cast<double>(j) + 5.0 / 12.0;
  const double low = (1.0 - g) * y[j - 1u] + g * y[j];
  arma::uword k = n - j + 1u;
  if (k > n) k = n;
  if (k < 2u) k = 2u;
  const double up = (1.0 - g) * y[k - 1u] + g * y[k - 2u];
  return up - low;
}

inline double idealf_iqr_sorted_std(const std::vector<double>& y) {
  const std::size_t n = y.size();
  if (n < 2u) return 0.0;

  const double j_d = std::floor(static_cast<double>(n) / 4.0 + 5.0 / 12.0);
  std::size_t j = static_cast<std::size_t>(std::max(1.0, j_d));
  if (j >= n) j = n - 1u;
  const double g = static_cast<double>(n) / 4.0 - static_cast<double>(j) + 5.0 / 12.0;
  const double low = (1.0 - g) * y[j - 1u] + g * y[j];

  std::size_t k = n - j + 1u;
  if (k > n) k = n;
  if (k < 2u) k = 2u;
  const double up = (1.0 - g) * y[k - 1u] + g * y[k - 2u];
  return up - low;
}

inline double idealf_iqr_select_std(std::vector<double>& y) {
  const std::size_t n = y.size();
  if (n < 2u) return 0.0;

  const double j_d = std::floor(static_cast<double>(n) / 4.0 + 5.0 / 12.0);
  std::size_t j = static_cast<std::size_t>(std::max(1.0, j_d));
  if (j >= n) j = n - 1u;
  const double g = static_cast<double>(n) / 4.0 - static_cast<double>(j) + 5.0 / 12.0;

  double low_lo = 0.0, low_hi = 0.0;
  order_stat_pair_std(y, j, low_lo, low_hi);
  const double low = (1.0 - g) * low_lo + g * low_hi;

  std::size_t k = n - j + 1u;
  if (k > n) k = n;
  if (k < 2u) k = 2u;
  double up_lo = 0.0, up_hi = 0.0;
  order_stat_pair_std(y, k - 1u, up_lo, up_hi);
  const double up = (1.0 - g) * up_hi + g * up_lo;
  return up - low;
}

inline double robust_scale_fallback(const arma::vec& x, const double med) {
  double sc = raw_mad(x, med);
  if (sc > 0.0 && std::isfinite(sc)) return sc;

  arma::vec xs = arma::sort(x);
  const double q1 = matrixCorr_detail::quantile_utils::quantile_sorted(xs, 0.25);
  const double q3 = matrixCorr_detail::quantile_utils::quantile_sorted(xs, 0.75);
  sc = (q3 - q1) / 1.3489795003921634;
  if (sc > 0.0 && std::isfinite(sc)) return sc;

  arma::vec centered = x - med;
  const double d2 = arma::dot(centered, centered);
  if (d2 > 0.0) return std::sqrt(d2 / static_cast<double>(x.n_elem - 1u));
  return 0.0;
}

inline double pearson_vec_core(const arma::vec& x, const arma::vec& y) {
  if (x.n_elem != y.n_elem || x.n_elem < 2u) return arma::datum::nan;
  const double mx = arma::mean(x);
  const double my = arma::mean(y);
  arma::vec xc = x - mx;
  arma::vec yc = y - my;
  const double dx = arma::dot(xc, xc);
  const double dy = arma::dot(yc, yc);
  if (!(dx > 0.0) || !(dy > 0.0)) return arma::datum::nan;
  return clamp_corr(arma::dot(xc, yc) / std::sqrt(dx * dy));
}

inline double spearman_vec_core(const arma::vec& x, const arma::vec& y) {
  if (x.n_elem != y.n_elem || x.n_elem < 2u) return arma::datum::nan;
  arma::vec rx(x.n_elem), ry(y.n_elem);
  rank_vector(x, rx);
  rank_vector(y, ry);
  return pearson_vec_core(rx, ry);
}

inline double pbos_scalar(const arma::vec& x,
                          const double beta,
                          double& omega,
                          bool& ok) {
  const arma::uword n = x.n_elem;
  omega = arma::datum::nan;
  ok = false;
  if (n < 2u) return arma::datum::nan;

  const double med = median_sorted_vec(arma::sort(x));
  arma::vec temp = arma::sort(arma::abs(x - med));
  const double raw_idx = (1.0 - beta) * static_cast<double>(n);
  const double fuzz = 8.0 * std::numeric_limits<double>::epsilon() *
    std::max(1.0, std::abs(raw_idx));
  arma::uword idx = static_cast<arma::uword>(std::floor(raw_idx + fuzz));
  if (idx < 1u) idx = 1u;
  if (idx > n) idx = n;

  omega = temp[idx - 1u];
  if (!(omega > 0.0) || !std::isfinite(omega)) return arma::datum::nan;

  arma::vec psi = (x - med) / omega;
  arma::uword i1 = 0u, i2 = 0u;
  double sx = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    if (psi[i] < -1.0) {
      ++i1;
    } else if (psi[i] > 1.0) {
      ++i2;
    } else {
      sx += x[i];
    }
  }

  const arma::uword denom_n = n - i1 - i2;
  if (denom_n == 0u) return arma::datum::nan;

  const std::ptrdiff_t tail_balance =
    static_cast<std::ptrdiff_t>(i2) - static_cast<std::ptrdiff_t>(i1);
  ok = true;
  return (sx + omega * static_cast<double>(tail_balance)) /
    static_cast<double>(denom_n);
}

inline void fill_pb_bent_scores(const arma::vec& x,
                                arma::vec& a,
                                const double beta,
                                bool& ok) {
  const arma::uword n = x.n_elem;
  if (a.n_elem != n) a.set_size(n);
  a.fill(arma::datum::nan);
  ok = false;
  if (n < 2u) return;

  double omega = arma::datum::nan;
  bool pbos_ok = false;
  const double theta = pbos_scalar(x, beta, omega, pbos_ok);
  if (!pbos_ok || !(omega > 0.0) || !std::isfinite(theta)) return;

  a = (x - theta) / omega;
  for (arma::uword i = 0; i < n; ++i) {
    if (a[i] <= -1.0) {
      a[i] = -1.0;
    } else if (a[i] >= 1.0) {
      a[i] = 1.0;
    }
  }
  ok = true;
}

inline double pbcor_pair_complete_core(const arma::vec& x,
                                       const arma::vec& y,
                                       const double beta) {
  if (x.n_elem != y.n_elem || x.n_elem < 2u) return arma::datum::nan;

  arma::vec a(x.n_elem), b(y.n_elem);
  bool okx = false, oky = false;
  fill_pb_bent_scores(x, a, beta, okx);
  fill_pb_bent_scores(y, b, beta, oky);
  if (!okx || !oky) return arma::datum::nan;

  const double da = arma::dot(a, a);
  const double db = arma::dot(b, b);
  if (!(da > 0.0) || !(db > 0.0)) return arma::datum::nan;
  return clamp_corr(arma::dot(a, b) / std::sqrt(da * db));
}

inline void standardise_win_column(const arma::vec& x,
                                   arma::vec& z,
                                   const double tr,
                                   bool& ok) {
  const arma::uword n = x.n_elem;
  if (z.n_elem != n) z.set_size(n);
  z.zeros();
  ok = false;
  if (n < 2u) {
    z.fill(arma::datum::nan);
    return;
  }

  arma::vec xs = arma::sort(x);
  arma::uword g = static_cast<arma::uword>(std::floor(tr * static_cast<double>(n)));
  arma::uword ibot = g;
  arma::uword itop = n - g - 1u;
  if (ibot >= n || itop >= n || ibot > itop) {
    z.fill(arma::datum::nan);
    return;
  }

  const double xbot = xs[ibot];
  const double xtop = xs[itop];
  arma::vec xw = x;
  for (arma::uword i = 0; i < n; ++i) {
    if (xw[i] <= xbot) xw[i] = xbot;
    else if (xw[i] >= xtop) xw[i] = xtop;
  }

  const double mu = arma::mean(xw);
  arma::vec centered = xw - mu;
  const double d2 = arma::dot(centered, centered);
  if (!(d2 > 0.0)) {
    z.fill(arma::datum::nan);
    return;
  }
  z = centered / std::sqrt(d2);
  ok = true;
}

inline void prepare_skip_detection_column(const arma::vec& x,
                                          arma::vec& z,
                                          const bool stand,
                                          double& center) {
  const arma::uword n = x.n_elem;
  if (z.n_elem != n) z.set_size(n);
  if (!stand) {
    z = x;
    center = median_sorted_vec(arma::sort(x));
    return;
  }

  const double med = median_sorted_vec(arma::sort(x));
  const double s = robust_scale_fallback(x, med);
  if (s > 0.0 && std::isfinite(s)) z = (x - med) / s;
  else z.zeros();
  center = 0.0;
}

inline void prepare_skip_detection_matrix(const arma::mat& X,
                                          const bool stand,
                                          arma::mat& Z,
                                          arma::vec& centers) {
  Z.set_size(X.n_rows, X.n_cols);
  centers.set_size(X.n_cols);

  for (arma::uword j = 0; j < X.n_cols; ++j) {
    arma::vec zcol;
    double center = 0.0;
    prepare_skip_detection_column(X.col(j), zcol, stand, center);
    Z.col(j) = zcol;
    centers[j] = center;
  }
}

inline bool extract_overlap_pair(const arma::uvec& idx_j,
                                 const arma::uvec& idx_k,
                                 const double* colj_ptr,
                                 const double* colk_ptr,
                                 const int min_n,
                                 std::vector<arma::uword>& overlap_idx,
                                 arma::vec& xbuf,
                                 arma::vec& ybuf) {
  const std::size_t possible = std::min(idx_j.n_elem, idx_k.n_elem);
  if (possible < static_cast<std::size_t>(min_n)) return false;

  overlap_idx.clear();
  overlap_idx.reserve(possible);

  arma::uword ia = 0, ib = 0;
  while (ia < idx_j.n_elem && ib < idx_k.n_elem) {
    const arma::uword a_val = idx_j[ia];
    const arma::uword b_val = idx_k[ib];
    if (a_val == b_val) {
      overlap_idx.push_back(a_val);
      ++ia;
      ++ib;
    } else if (a_val < b_val) {
      ++ia;
    } else {
      ++ib;
    }
  }
  if (overlap_idx.size() < static_cast<std::size_t>(min_n)) return false;

  xbuf.set_size(overlap_idx.size());
  ybuf.set_size(overlap_idx.size());
  for (std::size_t t = 0; t < overlap_idx.size(); ++t) {
    const arma::uword row = overlap_idx[t];
    xbuf[t] = colj_ptr[row];
    ybuf[t] = colk_ptr[row];
  }
  return true;
}

struct SkipMaskEntry {
  int j = 0;
  int k = 0;
  std::vector<int> skipped_rows;
};

struct SkipBootstrapResult {
  double lwr = arma::datum::nan;
  double upr = arma::datum::nan;
  double p_value = arma::datum::nan;
};

struct SkipMultipleResult {
  arma::mat p_adjusted;
  arma::mat reject;
  double critical_p_value = arma::datum::nan;
};

inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

inline double harrell_davis_quantile_sorted_std(const std::vector<double>& x, const double q) {
  const std::size_t n = x.size();
  if (n == 0u) return arma::datum::nan;
  if (q <= 0.0) return x.front();
  if (q >= 1.0) return x.back();
  if (n == 1u) return x.front();

  const double a = (static_cast<double>(n) + 1.0) * q;
  const double b = (static_cast<double>(n) + 1.0) * (1.0 - q);
  double out = 0.0;
  for (std::size_t i = 1u; i <= n; ++i) {
    const double lo = static_cast<double>(i - 1u) / static_cast<double>(n);
    const double hi = static_cast<double>(i) / static_cast<double>(n);
    const double w = R::pbeta(hi, a, b, 1, 0) - R::pbeta(lo, a, b, 1, 0);
    out += w * x[i - 1u];
  }
  return out;
}

inline double skipcor_pair_core_ptr(const double* x_orig,
                                    const double* y_orig,
                                    const double* x_det,
                                    const double* y_det,
                                    const arma::uword n,
                                    const double center_x,
                                    const double center_y,
                                    const int method_int,
                                    const bool use_mad,
                                    const double gval,
                                    const int min_n,
                                    arma::vec& bx,
                                    arma::vec& by,
                                    arma::vec& bot,
                                    std::vector<unsigned char>& outlier,
                                    std::vector<arma::uword>& keep_idx,
                                    std::vector<double>& dist_buf,
                                    std::vector<double>& work_buf,
                                    arma::vec& xkeep,
                                    arma::vec& ykeep,
                                    int* skipped_count = nullptr,
                                    std::vector<int>* skipped_rows = nullptr,
                                    const std::vector<arma::uword>* row_index = nullptr) {
  if (skipped_count != nullptr) *skipped_count = NA_INTEGER;
  if (static_cast<int>(n) < min_n) return arma::datum::nan;

  bx.set_size(n);
  by.set_size(n);
  bot.set_size(n);
  outlier.assign(static_cast<std::size_t>(n), 0u);
  dist_buf.resize(static_cast<std::size_t>(n));
  work_buf.resize(static_cast<std::size_t>(n));

  for (arma::uword i = 0; i < n; ++i) {
    const double dx = x_det[i] - center_x;
    const double dy = y_det[i] - center_y;
    bx[i] = dx;
    by[i] = dy;
    bot[i] = dx * dx + dy * dy;
  }

  for (arma::uword i = 0; i < n; ++i) {
    if (!(bot[i] > 0.0)) continue;
    const double bxi = bx[i];
    const double byi = by[i];
    const double denom = std::sqrt(bot[i]);
    for (arma::uword j = 0; j < n; ++j) {
      const double dij = std::abs(bx[j] * bxi + by[j] * byi) / denom;
      dist_buf[static_cast<std::size_t>(j)] = dij;
      work_buf[static_cast<std::size_t>(j)] = dij;
    }

    const double med = median_select_std(work_buf);
    double spread = 0.0;
    if (use_mad) {
      for (arma::uword j = 0; j < n; ++j) {
        work_buf[static_cast<std::size_t>(j)] = std::abs(dist_buf[static_cast<std::size_t>(j)] - med);
      }
      spread = median_select_std(work_buf);
    } else {
      spread = idealf_iqr_select_std(work_buf);
    }
    if (!(spread > 0.0) || !std::isfinite(spread)) continue;
    const double thresh = med + gval * spread;
    for (arma::uword j = 0; j < n; ++j) {
      if (dist_buf[static_cast<std::size_t>(j)] > thresh) outlier[static_cast<std::size_t>(j)] = 1u;
    }
  }

  if (skipped_rows != nullptr) {
    skipped_rows->clear();
    for (arma::uword i = 0; i < n; ++i) {
      if (outlier[static_cast<std::size_t>(i)] == 0u) continue;
      const int row = (row_index == nullptr)
        ? static_cast<int>(i + 1u)
        : static_cast<int>((*row_index)[static_cast<std::size_t>(i)] + 1u);
      skipped_rows->push_back(row);
    }
  }

  keep_idx.clear();
  keep_idx.reserve(static_cast<std::size_t>(n));
  for (arma::uword i = 0; i < n; ++i) {
    if (outlier[static_cast<std::size_t>(i)] == 0u) keep_idx.push_back(i);
  }
  if (skipped_count != nullptr) {
    *skipped_count = static_cast<int>(n - static_cast<arma::uword>(keep_idx.size()));
  }
  if (keep_idx.size() < static_cast<std::size_t>(min_n)) return arma::datum::nan;

  if (method_int == 1) {
    xkeep.set_size(keep_idx.size());
    ykeep.set_size(keep_idx.size());
    for (std::size_t i = 0; i < keep_idx.size(); ++i) {
      const arma::uword idx = keep_idx[i];
      xkeep[static_cast<arma::uword>(i)] = x_orig[idx];
      ykeep[static_cast<arma::uword>(i)] = y_orig[idx];
    }
    return spearman_vec_core(xkeep, ykeep);
  }

  double sum_x = 0.0, sum_y = 0.0;
  for (const arma::uword idx : keep_idx) {
    sum_x += x_orig[idx];
    sum_y += y_orig[idx];
  }
  const double inv_n = 1.0 / static_cast<double>(keep_idx.size());
  const double mean_x = sum_x * inv_n;
  const double mean_y = sum_y * inv_n;

  double sxy = 0.0, sxx = 0.0, syy = 0.0;
  for (const arma::uword idx : keep_idx) {
    const double dx = x_orig[idx] - mean_x;
    const double dy = y_orig[idx] - mean_y;
    sxy += dx * dy;
    sxx += dx * dx;
    syy += dy * dy;
  }
  if (!(sxx > 0.0) || !(syy > 0.0)) return arma::datum::nan;
  return clamp_corr(sxy / std::sqrt(sxx * syy));
}

template <class RNG>
inline SkipBootstrapResult skipcor_pair_bootstrap_b2_ptr(const double* x_orig,
                                                         const double* y_orig,
                                                         const arma::uword n,
                                                         const int method_int,
                                                         const bool stand,
                                                         const bool use_mad,
                                                         const double gval,
                                                         const int min_n,
                                                         const int n_boot,
                                                         const double conf_level,
                                                         RNG& rng,
                                                         arma::vec& xboot,
                                                         arma::vec& yboot,
                                                         arma::vec& xdet_boot,
                                                         arma::vec& ydet_boot,
                                                         arma::vec& bx,
                                                         arma::vec& by,
                                                         arma::vec& bot,
                                                         std::vector<unsigned char>& outlier,
                                                         std::vector<arma::uword>& keep_idx,
                                                         std::vector<double>& dist_buf,
                                                         std::vector<double>& work_buf,
                                                         arma::vec& xkeep,
                                                         arma::vec& ykeep) {
  SkipBootstrapResult out;
  if (n_boot <= 1 || !(conf_level > 0.0 && conf_level < 1.0) || static_cast<int>(n) < min_n) {
    return out;
  }

  std::uniform_int_distribution<arma::uword> draw(0u, n - 1u);
  std::vector<double> boots;
  boots.reserve(static_cast<std::size_t>(n_boot));

  xboot.set_size(n);
  yboot.set_size(n);

  for (int b = 0; b < n_boot; ++b) {
    for (arma::uword i = 0; i < n; ++i) {
      const arma::uword idx = draw(rng);
      xboot[i] = x_orig[idx];
      yboot[i] = y_orig[idx];
    }

    double center_x = 0.0;
    double center_y = 0.0;
    prepare_skip_detection_column(xboot, xdet_boot, stand, center_x);
    prepare_skip_detection_column(yboot, ydet_boot, stand, center_y);
    const double val = skipcor_pair_core_ptr(
      xboot.memptr(),
      yboot.memptr(),
      xdet_boot.memptr(),
      ydet_boot.memptr(),
      n,
      center_x,
      center_y,
      method_int,
      use_mad,
      gval,
      min_n,
      bx,
      by,
      bot,
      outlier,
      keep_idx,
      dist_buf,
      work_buf,
      xkeep,
      ykeep
    );
    if (!std::isfinite(val)) {
      return out;
    }
    boots.push_back(val);
  }

  std::sort(boots.begin(), boots.end());
  const double alpha = 1.0 - conf_level;
  int l = static_cast<int>(std::llround(alpha * static_cast<double>(n_boot) / 2.0));
  if (l < 0) l = 0;
  if (l > n_boot - 1) l = n_boot - 1;
  int u = n_boot - l - 1;
  if (u < 0) u = 0;
  if (u > n_boot - 1) u = n_boot - 1;
  if (u < l) u = l;

  out.lwr = boots[static_cast<std::size_t>(l)];
  out.upr = boots[static_cast<std::size_t>(u)];

  int less = 0;
  int equal = 0;
  for (const double val : boots) {
    if (val < 0.0) {
      ++less;
    } else if (val == 0.0) {
      ++equal;
    }
  }
  const double q = (static_cast<double>(less) + 0.5 * static_cast<double>(equal)) /
    static_cast<double>(n_boot);
  out.p_value = std::min(1.0, 2.0 * std::min(q, 1.0 - q));
  return out;
}

inline arma::mat skipcor_bootstrap_pvalue_matrix_complete(const arma::mat& X,
                                                          const int method_int,
                                                          const bool stand,
                                                          const bool use_mad,
                                                          const double gval,
                                                          const int min_n,
                                                          const int n_boot,
                                                          const int n_threads,
                                                          const std::uint64_t seed_base) {
  const std::size_t n = X.n_rows;
  const std::size_t p = X.n_cols;
  arma::mat p_mat(p, p, fill::none);
  p_mat.fill(arma::datum::nan);

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    static thread_local arma::vec bx;
    static thread_local arma::vec by;
    static thread_local arma::vec bot;
    static thread_local arma::vec xkeep;
    static thread_local arma::vec ykeep;
    static thread_local std::vector<unsigned char> outlier;
    static thread_local std::vector<arma::uword> keep_idx;
    static thread_local std::vector<double> dist_buf;
    static thread_local std::vector<double> work_buf;
    static thread_local arma::vec xboot;
    static thread_local arma::vec yboot;
    static thread_local arma::vec xdet_boot;
    static thread_local arma::vec ydet_boot;

    const double* xorig_ptr = X.colptr(static_cast<arma::uword>(j));
    for (std::size_t k = static_cast<std::size_t>(j + 1u); k < p; ++k) {
      const std::uint64_t pair_seed = splitmix64(
        seed_base ^
        (static_cast<std::uint64_t>(j + 1u) << 32) ^
        static_cast<std::uint64_t>(k + 1u)
      );
      std::mt19937_64 rng(pair_seed);
      const SkipBootstrapResult infer = skipcor_pair_bootstrap_b2_ptr(
        xorig_ptr,
        X.colptr(static_cast<arma::uword>(k)),
        static_cast<arma::uword>(n),
        method_int,
        stand,
        use_mad,
        gval,
        min_n,
        n_boot,
        0.95,
        rng,
        xboot,
        yboot,
        xdet_boot,
        ydet_boot,
        bx,
        by,
        bot,
        outlier,
        keep_idx,
        dist_buf,
        work_buf,
        xkeep,
        ykeep
      );
      p_mat(j, k) = infer.p_value;
      p_mat(k, j) = infer.p_value;
    }
  }

  for (std::size_t j = 0; j < p; ++j) p_mat(j, j) = 0.0;
  return p_mat;
}

inline SkipMultipleResult skipcor_hochberg_inference(const arma::mat& p_mat,
                                                     const double alpha) {
  const std::size_t p = p_mat.n_cols;
  SkipMultipleResult out;
  out.p_adjusted = arma::mat(p, p, fill::none);
  out.reject = arma::mat(p, p, fill::zeros);
  out.p_adjusted.fill(arma::datum::nan);

  std::vector<std::pair<double, std::pair<std::size_t, std::size_t>>> vals;
  vals.reserve((p * (p - 1u)) / 2u);
  for (std::size_t j = 0; j < p; ++j) {
    for (std::size_t k = j + 1u; k < p; ++k) {
      const double pv = p_mat(j, k);
      if (std::isfinite(pv)) vals.push_back({pv, {j, k}});
    }
  }

  if (vals.empty()) {
    for (std::size_t j = 0; j < p; ++j) {
      out.p_adjusted(j, j) = 0.0;
      out.reject(j, j) = 0.0;
    }
    return out;
  }

  std::sort(vals.begin(), vals.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

  std::vector<double> adj_desc(vals.size(), 1.0);
  double running = 1.0;
  for (std::size_t i = 0; i < vals.size(); ++i) {
    const double candidate = std::min(1.0, static_cast<double>(i + 1u) * vals[i].first);
    running = std::min(running, candidate);
    adj_desc[i] = running;
  }

  for (std::size_t i = 0; i < vals.size(); ++i) {
    const std::size_t j = vals[i].second.first;
    const std::size_t k = vals[i].second.second;
    out.p_adjusted(j, k) = adj_desc[i];
    out.p_adjusted(k, j) = adj_desc[i];
    const double rej = (adj_desc[i] <= alpha) ? 1.0 : 0.0;
    out.reject(j, k) = rej;
    out.reject(k, j) = rej;
  }
  for (std::size_t j = 0; j < p; ++j) {
    out.p_adjusted(j, j) = 0.0;
    out.reject(j, j) = 0.0;
  }
  return out;
}

inline SkipMultipleResult skipcor_ecp_inference_complete(const arma::mat& observed_p,
                                                         const arma::uword n,
                                                         const arma::uword p,
                                                         const int method_int,
                                                         const bool stand,
                                                         const bool use_mad,
                                                         const double gval,
                                                         const int min_n,
                                                         const int n_boot,
                                                         const int n_mc,
                                                         const int n_threads,
                                                         const double alpha,
                                                         const std::uint64_t seed_base) {
  SkipMultipleResult out;
  out.reject = arma::mat(p, p, fill::zeros);
  out.p_adjusted = arma::mat();

  std::vector<double> min_pvals;
  min_pvals.reserve(static_cast<std::size_t>(n_mc));
  std::normal_distribution<double> norm01(0.0, 1.0);

  for (int d = 0; d < n_mc; ++d) {
    std::mt19937_64 rng(splitmix64(seed_base ^ static_cast<std::uint64_t>(d + 1)));
    arma::mat X0(n, p, fill::none);
    for (arma::uword col = 0; col < p; ++col) {
      for (arma::uword row = 0; row < n; ++row) {
        X0(row, col) = norm01(rng);
      }
    }

    const arma::mat p0 = skipcor_bootstrap_pvalue_matrix_complete(
      X0,
      method_int,
      stand,
      use_mad,
      gval,
      min_n,
      n_boot,
      n_threads,
      splitmix64(seed_base ^ (static_cast<std::uint64_t>(d + 1u) << 16))
    );

    double minp = arma::datum::inf;
    for (arma::uword j = 0; j < p; ++j) {
      for (arma::uword k = j + 1u; k < p; ++k) {
        const double pv = p0(j, k);
        if (std::isfinite(pv) && pv < minp) minp = pv;
      }
    }
    if (std::isfinite(minp)) min_pvals.push_back(minp);
  }

  if (!min_pvals.empty()) {
    std::sort(min_pvals.begin(), min_pvals.end());
    out.critical_p_value = harrell_davis_quantile_sorted_std(min_pvals, alpha);
  }

  if (std::isfinite(out.critical_p_value)) {
    for (arma::uword j = 0; j < p; ++j) {
      for (arma::uword k = j + 1u; k < p; ++k) {
        const double rej = (std::isfinite(observed_p(j, k)) && observed_p(j, k) <= out.critical_p_value) ? 1.0 : 0.0;
        out.reject(j, k) = rej;
        out.reject(k, j) = rej;
      }
    }
  }
  for (arma::uword j = 0; j < p; ++j) out.reject(j, j) = 0.0;
  return out;
}

} // namespace

// [[Rcpp::export]]
arma::mat pbcor_matrix_cpp(const arma::mat& X,
                           const double beta = 0.2,
                           const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");
  if (!X.is_finite()) stop("X contains NA/NaN/Inf; please handle missingness upstream.");

  arma::mat A(n, p, fill::zeros);
  arma::vec col_ss(p, fill::zeros);
  std::vector<unsigned char> col_valid(p, 0u);

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    bool ok = false;
    arma::vec acol(n, fill::zeros);
    fill_pb_bent_scores(X.col(static_cast<arma::uword>(j)), acol, beta, ok);
    if (ok) {
      const double ss = arma::dot(acol, acol);
      if (ss > 0.0 && std::isfinite(ss)) {
        A.col(static_cast<arma::uword>(j)) = acol / std::sqrt(ss);
        col_ss[static_cast<arma::uword>(j)] = ss;
        col_valid[static_cast<std::size_t>(j)] = 1u;
      }
    }
  }

  arma::mat R = arma::symmatu(A.t() * A);
  R.transform([](double val) { return clamp_corr(val); });

  for (std::size_t j = 0; j < p; ++j) {
    if (col_valid[j] == 0u) {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
      R(j, j) = 1.0;
    } else {
      R(j, j) = 1.0;
    }
  }
  return R;
}

// [[Rcpp::export]]
arma::mat pbcor_matrix_pairwise_cpp(const arma::mat& X,
                                    const double beta = 0.2,
                                    const int min_n = 5,
                                    const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");

  arma::mat R(p, p, fill::none);
  R.fill(arma::datum::nan);

  std::vector<arma::uvec> finite_idx(p);
  for (std::size_t j = 0; j < p; ++j) finite_idx[j] = arma::find_finite(X.col(j));

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    const arma::uvec& idx_j = finite_idx[static_cast<std::size_t>(j)];
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {
      if (k == static_cast<std::size_t>(j)) continue;
      const arma::uvec& idx_k = finite_idx[k];
      static thread_local std::vector<arma::uword> overlap_idx;
      static thread_local arma::vec xbuf;
      static thread_local arma::vec ybuf;
      const double* colj_ptr = X.colptr(static_cast<arma::uword>(j));
      const double* colk_ptr = X.colptr(static_cast<arma::uword>(k));
      if (!extract_overlap_pair(idx_j, idx_k, colj_ptr, colk_ptr, min_n, overlap_idx, xbuf, ybuf)) continue;

      const double val = pbcor_pair_complete_core(xbuf, ybuf, beta);
      R(j, k) = val;
      R(k, j) = val;
    }
  }

  for (std::size_t j = 0; j < p; ++j) {
    if (finite_idx[j].n_elem < 2u) {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
      R(j, j) = 1.0;
    } else {
      R(j, j) = 1.0;
    }
  }
  return R;
}

// [[Rcpp::export]]
arma::mat wincor_matrix_cpp(const arma::mat& X,
                            const double tr = 0.2,
                            const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");
  if (!X.is_finite()) stop("X contains NA/NaN/Inf; please handle missingness upstream.");

  arma::mat Z(n, p, fill::zeros);
  std::vector<unsigned char> col_valid(p, 0u);

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(static)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    bool ok = false;
    arma::vec zcol(n, fill::zeros);
    standardise_win_column(X.col(static_cast<arma::uword>(j)), zcol, tr, ok);
    if (ok) {
      Z.col(static_cast<arma::uword>(j)) = zcol;
      col_valid[static_cast<std::size_t>(j)] = 1u;
    }
  }

  arma::mat R = arma::symmatu(Z.t() * Z);
  R.transform([](double val) { return clamp_corr(val); });

  for (std::size_t j = 0; j < p; ++j) {
    if (col_valid[j] == 0u) {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
      R(j, j) = 1.0;
    } else {
      R(j, j) = 1.0;
    }
  }
  return R;
}

// [[Rcpp::export]]
arma::mat wincor_matrix_pairwise_cpp(const arma::mat& X,
                                     const double tr = 0.2,
                                     const int min_n = 5,
                                     const int n_threads = 1) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");

  arma::mat R(p, p, fill::none);
  R.fill(arma::datum::nan);

  std::vector<arma::uvec> finite_idx(p);
  for (std::size_t j = 0; j < p; ++j) finite_idx[j] = arma::find_finite(X.col(j));

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
    const arma::uvec& idx_j = finite_idx[static_cast<std::size_t>(j)];
    for (std::size_t k = static_cast<std::size_t>(j); k < p; ++k) {
      if (k == static_cast<std::size_t>(j)) continue;
      const arma::uvec& idx_k = finite_idx[k];
      static thread_local std::vector<arma::uword> overlap_idx;
      static thread_local arma::vec xbuf;
      static thread_local arma::vec ybuf;
      const double* colj_ptr = X.colptr(static_cast<arma::uword>(j));
      const double* colk_ptr = X.colptr(static_cast<arma::uword>(k));
      if (!extract_overlap_pair(idx_j, idx_k, colj_ptr, colk_ptr, min_n, overlap_idx, xbuf, ybuf)) continue;

      static thread_local arma::vec zj;
      static thread_local arma::vec zk;
      bool okj = false, okk = false;
      standardise_win_column(xbuf, zj, tr, okj);
      standardise_win_column(ybuf, zk, tr, okk);
      double val = arma::datum::nan;
      if (okj && okk) val = clamp_corr(arma::dot(zj, zk));
      R(j, k) = val;
      R(k, j) = val;
    }
  }

  for (std::size_t j = 0; j < p; ++j) {
    if (finite_idx[j].n_elem < 2u) {
      R.row(j).fill(arma::datum::nan);
      R.col(j).fill(arma::datum::nan);
      R(j, j) = 1.0;
    } else {
      R(j, j) = 1.0;
    }
  }
  return R;
}

// [[Rcpp::export]]
Rcpp::List skipcor_matrix_cpp(const arma::mat& X,
                              const int method_int = 0,
                              const bool stand = true,
                              const bool use_mad = false,
                              const double gval = 2.717803,
                              const int min_n = 5,
                              const int n_threads = 1,
                              const bool return_masks = false,
                              const bool return_inference = false,
                              const double conf_level = 0.95,
                              const int n_boot = 2000,
                              const int seed = 0,
                              const int multiple_method_int = 0,
                              const double fwe_level = 0.05,
                              const int n_mc = 1000) {
  const std::size_t n = X.n_rows, p = X.n_cols;
  if (p == 0 || n < 2) stop("X must have >=2 rows and >=1 column.");
  if (return_inference && !X.is_finite()) {
    stop("Bootstrap inference for skipped correlation currently requires complete finite data.");
  }

  arma::mat R(p, p, fill::none);
  R.fill(arma::datum::nan);
  arma::mat n_complete_mat(p, p, fill::none);
  arma::mat skipped_n_mat(p, p, fill::none);
  arma::mat skipped_prop_mat(p, p, fill::none);
  n_complete_mat.fill(arma::datum::nan);
  skipped_n_mat.fill(arma::datum::nan);
  skipped_prop_mat.fill(arma::datum::nan);
  arma::mat lwr_mat;
  arma::mat upr_mat;
  arma::mat p_mat;
  arma::mat p_adj_mat;
  arma::mat reject_mat;
  double critical_p_value = arma::datum::nan;
  if (return_inference) {
    lwr_mat.set_size(p, p);
    upr_mat.set_size(p, p);
    p_mat.set_size(p, p);
    p_adj_mat.set_size(p, p);
    reject_mat.set_size(p, p);
    lwr_mat.fill(arma::datum::nan);
    upr_mat.fill(arma::datum::nan);
    p_mat.fill(arma::datum::nan);
    p_adj_mat.fill(arma::datum::nan);
    reject_mat.fill(0.0);
  }
  std::vector<std::vector<SkipMaskEntry>> mask_entries;
  if (return_masks) mask_entries.resize(p);

  if (X.is_finite()) {
    arma::mat Zdet;
    arma::vec centers;
    prepare_skip_detection_matrix(X, stand, Zdet, centers);

#ifdef _OPENMP
    omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
    for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
      static thread_local arma::vec bx;
      static thread_local arma::vec by;
      static thread_local arma::vec bot;
      static thread_local arma::vec xkeep;
      static thread_local arma::vec ykeep;
      static thread_local std::vector<unsigned char> outlier;
      static thread_local std::vector<arma::uword> keep_idx;
      static thread_local std::vector<double> dist_buf;
      static thread_local std::vector<double> work_buf;
      static thread_local std::vector<int> skipped_rows;
      static thread_local arma::vec xboot;
      static thread_local arma::vec yboot;
      static thread_local arma::vec xdet_boot;
      static thread_local arma::vec ydet_boot;

      const double* xorig_ptr = X.colptr(static_cast<arma::uword>(j));
      const double* xdet_ptr = Zdet.colptr(static_cast<arma::uword>(j));
      for (std::size_t k = static_cast<std::size_t>(j + 1); k < p; ++k) {
        int skipped_count = NA_INTEGER;
        const double val = skipcor_pair_core_ptr(
          xorig_ptr,
          X.colptr(static_cast<arma::uword>(k)),
          xdet_ptr,
          Zdet.colptr(static_cast<arma::uword>(k)),
          static_cast<arma::uword>(n),
          centers[static_cast<arma::uword>(j)],
          centers[static_cast<arma::uword>(k)],
          method_int,
          use_mad,
          gval,
          min_n,
          bx,
          by,
          bot,
          outlier,
          keep_idx,
          dist_buf,
          work_buf,
          xkeep,
          ykeep,
          &skipped_count,
          (return_masks ? &skipped_rows : nullptr),
          nullptr
        );
        R(j, k) = val;
        R(k, j) = val;
        n_complete_mat(j, k) = static_cast<double>(n);
        n_complete_mat(k, j) = static_cast<double>(n);
        if (skipped_count != NA_INTEGER) {
          skipped_n_mat(j, k) = static_cast<double>(skipped_count);
          skipped_n_mat(k, j) = static_cast<double>(skipped_count);
          skipped_prop_mat(j, k) = static_cast<double>(skipped_count) / static_cast<double>(n);
          skipped_prop_mat(k, j) = skipped_prop_mat(j, k);
        }
        if (return_masks && !skipped_rows.empty()) {
          SkipMaskEntry entry;
          entry.j = static_cast<int>(j + 1);
          entry.k = static_cast<int>(k + 1);
          entry.skipped_rows = skipped_rows;
          mask_entries[static_cast<std::size_t>(j)].push_back(std::move(entry));
        }
        if (return_inference && std::isfinite(val)) {
          const std::uint64_t pair_seed = splitmix64(
            static_cast<std::uint64_t>(seed <= 0 ? 1 : seed) ^
            (static_cast<std::uint64_t>(j + 1u) << 32) ^
            static_cast<std::uint64_t>(k + 1u)
          );
          std::mt19937_64 rng(pair_seed);
          const SkipBootstrapResult infer = skipcor_pair_bootstrap_b2_ptr(
            xorig_ptr,
            X.colptr(static_cast<arma::uword>(k)),
            static_cast<arma::uword>(n),
            method_int,
            stand,
            use_mad,
            gval,
            min_n,
            n_boot,
            conf_level,
            rng,
            xboot,
            yboot,
            xdet_boot,
            ydet_boot,
            bx,
            by,
            bot,
            outlier,
            keep_idx,
            dist_buf,
            work_buf,
            xkeep,
            ykeep
          );
          lwr_mat(j, k) = infer.lwr;
          lwr_mat(k, j) = infer.lwr;
          upr_mat(j, k) = infer.upr;
          upr_mat(k, j) = infer.upr;
          p_mat(j, k) = infer.p_value;
          p_mat(k, j) = infer.p_value;
        }
      }
    }

    for (std::size_t j = 0; j < p; ++j) {
      R(j, j) = 1.0;
      n_complete_mat(j, j) = static_cast<double>(n);
      skipped_n_mat(j, j) = 0.0;
      skipped_prop_mat(j, j) = 0.0;
      if (return_inference) {
        lwr_mat(j, j) = 1.0;
        upr_mat(j, j) = 1.0;
        p_mat(j, j) = 0.0;
        p_adj_mat(j, j) = 0.0;
        reject_mat(j, j) = 0.0;
      }
    }

    if (return_inference && multiple_method_int == 1) {
      const SkipMultipleResult mult = skipcor_hochberg_inference(p_mat, fwe_level);
      p_adj_mat = mult.p_adjusted;
      reject_mat = mult.reject;
    } else if (return_inference && multiple_method_int == 2) {
      const SkipMultipleResult mult = skipcor_ecp_inference_complete(
        p_mat,
        static_cast<arma::uword>(n),
        static_cast<arma::uword>(p),
        method_int,
        stand,
        use_mad,
        gval,
        min_n,
        n_boot,
        n_mc,
        n_threads,
        fwe_level,
        splitmix64(static_cast<std::uint64_t>(seed <= 0 ? 1 : seed) ^ 0xA5A5A5A5ULL)
      );
      critical_p_value = mult.critical_p_value;
      reject_mat = mult.reject;
    }
  } else {
    std::vector<arma::uvec> finite_idx(p);
    for (std::size_t j = 0; j < p; ++j) finite_idx[j] = arma::find_finite(X.col(j));

#ifdef _OPENMP
    omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(dynamic)
#endif
    for (std::ptrdiff_t j = 0; j < static_cast<std::ptrdiff_t>(p); ++j) {
      const arma::uvec& idx_j = finite_idx[static_cast<std::size_t>(j)];
      for (std::size_t k = static_cast<std::size_t>(j + 1); k < p; ++k) {
        const arma::uvec& idx_k = finite_idx[k];
        static thread_local std::vector<arma::uword> overlap_idx;
        static thread_local arma::vec xbuf;
        static thread_local arma::vec ybuf;
        const double* colj_ptr = X.colptr(static_cast<arma::uword>(j));
        const double* colk_ptr = X.colptr(static_cast<arma::uword>(k));
        if (!extract_overlap_pair(idx_j, idx_k, colj_ptr, colk_ptr, min_n, overlap_idx, xbuf, ybuf)) continue;

        static thread_local arma::vec xdet;
        static thread_local arma::vec ydet;
        static thread_local arma::vec bx;
        static thread_local arma::vec by;
        static thread_local arma::vec bot;
        static thread_local arma::vec xkeep;
        static thread_local arma::vec ykeep;
        static thread_local std::vector<unsigned char> outlier;
        static thread_local std::vector<arma::uword> keep_idx;
        static thread_local std::vector<double> dist_buf;
        static thread_local std::vector<double> work_buf;
        static thread_local std::vector<int> skipped_rows;
        double center_x = 0.0, center_y = 0.0;
        prepare_skip_detection_column(xbuf, xdet, stand, center_x);
        prepare_skip_detection_column(ybuf, ydet, stand, center_y);
        int skipped_count = NA_INTEGER;
        const double val = skipcor_pair_core_ptr(
          xbuf.memptr(),
          ybuf.memptr(),
          xdet.memptr(),
          ydet.memptr(),
          xbuf.n_elem,
          center_x,
          center_y,
          method_int,
          use_mad,
          gval,
          min_n,
          bx,
          by,
          bot,
          outlier,
          keep_idx,
          dist_buf,
          work_buf,
          xkeep,
          ykeep,
          &skipped_count,
          (return_masks ? &skipped_rows : nullptr),
          (return_masks ? &overlap_idx : nullptr)
        );
        R(j, k) = val;
        R(k, j) = val;
        n_complete_mat(j, k) = static_cast<double>(xbuf.n_elem);
        n_complete_mat(k, j) = static_cast<double>(xbuf.n_elem);
        if (skipped_count != NA_INTEGER && xbuf.n_elem > 0u) {
          skipped_n_mat(j, k) = static_cast<double>(skipped_count);
          skipped_n_mat(k, j) = static_cast<double>(skipped_count);
          skipped_prop_mat(j, k) = static_cast<double>(skipped_count) / static_cast<double>(xbuf.n_elem);
          skipped_prop_mat(k, j) = skipped_prop_mat(j, k);
        }
        if (return_masks && !skipped_rows.empty()) {
          SkipMaskEntry entry;
          entry.j = static_cast<int>(j + 1);
          entry.k = static_cast<int>(k + 1);
          entry.skipped_rows = skipped_rows;
          mask_entries[static_cast<std::size_t>(j)].push_back(std::move(entry));
        }
      }
    }

    for (std::size_t j = 0; j < p; ++j) {
      if (finite_idx[j].n_elem < 2u) {
        R.row(j).fill(arma::datum::nan);
        R.col(j).fill(arma::datum::nan);
        n_complete_mat.row(j).fill(arma::datum::nan);
        n_complete_mat.col(j).fill(arma::datum::nan);
        skipped_n_mat.row(j).fill(arma::datum::nan);
        skipped_n_mat.col(j).fill(arma::datum::nan);
        skipped_prop_mat.row(j).fill(arma::datum::nan);
        skipped_prop_mat.col(j).fill(arma::datum::nan);
        R(j, j) = 1.0;
      } else {
        R(j, j) = 1.0;
        n_complete_mat(j, j) = static_cast<double>(finite_idx[j].n_elem);
        skipped_n_mat(j, j) = 0.0;
        skipped_prop_mat(j, j) = 0.0;
      }
    }
  }

  Rcpp::IntegerVector pair_i;
  Rcpp::IntegerVector pair_j;
  Rcpp::List skipped;
  if (return_masks) {
    std::size_t n_entries = 0u;
    for (const auto& per_j : mask_entries) n_entries += per_j.size();
    pair_i = Rcpp::IntegerVector(static_cast<R_xlen_t>(n_entries));
    pair_j = Rcpp::IntegerVector(static_cast<R_xlen_t>(n_entries));
    skipped = Rcpp::List(static_cast<R_xlen_t>(n_entries));

    std::size_t pos = 0u;
    for (const auto& per_j : mask_entries) {
      for (const auto& entry : per_j) {
        pair_i[static_cast<R_xlen_t>(pos)] = entry.j;
        pair_j[static_cast<R_xlen_t>(pos)] = entry.k;
        skipped[static_cast<R_xlen_t>(pos)] = Rcpp::wrap(entry.skipped_rows);
        ++pos;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::_["cor"] = R,
    Rcpp::_["pair_i"] = pair_i,
    Rcpp::_["pair_j"] = pair_j,
    Rcpp::_["skipped_rows"] = skipped,
    Rcpp::_["n_complete"] = n_complete_mat,
    Rcpp::_["skipped_n"] = skipped_n_mat,
    Rcpp::_["skipped_prop"] = skipped_prop_mat,
    Rcpp::_["lwr"] = lwr_mat,
    Rcpp::_["upr"] = upr_mat,
    Rcpp::_["p_value"] = p_mat,
    Rcpp::_["p_value_adjusted"] = p_adj_mat,
    Rcpp::_["reject"] = reject_mat,
    Rcpp::_["critical_p_value"] = critical_p_value
  );
}
