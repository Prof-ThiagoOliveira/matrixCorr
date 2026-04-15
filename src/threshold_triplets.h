#ifndef MATRIXCORR_THRESHOLD_TRIPLETS_H
#define MATRIXCORR_THRESHOLD_TRIPLETS_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace matrixCorr_detail {
namespace threshold_triplets {

struct TripletBuffer {
  std::vector<int> i;
  std::vector<int> j;
  std::vector<double> x;
};

inline bool retain_value(const double value, const double threshold) {
  return std::isfinite(value) && std::abs(value) >= threshold;
}

template <class BlockComputer>
TripletBuffer collect_upper_triplets(const std::size_t p,
                                     const std::size_t block_size,
                                     const bool include_diag,
                                     const double threshold,
                                     BlockComputer&& compute_block) {
  TripletBuffer out;
  if (p == 0u) return out;
  const std::size_t bs = std::max<std::size_t>(1u, block_size);

  for (std::size_t j0 = 0u; j0 < p; j0 += bs) {
    const std::size_t j1 = std::min<std::size_t>(p, j0 + bs);
    for (std::size_t k0 = j0; k0 < p; k0 += bs) {
      const std::size_t k1 = std::min<std::size_t>(p, k0 + bs);
      const arma::mat block = compute_block(j0, j1, k0, k1);
      const std::size_t bj = j1 - j0;
      const std::size_t bk = k1 - k0;

      if (j0 == k0) {
        for (std::size_t r = 0u; r < bj; ++r) {
          const std::size_t c_start = include_diag ? r : (r + 1u);
          for (std::size_t c = c_start; c < bk; ++c) {
            const double value = block(
              static_cast<arma::uword>(r),
              static_cast<arma::uword>(c)
            );
            if (!retain_value(value, threshold)) continue;
            out.i.push_back(static_cast<int>(j0 + r + 1u));
            out.j.push_back(static_cast<int>(k0 + c + 1u));
            out.x.push_back(value);
          }
        }
      } else {
        for (std::size_t r = 0u; r < bj; ++r) {
          for (std::size_t c = 0u; c < bk; ++c) {
            const double value = block(
              static_cast<arma::uword>(r),
              static_cast<arma::uword>(c)
            );
            if (!retain_value(value, threshold)) continue;
            out.i.push_back(static_cast<int>(j0 + r + 1u));
            out.j.push_back(static_cast<int>(k0 + c + 1u));
            out.x.push_back(value);
          }
        }
      }
    }
  }

  return out;
}

inline Rcpp::List as_list(const TripletBuffer& out) {
  return Rcpp::List::create(
    Rcpp::_["i"] = Rcpp::wrap(out.i),
    Rcpp::_["j"] = Rcpp::wrap(out.j),
    Rcpp::_["x"] = Rcpp::wrap(out.x)
  );
}

} // namespace threshold_triplets
} // namespace matrixCorr_detail

#endif // MATRIXCORR_THRESHOLD_TRIPLETS_H
