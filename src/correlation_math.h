#pragma once

#include <R_ext/Arith.h>
#include <cmath>
#include <limits>

namespace matrixCorr {
namespace correlation_math {

inline double clamp_corr_na(double x) noexcept {
  if (!std::isfinite(x)) return NA_REAL;
  if (x > 1.0) return 1.0;
  if (x < -1.0) return -1.0;
  return x;
}

inline double clamp_corr_nan(double x) noexcept {
  if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
  if (x > 1.0) return 1.0;
  if (x < -1.0) return -1.0;
  return x;
}

inline double fisher_z(double r) noexcept {
  const double one_minus = std::nextafter(1.0, 0.0);
  if (r > one_minus) r = one_minus;
  if (r < -one_minus) r = -one_minus;
  return std::atanh(r);
}

inline double fisher_z_inv(double z) noexcept {
  double r = std::tanh(z);
  const double one_minus = std::nextafter(1.0, 0.0);
  if (r > one_minus) r = one_minus;
  if (r < -one_minus) r = -one_minus;
  return r;
}

} // namespace correlation_math
} // namespace matrixCorr
