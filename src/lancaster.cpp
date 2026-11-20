// Pairwise Lancaster estimates and inference (Holzmann & Klar, 2025).
#include "matrixCorr_detail.h"
#include "correlation_math.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <random>
#include <vector>

namespace {
using Values = std::vector<double>;
using Pair = std::pair<int, int>;

// Include half ranks for bootstrap ties. Reflect the lower tail exactly so
// equal absolute scores remain equal when their squares are correlated.
Values normal_scores(int n) {
  Values table(2 * static_cast<size_t>(n) + 1);
  for (size_t k = 2; k <= static_cast<size_t>(n) + 1; ++k) {
    double z = R::qnorm(0.5 * k / (n + 1.0), 0, 1, true, false);
    table[k] = z;
    table[table.size() + 1 - k] = -z;
  }
  return table;
}

Values rank_transform(const Values& raw, const Values& table) {
  std::vector<size_t> order(raw.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b) { return raw[a] < raw[b]; });
  Values scores(raw.size());
  for (size_t i = 0; i < order.size();) {
    size_t end = i + 1;
    while (end < order.size() && raw[order[end]] == raw[order[i]]) ++end;
    const double z = table[i + 1 + end];
    for (; i < end; ++i) scores[order[i]] = z;
  }
  return scores;
}

// Shifting before summation protects tightly clustered observations. Binary
// rescaling protects the subtraction when the data span the double range.
bool center_and_normalize(Values& values) {
  const auto range = std::minmax_element(values.begin(), values.end());
  if (*range.first == *range.second) return false;
  const int exponent = std::ilogb(std::max(std::abs(*range.first), std::abs(*range.second)));
  const double origin = std::scalbn(values.front(), -exponent);
  long double sum = 0;
  for (double& x : values) {
    x = std::scalbn(x, -exponent) - origin;
    sum += x;
  }
  const long double mean = sum / values.size();
  long double ss = 0;
  for (double& x : values) {
    x = static_cast<double>(x - mean);
    ss += static_cast<long double>(x) * x;
  }
  const double norm = std::sqrt(static_cast<double>(ss));
  for (double& x : values) x /= norm;
  return true;
}

struct Scores {
  Values linear, quadratic;
  bool valid = false;
};

Scores make_scores(const Values& raw, const Values& quantiles) {
  Scores scores;
  if (raw.size() < 3) return scores;
  scores.linear = quantiles.empty() ? raw : rank_transform(raw, quantiles);
  // For tied ranks the normal scores can have nonzero mean. Their squares
  // must be taken before centering; moment scores must be centered first.
  if (!quantiles.empty()) scores.quadratic = scores.linear;
  if (!center_and_normalize(scores.linear)) return scores;
  if (quantiles.empty()) scores.quadratic = scores.linear;
  for (double& x : scores.quadratic) x *= x;
  scores.valid = center_and_normalize(scores.quadratic);
  return scores;
}

struct Components {
  double first, second;
  double estimate() const { return std::max(std::abs(first), std::abs(second)); }
};

Components components(const Scores& a, const Scores& b, const std::vector<int>* perm = nullptr) {
  long double first = 0, second = 0;
  for (size_t i = 0; i < a.linear.size(); ++i) {
    const size_t k = perm ? (*perm)[i] : i;
    first += static_cast<long double>(a.linear[i]) * b.linear[k];
    second += static_cast<long double>(a.quadratic[i]) * b.quadratic[k];
  }
  return {matrixCorr::correlation_math::clamp_corr_na(first),
          matrixCorr::correlation_math::clamp_corr_na(second)};
}

std::mt19937 pair_rng(int seed, int j, int k, int stream) {
  std::seed_seq sequence{seed, j, k, stream};
  return std::mt19937(sequence);
}

struct Interval {
  double lower = NA_REAL, upper = NA_REAL;
  int count = 0;
};

// Lemma 1: P(max(U,V) <= q) is the bivariate normal CDF at (q,q).
// Evaluate the complementary probability in the lower tails for precision.
double max_normal_quantile(double sx, double sy, double rho, double critical,
                           double tail) {
  const double scale = std::max(sx, sy);
  double low = scale * critical;
  if (std::min(sx, sy) == 0 || rho >= 1) return low;
  double high = scale * (critical + 2);
  for (int i = 0; i < 45; ++i) {
    const double q = 0.5 * (low + high);
    const double x = -q / sx, y = -q / sy;
    const double upper = matrixCorr_detail::norm1::Phi(x) +
      matrixCorr_detail::norm1::Phi(y) - matrixCorr_detail::bvn_fast::Phi2(x, y, rho);
    if (upper > tail) low = q;
    else high = q;
  }
  return high;
}

Interval bootstrap_ci(const Values& x, const Values& y, const Values& quantiles,
                      Components obs, int draws, double critical, double tail,
                      bool conservative, std::mt19937 rng) {
  Interval ci;
  const bool first = std::abs(obs.first) >= std::abs(obs.second);
  std::uniform_int_distribution<int> sample(0, static_cast<int>(x.size()) - 1);
  Values bx(x.size()), by(x.size());
  double mean1 = 0, mean2 = 0, ss1 = 0, ss2 = 0, cross = 0;
  for (int b = 0; b < draws; ++b) {
    for (size_t i = 0; i < x.size(); ++i) {
      const int index = sample(rng);
      bx[i] = x[index];
      by[i] = y[index];
    }
    const Scores a = make_scores(bx, quantiles), c = make_scores(by, quantiles);
    if (!a.valid || !c.valid) continue;
    const Components replicate = components(a, c);
    const double d1 = replicate.first - mean1, d2 = replicate.second - mean2;
    ++ci.count;
    mean1 += d1 / ci.count;
    mean2 += d2 / ci.count;
    ss1 += d1 * (replicate.first - mean1);
    ss2 += d2 * (replicate.second - mean2);
    cross += d1 * (replicate.second - mean2);
  }
  if (ci.count >= 2) {
    const double sx = std::sqrt(ss1 / (ci.count - 1));
    const double sy = std::sqrt(ss2 / (ci.count - 1));
    const double width = critical * (first ? sx : sy);
    double left = width;
    if (conservative) {
      double rho = sx > 0 && sy > 0 ? cross / std::sqrt(ss1 * ss2) : 0;
      if (obs.first * obs.second < 0) rho = -rho;
      left = max_normal_quantile(sx, sy, rho, critical, tail);
    }
    ci.lower = std::max(0.0, obs.estimate() - left);
    ci.upper = std::min(1.0, obs.estimate() + width);
  }
  return ci;
}

double permutation_p(const Scores& x, const Scores& y, double observed,
                     int draws, std::mt19937 rng) {
  std::vector<int> perm(x.linear.size());
  std::iota(perm.begin(), perm.end(), 0);
  int exceed = 0;
  const double cutoff = observed - 32 * std::numeric_limits<double>::epsilon();
  for (int b = 0; b < draws; ++b) {
    std::shuffle(perm.begin(), perm.end(), rng);
    if (components(x, y, &perm).estimate() >= cutoff) ++exceed;
  }
  return (exceed + 1.0) / (draws + 1.0);
}
} // namespace

// [[Rcpp::export]]
Rcpp::List lancaster_matrix_cpp(Rcpp::NumericMatrix X, bool rank_based, bool ci,
                                bool p_value, double conf_level, int n_boot,
                                int n_perm, int seed, int n_threads, bool conservative) {
  const int n = X.nrow(), p = X.ncol();
  const double critical = ci ? R::qnorm((1 - conf_level) / 2, 0, 1, false, false) : 0;
  auto matrix = [p](bool needed) {
    Rcpp::NumericMatrix out(needed ? p : 0, needed ? p : 0);
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  };
  auto est = matrix(true), first = matrix(true), second = matrix(true);
  auto lower = matrix(ci), upper = matrix(ci), prob = matrix(p_value);
  Rcpp::IntegerMatrix count(p, p), successes(ci ? p : 0, ci ? p : 0);
  const double* data = REAL(X);
  auto column = [&](int j) { return data + static_cast<size_t>(n) * j; };

  std::vector<bool> complete(p);
  for (int j = 0; j < p; ++j) {
    complete[j] = std::all_of(column(j), column(j) + n,
                             [](double x) { return std::isfinite(x); });
  }
  std::map<int, std::vector<Pair>> groups;
  for (int j = 0; j < p; ++j) {
    for (int k = j; k < p; ++k) {
      int size = n;
      if (!complete[j] || !complete[k]) {
        size = 0;
        for (int i = 0; i < n; ++i)
          size += std::isfinite(column(j)[i]) && std::isfinite(column(k)[i]);
      }
      count(j,k) = count(k,j) = size;
      if (size >= 3) groups[size].emplace_back(j,k);
    }
  }

  // Group by pair count so quantiles are computed once on the main thread
  // and only one table is live. Complete columns also share their transforms.
  for (const auto& group : groups) {
    Rcpp::checkUserInterrupt();
    const int size = group.first;
    const Values quantiles = rank_based ? normal_scores(size) : Values{};
    std::vector<Scores> cached(p);
    if (size == n) {
      for (int j = 0; j < p; ++j)
        if (complete[j]) cached[j] = make_scores(Values(column(j), column(j) + n), quantiles);
    }
    const auto& pairs = group.second;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads) if(n_threads > 1)
#endif
    for (std::ptrdiff_t t = 0; t < static_cast<std::ptrdiff_t>(pairs.size()); ++t) {
      const auto [j,k] = pairs[t];
      Values x, y;
      if (size != n || (ci && j != k)) {
        x.reserve(size);
        y.reserve(size);
        for (int i = 0; i < n; ++i) {
          if (std::isfinite(column(j)[i]) && std::isfinite(column(k)[i])) {
            x.push_back(column(j)[i]);
            y.push_back(column(k)[i]);
          }
        }
      }
      Scores a, b;
      if (size != n) {
        a = make_scores(x, quantiles);
        b = make_scores(y, quantiles);
      }
      const Scores& sx = size == n ? cached[j] : a;
      const Scores& sy = size == n ? cached[k] : b;
      if (!sx.valid || !sy.valid) continue;
      const Components obs = j == k ? Components{1, 1} : components(sx, sy);
      est(j,k) = est(k,j) = obs.estimate();
      first(j,k) = first(k,j) = obs.first;
      second(j,k) = second(k,j) = obs.second;
      if (j == k) {
        if (ci) lower(j,j) = upper(j,j) = 1;
        continue;
      }
      if (ci) {
        const auto interval = bootstrap_ci(x, y, quantiles, obs, n_boot, critical,
                                           (1 - conf_level) / 2, conservative, pair_rng(seed, j, k, 1));
        lower(j,k) = lower(k,j) = interval.lower;
        upper(j,k) = upper(k,j) = interval.upper;
        successes(j,k) = successes(k,j) = interval.count;
      }
      if (p_value) prob(j,k) = prob(k,j) = permutation_p(
        sx, sy, obs.estimate(), n_perm, pair_rng(seed, j, k, 2));
    }
  }
  Rcpp::List out = Rcpp::List::create(Rcpp::_["est"] = est, Rcpp::_["first"] = first,
                                      Rcpp::_["second"] = second, Rcpp::_["n_complete"] = count);
  if (ci) {
    out["lwr"] = lower;
    out["upr"] = upper;
    out["n_boot_success"] = successes;
  }
  if (p_value) out["p_value"] = prob;
  return out;
}
