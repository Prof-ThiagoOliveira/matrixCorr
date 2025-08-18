// Thiago de Paula Oliveira
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <vector>
using namespace Rcpp;

// ---------- tau-a helpers (merge-sort inversion count) ----------
static long long merge_count(IntegerVector& a, IntegerVector& tmp,
                             int L, int M, int R){
  int i = L, j = M, k = L;
  long long inv = 0;
  while (i < M && j <= R){
    if (a[i] <= a[j]) tmp[k++] = a[i++];
    else { tmp[k++] = a[j++]; inv += (M - i); }
  }
  while (i < M)    tmp[k++] = a[i++];
  while (j <= R)   tmp[k++] = a[j++];
  for (int t = L; t <= R; ++t) a[t] = tmp[t];
  return inv;
}

static long long sort_count(IntegerVector& a, IntegerVector& tmp,
                            int L, int R){
  if (R - L < 1) return 0LL;
  int M = L + (R - L) / 2;
  long long inv = 0;
  inv += sort_count(a, tmp, L, M);
  inv += sort_count(a, tmp, M + 1, R);
  inv += merge_count(a, tmp, L, M + 1, R);
  return inv;
}

// ---------- tau-b helpers (Fenwick tree) ----------
struct Fenwick {
  std::vector<long long> t;
  int n;
  Fenwick(int n): t(n + 1, 0), n(n) {}
  void add(int i, long long v){ for (; i <= n; i += i & -i) t[i] += v; }
  long long sum(int i) const { long long s = 0; for (; i > 0; i -= i & -i) s += t[i]; return s; }
};

// [[Rcpp::export]]
double kendall_tau_auto_cpp(NumericVector x, NumericVector y, double scale = 1e8){
  const int n = x.size();
  if (n != y.size() || n < 2) return NA_REAL;

  // Discretise to integers for stable equality checks
  std::vector<long long> xi(n), yi(n);
  for (int i = 0; i < n; ++i){
    xi[i] = (long long) std::floor(x[i] * scale);
    yi[i] = (long long) std::floor(y[i] * scale);
  }

  // Detect ties on each margin
  auto has_ties = [](const std::vector<long long>& v){
    std::vector<long long> s = v;
    std::sort(s.begin(), s.end());
    return std::adjacent_find(s.begin(), s.end()) != s.end();
  };
  const bool ties_x = has_ties(xi);
  const bool ties_y = has_ties(yi);

  const double n0 = (double)n * (double)(n - 1) / 2.0;

  // --------- No ties anywhere -> tau-a (fast path) ----------
  if (!ties_x && !ties_y){
    // order indices by x
    std::vector<int> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return xi[a] < xi[b]; });

    // y in x-sorted order
    IntegerVector y_ord(n), tmp(n);
    for (int i = 0; i < n; ++i) y_ord[i] = (int) yi[ord[i]];

    long long discord = sort_count(y_ord, tmp, 0, n - 1);
    return (n0 - 2.0 * (double)discord) / n0;  // tau-a
  }

  // --------- Ties present -> tau-b ----------
  // Order by x then y (stable within equal x by y)
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b){
    if (xi[a] != xi[b]) return xi[a] < xi[b];
    return yi[a] < yi[b];
  });

  // T_x (pairs tied on x)
  long long T_x = 0;
  for (int i = 0; i < n; ){
    int j = i + 1;
    while (j < n && xi[ord[j]] == xi[ord[i]]) ++j;
    long long g = j - i; T_x += g * (g - 1) / 2;
    i = j;
  }

  // T_y (pairs tied on y)
  std::vector<long long> y_sorted = yi;
  std::sort(y_sorted.begin(), y_sorted.end());
  long long T_y = 0;
  for (int i = 0; i < n; ){
    int j = i + 1;
    while (j < n && y_sorted[j] == y_sorted[i]) ++j;
    long long g = j - i; T_y += g * (g - 1) / 2;
    i = j;
  }

  // Coordinate-compress y
  std::vector<long long> y_unique = y_sorted;
  y_unique.erase(std::unique(y_unique.begin(), y_unique.end()), y_unique.end());
  auto rank_of = [&](long long v){
    return (int)(std::lower_bound(y_unique.begin(), y_unique.end(), v) - y_unique.begin()) + 1; // 1-based
  };

  // S = C - D via BIT; compare each x-group only to prior items (strict < and >)
  Fenwick bit((int) y_unique.size());
  long long processed = 0, S = 0;

  for (int i = 0; i < n; ){
    int j = i + 1;
    while (j < n && xi[ord[j]] == xi[ord[i]]) ++j;

    // Query against prior items
    for (int k = i; k < j; ++k){
      int idx = ord[k];
      int r = rank_of(yi[idx]);
      long long less    = bit.sum(r - 1);      // y_prev < y
      long long leq     = bit.sum(r);          // y_prev <= y
      long long greater = processed - leq;     // y_prev > y
      S += (less - greater);                   // +1 for concordant, -1 for discordant
    }
    // Insert current group
    for (int k = i; k < j; ++k){
      int r = rank_of(yi[ord[k]]);
      bit.add(r, 1);
    }
    processed += (j - i);
    i = j;
  }

  const double denom_x = n0 - (double)T_x;
  const double denom_y = n0 - (double)T_y;
  if (denom_x <= 0.0 || denom_y <= 0.0) return NA_REAL;

  const double denom = std::sqrt(denom_x * denom_y);
  return (double)S / denom;  // tau-b
}

// [[Rcpp::export]]
double kendall_tau_a_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y, double scale = 1e8) {
  const int n = x.size();
  if (n != y.size() || n < 2) return NA_REAL;

  // order by x
  std::vector<long long> xi(n), yi(n);
  for (int i = 0; i < n; ++i) {
    xi[i] = (long long) std::floor(x[i] * scale);
    yi[i] = (long long) std::floor(y[i] * scale);
  }
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b){ return xi[a] < xi[b]; });

  // y in x-sorted order
  Rcpp::IntegerVector y_ord(n), tmp(n);
  for (int i = 0; i < n; ++i) y_ord[i] = (int) yi[ord[i]];

  long long discord = sort_count(y_ord, tmp, 0, n - 1);
  const double n0 = (double)n * (double)(n - 1) / 2.0;
  return (n0 - 2.0 * (double)discord) / n0;
}

// [[Rcpp::export]]
double kendall_tau_b_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y, double scale = 1e8) {
  const int n = x.size();
  if (n != y.size() || n < 2) return NA_REAL;

  // discretise
  std::vector<long long> xi(n), yi(n);
  for (int i = 0; i < n; ++i) {
    xi[i] = (long long) std::floor(x[i] * scale);
    yi[i] = (long long) std::floor(y[i] * scale);
  }

  // order by x then y
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b){
    if (xi[a] != xi[b]) return xi[a] < xi[b];
    return yi[a] < yi[b];
  });

  // T_x
  long long T_x = 0;
  for (int i = 0; i < n; ) {
    int j = i + 1;
    while (j < n && xi[ord[j]] == xi[ord[i]]) ++j;
    long long g = j - i; T_x += g * (g - 1) / 2;
    i = j;
  }

  // T_y
  std::vector<long long> y_sorted = yi;
  std::sort(y_sorted.begin(), y_sorted.end());
  long long T_y = 0;
  for (int i = 0; i < n; ) {
    int j = i + 1;
    while (j < n && y_sorted[j] == y_sorted[i]) ++j;
    long long g = j - i; T_y += g * (g - 1) / 2;
    i = j;
  }

  // coordinate-compress y
  std::vector<long long> y_unique = y_sorted;
  y_unique.erase(std::unique(y_unique.begin(), y_unique.end()), y_unique.end());
  auto rank_of = [&](long long v){
    return (int)(std::lower_bound(y_unique.begin(), y_unique.end(), v) - y_unique.begin()) + 1;
  };

  // S = C - D (exclude ties in numerator via strict < and >)
  Fenwick bit((int) y_unique.size());
  long long processed = 0, S = 0;
  for (int i = 0; i < n; ) {
    int j = i + 1;
    while (j < n && xi[ord[j]] == xi[ord[i]]) ++j;

    for (int k = i; k < j; ++k) {
      int idx = ord[k];
      int r = rank_of(yi[idx]);
      long long less    = bit.sum(r - 1);
      long long leq     = bit.sum(r);
      long long greater = processed - leq;
      S += (less - greater);
    }
    for (int k = i; k < j; ++k) {
      int r = rank_of(yi[ord[k]]);
      bit.add(r, 1);
    }
    processed += (j - i);
    i = j;
  }

  const double n0      = (double)n * (double)(n - 1) / 2.0;
  const double denom_x = n0 - (double)T_x;
  const double denom_y = n0 - (double)T_y;
  if (denom_x <= 0.0 || denom_y <= 0.0) return NA_REAL;

  return (double)S / std::sqrt(denom_x * denom_y);
}

// [[Rcpp::export]]
NumericMatrix kendall_matrix_cpp(NumericMatrix mat) {
  int n = mat.ncol();
  NumericMatrix result(n, n);

  for (int i = 0; i < n; ++i) {
    result(i, i) = 1.0;
    for (int j = i + 1; j < n; ++j) {
      double tau = kendall_tau_auto_cpp(mat(_, i), mat(_, j));
      result(i, j) = tau;
      result(j, i) = tau;
    }
  }

  return result;
}
