// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include "matrixCorr_omp.h"
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

namespace {

inline bool is_diagonal_matrix(const arma::mat& M) {
  const arma::uword nr = M.n_rows;
  const arma::uword nc = M.n_cols;
  if (nr != nc) return false;
  for (arma::uword j = 0; j < nc; ++j) {
    for (arma::uword i = 0; i < nr; ++i) {
      if (i != j && M(i, j) != 0.0) return false;
    }
  }
  return true;
}

inline void fill_abs_delta_diff(const arma::mat& A,
                                const int row_a,
                                const arma::mat& B,
                                const int row_b,
                                const int ntime,
                                const bool delta_zero,
                                const bool delta_one,
                                const bool delta_two,
                                const double delta,
                                std::vector<double>& out) {
  for (int t = 0; t < ntime; ++t) {
    const double d = std::fabs(A(row_a, t) - B(row_b, t));
    if (delta_zero) {
      out[static_cast<std::size_t>(t)] = (d != 0.0) ? 1.0 : 0.0;
    } else if (delta_one) {
      out[static_cast<std::size_t>(t)] = d;
    } else if (delta_two) {
      out[static_cast<std::size_t>(t)] = d * d;
    } else {
      out[static_cast<std::size_t>(t)] = std::pow(d, delta);
    }
  }
}

inline double qform_diag(const std::vector<double>& v,
                         const std::vector<double>& diag_w) {
  const std::size_t n = v.size();
  double acc = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double vi = v[i];
    acc += diag_w[i] * vi * vi;
  }
  return acc;
}

inline double qform_dense(const std::vector<double>& v,
                          const arma::mat& D,
                          const int ntime) {
  double acc = 0.0;
  for (int r = 0; r < ntime; ++r) {
    const double vr = v[static_cast<std::size_t>(r)];
    double row_acc = 0.0;
    for (int c = 0; c < ntime; ++c) {
      row_acc += D(r, c) * v[static_cast<std::size_t>(c)];
    }
    acc += vr * row_acc;
  }
  return acc;
}

} // namespace

// [[Rcpp::export]]
List cccUst_rcpp(NumericVector y_vec,
                 IntegerVector met_vec,
                 IntegerVector time_vec,
                 IntegerVector subj_vec,
                 int nmet0,
                 int nmet1,
                 int ntime,
                 int ns,
                 NumericMatrix Dmat,
                 double delta,
                 double cl) {

  arma::mat D = Rcpp::as<arma::mat>(Dmat); // Convert once
  const bool D_is_diag = is_diagonal_matrix(D);
  std::vector<double> D_diag;
  if (D_is_diag) {
    D_diag.resize(static_cast<std::size_t>(ntime));
    for (int t = 0; t < ntime; ++t) {
      D_diag[static_cast<std::size_t>(t)] = D(t, t);
    }
  }

  // Reshape into ns x ntime matrices
  if (y_vec.size() != met_vec.size() ||
      y_vec.size() != time_vec.size() ||
      y_vec.size() != subj_vec.size()) {
    Rcpp::stop("Input vectors must have the same length");
  }

  if (ns < 2) {
    Rcpp::stop("At least two subjects are required to compute repeated-measures CCC");
  }

  arma::mat X(ns, ntime, arma::fill::zeros);
  arma::mat Y(ns, ntime, arma::fill::zeros);
  arma::umat seenX(ns, ntime, arma::fill::zeros);
  arma::umat seenY(ns, ntime, arma::fill::zeros);

  for (int i = 0; i < y_vec.size(); ++i) {
    int mi = met_vec[i];
    int ti = time_vec[i];
    int si = subj_vec[i];
    if (ti < 0 || ti >= ntime) {
      Rcpp::stop("time index out of bounds");
    }
    if (si < 0 || si >= ns) {
      Rcpp::stop("subject index out of bounds");
    }
    if (mi == nmet0) {
      if (seenX(si, ti)) Rcpp::stop("duplicate observations for method 0");
      X(si, ti) = y_vec[i];
      seenX(si, ti) = 1u;
    } else if (mi == nmet1) {
      if (seenY(si, ti)) Rcpp::stop("duplicate observations for method 1");
      Y(si, ti) = y_vec[i];
      seenY(si, ti) = 1u;
    }
  }

  if (arma::accu(seenX) != static_cast<unsigned int>(ns * ntime) ||
      arma::accu(seenY) != static_cast<unsigned int>(ns * ntime)) {
    Rcpp::stop("Missing observations detected for at least one subject/time cell");
  }

  const bool delta_zero = (delta == 0.0);
  const bool delta_one = (!delta_zero && delta == 1.0);
  const bool delta_two = (!delta_zero && delta == 2.0);

  arma::vec within_quad(ns, arma::fill::zeros);
  std::vector<double> diff_within(static_cast<std::size_t>(ntime));
  for (int i = 0; i < ns; ++i) {
    fill_abs_delta_diff(
      X, i, Y, i, ntime,
      delta_zero, delta_one, delta_two, delta,
      diff_within
    );
    within_quad[i] = D_is_diag ? qform_diag(diff_within, D_diag) : qform_dense(diff_within, D, ntime);
  }

  arma::vec phi1_sum(ns, arma::fill::zeros);
  arma::vec phi2_sum(ns, arma::fill::zeros);
  double U = 0.0, V = 0.0;

  // OpenMP ordered-pair loop
#ifdef _OPENMP
#pragma omp parallel for reduction(+:U,V)
#endif
  for (int i = 0; i < ns; ++i) {
    std::vector<double> d3(static_cast<std::size_t>(ntime));
    std::vector<double> d4(static_cast<std::size_t>(ntime));

    const double within_i = within_quad[i];
    double phi1_acc = 0.0;
    double phi2_acc = 0.0;

    for (int j = 0; j < ns; ++j) {
      if (i == j) continue;

      const double phi1 = 0.5 * (within_i + within_quad[j]);

      fill_abs_delta_diff(
        X, i, Y, j, ntime,
        delta_zero, delta_one, delta_two, delta,
        d3
      );
      fill_abs_delta_diff(
        X, j, Y, i, ntime,
        delta_zero, delta_one, delta_two, delta,
        d4
      );

      const double phi2 = D_is_diag
        ? 0.5 * (qform_diag(d3, D_diag) + qform_diag(d4, D_diag))
        : 0.5 * (qform_dense(d3, D, ntime) + qform_dense(d4, D, ntime));

      phi1_acc += phi1;
      phi2_acc += phi2;
      U += phi1;
      V += phi2;
    }

    phi1_sum[i] = phi1_acc / (ns - 1);
    phi2_sum[i] = phi2_acc / (ns - 1);
  }

  // Mean U and V
  double denom = ns * (ns - 1);
  U /= denom;
  V /= denom;

  double CCC = ((ns - 1) * (V - U)) / (U + (ns - 1) * V);

  // Variance and CI
  double s11 = 0.0;
  double s12 = 0.0;
  double s22 = 0.0;
  for (int i = 0; i < ns; ++i) {
    const double d1 = phi1_sum[i] - U;
    const double d2 = phi2_sum[i] - V;
    s11 += d1 * d1;
    s12 += d1 * d2;
    s22 += d2 * d2;
  }

  arma::mat Saux(2, 2, arma::fill::zeros);
  Saux(0, 0) = s11;
  Saux(0, 1) = s12;
  Saux(1, 0) = s12;
  Saux(1, 1) = s22;

  arma::mat C; C.eye(2, 2); C(1, 1) = 2;
  arma::mat Smat = C * (Saux / (ns * ns)) * C;

  arma::rowvec dev(2);
  dev[0] = (-ns * (ns - 1) * V) / std::pow(U + (ns - 1) * V, 2);
  dev[1] = (ns * (ns - 1) * U) / std::pow(U + (ns - 1) * V, 2);

  double VarCCC = arma::as_scalar(dev * Smat * dev.t());

  // Fisher Z-based confidence interval
  double alpha = 1.0 - cl;
  double z = 0.5 * std::log((1 + CCC) / (1 - CCC));
  double VarZ = VarCCC / std::pow(1 - CCC * CCC, 2);
  double seZ = std::sqrt(VarZ);

  double zcrit = R::qnorm(1.0 - alpha / 2.0, 0.0, 1.0, 1, 0);
  double z_low = z - zcrit * seZ;
  double z_upp = z + zcrit * seZ;

  double ci_low = (std::exp(2 * z_low) - 1) / (std::exp(2 * z_low) + 1);
  double ci_upp = (std::exp(2 * z_upp) - 1) / (std::exp(2 * z_upp) + 1);

  return List::create(
    Named("CCC")     = CCC,
    Named("LL_CI")   = ci_low,
    Named("UL_CI")   = ci_upp,
    Named("SE_CCC")  = std::sqrt(VarCCC),
    Named("Z")       = z,
    Named("SE_Z")    = seZ
  );
}

// [[Rcpp::export]]
void set_omp_threads(const int n) {
#ifdef _OPENMP
  omp_set_num_threads(n);
#endif
}

// [[Rcpp::export]]
int get_omp_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}
