// Thiago de Paula Oliveira
#include <RcppArmadillo.h>
#include "matrixCorr_omp.h"
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

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

  arma::mat D = Rcpp::as<arma::mat>(Dmat);  // Convert once

  // Reshape into ns × ntime matrices
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
  arma::vec within_quad(ns, arma::fill::zeros);
  for (int i = 0; i < ns; ++i) {
    arma::rowvec d = arma::abs(X.row(i) - Y.row(i));
    if (delta_zero) {
      d = arma::conv_to<arma::rowvec>::from(d != 0.0);
    } else {
      d = arma::pow(d, delta);
    }
    within_quad[i] = arma::as_scalar(d * D * d.t());
  }

  arma::vec phi1_sum(ns, arma::fill::zeros);
  arma::vec phi2_sum(ns, arma::fill::zeros);
  double U = 0.0, V = 0.0;
  
  // ===== OpenMP ordered-pair loop =====
#ifdef _OPENMP
#pragma omp parallel for reduction(+:U,V)
#endif
  for (int i = 0; i < ns; ++i) {
    const arma::rowvec Xi = X.row(i);
    const arma::rowvec Yi = Y.row(i);
    const double within_i = within_quad[i];
    double phi1_acc = 0.0;
    double phi2_acc = 0.0;

    for (int j = 0; j < ns; ++j) {
      if (i == j) continue;
      const arma::rowvec Xj = X.row(j);
      const arma::rowvec Yj = Y.row(j);

      const double phi1 = 0.5 * (within_i + within_quad[j]);

      arma::rowvec d3 = arma::abs(Xi - Yj);
      arma::rowvec d4 = arma::abs(Xj - Yi);

      double phi2;
      if (delta_zero) {
        d3 = arma::conv_to<arma::rowvec>::from(d3 != 0.0);
        d4 = arma::conv_to<arma::rowvec>::from(d4 != 0.0);
        phi2 = 0.5 * (arma::as_scalar(d3 * D * d3.t()) +
          arma::as_scalar(d4 * D * d4.t()));
      } else {
        arma::rowvec v3 = arma::pow(d3, delta);
        arma::rowvec v4 = arma::pow(d4, delta);
        phi2 = 0.5 * (arma::as_scalar(v3 * D * v3.t()) +
          arma::as_scalar(v4 * D * v4.t()));
      }

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

  // ===== Variance & CI =====
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

