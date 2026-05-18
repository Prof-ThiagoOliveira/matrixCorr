// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

struct SubjectPoissonBlock {
  double y1 = 0.0;
  double y2 = 0.0;
  double n1 = 0.0;
  double n2 = 0.0;
  double log_factorial = 0.0;
};

static std::vector<SubjectPoissonBlock> aggregate_subject_blocks(
    const Rcpp::NumericVector& y,
    const Rcpp::IntegerVector& subject,
    const Rcpp::IntegerVector& method_code,
    const int n_subjects
) {
  const int n = y.size();
  std::vector<SubjectPoissonBlock> blocks(n_subjects);

  for (int r = 0; r < n; ++r) {
    if (!R_finite(y[r]) || y[r] < 0.0 ||
        IntegerVector::is_na(subject[r]) ||
        IntegerVector::is_na(method_code[r]) ||
        subject[r] < 1 || subject[r] > n_subjects ||
        method_code[r] < 1 || method_code[r] > 2) {
      Rcpp::stop("Invalid Poisson GLMM likelihood inputs.");
    }

    SubjectPoissonBlock& block = blocks[subject[r] - 1];
    if (method_code[r] == 1) {
      block.y1 += y[r];
      block.n1 += 1.0;
    } else {
      block.y2 += y[r];
      block.n2 += 1.0;
    }
    block.log_factorial += std::lgamma(y[r] + 1.0);
  }

  return blocks;
}

static std::vector<SubjectPoissonBlock> subject_blocks_from_vectors(
    const Rcpp::NumericVector& y1,
    const Rcpp::NumericVector& y2,
    const Rcpp::NumericVector& n1,
    const Rcpp::NumericVector& n2,
    const Rcpp::NumericVector& log_factorial
) {
  const int n_subjects = y1.size();
  if (y2.size() != n_subjects ||
      n1.size() != n_subjects ||
      n2.size() != n_subjects ||
      log_factorial.size() != n_subjects) {
    Rcpp::stop("Invalid Poisson GLMM likelihood inputs.");
  }

  std::vector<SubjectPoissonBlock> blocks(n_subjects);
  for (int s = 0; s < n_subjects; ++s) {
    if (!R_finite(y1[s]) || !R_finite(y2[s]) ||
        !R_finite(n1[s]) || !R_finite(n2[s]) ||
        !R_finite(log_factorial[s])) {
      Rcpp::stop("Invalid Poisson GLMM likelihood inputs.");
    }
    blocks[s].y1 = y1[s];
    blocks[s].y2 = y2[s];
    blocks[s].n1 = n1[s];
    blocks[s].n2 = n2[s];
    blocks[s].log_factorial = log_factorial[s];
  }
  return blocks;
}

static double log_sum_exp(const std::vector<double>& values) {
  double max_log = values[0];
  for (std::size_t i = 1; i < values.size(); ++i) {
    if (values[i] > max_log) max_log = values[i];
  }

  double sum_exp = 0.0;
  for (double value : values) {
    sum_exp += std::exp(value - max_log);
  }

  if (!std::isfinite(sum_exp) || sum_exp <= 0.0) return R_NegInf;
  return max_log + std::log(sum_exp);
}

static double ccc_glmm_poisson_ghq_nll_blocks_impl(
    Rcpp::NumericVector par,
    const std::vector<SubjectPoissonBlock>& blocks,
    bool include_subject_method,
    Rcpp::NumericVector gh_nodes,
    Rcpp::NumericVector gh_weights
) {
  const int expected_par = include_subject_method ? 4 : 3;
  if (par.size() != expected_par) {
    Rcpp::stop("Invalid Poisson GLMM likelihood inputs.");
  }
  if (gh_nodes.size() != gh_weights.size() || gh_nodes.size() < 3) {
    Rcpp::stop("Gauss-Hermite nodes and weights must have matching length >= 3.");
  }

  const double beta0 = par[0];
  const double beta_method = par[1];
  const double sigma_subject = std::exp(par[2]);
  const double sigma_subject_method = include_subject_method ? std::exp(par[3]) : 0.0;

  if (!std::isfinite(beta0) ||
      !std::isfinite(beta_method) ||
      !std::isfinite(sigma_subject) ||
      sigma_subject <= 0.0 ||
      (include_subject_method &&
       (!std::isfinite(sigma_subject_method) || sigma_subject_method <= 0.0))) {
    return R_PosInf;
  }

  const int Q = gh_nodes.size();
  const int n_subjects = blocks.size();
  const double log_sqrt_pi = 0.5 * std::log(M_PI);
  std::vector<double> log_weights(Q);
  for (int q = 0; q < Q; ++q) {
    log_weights[q] = std::log(gh_weights[q]);
  }
  double total_loglik = 0.0;

  if (!include_subject_method) {
    std::vector<double> log_terms(Q);

    for (int s = 0; s < n_subjects; ++s) {
      const SubjectPoissonBlock& block = blocks[s];

      for (int q = 0; q < Q; ++q) {
        const double alpha = std::sqrt(2.0) * sigma_subject * gh_nodes[q];
        const double eta1 = beta0 + alpha;
        const double eta2 = beta0 + beta_method + alpha;
        if (eta1 > 700.0 || eta2 > 700.0) return R_PosInf;

        log_terms[q] =
            log_weights[q] - log_sqrt_pi +
            block.y1 * eta1 - block.n1 * std::exp(eta1) +
            block.y2 * eta2 - block.n2 * std::exp(eta2) -
            block.log_factorial;
      }

      const double subject_loglik = log_sum_exp(log_terms);
      if (!std::isfinite(subject_loglik)) return R_PosInf;
      total_loglik += subject_loglik;
    }

    return -total_loglik;
  }

  std::vector<double> log_terms(static_cast<std::size_t>(Q) * Q * Q);
  std::vector<double> tensor_log_weights(static_cast<std::size_t>(Q) * Q * Q);
  {
    std::size_t pos = 0;
    for (int qa = 0; qa < Q; ++qa) {
      for (int q1 = 0; q1 < Q; ++q1) {
        for (int q2 = 0; q2 < Q; ++q2) {
          tensor_log_weights[pos++] =
              log_weights[qa] + log_weights[q1] + log_weights[q2] -
              3.0 * log_sqrt_pi;
        }
      }
    }
  }

  for (int s = 0; s < n_subjects; ++s) {
    const SubjectPoissonBlock& block = blocks[s];
    std::size_t pos = 0;

    for (int qa = 0; qa < Q; ++qa) {
      const double alpha = std::sqrt(2.0) * sigma_subject * gh_nodes[qa];

      for (int q1 = 0; q1 < Q; ++q1) {
        const double gamma1 = std::sqrt(2.0) * sigma_subject_method * gh_nodes[q1];

        for (int q2 = 0; q2 < Q; ++q2) {
          const double gamma2 = std::sqrt(2.0) * sigma_subject_method * gh_nodes[q2];
          const double eta1 = beta0 + alpha + gamma1;
          const double eta2 = beta0 + beta_method + alpha + gamma2;
          if (eta1 > 700.0 || eta2 > 700.0) return R_PosInf;

          const std::size_t idx = pos++;
          log_terms[idx] =
              tensor_log_weights[idx] +
              block.y1 * eta1 - block.n1 * std::exp(eta1) +
              block.y2 * eta2 - block.n2 * std::exp(eta2) -
              block.log_factorial;
        }
      }
    }

    const double subject_loglik = log_sum_exp(log_terms);
    if (!std::isfinite(subject_loglik)) return R_PosInf;
    total_loglik += subject_loglik;
  }

  return -total_loglik;
}

// [[Rcpp::export]]
Rcpp::List ccc_glmm_poisson_prepare_blocks_cpp(
    Rcpp::NumericVector y,
    Rcpp::IntegerVector subject,
    Rcpp::IntegerVector method_code,
    int n_subjects
) {
  const int n = y.size();
  if (subject.size() != n || method_code.size() != n) {
    Rcpp::stop("Invalid Poisson GLMM likelihood inputs.");
  }
  if (n_subjects < 1) {
    Rcpp::stop("n_subjects must be positive.");
  }

  const std::vector<SubjectPoissonBlock> blocks =
      aggregate_subject_blocks(y, subject, method_code, n_subjects);
  Rcpp::NumericVector y1(n_subjects);
  Rcpp::NumericVector y2(n_subjects);
  Rcpp::NumericVector n1(n_subjects);
  Rcpp::NumericVector n2(n_subjects);
  Rcpp::NumericVector log_factorial(n_subjects);
  for (int s = 0; s < n_subjects; ++s) {
    y1[s] = blocks[s].y1;
    y2[s] = blocks[s].y2;
    n1[s] = blocks[s].n1;
    n2[s] = blocks[s].n2;
    log_factorial[s] = blocks[s].log_factorial;
  }

  return Rcpp::List::create(
    Rcpp::Named("y1") = y1,
    Rcpp::Named("y2") = y2,
    Rcpp::Named("n1") = n1,
    Rcpp::Named("n2") = n2,
    Rcpp::Named("log_factorial") = log_factorial
  );
}

// [[Rcpp::export]]
double ccc_glmm_poisson_ghq_nll_blocks_cpp(
    Rcpp::NumericVector par,
    Rcpp::NumericVector y1,
    Rcpp::NumericVector y2,
    Rcpp::NumericVector n1,
    Rcpp::NumericVector n2,
    Rcpp::NumericVector log_factorial,
    bool include_subject_method,
    Rcpp::NumericVector gh_nodes,
    Rcpp::NumericVector gh_weights
) {
  const std::vector<SubjectPoissonBlock> blocks =
      subject_blocks_from_vectors(y1, y2, n1, n2, log_factorial);
  return ccc_glmm_poisson_ghq_nll_blocks_impl(
    par,
    blocks,
    include_subject_method,
    gh_nodes,
    gh_weights
  );
}

// [[Rcpp::export]]
double ccc_glmm_poisson_ghq_nll_cpp(
    Rcpp::NumericVector par,
    Rcpp::NumericVector y,
    Rcpp::IntegerVector subject,
    Rcpp::IntegerVector method_code,
    int n_subjects,
    bool include_subject_method,
    Rcpp::NumericVector gh_nodes,
    Rcpp::NumericVector gh_weights
) {
  const int n = y.size();
  const int expected_par = include_subject_method ? 4 : 3;
  if (par.size() != expected_par || subject.size() != n || method_code.size() != n) {
    Rcpp::stop("Invalid Poisson GLMM likelihood inputs.");
  }
  if (n_subjects < 1) {
    Rcpp::stop("n_subjects must be positive.");
  }

  const std::vector<SubjectPoissonBlock> blocks =
      aggregate_subject_blocks(y, subject, method_code, n_subjects);
  return ccc_glmm_poisson_ghq_nll_blocks_impl(
    par,
    blocks,
    include_subject_method,
    gh_nodes,
    gh_weights
  );
}
