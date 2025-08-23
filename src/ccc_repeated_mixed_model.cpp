// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

//---------------- helpers ----------------//
static inline mat solve_sympd_safe(const mat& A, const mat& B) {
  mat Z;
  if (solve(Z, A, B, solve_opts::likely_sympd)) return Z;
  mat Aj = A;
  double jitter = 1e-8;
  if (A.n_rows > 0) {
    double tr = trace(A);
    if (std::isfinite(tr) && tr > 0.0) jitter = std::max(1e-12, 1e-8 * tr / A.n_rows);
  }
  Aj.diag() += jitter;
  if (solve(Z, Aj, B, solve_opts::likely_sympd)) return Z;
  return pinv(A) * B;
}

static inline bool inv_sympd_safe(mat& out, const mat& A) {
  if (inv_sympd(out, A)) return true;
  mat Aj = A;
  double jitter = 1e-8;
  if (A.n_rows > 0) {
    double tr = trace(A);
    if (std::isfinite(tr) && tr > 0.0) jitter = std::max(1e-12, 1e-8 * tr / A.n_rows);
  }
  Aj.diag() += jitter;
  if (inv_sympd(out, Aj)) return true;
  out = pinv(A);
  return out.is_finite();
}

// unbiased sample variance (n-1 in denom); returns 0 if <2
static inline double sample_var(const arma::vec& v) {
  const arma::uword n = v.n_elem;
  if (n < 2u) return 0.0;
  const double mu = arma::mean(v);
  double acc = 0.0;
  for (arma::uword i=0;i<n;++i) { double d = v[i]-mu; acc += d*d; }
  return acc / (double)(n-1);
}

//---------------- reindex subjects ----------------//
static inline void reindex_subject(const IntegerVector& subject,
                                   std::vector<int>& subj_idx,
                                   int& m_out) {
  subj_idx.resize(subject.size());
  std::unordered_map<int,int> map;
  map.reserve(subject.size());
  int next = 0;
  for (int i = 0; i < subject.size(); ++i) {
    int s = subject[i];
    if (s == NA_INTEGER) stop("subject contains NA");
    auto it = map.find(s);
    if (it == map.end()) { map.emplace(s,next); subj_idx[i]=next; ++next; }
    else subj_idx[i] = it->second;
  }
  m_out = next;
}

//---------------- per-subject indexing ----------------//
struct BySubj {
  std::vector< std::vector<int> > rows; // row ids per subject (global ids)
  std::vector< std::vector<int> > met;  // 0..nm-1 (or -1 for NA)
  std::vector< std::vector<int> > tim;  // 0..nt-1 (or -1 for NA)
};

static inline BySubj index_by_subject(const std::vector<int>& subj_idx,
                                      const IntegerVector& method,
                                      const IntegerVector& time,
                                      int m) {
  BySubj S;
  S.rows.assign(m, {});
  S.met.assign(m, {});
  S.tim.assign(m, {});
  const int n = subj_idx.size();
  for (int i = 0; i < n; ++i) {
    int j = subj_idx[i];
    S.rows[j].push_back(i);
    if (method.size() > 0) {
      int v = method[i];
      S.met[j].push_back(v == NA_INTEGER ? -1 : (v - 1));
    }
    if (time.size() > 0) {
      int v = time[i];
      S.tim[j].push_back(v == NA_INTEGER ? -1 : (v - 1));
    }
  }
  return S;
}

//---------------- build U'U for one subject ----------------//
static inline void make_UtU(const std::vector<int>& met_i,
                            const std::vector<int>& tim_i,
                            int n_i, int nm, int nt,
                            mat& UtU) {
  const int r = 1 + (nm>0?nm:0) + (nt>0?nt:0);
  UtU.zeros(r,r);
  UtU(0,0) = n_i;

  if (nm > 0) {
    vec cm(nm, fill::zeros);
    for (int v : met_i) if (v >= 0) cm[v] += 1.0;
    for (int l=0; l<nm; ++l) {
      UtU(0,1+l) = UtU(1+l,0) = cm[l];
      UtU(1+l,1+l) = cm[l];
    }
    if (nt > 0) {
      vec ct(nt, fill::zeros);
      for (int v : tim_i) if (v >= 0) ct[v] += 1.0;
      for (int t=0; t<nt; ++t) {
        int jt = 1 + nm + t;
        UtU(0,jt) = UtU(jt,0) = ct[t];
        UtU(jt,jt) = ct[t];
      }
      mat cmt(nm, nt, fill::zeros);
      const int n_i_rows = (int)met_i.size();
      for (int k=0; k<n_i_rows; ++k) {
        int l = met_i[k], t = tim_i[k];
        if (l>=0 && t>=0) cmt(l,t) += 1.0;
      }
      for (int l=0; l<nm; ++l) for (int t=0; t<nt; ++t) {
        int il = 1 + l;
        int jt = 1 + nm + t;
        UtU(il,jt) = UtU(jt,il) = cmt(l,t);
      }
    }
  } else if (nt > 0) {
    vec ct(nt, fill::zeros);
    for (int v : tim_i) if (v >= 0) ct[v] += 1.0;
    for (int t=0; t<nt; ++t) {
      int jt = 1 + t;
      UtU(0,jt) = UtU(jt,0) = ct[t];
      UtU(jt,jt) = ct[t];
    }
  }
}

//---------------- accumulate U^T * v for one subject ----------------//
template <typename Accessor>
static inline void accum_Ut_vec(const std::vector<int>& rows_i,
                                const std::vector<int>& met_i,
                                const std::vector<int>& tim_i,
                                int nm, int nt,
                                Accessor v, vec& Utv) {
  Utv.zeros();

  // intercept
  double s0 = 0.0;
  for (int ridx : rows_i) s0 += v(ridx);
  Utv[0] = s0;

  // method sums
  if (nm > 0) {
    for (int l = 0; l < nm; ++l) {
      double sm = 0.0;
      for (size_t k = 0; k < rows_i.size(); ++k) {
        if (met_i[k] == l) sm += v(rows_i[k]);
      }
      Utv[1 + l] = sm;
    }
  }
  // time sums
  if (nt > 0) {
    for (int t = 0; t < nt; ++t) {
      double st = 0.0;
      for (size_t k = 0; k < rows_i.size(); ++k) {
        if (tim_i[k] == t) st += v(rows_i[k]);
      }
      int jt = 1 + (nm > 0 ? nm : 0) + t;
      Utv[jt] = st;
    }
  }

  const int r_expected = 1 + (nm > 0 ? nm : 0) + (nt > 0 ? nt : 0);
  if ((int)Utv.n_rows != r_expected)
    Rcpp::stop("accum_Ut_vec: Utv length mismatch (got %d, expected %d)",
               (int)Utv.n_rows, r_expected);
}

//---------------- Ua builder (fixed) ----------------//
static inline void add_U_times(const std::vector<int>& rows_i,
                               const std::vector<int>& met_i,
                               const std::vector<int>& tim_i,
                               int nm, int nt,
                               const vec& a, vec& out /* length n_i */) {
  const int n_i = (int)rows_i.size();
  const double a0 = a[0];

  std::vector<double> am(nm>0?nm:0), at(nt>0?nt:0);
  if (nm>0) for (int l=0;l<nm;++l) am[l] = a[1+l];
  if (nt>0) for (int t=0;t<nt;++t) at[t] = a[1 + (nm>0?nm:0) + t];

  for (int k=0; k<n_i; ++k) {
    double val = a0;
    if (nm>0 && met_i[k] >= 0) val += am[ met_i[k] ];
    if (nt>0 && tim_i[k] >= 0) val += at[ tim_i[k] ];
    out[k] += val;
  }
}

//---------------- caches ----------------//
struct Cache {
  arma::mat UtU;   // r x r
  arma::vec Uty;   // r
  arma::mat Utx;   // r x p
  arma::vec Xty;   // p
  arma::mat XtX;   // p x p
  int n_i;
};

//---------------- main entry ----------------//
// [[Rcpp::export]]
Rcpp::List ccc_vc_cpp(
    Rcpp::NumericMatrix Xr,
    Rcpp::NumericVector yr,
    Rcpp::IntegerVector subject,
    Rcpp::IntegerVector method,
    Rcpp::IntegerVector time,
    int nm, int nt,
    int max_iter = 200,
    double tol = 1e-6,
    double conf_level = 0.95,
    Rcpp::Nullable<Rcpp::NumericMatrix> Lr = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> auxDr = R_NilValue
) {
  const int n = yr.size();
  if (Xr.nrow() != n) stop("nrow(X) must match length(y)");
  if (subject.size() != n) stop("length(subject) mismatch");
  if (method.size()  && method.size()!=n)  stop("length(method) mismatch");
  if (time.size()    && time.size()!=n)    stop("length(time) mismatch");

  mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
  vec y(yr.begin(), yr.size(), false);
  const int p = X.n_cols;

  // subjects
  std::vector<int> subj_idx; int m = 0;
  reindex_subject(subject, subj_idx, m);
  BySubj S = index_by_subject(subj_idx, method, time, m);

  // EM init
  double sa  = 1.0;
  double sab = (nm>0 ? 0.5 : 0.0);
  double sag = (nt>0 ? 0.5 : 0.0);
  double se  = 1.0;
  const double eps = 1e-10;

  vec beta(p, fill::zeros);

  // workspace sizes
  const int r = 1 + (nm>0?nm:0) + (nt>0?nt:0);

  // -------------------- precompute per-subject invariants --------------------
  std::vector<Cache> C(m);
  for (int i=0; i<m; ++i) {
    const auto& rows = S.rows[i];
    const auto& met  = S.met[i];
    const auto& tim  = S.tim[i];
    const int   n_i  = (int)rows.size();
    C[i].n_i = n_i;

    // UtU
    C[i].UtU.set_size(r,r);
    make_UtU(met, tim, n_i, nm, nt, C[i].UtU);

    // U^T y
    C[i].Uty.set_size(r);
    accum_Ut_vec(rows, met, tim, nm, nt, [&](int idx){ return y[idx]; }, C[i].Uty);

    // U^T X
    C[i].Utx.set_size(r, p);
    for (int k=0;k<p;++k) {
      vec tmp(r, fill::zeros);
      accum_Ut_vec(rows, met, tim, nm, nt, [&](int idx){ return X(idx,k); }, tmp);
      C[i].Utx.col(k) = tmp;
    }

    // X^T y and X^T X (symmetric)
    C[i].Xty.set_size(p);
    C[i].XtX.zeros(p,p);
    for (int k=0; k<p; ++k) {
      double sxy = 0.0;
      for (int idx : rows) sxy += X(idx,k) * y[idx];
      C[i].Xty[k] = sxy;
      for (int l=k; l<p; ++l) {
        double sxx = 0.0;
        for (int idx : rows) sxx += X(idx,k) * X(idx,l);
        C[i].XtX(k,l) = sxx;
        if (l!=k) C[i].XtX(l,k) = sxx;
      }
    }
  }

  // for delta-method (filled while running)
  std::vector<double> sa_term(m, 0.0);
  std::vector<double> sab_term(m, 0.0);
  std::vector<double> sag_term(m, 0.0);
  std::vector<double> se_term(m, 0.0);
  std::vector<int>    n_per_subj(m, 0);

  // ---------------------------- EM iterations -------------------------------
  for (int iter=0; iter<max_iter; ++iter) {

    // (1) Assemble XtViX, XtViy using rank-r corrections (parallel)
    mat XtViX(p,p,fill::zeros);
    vec XtViy(p, fill::zeros);
    const double inv_se = 1.0 / std::max(se, eps);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    std::vector<mat> XtViX_tls(nthreads, mat(p,p, fill::zeros));
    std::vector<vec> XtViy_tls(nthreads, vec(p, fill::zeros));

#pragma omp parallel
{
  int tid = omp_get_thread_num();
  mat& XtViX_loc = XtViX_tls[tid];
  vec& XtViy_loc = XtViy_tls[tid];

  // thread-local temporaries
  mat M(r,r), A(r, 1+p), Z(r, 1+p), TX;
  vec Z_y;

#pragma omp for schedule(static)
  for (int i=0; i<m; ++i) {
    const Cache& Ci = C[i];

    // M = G^{-1} + inv_se * UtU
    M.zeros();
    M(0,0) = 1.0 / std::max(sa, eps);
    int off = 1;
    if (nm>0) { for (int l=0;l<nm;++l) M(off+l, off+l) = 1.0 / std::max(sab, eps); off += nm; }
    if (nt>0) { for (int t=0;t<nt;++t) M(off+t, off+t) = 1.0 / std::max(sag, eps); }
    M += inv_se * Ci.UtU;

    // A = [U^T y, U^T X]
    A.col(0)      = Ci.Uty;
    A.cols(1,p)   = Ci.Utx;

    // Z = M^{-1} * (A / se)
    Z = solve_sympd_safe(M, A * inv_se);

    TX   = A.cols(1,p);
    Z_y  = Z.col(0);
    mat Z_X = Z.cols(1,p);

    for (int k=0; k<p; ++k) {
      XtViy_loc[k] += inv_se * (Ci.Xty[k] - dot(TX.col(k), Z_y));
      for (int l=k; l<p; ++l) {
        double val = inv_se * (Ci.XtX(k,l) - dot(TX.col(k), Z_X.col(l)));
        XtViX_loc(k,l) += val;
        if (l!=k) XtViX_loc(l,k) += val;
      }
    }
  } // for i
}   // parallel

// reduce tls -> global
for (int t=0; t<(int)XtViX_tls.size(); ++t) {
  XtViX += XtViX_tls[t];
  XtViy += XtViy_tls[t];
}
#else
// serial fallback
mat M(r,r), A(r, 1+p), Z(r, 1+p);
for (int i=0; i<m; ++i) {
  const Cache& Ci = C[i];
  M.zeros();
  M(0,0) = 1.0 / std::max(sa, eps);
  int off = 1;
  if (nm>0) { for (int l=0;l<nm;++l) M(off+l, off+l) = 1.0 / std::max(sab, eps); off += nm; }
  if (nt>0) { for (int t=0;t<nt;++t) M(off+t, off+t) = 1.0 / std::max(sag, eps); }
  M += inv_se * Ci.UtU;

  A.col(0)    = Ci.Uty;
  A.cols(1,p) = Ci.Utx;
  Z = solve_sympd_safe(M, A * inv_se);

  mat TX  = A.cols(1,p);
  vec Z_y = Z.col(0);
  mat Z_X = Z.cols(1,p);

  for (int k=0; k<p; ++k) {
    XtViy[k] += inv_se * (Ci.Xty[k] - dot(TX.col(k), Z_y));
    for (int l=k; l<p; ++l) {
      double val = inv_se * (Ci.XtX(k,l) - dot(TX.col(k), Z_X.col(l)));
      XtViX(k,l) += val; if (l!=k) XtViX(l,k) += val;
    }
  }
}
#endif

// (2) GLS beta
{
  vec btmp;
  if (!solve(btmp, XtViX, XtViy, solve_opts::likely_sympd)) {
    mat XtViXj = XtViX;
    // slightly stronger ridge helps avoid pinv on ill-conditioned designs
    XtViXj.diag() += 1e-8;
    if (!solve(btmp, XtViXj, XtViy, solve_opts::likely_sympd))
      btmp = pinv(XtViX) * XtViy;
  }
  beta = btmp;
}

// (3) M-step-ish (parallel): update sa, sab, sag, se pieces
double sa_acc = 0.0, sab_acc = 0.0, sag_acc = 0.0;
double se_sumsq = 0.0, se_trace = 0.0;

// global residuals once (BLAS gemv)
vec r_global = y - X * beta;

#ifdef _OPENMP
int nthreads2 = omp_get_max_threads();
std::vector<double> sa_tls(nthreads2,0.0), sab_tls(nthreads2,0.0),
sag_tls(nthreads2,0.0), ss_tls(nthreads2,0.0),
tr_tls(nthreads2,0.0);

#pragma omp parallel
{
  int tid = omp_get_thread_num();
  // thread-local temporaries
  mat M(r,r), Minv;
  vec Utr(r), b_i;
#pragma omp for schedule(static)
  for (int i=0; i<m; ++i) {
    const Cache& Ci = C[i];
    const auto& rows_i = S.rows[i];
    const auto& met_i  = S.met[i];
    const auto& tim_i  = S.tim[i];
    const int n_i = Ci.n_i;

    // M = G^{-1} + (1/se) * UtU
    M.zeros();
    M(0,0) = 1.0 / std::max(sa, eps);
    int off = 1;
    if (nm>0) { for (int l=0;l<nm;++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
    if (nt>0) { for (int t=0;t<nt;++t) M(off+t, off+t) = 1.0/std::max(sag,eps); }
    M += (1.0/std::max(se,eps)) * Ci.UtU;

    // U^T r_i = U^T y_i - U^T X_i * beta
    Utr = Ci.Uty - Ci.Utx * beta;

    // b_i and Var(b_i|y)
    b_i  = solve_sympd_safe(M, Utr / std::max(se,eps));
    inv_sympd_safe(Minv, M);

    // r_i from r_global
    vec r_i(n_i);
    for (int t=0; t<n_i; ++t) r_i[t] = r_global[ rows_i[t] ];

    // residual part: r_i - U b_i
    vec Ub(n_i, fill::zeros);
    add_U_times(rows_i, met_i, tim_i, nm, nt, b_i, Ub);
    double ss = 0.0; for (int t=0;t<n_i;++t) { double e = r_i[t] - Ub[t]; ss += e*e; }
    ss_tls[tid] += ss;
    tr_tls[tid] += trace(Minv * Ci.UtU);

    // accumulate E[b^2] (means)
    sa_tls[tid] += b_i[0]*b_i[0] + Minv(0,0);
    int pos = 1;
    if (nm>0) { for (int l=0;l<nm;++l) sab_tls[tid] += b_i[pos+l]*b_i[pos+l] + Minv(pos+l, pos+l); pos += nm; }
    if (nt>0) { for (int t=0;t<nt;++t)  sag_tls[tid] += b_i[pos+t]*b_i[pos+t] + Minv(pos+t, pos+t); }

    // store per-subject contributions for delta method
    n_per_subj[i] = n_i;

    sa_term[i] = b_i[0]*b_i[0] + Minv(0,0);

    int pos_cov = 1;
    if (nm > 0) {
      double acc_m = 0.0;
      for (int l = 0; l < nm; ++l)
        acc_m += b_i[pos_cov + l]*b_i[pos_cov + l] + Minv(pos_cov + l, pos_cov + l);
      sab_term[i] = acc_m / (double)nm;
      pos_cov += nm;
    } else {
      sab_term[i] = 0.0;
    }

    if (nt > 0) {
      double acc_t = 0.0;
      for (int t = 0; t < nt; ++t)
        acc_t += b_i[pos_cov + t]*b_i[pos_cov + t] + Minv(pos_cov + t, pos_cov + t);
      sag_term[i] = acc_t / (double)nt;
    } else {
      sag_term[i] = 0.0;
    }

    double Si = ss + trace(Minv * Ci.UtU);
    se_term[i] = Si / (double)n_i;
  }
} // parallel

// reductions
for (int t=0; t<nthreads2; ++t) {
  sa_acc  += sa_tls[t];
  sab_acc += sab_tls[t];
  sag_acc += sag_tls[t];
  se_sumsq+= ss_tls[t];
  se_trace+= tr_tls[t];
}
#else
// serial fallback
mat M(r,r), Minv;
vec Utr(r), b_i;
for (int i=0; i<m; ++i) {
  const Cache& Ci = C[i];
  const auto& rows_i = S.rows[i];
  const auto& met_i  = S.met[i];
  const auto& tim_i  = S.tim[i];
  const int n_i = Ci.n_i;

  M.zeros();
  M(0,0) = 1.0 / std::max(sa, eps);
  int off = 1;
  if (nm>0) { for (int l=0;l<nm;++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
  if (nt>0) { for (int t=0;t<nt;++t) M(off+t, off+t) = 1.0/std::max(sag,eps); }
  M += (1.0/std::max(se,eps)) * Ci.UtU;

  Utr = Ci.Uty - Ci.Utx * beta;
  b_i = solve_sympd_safe(M, Utr / std::max(se,eps));
  inv_sympd_safe(Minv, M);

  vec r_i(n_i);
  for (int t=0; t<n_i; ++t) r_i[t] = r_global[ rows_i[t] ];

  vec Ub(n_i, fill::zeros);
  add_U_times(rows_i, met_i, tim_i, nm, nt, b_i, Ub);
  double ss = 0.0; for (int t=0;t<n_i;++t) { double e = r_i[t] - Ub[t]; ss += e*e; }
  se_sumsq += ss;
  se_trace += trace(Minv * Ci.UtU);

  sa_acc += b_i[0]*b_i[0] + Minv(0,0);
  int pos = 1;
  if (nm>0) { for (int l=0;l<nm;++l) sab_acc += b_i[pos+l]*b_i[pos+l] + Minv(pos+l, pos+l); pos += nm; }
  if (nt>0) { for (int t=0;t<nt;++t)  sag_acc += b_i[pos+t]*b_i[pos+t] + Minv(pos+t, pos+t); }

  n_per_subj[i] = n_i;

  sa_term[i] = b_i[0]*b_i[0] + Minv(0,0);

  int pos_cov = 1;
  if (nm > 0) {
    double acc_m = 0.0;
    for (int l = 0; l < nm; ++l)
      acc_m += b_i[pos_cov + l]*b_i[pos_cov + l] + Minv(pos_cov + l, pos_cov + l);
    sab_term[i] = acc_m / (double)nm;
    pos_cov += nm;
  } else {
    sab_term[i] = 0.0;
  }

  if (nt > 0) {
    double acc_t = 0.0;
    for (int t = 0; t < nt; ++t)
      acc_t += b_i[pos_cov + t]*b_i[pos_cov + t] + Minv(pos_cov + t, pos_cov + t);
    sag_term[i] = acc_t / (double)nt;
  } else {
    sag_term[i] = 0.0;
  }

  double Si = ss + trace(Minv * Ci.UtU);
  se_term[i] = Si / (double)n_i;
}
#endif

double sa_new  = std::max(sa_acc  / (double)m, eps);
double sab_new = (nm>0 ? std::max(sab_acc / (double)(m*nm), eps) : 0.0);
double sag_new = (nt>0 ? std::max(sag_acc / (double)(m*nt), eps) : 0.0);
double se_new  = std::max( (se_sumsq + se_trace) / (double)n, eps );

double diff = std::fabs(sa_new-sa) + std::fabs(se_new-se)
  + std::fabs(sab_new-sab) + std::fabs(sag_new-sag);
sa=sa_new; se=se_new; sab=sab_new; sag=sag_new;
if (diff < tol) break;
  }

  // ---------------- VarFix using caches (parallel) ----------------
  mat XtViX(p,p,fill::zeros);
  const double inv_se_final = 1.0 / std::max(se, eps);

#ifdef _OPENMP
  int nthreads3 = omp_get_max_threads();
  std::vector<mat> XtViX_tls2(nthreads3, mat(p,p, fill::zeros));
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  mat& XtViX_loc = XtViX_tls2[tid];

  mat M(r,r), Zx(r,p);

#pragma omp for schedule(static)
  for (int i=0; i<m; ++i) {
    const Cache& Ci = C[i];
    M.zeros();
    M(0,0) = 1.0 / std::max(sa, eps);
    int off = 1;
    if (nm>0) { for (int l=0;l<nm;++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
    if (nt>0) { for (int t=0;t<nt;++t) M(off+t, off+t) = 1.0/std::max(sag,eps); }
    M += inv_se_final * Ci.UtU;

    Zx = solve_sympd_safe(M, Ci.Utx * inv_se_final);

    for (int k=0;k<p;++k) for (int l=k;l<p;++l) {
      double val = inv_se_final * (Ci.XtX(k,l) - dot(Ci.Utx.col(k), Zx.col(l)));
      XtViX_loc(k,l) += val; if (l!=k) XtViX_loc(l,k) += val;
    }
  }
}
for (int t=0; t<nthreads3; ++t) XtViX += XtViX_tls2[t];
#else
{
  mat M(r,r), Zx(r,p);
  for (int i=0; i<m; ++i) {
    const Cache& Ci = C[i];
    M.zeros();
    M(0,0) = 1.0 / std::max(sa, eps);
    int off = 1;
    if (nm>0) { for (int l=0;l<nm;++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
    if (nt>0) { for (int t=0;t<nt;++t) M(off+t, off+t) = 1.0/std::max(sag,eps); }
    M += inv_se_final * Ci.UtU;

    Zx = solve_sympd_safe(M, Ci.Utx * inv_se_final);

    for (int k=0;k<p;++k) for (int l=k;l<p;++l) {
      double val = inv_se_final * (Ci.XtX(k,l) - dot(Ci.Utx.col(k), Zx.col(l)));
      XtViX(k,l) += val; if (l!=k) XtViX(l,k) += val;
    }
  }
}
#endif

mat VarFix;
if (!inv_sympd_safe(VarFix, XtViX)) stop("Failed to invert XtViX.");

// ---------------- SB and CCC ----------------
double SB = 0.0, varSB = 0.0;
if (nm > 0 && Lr.isNotNull() && auxDr.isNotNull()) {
  NumericMatrix Lrm = as<NumericMatrix>(Lr);
  NumericMatrix Drm = as<NumericMatrix>(auxDr);
  mat L(Lrm.begin(), X.n_cols, Lrm.ncol(), false);
  mat auxD(Drm.begin(), Drm.nrow(), Drm.ncol(), false);
  const double den = (double)nm * (double)(nm-1) * (double)std::max(nt,1);

  vec difmed = L.t() * beta;
  mat Afix   = L * auxD * L.t();
  double num = as_scalar(difmed.t() * auxD * difmed) - trace(Afix * VarFix);
  SB = std::max(num / den, 0.0);

  mat AV = Afix * VarFix;
  double term1 = 2.0 * trace(AV * AV);
  double term2 = 4.0 * as_scalar(beta.t() * Afix * VarFix * Afix * beta);
  varSB = std::max((term1 + term2) / (den * den), 0.0);
}

const double ccc = (sa + sag) / (sa + sab + sag + SB + se);

// ---------------- delta-method SE & CI for CCC ----------------
double Nnum = sa + sag;
double Dden = sa + sab + sag + SB + se;
if (Dden < 1e-14) Dden = 1e-14;

const double d_sa  = (sab + SB + se) / (Dden * Dden);
const double d_sab = -Nnum / (Dden * Dden);
const double d_sag = (sab + SB + se) / (Dden * Dden);
const double d_se  = -Nnum / (Dden * Dden);
const double d_SB  = -Nnum / (Dden * Dden);

arma::mat Z;
if (nm > 0 && nt > 0) {
  Z.set_size(m, 3);
  for (int i = 0; i < m; ++i) { Z(i,0)=sa_term[i]; Z(i,1)=sab_term[i]; Z(i,2)=sag_term[i]; }
} else if (nm > 0 && nt == 0) {
  Z.set_size(m, 2);
  for (int i = 0; i < m; ++i) { Z(i,0)=sa_term[i]; Z(i,1)=sab_term[i]; }
} else if (nm == 0 && nt > 0) {
  Z.set_size(m, 2);
  for (int i = 0; i < m; ++i) { Z(i,0)=sa_term[i]; Z(i,1)=sag_term[i]; }
} else {
  Z.set_size(m, 1);
  for (int i = 0; i < m; ++i) { Z(i,0)=sa_term[i]; }
}
arma::mat Sigma_vc = arma::cov(Z);
if (Z.n_cols == 1u) Sigma_vc(0,0) = sample_var(Z.col(0));
Sigma_vc /= std::max(1.0, (double)m);

arma::vec se_vec(m);
double n_total = 0.0;
for (int i = 0; i < m; ++i) n_total += (double)C[i].n_i;
double w2sum = 0.0;
for (int i = 0; i < m; ++i) {
  se_vec[i] = se_term[i];
  const double wi = ((double)C[i].n_i) / std::max(1.0, n_total);
  w2sum += wi * wi;
}
double var_sehat = sample_var(se_vec) * w2sum;

double var_ccc = 0.0;
if (nm > 0 && nt > 0) {
  arma::vec g3(3); g3[0]=d_sa; g3[1]=d_sab; g3[2]=d_sag;
  var_ccc += arma::as_scalar(g3.t() * Sigma_vc * g3);
} else if (nm > 0 && nt == 0) {
  arma::vec g2(2); g2[0]=d_sa; g2[1]=d_sab;
  var_ccc += arma::as_scalar(g2.t() * Sigma_vc * g2);
} else if (nm == 0 && nt > 0) {
  arma::vec g2(2); g2[0]=d_sa; g2[1]=d_sag;
  var_ccc += arma::as_scalar(g2.t() * Sigma_vc * g2);
} else {
  var_ccc += d_sa * d_sa * Sigma_vc(0,0);
}
var_ccc += d_se * d_se * var_sehat;
var_ccc += d_SB * d_SB * varSB;

double se_ccc = std::sqrt(std::max(0.0, var_ccc));

const double alpha = 1.0 - std::min(std::max(conf_level, 0.0), 1.0);
const double z = R::qnorm(1.0 - 0.5 * alpha, 0.0, 1.0, 1, 0);
double lwr = std::min(1.0, std::max(0.0, ccc - z * se_ccc));
double upr = std::min(1.0, std::max(0.0, ccc + z * se_ccc));

return Rcpp::List::create(
  _["sigma2_subject"]        = sa,
  _["sigma2_subject_method"] = sab,
  _["sigma2_subject_time"]   = sag,
  _["sigma2_error"]          = se,
  _["SB"]                    = SB,
  _["beta"]                  = beta,
  _["varFix"]                = VarFix,
  _["ccc"]                   = ccc,
  _["lwr"]                   = lwr,
  _["upr"]                   = upr,
  _["se_ccc"]                = se_ccc,
  _["conf_level"]            = conf_level
);
}
