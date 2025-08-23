// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

//---------------- helpers ----------------//
static inline mat solve_sympd_safe(const mat& A, const mat& B) {
  mat Z;
  // first try: symmetric PD solve (fast path)
  if (solve(Z, A, B, solve_opts::likely_sympd)) return Z;

  // second try: tiny ridge on the diagonal
  mat Aj = A;
  double jitter = 1e-8;
  if (A.n_rows > 0) {
    double tr = trace(A);
    if (std::isfinite(tr) && tr > 0.0) jitter = std::max(1e-12, 1e-8 * tr / A.n_rows);
  }
  Aj.diag() += jitter;
  if (solve(Z, Aj, B, solve_opts::likely_sympd)) return Z;

  // last resort: pseudo-inverse (robust but slower)
  return pinv(A) * B;
}

static inline bool inv_sympd_safe(mat& out, const mat& A) {
  // try fast SPD inverse
  if (inv_sympd(out, A)) return true;
  // ridge + SPD inverse
  mat Aj = A;
  double jitter = 1e-8;
  if (A.n_rows > 0) {
    double tr = trace(A);
    if (std::isfinite(tr) && tr > 0.0) jitter = std::max(1e-12, 1e-8 * tr / A.n_rows);
  }
  Aj.diag() += jitter;
  if (inv_sympd(out, Aj)) return true;
  // pseudo-inverse
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
      S.met[j].push_back(v == NA_INTEGER ? -1 : (v - 1));   // NA -> -1
    }
    if (time.size() > 0) {
      int v = time[i];
      S.tim[j].push_back(v == NA_INTEGER ? -1 : (v - 1));   // NA -> -1
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
  const int r = Utv.n_rows;
  Utv.zeros();

  // intercept
  double s0 = 0.0;
  for (int ridx : rows_i) s0 += v(ridx);
  Utv[0] = s0;

  // method sums
  if (nm > 0) {
    for (int l=0; l<nm; ++l) {
      double sm = 0.0;
      for (size_t k=0; k<rows_i.size(); ++k)
        if (met_i[k] == l) sm += v(rows_i[k]);
        Utv[1+l] = sm;
    }
  }
  // time sums
  if (nt > 0) {
    for (int t=0; t<nt; ++t) {
      double st = 0.0;
      for (size_t k=0; k<rows_i.size(); ++k)
        if (tim_i[k] == t) st += v(rows_i[k]);
        int jt = 1 + (nm>0?nm:0) + t;
        Utv[jt] = st;
    }
  }
}

//---------------- Ua builder (fixed) ----------------//
static inline void add_U_times(const std::vector<int>& rows_i,
                               const std::vector<int>& met_i,
                               const std::vector<int>& tim_i,
                               int nm, int nt,
                               const vec& a, vec& out /* length n_i */) {
  const int n_i = (int)rows_i.size();
  const double a0 = a[0];

  // copy coefficients for quick access
  std::vector<double> am(nm>0?nm:0), at(nt>0?nt:0);
  if (nm>0) for (int l=0;l<nm;++l) am[l] = a[1+l];
  if (nt>0) for (int t=0;t<nt;++t) at[t] = a[1 + (nm>0?nm:0) + t];

  for (int k=0; k<n_i; ++k) {
    double val = a0;
    if (nm>0 && met_i[k] >= 0) val += am[ met_i[k] ];
    if (nt>0 && tim_i[k] >= 0) val += at[ tim_i[k] ];
    out[k] += val;                 // <<<<<< LOCAL index (fixed)
  }
}

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

  // workspace
  const int r = 1 + (nm>0?nm:0) + (nt>0?nt:0);
  mat M(r,r), UtU(r,r), A(r, 1 + p);
  vec Ut_y(r);
  mat XtViX(p,p,fill::zeros);
  vec XtViy(p, fill::zeros);

  for (int iter=0; iter<max_iter; ++iter) {
    XtViX.zeros(); XtViy.zeros();
    const double inv_se = 1.0 / std::max(se, eps);

    for (int i=0; i<m; ++i) {
      const std::vector<int>& rows_i = S.rows[i];
      const std::vector<int>& met_i  = S.met[i];
      const std::vector<int>& tim_i  = S.tim[i];

      make_UtU(met_i, tim_i, (int)rows_i.size(), nm, nt, UtU);

      M.zeros(r,r);
      M(0,0) = 1.0 / std::max(sa, eps);
      int off = 1;
      if (nm>0) { for (int l=0;l<nm;++l) M(off+l, off+l) = 1.0 / std::max(sab, eps); off += nm; }
      if (nt>0) { for (int t=0;t<nt;++t) M(off+t, off+t) = 1.0 / std::max(sag, eps); }
      M += inv_se * UtU;

      // A = [U^T y, U^T X]
      A.set_size(r, 1 + p);
      accum_Ut_vec(rows_i, met_i, tim_i, nm, nt,
                   [&](int ridx){ return y[ridx]; }, Ut_y);
      A.col(0) = Ut_y;
      for (int k=0;k<p;++k) {
        vec UtXk(r, fill::zeros);
        accum_Ut_vec(rows_i, met_i, tim_i, nm, nt,
                     [&](int ridx){ return X(ridx,k); }, UtXk);
        A.col(1+k) = UtXk;
      }

      // Z = M^{-1} * (A / se)
      mat Z = solve_sympd_safe(M, A * inv_se);

      mat TX = A.cols(1, p); // r x p
      vec Z_y = Z.col(0);
      mat Z_X = Z.cols(1, p);

      // naive sums and rank-k corrections
      for (int k=0; k<p; ++k) {
        double sum_xy = 0.0;
        for (int idx : rows_i) sum_xy += X(idx,k) * y[idx];
        XtViy[k] += (sum_xy * inv_se) - inv_se * dot(TX.col(k), Z_y);

        for (int l=k; l<p; ++l) {
          double sum_xx = 0.0;
          for (int idx : rows_i) sum_xx += X(idx,k) * X(idx,l);
          double corr = dot(TX.col(k), Z_X.col(l));
          double val  = (sum_xx * inv_se) - inv_se * corr;
          XtViX(k,l) += val;
          if (l!=k) XtViX(l,k) += val;
        }
      }
    }

    // GLS beta
    {
      vec btmp;
      if (!solve(btmp, XtViX, XtViy, solve_opts::likely_sympd)) {
        // ridge if needed
        mat XtViXj = XtViX;
        XtViXj.diag() += 1e-8;
        if (!solve(btmp, XtViXj, XtViy, solve_opts::likely_sympd))
          btmp = pinv(XtViX) * XtViy;
      }
      beta = btmp;
    }

    // --- M-step-ish ---
    double sa_acc = 0.0, sab_acc = 0.0, sag_acc = 0.0;
    double se_sumsq = 0.0, se_trace = 0.0;

    for (int i=0; i<m; ++i) {
      const std::vector<int>& rows_i = S.rows[i];
      const std::vector<int>& met_i  = S.met[i];
      const std::vector<int>& tim_i  = S.tim[i];
      const int n_i = (int)rows_i.size();

      // r_i
      vec r_i(n_i);
      for (int t=0; t<n_i; ++t) {
        int idx = rows_i[t];
        r_i[t] = y[idx] - dot(X.row(idx), beta.t());
      }

      // UtU, M
      make_UtU(met_i, tim_i, n_i, nm, nt, UtU);
      M.zeros(r,r);
      M(0,0) = 1.0 / std::max(sa, eps);
      int off = 1;
      if (nm>0) { for (int l=0;l<nm;++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
      if (nt>0) { for (int t=0;t<nt;++t) M(off+t, off+t) = 1.0/std::max(sag,eps); }
      M += (1.0/std::max(se,eps)) * UtU;

      // U^T r_i
      vec Utr(r, fill::zeros);
      // intercept
      Utr[0] = accu(r_i);
      if (nm>0) {
        for (int l=0;l<nm;++l) {
          double sm = 0.0; for (int t=0;t<n_i;++t) if (met_i[t]==l) sm += r_i[t];
          Utr[1+l] = sm;
        }
      }
      if (nt>0) {
        for (int tt=0; tt<nt; ++tt) {
          double st = 0.0; for (int t=0; t<n_i; ++t) if (tim_i[t]==tt) st += r_i[t];
          int jt = 1 + (nm>0?nm:0) + tt;
          Utr[jt] = st;
        }
      }

      // b_i and Var(b_i|y)
      vec b_i = solve_sympd_safe(M, Utr / std::max(se,eps));
      mat Minv;
      inv_sympd_safe(Minv, M);

      // residual part: r_i - U b_i
      vec Ub(n_i, fill::zeros);
      add_U_times(rows_i, met_i, tim_i, nm, nt, b_i, Ub);
      double ss = 0.0; for (int t=0;t<n_i;++t) { double e = r_i[t] - Ub[t]; ss += e*e; }
      se_sumsq += ss;
      se_trace += trace(Minv * UtU);

      // accumulate E[b^2]
      sa_acc += b_i[0]*b_i[0] + Minv(0,0);
      int pos = 1;
      if (nm>0) { for (int l=0;l<nm;++l) sab_acc += b_i[pos+l]*b_i[pos+l] + Minv(pos+l, pos+l); pos += nm; }
      if (nt>0) { for (int t=0;t<nt;++t)  sag_acc += b_i[pos+t]*b_i[pos+t] + Minv(pos+t, pos+t); }
    }

    double sa_new  = std::max(sa_acc  / (double)m, eps);
    double sab_new = (nm>0 ? std::max(sab_acc / (double)(m*nm), eps) : 0.0);
    double sag_new = (nt>0 ? std::max(sag_acc / (double)(m*nt), eps) : 0.0);
    double se_new  = std::max( (se_sumsq + se_trace) / (double)n, eps );

    double diff = std::fabs(sa_new-sa) + std::fabs(se_new-se)
      + std::fabs(sab_new-sab) + std::fabs(sag_new-sag);
    sa=sa_new; se=se_new; sab=sab_new; sag=sag_new;
    if (diff < tol) break;
  }

  // VarFix
  XtViX.zeros();
  const double inv_se_final = 1.0 / std::max(se, eps);
  for (int i=0; i<m; ++i) {
    const std::vector<int>& rows_i = S.rows[i];
    const std::vector<int>& met_i  = S.met[i];
    const std::vector<int>& tim_i  = S.tim[i];

    make_UtU(met_i, tim_i, (int)rows_i.size(), nm, nt, UtU);
    M.zeros(r,r);
    M(0,0) = 1.0 / std::max(sa, eps);
    int off = 1;
    if (nm>0) { for (int l=0;l<nm;++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
    if (nt>0) { for (int t=0;t<nt;++t) M(off+t, off+t) = 1.0/std::max(sag,eps); }
    M += inv_se_final * UtU;

    mat A2(r, p);
    for (int k=0;k<p;++k) {
      vec UtXk(r, fill::zeros);
      double s0 = 0.0; for (int idx: rows_i) s0 += X(idx,k);
      UtXk[0] = s0;
      if (nm>0) {
        for (int l=0;l<nm;++l) {
          double sm=0.0; for (size_t t=0;t<rows_i.size();++t) if (met_i[t]==l) sm += X(rows_i[t],k);
          UtXk[1+l]=sm;
        }
      }
      if (nt>0) {
        for (int tt=0; tt<nt; ++tt) {
          double st=0.0; for (size_t t=0;t<rows_i.size();++t) if (tim_i[t]==tt) st += X(rows_i[t],k);
          int jt = 1 + (nm>0?nm:0) + tt;
          UtXk[jt]=st;
        }
      }
      A2.col(k) = UtXk;
    }

    mat Zx = solve_sympd_safe(M, A2 * inv_se_final);
    for (int k=0;k<p;++k) for (int l=k;l<p;++l) {
      double sum_xx = 0.0; for (int idx: rows_i) sum_xx += X(idx,k)*X(idx,l);
      double corr   = dot(A2.col(k), Zx.col(l));
      double val    = (sum_xx * inv_se_final) - inv_se_final * corr;
      XtViX(k,l) += val; if (l!=k) XtViX(l,k) += val;
    }
  }

  mat VarFix;
  if (!inv_sympd_safe(VarFix, XtViX)) stop("Failed to invert XtViX.");

  // SB and CCC
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

  return Rcpp::List::create(
    _["sigma2_subject"]        = sa,
    _["sigma2_subject_method"] = sab,
    _["sigma2_subject_time"]   = sag,
    _["sigma2_error"]          = se,
    _["SB"]                    = SB,
    _["beta"]                  = beta,
    _["varFix"]                = VarFix,
    _["ccc"]                   = ccc,
    _["lwr"]                   = R_NilValue,
    _["upr"]                   = R_NilValue,
    _["se_ccc"]                = R_NilValue,
    _["conf_level"]            = conf_level
  );
}
