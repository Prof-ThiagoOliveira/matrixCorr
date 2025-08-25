// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <algorithm>   // stable_sort
#include <numeric>     // iota
#include <cmath>       // std::isfinite, fabs
#include <cstdlib>     // setenv
#ifdef _OPENMP
#include <omp.h>
#if defined(__unix__) || defined(__APPLE__)
#include <dlfcn.h> // dlsym for lazy binding to BLAS controls
#endif
#endif
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
namespace detail_blas_guard {

struct BLASThreadGuard {
  // saved values (best-effort)
  int saved_openblas = -1;
  int saved_mkl      = -1;

  // function pointer types
  using OB_get_t = int  (*)();
  using OB_set_t = void (*)(int);
  using MKL_get_t = int  (*)();
  using MKL_set_t = void (*)(int);
  using MKL_dyn_t = void (*)(int);

  // resolved symbols (may be null)
  OB_get_t  ob_get  = nullptr;
  OB_set_t  ob_set  = nullptr;             // openblas_set_num_threads
  OB_set_t  ob_set_local = nullptr;        // openblas_set_num_threads_local
  MKL_get_t mkl_get = nullptr;             // mkl_get_max_threads / MKL_Get_Max_Threads
  MKL_set_t mkl_set_local = nullptr;       // mkl_set_num_threads_local / MKL_Set_Num_Threads_Local
  MKL_dyn_t mkl_set_dynamic = nullptr;     // mkl_set_dynamic

  // ctor sets BLAS threads to 1 (best-effort)
  explicit BLASThreadGuard(int one = 1) {
#if defined(__APPLE__)
    setenv("VECLIB_MAXIMUM_THREADS", "1", 1);
    setenv("ACCELERATE_MAX_THREADS", "1", 1);
#endif

#if defined(__unix__) || defined(__APPLE__)
    // Resolve OpenBLAS symbols if present
    ob_get       = reinterpret_cast<OB_get_t>( dlsym(RTLD_DEFAULT, "openblas_get_num_threads") );
    ob_set       = reinterpret_cast<OB_set_t>( dlsym(RTLD_DEFAULT, "openblas_set_num_threads") );
    ob_set_local = reinterpret_cast<OB_set_t>( dlsym(RTLD_DEFAULT, "openblas_set_num_threads_local") );

    // Resolve MKL symbols (names differ slightly across versions)
    mkl_get = reinterpret_cast<MKL_get_t>( dlsym(RTLD_DEFAULT, "mkl_get_max_threads") );
    if (!mkl_get)
      mkl_get = reinterpret_cast<MKL_get_t>( dlsym(RTLD_DEFAULT, "MKL_Get_Max_Threads") );

    mkl_set_local = reinterpret_cast<MKL_set_t>( dlsym(RTLD_DEFAULT, "mkl_set_num_threads_local") );
    if (!mkl_set_local)
      mkl_set_local = reinterpret_cast<MKL_set_t>( dlsym(RTLD_DEFAULT, "MKL_Set_Num_Threads_Local") );

    mkl_set_dynamic = reinterpret_cast<MKL_dyn_t>( dlsym(RTLD_DEFAULT, "mkl_set_dynamic") );
#endif

    // Save current settings where possible
    if (ob_get)   saved_openblas = ob_get();
    if (mkl_get)  saved_mkl      = mkl_get();

    // Force 1 thread (prefer local setters if present)
    if (ob_set_local)      ob_set_local(one);
    else if (ob_set)       ob_set(one);

    if (mkl_set_local)     mkl_set_local(one);
    if (mkl_set_dynamic)   mkl_set_dynamic(0); // disable MKL dynamic teams
  }

  // dtor restores prior settings (best-effort)
  ~BLASThreadGuard() {
    if (saved_openblas > 0) {
      if (ob_set_local) ob_set_local(saved_openblas);
      else if (ob_set)  ob_set(saved_openblas);
    }
    if (saved_mkl > 0 && mkl_set_local) mkl_set_local(saved_mkl);
    if (mkl_set_dynamic) mkl_set_dynamic(1);
  }
};

inline void harden_omp_runtime_once() {
  // Ensure stable teams: no nested, no dynamic resizing
  omp_set_dynamic(0);
#if defined(_OPENMP) && (_OPENMP >= 201307) // OpenMP 4.0+
  omp_set_max_active_levels(1);
#else
  omp_set_nested(0);
#endif
}

} // namespace detail_blas_guard
#endif // _OPENMP

//---------------- helpers ----------------//
static inline mat solve_sympd_safe(const mat& A, const mat& B) {
  // Try SPD inverse first; if it fails, add an adaptive ridge; else pinv.
  mat Ai;
  if (inv_sympd(Ai, A)) return Ai * B;

  // adaptive ridge based on average diag scale
  double base = 1.0;
  if (A.n_rows > 0) {
    double tr = trace(A);
    if (std::isfinite(tr) && tr > 0.0) base = std::max(1.0, tr / A.n_rows);
  }
  double lam = std::max(1e-12, 1e-8 * base);

  for (int k = 0; k < 6; ++k) {
    mat Aj = A;
    Aj.diag() += lam;
    if (inv_sympd(Ai, Aj)) return Ai * B;
    lam *= 10.0;
  }
  // last resort
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

// AR(1) variance reduction factor for the mean of T equally spaced points
static inline double ar1_kappa_T(double rho, int T) {
  if (T <= 1) return 1.0;                    // mean of 1 point has no reduction
  const double r = rho;
  double acc = (double)T;                    // the 'T' term
  double rpow = r;
  for (int k = 1; k <= T-1; ++k) {           // add 2*(T-k)*rho^k
    acc += 2.0 * (double)(T - k) * rpow;
    rpow *= r;
  }
  double TT = (double)T * (double)T;
  double kappa = acc / TT;
  // guard against tiny negatives from rounding (rho near ±1)
  return std::max(kappa, 1e-12);
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

// ---------------- build AR(1) precision C^{-1} for a subject -------------
static inline void make_ar1_Cinv(const std::vector<int>& tim_i, double rho, arma::mat& Cinv) {
  const int n_i = (int)tim_i.size();
  Cinv.zeros(n_i, n_i);
  if (n_i == 0) return;

  const double r2 = rho * rho;
  const double denom = std::max(1.0 - r2, 1e-12);

  int s = 0;
  while (s < n_i) {
    if (tim_i[s] < 0) { Cinv(s,s) += 1.0; ++s; continue; }
    int e = s;
    while (e + 1 < n_i && tim_i[e+1] >= 0) ++e;

    const int L = e - s + 1;
    if (L == 1) {
      Cinv(s,s) += 1.0;
    } else {
      Cinv(s,s)         += 1.0 / denom;
      Cinv(e,e)         += 1.0 / denom;
      for (int t = s + 1; t <= e - 1; ++t)
        Cinv(t,t) += (1.0 + r2) / denom;
      for (int t = s; t <= e - 1; ++t) {
        Cinv(t, t+1) += -rho / denom;
        Cinv(t+1, t) += -rho / denom;
      }
    }
    s = e + 1;
  }

  // light jitter to keep SPD numerically
  Cinv.diag() += 1e-10;
}


// --------------- build base U_i matrix explicitly (n_i x r_base) ---------
// Columns: [ intercept | method dummies (nm) | time dummies (nt) ]
static inline void build_U_base_matrix(const std::vector<int>& met_i,
                                       const std::vector<int>& tim_i,
                                       int nm, int nt,
                                       arma::mat& U) {
  const int n_i = (int)tim_i.size();            // <<-- was met_i.size()
  const int r_base = 1 + (nm>0?nm:0) + (nt>0?nt:0);
  U.zeros(n_i, r_base);
  for (int t=0; t<n_i; ++t) U(t,0) = 1.0;
  int col = 1;
  if (nm > 0) {
    for (int t=0; t<n_i; ++t) {
      int l = met_i[t];
      if (l >= 0) U(t, col + l) = 1.0;
    }
    col += nm;
  }
  if (nt > 0) {
    for (int t=0; t<n_i; ++t) {
      int tt = tim_i[t];
      if (tt >= 0) U(t, col + tt) = 1.0;
    }
  }
}


// --------------- utility to extract per-subject blocks -------------------
static inline void rows_take_to(const arma::mat& M, const std::vector<int>& rows, arma::mat& out) {
  const int n_i = (int)rows.size();
  const int q = M.n_cols;
  out.set_size(n_i, q);
  for (int k=0; k<n_i; ++k) out.row(k) = M.row(rows[k]);
}

// --------------- add U_extra * a_extra to Ub (length n_i) ----------------
static inline void add_Uextra_times(const arma::mat& Uextra, const arma::vec& a_extra, arma::vec& out) {
  if (Uextra.n_rows == 0 || Uextra.n_cols == 0) return;
  out += Uextra * a_extra;
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
    Rcpp::Nullable<Rcpp::NumericMatrix> auxDr = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> Zr = R_NilValue,
    bool use_ar1 = false,
    double ar1_rho = 0.0
) {

#ifdef _OPENMP
  // Harden OMP runtime and prevent BLAS conflicts for the whole call.
  detail_blas_guard::harden_omp_runtime_once();
  detail_blas_guard::BLASThreadGuard _guard_one_thread_blas(1);
#endif

  const int n = yr.size();
  if (Xr.nrow() != n) stop("nrow(X) must match length(y)");
  if (subject.size() != n) stop("length(subject) mismatch");
  if (method.size()  && method.size()!=n)  stop("length(method) mismatch");
  if (time.size()    && time.size()!=n)    stop("length(time) mismatch");

  arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
  arma::vec y(yr.begin(), yr.size(), false);
  const int p = X.n_cols;

  // ---------------- optional extra random-effect design Zr ---------
  arma::mat Z; // n x qZ (may be empty)
  int qZ = 0;
  bool has_extra = false;
  if (Zr.isNotNull()) {
    Rcpp::NumericMatrix Zrm = Zr.get();
    Z = arma::mat(Zrm.begin(), Zrm.nrow(), Zrm.ncol(), false);
    if ((int)Z.n_rows != n) stop("Zr must have n rows");
    qZ = (int)Z.n_cols;
    has_extra = (qZ > 0);
  }
  if (use_ar1) {
    if (!std::isfinite(ar1_rho))
      Rcpp::stop("ar1_rho is NA/NaN; pass a finite value in (-0.999, 0.999) or disable ar1.");
    if (std::fabs(ar1_rho) >= 0.999)
      Rcpp::stop("ar1_rho must be in (-0.999, 0.999).");
    if (nt == 0)
      Rcpp::warning("use_ar1=TRUE but nt==0; AR(1) will be ignored.");
  }

  // subjects
  std::vector<int> subj_idx; int m = 0;
  reindex_subject(subject, subj_idx, m);
  BySubj S = index_by_subject(subj_idx, method, time, m);

  // EM init
  double sa  = 1.0;
  double sab = (nm>0 ? 0.5 : 0.0);
  double sag = (nt>0 ? 0.5 : 0.0);
  double se  = 1.0;
  double tau2 = 0.5; // only used when has_extra
  const double eps = 1e-10;

  arma::vec beta(p, arma::fill::zeros);

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

    C[i].UtU.set_size(r,r);
    make_UtU(met, tim, n_i, nm, nt, C[i].UtU);

    C[i].Uty.set_size(r);
    accum_Ut_vec(rows, met, tim, nm, nt, [&](int idx){ return y[idx]; }, C[i].Uty);

    C[i].Utx.set_size(r, p);
    for (int k=0;k<p;++k) {
      arma::vec tmp(r, arma::fill::zeros);
      accum_Ut_vec(rows, met, tim, nm, nt, [&](int idx){ return X(idx,k); }, tmp);
      C[i].Utx.col(k) = tmp;
    }

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

  // ---------------------------- EM iterations -------------------------------
  for (int iter=0; iter<max_iter; ++iter) {

    // (1) Assemble XtViX, XtViy
    arma::mat XtViX(p,p,arma::fill::zeros);
    arma::vec XtViy(p, arma::fill::zeros);

    if (use_ar1 || has_extra) {
      // -------- generic Woodbury with R_i^{-1} and extra U --------
#ifdef _OPENMP
      int nthreads = omp_get_max_threads();
      std::vector<arma::mat> XtViX_tls(nthreads, arma::mat(p,p, arma::fill::zeros));
      std::vector<arma::vec> XtViy_tls(nthreads, arma::vec(p, arma::fill::zeros));

#pragma omp parallel
{
  int tid = omp_get_thread_num();
  arma::mat& XtViX_loc = XtViX_tls[tid];
  arma::vec& XtViy_loc = XtViy_tls[tid];

  arma::mat Ubase, Ueff, Zi, Cinv, Rinv, M, A, Zsol, X_i, RinvX, RinvU, S_ux, XTRinvX;
  arma::vec y_i, S_uy, XTRinvY;

#pragma omp for schedule(static)
  for (int i=0; i<m; ++i) {
    const auto& rows_i = S.rows[i];
    const auto& met_i  = S.met[i];
    const auto& tim_i  = S.tim[i];
    const int n_i = (int)rows_i.size();
    if (n_i == 0) continue;

    std::vector<int> ord(n_i);
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){
      int ta = tim_i[a], tb = tim_i[b];
      if (ta < 0 && tb < 0) return a < b; // keep original order for NA/NA
      if (ta < 0) return false;           // NA goes after any valid time
      if (tb < 0) return true;
      return ta < tb;                      // ascending time
    });

    X_i.set_size(n_i, p);
    y_i.set_size(n_i);
    for (int k = 0; k < n_i; ++k) {
      int g = rows_i[ ord[k] ];
      X_i.row(k) = X.row(g);
      y_i[k]     = y[g];
    }

    std::vector<int> tim_ord(n_i);
    std::vector<int> met_ord(n_i, -1);
    for (int k = 0; k < n_i; ++k) {
      const int ok = ord[k];
      tim_ord[k] = tim_i[ ok ];
      if (nm > 0) met_ord[k] = met_i[ ok ];
    }

#ifndef NDEBUG
    if ((int)met_ord.size() != n_i) Rcpp::stop("met_ord size mismatch");
    if ((int)tim_ord.size() != n_i) Rcpp::stop("tim_ord size mismatch");
#endif

    const int r_base = 1 + (nm>0?nm:0) + (nt>0?nt:0);
    build_U_base_matrix(met_ord, tim_ord, nm, nt, Ubase);

    if (has_extra) {
      rows_take_to(Z, rows_i, Zi);
      arma::uvec ord_u = arma::conv_to<arma::uvec>::from(ord);
      Zi = Zi.rows(ord_u);
    } else {
      Zi.reset();
    }

    const int r_eff = r_base + (has_extra ? qZ : 0);
    Ueff.set_size(n_i, r_eff);
    Ueff.cols(0, r_base-1) = Ubase;
    if (has_extra) Ueff.cols(r_base, r_eff-1) = Zi;

    Cinv.zeros(n_i, n_i);
    if (use_ar1 && nt > 0) make_ar1_Cinv(tim_ord, ar1_rho, Cinv);
    else Cinv.eye(n_i, n_i);
    Rinv = (1.0 / std::max(se, eps)) * Cinv;

    RinvX   = Rinv * X_i;
    XTRinvX = X_i.t() * RinvX;
    XTRinvY = X_i.t() * (Rinv * y_i);

    RinvU = Rinv * Ueff;
    S_ux  = Ueff.t() * RinvX;           // r_eff x p
    S_uy  = Ueff.t() * (Rinv * y_i);    // r_eff

    M.zeros(r_eff, r_eff);
    M(0,0) = 1.0 / std::max(sa,  eps);
    int off = 1;
    if (nm>0) { for (int l=0; l<nm; ++l) M(off+l, off+l) = 1.0 / std::max(sab, eps); off += nm; }
    if (nt>0) { for (int t=0; t<nt; ++t) M(off+t, off+t) = 1.0 / std::max(sag, eps); off += nt; }
    if (has_extra) { for (int j=0; j<qZ; ++j) M(off+j, off+j) = 1.0 / std::max(tau2, eps); }

    M += Ueff.t() * RinvU;

    A.set_size(r_eff, 1 + p);
    A.col(0) = S_uy;
    A.cols(1, p) = S_ux;
    Zsol = solve_sympd_safe(M, A);

    arma::vec Z_y = Zsol.col(0);
    arma::mat Z_X = Zsol.cols(1, p);

    for (int k=0; k<p; ++k)
      XtViy_loc[k] += XTRinvY[k] - dot(S_ux.col(k), Z_y);

    for (int k=0; k<p; ++k) for (int l=k; l<p; ++l) {
      double val = XTRinvX(k,l) - dot(S_ux.col(k), Z_X.col(l));
      XtViX_loc(k,l) += val; if (l!=k) XtViX_loc(l,k) += val;
    }
  }
}
for (int t=0; t<nthreads; ++t) { XtViX += XtViX_tls[t]; XtViy += XtViy_tls[t]; }
#else
// serial
arma::mat Ubase, Ueff, Zi, Cinv, Rinv, M, A, Zsol, X_i, RinvX, RinvU, S_ux, XTRinvX;
arma::vec y_i, S_uy, XTRinvY;

for (int i=0; i<m; ++i) {
  const auto& rows_i = S.rows[i];
  const auto& met_i  = S.met[i];
  const auto& tim_i  = S.tim[i];
  const int n_i = (int)rows_i.size();
  if (n_i == 0) continue;

  std::vector<int> ord(n_i);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){
    int ta = tim_i[a], tb = tim_i[b];
    if (ta < 0 && tb < 0) return a < b;
    if (ta < 0) return false;
    if (tb < 0) return true;
    return ta < tb;
  });

  X_i.set_size(n_i, p);
  y_i.set_size(n_i);
  for (int k = 0; k < n_i; ++k) {
    int g = rows_i[ ord[k] ];
    X_i.row(k) = X.row(g);
    y_i[k]     = y[g];
  }

  std::vector<int> tim_ord(n_i);
  std::vector<int> met_ord(n_i, -1);
  for (int k = 0; k < n_i; ++k) {
    const int ok = ord[k];
    tim_ord[k] = tim_i[ ok ];
    if (nm > 0) met_ord[k] = met_i[ ok ];
  }

  const int r_base = 1 + (nm>0?nm:0) + (nt>0?nt:0);
  build_U_base_matrix(met_ord, tim_ord, nm, nt, Ubase);

  if (has_extra) {
    rows_take_to(Z, rows_i, Zi);
    arma::uvec ord_u = arma::conv_to<arma::uvec>::from(ord);
    Zi = Zi.rows(ord_u);
  } else {
    Zi.reset();
  }

  const int r_eff = r_base + (has_extra ? qZ : 0);
  Ueff.set_size(n_i, r_eff);
  Ueff.cols(0, r_base-1) = Ubase;
  if (has_extra) Ueff.cols(r_base, r_eff-1) = Zi;

  Cinv.zeros(n_i, n_i);
  if (use_ar1 && nt > 0) make_ar1_Cinv(tim_ord, ar1_rho, Cinv);
  else Cinv.eye(n_i, n_i);
  Rinv = (1.0 / std::max(se, eps)) * Cinv;

  RinvX   = Rinv * X_i;
  XTRinvX = X_i.t() * RinvX;
  XTRinvY = X_i.t() * (Rinv * y_i);

  RinvU = Rinv * Ueff;
  S_ux  = Ueff.t() * RinvX;
  S_uy  = Ueff.t() * (Rinv * y_i);

  M.zeros(r_eff, r_eff);
  M(0,0) = 1.0 / std::max(sa,  eps);
  int off = 1;
  if (nm>0) { for (int l=0; l<nm; ++l) M(off+l, off+l) = 1.0 / std::max(sab, eps); off += nm; }
  if (nt>0) { for (int t=0; t<nt; ++t) M(off+t, off+t) = 1.0 / std::max(sag, eps); off += nt; }
  if (has_extra) { for (int j=0; j<qZ; ++j) M(off+j, off+j) = 1.0 / std::max(tau2, eps); }
  M += Ueff.t() * RinvU;

  A.set_size(r_eff, 1 + p);
  A.col(0) = S_uy;
  A.cols(1, p) = S_ux;
  Zsol = solve_sympd_safe(M, A);

  arma::vec Z_y = Zsol.col(0);
  arma::mat Z_X = Zsol.cols(1, p);

  for (int k=0; k<p; ++k)
    XtViy[k] += XTRinvY[k] - dot(S_ux.col(k), Z_y);

  for (int k=0; k<p; ++k) for (int l=k; l<p; ++l) {
    double val = XTRinvX(k,l) - dot(S_ux.col(k), Z_X.col(l));
    XtViX(k,l) += val; if (l!=k) XtViX(l,k) += val;
  }
}
#endif

    } else {
      // ---------------- cached iid residuals ----------------
      const double inv_se = 1.0 / std::max(se, eps);
#ifdef _OPENMP
      int nthreads = omp_get_max_threads();
      std::vector<arma::mat> XtViX_tls(nthreads, arma::mat(p,p, arma::fill::zeros));
      std::vector<arma::vec> XtViy_tls(nthreads, arma::vec(p, arma::fill::zeros));

#pragma omp parallel
{
  int tid = omp_get_thread_num();
  arma::mat& XtViX_loc = XtViX_tls[tid];
  arma::vec& XtViy_loc = XtViy_tls[tid];
  arma::mat M(r,r), A(r, 1+p), Zsol(r, 1+p), TX;

#pragma omp for schedule(static)
  for (int i=0; i<m; ++i) {
    const Cache& Ci = C[i];
    M.zeros();
    M(0,0) = 1.0 / std::max(sa,  eps);
    int off = 1;
    if (nm>0) { for (int l=0; l<nm; ++l) M(off+l, off+l) = 1.0 / std::max(sab, eps); off += nm; }
    if (nt>0) { for (int t=0; t<nt; ++t) M(off+t, off+t) = 1.0 / std::max(sag, eps); }
    M += inv_se * Ci.UtU;

    A.col(0)    = Ci.Uty;
    A.cols(1,p) = Ci.Utx;
    Zsol = solve_sympd_safe(M, A * inv_se);

    TX = A.cols(1,p);
    arma::vec Z_y = Zsol.col(0);
    arma::mat Z_X = Zsol.cols(1,p);

    for (int k=0; k<p; ++k) {
      XtViy_loc[k] += inv_se * (Ci.Xty[k] - dot(TX.col(k), Z_y));
      for (int l=k; l<p; ++l) {
        double val = inv_se * (Ci.XtX(k,l) - dot(TX.col(k), Z_X.col(l)));
        XtViX_loc(k,l) += val;
        if (l!=k) XtViX_loc(l,k) += val;
      }
    }
  }
}
for (int t=0; t<nthreads; ++t) { XtViX += XtViX_tls[t]; XtViy += XtViy_tls[t]; }
#else
arma::mat M(r,r), A(r, 1+p), Zsol(r, 1+p);

for (int i=0; i<m; ++i) {
  const Cache& Ci = C[i];
  M.zeros();
  M(0,0) = 1.0 / std::max(sa,  eps);
  int off = 1;
  if (nm>0) { for (int l=0; l<nm; ++l) M(off+l, off+l) = 1.0 / std::max(sab, eps); off += nm; }
  if (nt>0) { for (int t=0; t<nt; ++t) M(off+t, off+t) = 1.0 / std::max(sag, eps); }
  M += inv_se * Ci.UtU;

  A.col(0)    = Ci.Uty;
  A.cols(1,p) = Ci.Utx;
  Zsol = solve_sympd_safe(M, A * inv_se);

  arma::mat TX  = A.cols(1,p);
  arma::vec Z_y = Zsol.col(0);
  arma::mat Z_X = Zsol.cols(1,p);

  for (int k=0; k<p; ++k) {
    XtViy[k] += inv_se * (Ci.Xty[k] - dot(TX.col(k), Z_y));
    for (int l=k; l<p; ++l) {
      double val = inv_se * (Ci.XtX(k,l) - dot(TX.col(k), Z_X.col(l)));
      XtViX(k,l) += val; if (l!=k) XtViX(l,k) += val;
    }
  }
}
#endif
    }

    // (2) GLS beta
    {
      arma::mat XtViX_inv;
      if (!inv_sympd_safe(XtViX_inv, XtViX)) {
        // adaptive ridge if needed
        arma::mat XtViXj = XtViX;
        double base = 1.0;
        if (XtViX.n_rows > 0) {
          double tr = arma::trace(XtViX);
          if (std::isfinite(tr) && tr > 0.0) base = std::max(1.0, tr / XtViX.n_rows);
        }
        double lam = std::max(1e-12, 1e-8 * base);
        bool ok = false;
        for (int k=0; k<6 && !ok; ++k) {
          XtViXj = XtViX;
          XtViXj.diag() += lam;
          ok = arma::inv_sympd(XtViX_inv, XtViXj);
          lam *= 10.0;
        }
        if (!ok) XtViX_inv = arma::pinv(XtViX);
      }
      beta = XtViX_inv * XtViy;
    }

    // (3) M-step: update sa, sab, sag, se (and tau2 if has_extra)
    double sa_acc = 0.0, sab_acc = 0.0, sag_acc = 0.0, tau2_acc = 0.0;
    double se_sumsq = 0.0, se_trace = 0.0;

    arma::vec r_global = y - X * beta;

    if (use_ar1 || has_extra) {
#ifdef _OPENMP
      int nthreads2 = omp_get_max_threads();
      std::vector<double> sa_tls(nthreads2,0.0), sab_tls(nthreads2,0.0),
      sag_tls(nthreads2,0.0), tau2_tls(nthreads2,0.0),
      ss_tls(nthreads2,0.0),   tr_tls(nthreads2,0.0);
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  arma::mat Ubase, Ueff, Zi, Cinv, Rinv, M, Minv, RinvU;
  arma::mat X_i;
  arma::vec y_i, Utr, b_i;

#pragma omp for schedule(static)
  for (int i=0; i<m; ++i) {
    const auto& rows_i = S.rows[i];
    const auto& met_i  = S.met[i];
    const auto& tim_i  = S.tim[i];
    const int n_i = (int)rows_i.size();
    if (n_i == 0) continue;

    std::vector<int> ord(n_i);
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){
      int ta = tim_i[a], tb = tim_i[b];
      if (ta < 0 && tb < 0) return a < b; // keep original order for NA/NA
      if (ta < 0) return false;           // NA goes after any valid time
      if (tb < 0) return true;
      return ta < tb;                      // ascending time
    });

    X_i.set_size(n_i, p);
    y_i.set_size(n_i);
    for (int k = 0; k < n_i; ++k) {
      int g = rows_i[ ord[k] ];
      X_i.row(k) = X.row(g);
      y_i[k]     = y[g];
    }

    std::vector<int> tim_ord(n_i);
    std::vector<int> met_ord(n_i, -1);
    for (int k = 0; k < n_i; ++k) {
      const int ok = ord[k];
      tim_ord[k] = tim_i[ ok ];
      if (nm > 0) met_ord[k] = met_i[ ok ];
    }

#ifndef NDEBUG
    if ((int)met_ord.size() != n_i) Rcpp::stop("met_ord size mismatch");
    if ((int)tim_ord.size() != n_i) Rcpp::stop("tim_ord size mismatch");
#endif

    const int r_base = 1 + (nm>0?nm:0) + (nt>0?nt:0);
    build_U_base_matrix(met_ord, tim_ord, nm, nt, Ubase);

    if (has_extra) {
      rows_take_to(Z, rows_i, Zi);
      arma::uvec ord_u = arma::conv_to<arma::uvec>::from(ord);
      Zi = Zi.rows(ord_u);
    } else {
      Zi.reset();
    }

    const int r_eff = r_base + (has_extra ? qZ : 0);
    Ueff.set_size(n_i, r_eff);
    Ueff.cols(0, r_base-1) = Ubase;
    if (has_extra) Ueff.cols(r_base, r_eff-1) = Zi;

    Cinv.zeros(n_i, n_i);
    if (use_ar1 && nt > 0) make_ar1_Cinv(tim_ord, ar1_rho, Cinv);
    else Cinv.eye(n_i, n_i);
    Rinv = (1.0 / std::max(se, eps)) * Cinv;

    M.zeros(r_eff, r_eff);
    M(0,0) = 1.0 / std::max(sa,  eps);
    int off = 1;
    if (nm>0) { for (int l=0; l<nm; ++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
    if (nt>0) { for (int t=0; t<nt; ++t) M(off+t, off+t) = 1.0/std::max(sag,eps); off += nt; }
    if (has_extra) { for (int j=0; j<qZ; ++j) M(off+j, off+j) = 1.0/std::max(tau2,eps); }

    RinvU = Rinv * Ueff;
    M += Ueff.t() * RinvU;

    arma::vec r_i(n_i);
    for (int t=0; t<n_i; ++t) r_i[t] = r_global[ rows_i[ ord[t] ] ];

    Utr = Ueff.t() * (Rinv * r_i);
    b_i = solve_sympd_safe(M, Utr);
    inv_sympd_safe(Minv, M);

    arma::vec e = r_i - Ueff * b_i;

    double quad = arma::as_scalar(e.t() * Rinv * e);
    double trce = arma::trace(Minv * (Ueff.t() * Rinv * Ueff));
    ss_tls[tid] += quad;
    tr_tls[tid] += trce;
    se_term[i] = (quad + trce) / std::max(1, n_i);
    sa_tls[tid] += b_i[0]*b_i[0] + Minv(0,0);
    int pos = 1;
    if (nm>0) { for (int l=0; l<nm; ++l) sab_tls[tid] += b_i[pos+l]*b_i[pos+l] + Minv(pos+l,pos+l); pos += nm; }
    if (nt>0) { for (int t=0; t<nt; ++t) sag_tls[tid] += b_i[pos+t]*b_i[pos+t] + Minv(pos+t,pos+t); pos += nt; }
    if (has_extra) { for (int j=0; j<qZ; ++j) tau2_tls[tid] += b_i[pos+j]*b_i[pos+j] + Minv(pos+j,pos+j); }

    // delta-method storage
    sa_term[i] = b_i[0]*b_i[0] + Minv(0,0);
    if (nm > 0) {
      double acc_m = 0.0;
      for (int l = 0; l < nm; ++l) acc_m += b_i[1 + l]*b_i[1 + l] + Minv(1 + l, 1 + l);
      sab_term[i] = acc_m / (double)nm;
    } else sab_term[i] = 0.0;
    if (nt > 0) {
      int base = 1 + (nm>0?nm:0);
      double acc_t = 0.0;
      for (int t = 0; t < nt; ++t) acc_t += b_i[base + t]*b_i[base + t] + Minv(base + t, base + t);
      sag_term[i] = acc_t / (double)nt;
    } else sag_term[i] = 0.0;
  }
}
for (int t=0; t<nthreads2; ++t) {
  sa_acc  += sa_tls[t];
  sab_acc += sab_tls[t];
  sag_acc += sag_tls[t];
  tau2_acc+= tau2_tls[t];
  se_sumsq+= ss_tls[t];
  se_trace+= tr_tls[t];
}
#else
// serial
arma::mat Ubase, Ueff, Zi, Cinv, Rinv, M, Minv, RinvU;
arma::mat X_i, RinvX, S_ux;
arma::vec y_i, Utr, b_i;

for (int i=0; i<m; ++i) {
  const auto& rows_i = S.rows[i];
  const auto& met_i  = S.met[i];
  const auto& tim_i  = S.tim[i];
  const int n_i = (int)rows_i.size();
  if (n_i == 0) continue;

  // order by time (NA last), stable within ties
  std::vector<int> ord(n_i);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){
    int ta = tim_i[a], tb = tim_i[b];
    if (ta < 0 && tb < 0) return a < b;
    if (ta < 0) return false;
    if (tb < 0) return true;
    return ta < tb;
  });

  X_i.set_size(n_i, p);
  y_i.set_size(n_i);
  for (int k = 0; k < n_i; ++k) {
    int g = rows_i[ ord[k] ];
    X_i.row(k) = X.row(g);
    y_i[k]     = y[g];
  }

  std::vector<int> tim_ord(n_i);
  std::vector<int> met_ord(n_i, -1);
  for (int k = 0; k < n_i; ++k) {
    const int ok = ord[k];
    tim_ord[k] = tim_i[ ok ];
    if (nm > 0) met_ord[k] = met_i[ ok ];
  }

#ifndef NDEBUG
  if ((int)met_ord.size() != n_i) Rcpp::stop("met_ord size mismatch");
  if ((int)tim_ord.size() != n_i) Rcpp::stop("tim_ord size mismatch");
#endif

  const int r_base = 1 + (nm>0?nm:0) + (nt>0?nt:0);
  build_U_base_matrix(met_ord, tim_ord, nm, nt, Ubase);

  if (has_extra) {
    rows_take_to(Z, rows_i, Zi);
    arma::uvec ord_u = arma::conv_to<arma::uvec>::from(ord);
    Zi = Zi.rows(ord_u);
  } else {
    Zi.reset();
  }

  const int r_eff = r_base + (has_extra ? qZ : 0);
  Ueff.set_size(n_i, r_eff);
  Ueff.cols(0, r_base-1) = Ubase;
  if (has_extra) Ueff.cols(r_base, r_eff-1) = Zi;

  Cinv.zeros(n_i, n_i);
  if (use_ar1 && nt > 0) make_ar1_Cinv(tim_ord, ar1_rho, Cinv);
  else Cinv.eye(n_i, n_i);
  Rinv = (1.0 / std::max(se, eps)) * Cinv;

  M.zeros(r_eff, r_eff);
  M(0,0) = 1.0 / std::max(sa,  eps);
  int off = 1;
  if (nm>0) { for (int l=0; l<nm; ++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
  if (nt>0) { for (int t=0; t<nt; ++t) M(off+t, off+t) = 1.0/std::max(sag,eps); off += nt; }
  if (has_extra) { for (int j=0; j<qZ; ++j) M(off+j, off+j) = 1.0/std::max(tau2,eps); }

  RinvU = Rinv * Ueff;
  M += Ueff.t() * RinvU;

  arma::vec r_i(n_i);
  for (int t=0; t<n_i; ++t) r_i[t] = r_global[ rows_i[ ord[t] ] ];

  Utr = Ueff.t() * (Rinv * r_i);
  b_i = solve_sympd_safe(M, Utr);
  inv_sympd_safe(Minv, M);

  arma::vec e = r_i - Ueff * b_i;
  double quad = arma::as_scalar(e.t() * Rinv * e);
  double trce = arma::trace(Minv * (Ueff.t() * Rinv * Ueff));
  se_sumsq += quad;
  se_trace += trce;
  se_term[i] = (quad + trce) / std::max(1, n_i);

  sa_acc += b_i[0]*b_i[0] + Minv(0,0);
  int pos = 1;
  if (nm>0) { for (int l=0; l<nm; ++l) sab_acc += b_i[pos+l]*b_i[pos+l] + Minv(pos+l,pos+l); pos += nm; }
  if (nt>0) { for (int t=0; t<nt; ++t) sag_acc += b_i[pos+t]*b_i[pos+t] + Minv(pos+t,pos+t); pos += nt; }
  if (has_extra) { for (int j=0; j<qZ; ++j) tau2_acc += b_i[pos+j]*b_i[pos+j] + Minv(pos+j,pos+j); }

  // delta-method storage
  sa_term[i] = b_i[0]*b_i[0] + Minv(0,0);
  if (nm > 0) {
    double acc_m = 0.0;
    for (int l = 0; l < nm; ++l) acc_m += b_i[1 + l]*b_i[1 + l] + Minv(1 + l, 1 + l);
    sab_term[i] = acc_m / (double)nm;
  } else sab_term[i] = 0.0;
  if (nt > 0) {
    int base = 1 + (nm>0?nm:0);
    double acc_t = 0.0;
    for (int t = 0; t < nt; ++t) acc_t += b_i[base + t]*b_i[base + t] + Minv(base + t, base + t);
    sag_term[i] = acc_t / (double)nt;
  } else sag_term[i] = 0.0;
}
#endif

    } else {
      // iid residuals path (no AR1, no extra Z)
#ifdef _OPENMP
      int nthreads2 = omp_get_max_threads();
      std::vector<double> sa_tls(nthreads2,0.0), sab_tls(nthreads2,0.0),
      sag_tls(nthreads2,0.0), ss_tls(nthreads2,0.0),
      tr_tls(nthreads2,0.0);
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  arma::mat M(r,r), Minv;
  arma::vec Utr(r), b_i;

#pragma omp for schedule(static)
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

    arma::vec r_i(n_i);
    for (int t=0; t<n_i; ++t) r_i[t] = r_global[ rows_i[t] ];

    arma::vec Ub(n_i, arma::fill::zeros);
    add_U_times(rows_i, met_i, tim_i, nm, nt, b_i, Ub);

    double ss = 0.0;
    for (int t=0; t<n_i; ++t) { double e = r_i[t] - Ub[t]; ss += e*e; }
    ss_tls[tid] += ss;
    double trce = arma::trace(Minv * Ci.UtU);
    tr_tls[tid] += trce;
    se_term[i] = (ss + trce) / std::max(1, n_i);
    sa_tls[tid] += b_i[0]*b_i[0] + Minv(0,0);
    int pos = 1;
    if (nm>0) { for (int l=0;l<nm;++l) sab_tls[tid] += b_i[pos+l]*b_i[pos+l] + Minv(pos+l, pos+l); pos += nm; }
    if (nt>0) { for (int t=0;t<nt;++t)  sag_tls[tid] += b_i[pos+t]*b_i[pos+t] + Minv(pos+t, pos+t); }
  }
}
for (int t=0; t<nthreads2; ++t) {
  sa_acc  += sa_tls[t];
  sab_acc += sab_tls[t];
  sag_acc += sag_tls[t];
  se_sumsq+= ss_tls[t];
  se_trace+= tr_tls[t];
}
#else
// serial
arma::mat M(r,r), Minv;
arma::vec Utr(r), b_i;
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

  arma::vec r_i(n_i);
  for (int t=0; t<n_i; ++t) r_i[t] = r_global[ rows_i[t] ];

  arma::vec Ub(n_i, arma::fill::zeros);
  add_U_times(rows_i, met_i, tim_i, nm, nt, b_i, Ub);

  double ss = 0.0; for (int t=0;t<n_i;++t) { double e = r_i[t] - Ub[t]; ss += e*e; }
  se_sumsq += ss;
  se_trace += arma::trace(Minv * Ci.UtU);

  sa_acc += b_i[0]*b_i[0] + Minv(0,0);
  int pos = 1;
  if (nm>0) { for (int l=0;l<nm;++l) sab_acc += b_i[pos+l]*b_i[pos+l] + Minv(pos+l, pos+l); pos += nm; }
  if (nt>0) { for (int t=0;t<nt;++t)  sag_acc += b_i[pos+t]*b_i[pos+t] + Minv(pos+t, pos+t); }
}
#endif
    }

    // ---- updates ----
    double sa_new   = std::max(sa_acc  / (double)m, eps);
    double sab_new  = (nm>0 ? std::max(sab_acc / (double)(m*nm), eps) : 0.0);
    double sag_new  = (nt>0 ? std::max(sag_acc / (double)(m*nt), eps) : 0.0);
    double tau2_new = (has_extra ? std::max(tau2_acc / (double)(m*std::max(qZ,1)), eps) : tau2);
    double se_new   = std::max( (se_sumsq + se_trace) / (double)n, eps );

    double diff = std::fabs(sa_new - sa)
      + std::fabs(sab_new - sab)
      + std::fabs(sag_new - sag)
      + (has_extra ? std::fabs(tau2_new - tau2) : 0.0)
      + std::fabs(se_new  - se);

      sa = sa_new; sab = sab_new; sag = sag_new; se = se_new; if (has_extra) tau2 = tau2_new;

      if (diff < tol) break;
  } // end EM iterations

  // ---------------- VarFix using caches ----------------
  arma::mat XtViX_final(p,p,arma::fill::zeros);

  if (use_ar1 || has_extra) {
#ifdef _OPENMP
    int nthreads3 = omp_get_max_threads();
    std::vector<arma::mat> XtViX_tls2(nthreads3, arma::mat(p,p, arma::fill::zeros));
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  arma::mat& XtViX_loc = XtViX_tls2[tid];

  arma::mat Ubase, Ueff, Zi, Cinv, Rinv, M, Zx, X_i, RinvU, S_ux, XTRinvX;

#pragma omp for schedule(static)
  for (int i=0; i<m; ++i) {
    const auto& rows_i = S.rows[i];
    const auto& met_i  = S.met[i];
    const auto& tim_i  = S.tim[i];
    const int n_i = (int)rows_i.size();
    if (n_i == 0) continue;

    std::vector<int> ord(n_i);
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){
      int ta = tim_i[a], tb = tim_i[b];
      if (ta < 0 && tb < 0) return a < b;
      if (ta < 0) return false;
      if (tb < 0) return true;
      return ta < tb;
    });

    X_i.set_size(n_i, p);
    for (int k = 0; k < n_i; ++k) {
      int g = rows_i[ ord[k] ];
      X_i.row(k) = X.row(g);
    }

    std::vector<int> tim_ord(n_i);
    std::vector<int> met_ord(n_i, -1);
    for (int k = 0; k < n_i; ++k) {
      const int ok = ord[k];
      tim_ord[k] = tim_i[ ok ];
      if (nm > 0) met_ord[k] = met_i[ ok ];
    }

    const int r_base = 1 + (nm>0?nm:0) + (nt>0?nt:0);
    build_U_base_matrix(met_ord, tim_ord, nm, nt, Ubase);

    if (has_extra) {
      rows_take_to(Z, rows_i, Zi);
      arma::uvec ord_u = arma::conv_to<arma::uvec>::from(ord);
      Zi = Zi.rows(ord_u);
    } else {
      Zi.reset();
    }

    const int r_eff = r_base + (has_extra ? qZ : 0);
    Ueff.set_size(n_i, r_eff);
    Ueff.cols(0, r_base-1) = Ubase;
    if (has_extra) Ueff.cols(r_base, r_eff-1) = Zi;

    Cinv.zeros(n_i, n_i);
    if (use_ar1 && nt > 0) make_ar1_Cinv(tim_ord, ar1_rho, Cinv);
    else Cinv.eye(n_i, n_i);
    Rinv = (1.0 / std::max(se, eps)) * Cinv;

    M.zeros(r_eff, r_eff);
    M(0,0) = 1.0 / std::max(sa, eps);
    int off = 1;
    if (nm>0) { for (int l=0; l<nm; ++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
    if (nt>0) { for (int t=0; t<nt; ++t) M(off+t, off+t) = 1.0/std::max(sag,eps); off += nt; }
    if (has_extra) { for (int j=0; j<qZ; ++j) M(off+j, off+j) = 1.0/std::max(tau2,eps); }

    RinvU = Rinv * Ueff;
    S_ux = Ueff.t() * (Rinv * X_i);
    M += Ueff.t() * RinvU;
    Zx = solve_sympd_safe(M, S_ux);

    XTRinvX = X_i.t() * (Rinv * X_i);
    for (int k=0;k<p;++k) for (int l=k;l<p;++l) {
      double val = XTRinvX(k,l) - dot(S_ux.col(k), Zx.col(l));
      XtViX_loc(k,l) += val; if (l!=k) XtViX_loc(l,k) += val;
    }
  }
}
for (int t=0; t<nthreads3; ++t) XtViX_final += XtViX_tls2[t];
#else
// serial
arma::mat Ubase, Ueff, Zi, Cinv, Rinv, M, Zx, X_i, RinvU, S_ux, XTRinvX;

for (int i=0; i<m; ++i) {
  const auto& rows_i = S.rows[i];
  const auto& met_i  = S.met[i];
  const auto& tim_i  = S.tim[i];
  const int n_i = (int)rows_i.size();
  if (n_i == 0) continue;

  std::vector<int> ord(n_i);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b){
    int ta = tim_i[a], tb = tim_i[b];
    if (ta < 0 && tb < 0) return a < b;
    if (ta < 0) return false;
    if (tb < 0) return true;
    return ta < tb;
  });

  X_i.set_size(n_i, p);
  for (int k = 0; k < n_i; ++k) {
    int g = rows_i[ ord[k] ];
    X_i.row(k) = X.row(g);
  }

  std::vector<int> tim_ord(n_i);
  std::vector<int> met_ord(n_i, -1);
  for (int k = 0; k < n_i; ++k) {
    const int ok = ord[k];
    tim_ord[k] = tim_i[ ok ];
    if (nm > 0) met_ord[k] = met_i[ ok ];
  }

  const int r_base = 1 + (nm>0?nm:0) + (nt>0?nt:0);
  build_U_base_matrix(met_ord, tim_ord, nm, nt, Ubase);

  if (has_extra) {
    rows_take_to(Z, rows_i, Zi);
    arma::uvec ord_u = arma::conv_to<arma::uvec>::from(ord);
    Zi = Zi.rows(ord_u);
  } else {
    Zi.reset();
  }

  const int r_eff = r_base + (has_extra ? qZ : 0);
  Ueff.set_size(n_i, r_eff);
  Ueff.cols(0, r_base-1) = Ubase;
  if (has_extra) Ueff.cols(r_base, r_eff-1) = Zi;

  Cinv.zeros(n_i, n_i);
  if (use_ar1 && nt > 0) make_ar1_Cinv(tim_ord, ar1_rho, Cinv);
  else Cinv.eye(n_i, n_i);
  Rinv = (1.0 / std::max(se, eps)) * Cinv;

  M.zeros(r_eff, r_eff);
  M(0,0) = 1.0 / std::max(sa, eps);
  int off = 1;
  if (nm>0) { for (int l=0; l<nm; ++l) M(off+l, off+l) = 1.0/std::max(sab,eps); off += nm; }
  if (nt>0) { for (int t=0; t<nt; ++t) M(off+t, off+t) = 1.0/std::max(sag,eps); off += nt; }
  if (has_extra) { for (int j=0; j<qZ; ++j) M(off+j, off+j) = 1.0/std::max(tau2,eps); }

  RinvU = Rinv * Ueff;
  S_ux = Ueff.t() * (Rinv * X_i);
  M += Ueff.t() * RinvU;
  Zx = solve_sympd_safe(M, S_ux);

  XTRinvX = X_i.t() * (Rinv * X_i);
  for (int k=0;k<p;++k) for (int l=k;l<p;++l) {
    double val = XTRinvX(k,l) - dot(S_ux.col(k), Zx.col(l));
    XtViX_final(k,l) += val; if (l!=k) XtViX_final(l,k) += val;
  }
}
#endif
  } else {
    // ---------------- OLD BRANCH (iid) ----------------
    const double inv_se_final = 1.0 / std::max(se, eps);
#ifdef _OPENMP
    int nthreads3 = omp_get_max_threads();
    std::vector<arma::mat> XtViX_tls2(nthreads3, arma::mat(p,p, arma::fill::zeros));
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  arma::mat& XtViX_loc = XtViX_tls2[tid];
  arma::mat M(r,r), Zx(r,p);
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
for (int t=0; t<nthreads3; ++t) XtViX_final += XtViX_tls2[t];
#else
{
  arma::mat M(r,r), Zx(r,p);
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
      XtViX_final(k,l) += val; if (l!=k) XtViX_final(l,k) += val;
    }
  }
}
#endif
  }

  arma::mat VarFix;
  if (!inv_sympd_safe(VarFix, XtViX_final)) Rcpp::stop("Failed to invert XtViX.");
  if (!arma::is_finite(VarFix)) Rcpp::stop("VarFix is not finite");

  // ---------------- SB and CCC ----------------
  double SB = 0.0, varSB = 0.0;
  if (nm > 0 && Lr.isNotNull() && auxDr.isNotNull()) {
    Rcpp::NumericMatrix Lrm = Rcpp::as<Rcpp::NumericMatrix>(Lr);
    Rcpp::NumericMatrix Drm = Rcpp::as<Rcpp::NumericMatrix>(auxDr);
    arma::mat L(Lrm.begin(), X.n_cols, Lrm.ncol(), false);
    arma::mat auxD(Drm.begin(), Drm.nrow(), Drm.ncol(), false);
    const double den = (double)nm * (double)(nm-1) * (double)std::max(nt,1);

    arma::vec difmed = L.t() * beta;
    arma::mat Afix   = L * auxD * L.t();
    double num = arma::as_scalar(difmed.t() * auxD * difmed) - arma::trace(Afix * VarFix);
    SB = std::max(num / den, 0.0);

    if (!std::isfinite(SB))    SB    = 0.0;
    if (!std::isfinite(varSB) || varSB < 0.0) varSB = 0.0;

    arma::mat AV = Afix * VarFix;
    double term1 = 2.0 * arma::trace(AV * AV);
    double term2 = 4.0 * arma::as_scalar(beta.t() * Afix * VarFix * Afix * beta);
    varSB = std::max((term1 + term2) / (den * den), 0.0);
  }

  // ---- kappa factors for time-averaged CCC ----
  // for subject by time variance (sag)
  double kappa_g_bar = 1.0;
  // for residual variance (se)
  double kappa_e_bar = 1.0;
  int    units = 0;

  if (nt > 0) {
    kappa_g_bar = 0.0;
    kappa_e_bar = 0.0;

    for (int i = 0; i < m; ++i) {
      const auto& met_i = S.met[i];
      const auto& tim_i = S.tim[i];

      if (nm > 0) {
        // per-(subject,method) blocks
        for (int l = 0; l < nm; ++l) {
          // count observed times for this (i,l)
          int T = 0;
          for (size_t k = 0; k < tim_i.size(); ++k)
            if (met_i[k] == l && tim_i[k] >= 0) ++T;
            if (T <= 0) continue;

            kappa_g_bar += 1.0 / static_cast<double>(T);
            kappa_e_bar += use_ar1 ? ar1_kappa_T(ar1_rho, T)
              : 1.0 / static_cast<double>(T);
            ++units;
        }
      } else {
        // no method factor: one block per subject
        int T = 0;
        for (size_t k = 0; k < tim_i.size(); ++k)
          if (tim_i[k] >= 0) ++T;
          if (T <= 0) continue;

          kappa_g_bar += 1.0 / static_cast<double>(T);
          kappa_e_bar += use_ar1 ? ar1_kappa_T(ar1_rho, T)
            : 1.0 / static_cast<double>(T);
          ++units;
      }
    }

    if (units > 0) {
      kappa_g_bar /= static_cast<double>(units);
      kappa_e_bar /= static_cast<double>(units);
    } else {
      // no usable time info despite nt>0
      kappa_g_bar = 1.0;
      kappa_e_bar = 1.0;
    }

    // numerical guard (helps when |rho| ~ 1 and T is large)
    kappa_e_bar = std::min(1.0, std::max(kappa_e_bar, 1e-12));
  } else {
    // no time factor -> no averaging shrinkage (sag==0), residual unchanged
    kappa_g_bar = 0.0; // unused because sag==0
    kappa_e_bar = 1.0;
  }

  // Effective (averaged) time-varying components
  const double sag_bar = (nt > 0 ? kappa_g_bar * sag : 0.0);
  const double se_bar  = kappa_e_bar * se;

  // ---- CCC of the time-averaged measurement ----
  const double ccc = (sa + sag_bar) / (sa + sab + sag_bar + SB + se_bar);

  // ---------------- delta-method SE & CI for CCC ----------------
  double Nnum = sa + sag_bar;
  double Dden = sa + sab + sag_bar + SB + se_bar;
  if (Dden < 1e-14) Dden = 1e-14;

  const double d_sa  = (sab + SB + se_bar) / (Dden * Dden);
  const double d_sab = -Nnum / (Dden * Dden);
  const double d_sag =  kappa_g_bar * (sab + SB + se_bar) / (Dden * Dden);
  const double d_se  = -kappa_e_bar * Nnum / (Dden * Dden);
  const double d_SB  = -Nnum / (Dden * Dden);

  arma::mat Zdm;
  if (nm > 0 && nt > 0) {
    Zdm.set_size(m, 3);
    for (int i = 0; i < m; ++i) { Zdm(i,0)=sa_term[i]; Zdm(i,1)=sab_term[i]; Zdm(i,2)=sag_term[i]; }
  } else if (nm > 0 && nt == 0) {
    Zdm.set_size(m, 2);
    for (int i = 0; i < m; ++i) { Zdm(i,0)=sa_term[i]; Zdm(i,1)=sab_term[i]; }
  } else if (nm == 0 && nt > 0) {
    Zdm.set_size(m, 2);
    for (int i = 0; i < m; ++i) { Zdm(i,0)=sa_term[i]; Zdm(i,1)=sag_term[i]; }
  } else {
    Zdm.set_size(m, 1);
    for (int i = 0; i < m; ++i) { Zdm(i,0)=sa_term[i]; }
  }
  arma::mat Sigma_vc = arma::cov(Zdm);
  if (Zdm.n_cols == 1u) Sigma_vc(0,0) = sample_var(Zdm.col(0));
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
    _["sigma2_extra"]          = (has_extra ? tau2 : NA_REAL),
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
