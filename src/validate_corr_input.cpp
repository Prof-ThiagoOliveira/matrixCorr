// Thiago de Paula Oliveira
// validate_corr_input.cpp
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// ---- helpers ----
inline bool any_na_real_matrix(SEXP x, int nr, int nc){
  const double* pr = REAL(x);
  const R_xlen_t N = static_cast<R_xlen_t>(nr) * static_cast<R_xlen_t>(nc);
  for (R_xlen_t i = 0; i < N; ++i) if (ISNAN(pr[i])) return true;
  return false;
}
inline bool any_na_int_vec(const int* pi, int n){
  for (int i = 0; i < n; ++i) if (pi[i] == NA_INTEGER) return true;
  return false;
}
inline bool any_na_logi_vec(const int* pl, int n){      // LOGICALSXP uses int storage
  for (int i = 0; i < n; ++i) if (pl[i] == NA_LOGICAL) return true;
  return false;
}

//' Input validator for correlation
//'
//' @param data matrix or data.frame
//' @param check_na logical; if TRUE, scan and error on any NA/NaN
//' @return numeric (double) matrix with colnames preserved when available
//'
// [[Rcpp::export]]
 Rcpp::NumericMatrix validate_corr_input_cpp(SEXP data, bool check_na = true){
   // Accept matrix or data.frame
   const bool is_mat = Rf_isMatrix(data);
   const bool is_df  = Rf_inherits(data, "data.frame");
   if (!is_mat && !is_df) Rcpp::stop("Input must be a matrix or data frame.");

   // ---------------- matrix path ----------------
   if (is_mat){
     const int nr = Rf_nrows(data);
     const int nc = Rf_ncols(data);
     if (nc < 2) Rcpp::stop("At least two numeric columns are required.");
     if (nr < 2) Rcpp::stop("Each column must contain at least two values.");

     const SEXPTYPE t = TYPEOF(data);

     // double matrix
     if (t == REALSXP){
       if (check_na && any_na_real_matrix(data, nr, nc))
         Rcpp::stop("Missing values are not allowed.");
       return Rcpp::NumericMatrix(data); // shallow wrap, no copy
     }

     // integer matrix
     if (t == INTSXP || t == LGLSXP){
       Rcpp::IntegerMatrix Mi(data);
       if (check_na){
         if (t == INTSXP){
           for (int j = 0; j < nc; ++j) if (any_na_int_vec(&(Mi(0, j)), nr))
             Rcpp::stop("Missing values are not allowed.");
         } else { // LOGICALSXP
           for (int j = 0; j < nc; ++j) if (any_na_logi_vec(&(Mi(0, j)), nr))
             Rcpp::stop("Missing values are not allowed.");
         }
       }
       Rcpp::NumericMatrix M(nr, nc);
       for (int j = 0; j < nc; ++j)
         for (int i = 0; i < nr; ++i)
           M(i, j) = static_cast<double>(Mi(i, j));

       // preserve colnames
       Rcpp::List dn = Mi.attr("dimnames");
       if (dn.size() == 2){
         SEXP cn = dn[1];
         if (!Rf_isNull(cn)) M.attr("dimnames") = Rcpp::List::create(R_NilValue, cn);
       }
       return M;
     }

     Rcpp::stop("Matrix must be numeric (integer, logical, or double).");
   }

   // ---------------- data.frame path ----------------
   Rcpp::List df(data);
   const int p_all = df.size();
   if (p_all < 1) Rcpp::stop("No columns found.");

   // collect numeric columns (REALSXP/INTSXP/LOGICALSXP) and validate lengths
   std::vector<int> idx; idx.reserve(p_all);
   int n = -1;
   Rcpp::CharacterVector src_names = df.attr("names");

   for (int j = 0; j < p_all; ++j){
     SEXP col = df[j];
     const SEXPTYPE t = TYPEOF(col);
     if (t == REALSXP || t == INTSXP || t == LGLSXP){
       const int len = Rf_length(col);
       if (len < 2) Rcpp::stop("All columns must contain at least two values.");
       if (n == -1) n = len; else if (len != n)
         Rcpp::stop("All numeric columns must have the same length.");

       if (check_na){
         if (t == REALSXP){
           const double* pr = REAL(col);
           for (int i = 0; i < len; ++i) if (ISNAN(pr[i])) Rcpp::stop("Missing values are not allowed.");
         } else if (t == INTSXP){
           if (any_na_int_vec(INTEGER(col), len)) Rcpp::stop("Missing values are not allowed.");
         } else { // LOGICALSXP
           if (any_na_logi_vec(LOGICAL(col), len)) Rcpp::stop("Missing values are not allowed.");
         }
       }

       idx.push_back(j);
     }
   }

   if (idx.empty())  Rcpp::stop("No numeric columns found in the input.");
   if (idx.size() < 2) Rcpp::stop("At least two numeric columns are required.");

   // materialize double matrix
   const int k = static_cast<int>(idx.size());
   Rcpp::NumericMatrix M(n, k);
   Rcpp::CharacterVector cn(k);

   for (int jj = 0; jj < k; ++jj){
     const int j = idx[jj];
     SEXP col = df[j];
     const SEXPTYPE t = TYPEOF(col);

     if (t == REALSXP){
       Rcpp::NumericVector v(col);
       std::copy(v.begin(), v.end(), &(M(0, jj)));
     } else if (t == INTSXP){
       const int* pi = INTEGER(col);
       for (int i = 0; i < n; ++i) M(i, jj) = static_cast<double>(pi[i]);
     } else { // LOGICALSXP
       const int* pl = LOGICAL(col);
       for (int i = 0; i < n; ++i) M(i, jj) = static_cast<double>(pl[i] == NA_LOGICAL ? NA_REAL : pl[i]);
       if (check_na){
         for (int i = 0; i < n; ++i) if (ISNAN(M(i, jj))) Rcpp::stop("Missing values are not allowed.");
       }
     }

     if (src_names.size() > j) cn[jj] = src_names[j];
   }

   M.attr("dimnames") = Rcpp::List::create(R_NilValue, cn);
   return M;
 }
