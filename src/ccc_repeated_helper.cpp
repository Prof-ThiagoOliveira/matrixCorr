// Thiago de Paula Oliveira
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using arma::mat;
using arma::uword;

// find column index or -1
static inline int find_col(const std::unordered_map<std::string,int>& pos,
                           const std::string& key) {
  auto it = pos.find(key);
  return (it==pos.end()? -1 : it->second);
}

// try "a:b" then "b:a"
static inline int find_inter(const std::unordered_map<std::string,int>& pos,
                             const std::string& a, const std::string& b) {
  int id = find_col(pos, a + ":" + b);
  if (id >= 0) return id;
  return find_col(pos, b + ":" + a);
}

// [[Rcpp::export]]
Rcpp::List build_L_Dm_cpp(
  Rcpp::CharacterVector colnames_X,      // colnames(model.matrix(..., df_sub))
  Rcpp::Nullable<Rcpp::CharacterVector> rmet_name,   // e.g. "method" or NULL
  Rcpp::Nullable<Rcpp::CharacterVector> rtime_name,  // e.g. "time"   or NULL
  Rcpp::CharacterVector method_levels,   // levels(df_sub[[rmet]]) or character(0)
  Rcpp::CharacterVector time_levels,     // levels(df_sub[[rtime]]) or character(0)
  bool has_interaction,
  Rcpp::Nullable<Rcpp::NumericMatrix> Dmat_global = R_NilValue
) {
  const bool have_met  = rmet_name.isNotNull()  && method_levels.size() >= 1;
  const bool have_time = rtime_name.isNotNull() && time_levels.size()  >= 1;

  const std::string rmet  = have_met  ? as<std::string>(as<CharacterVector>(rmet_name)[0])  : std::string();
  const std::string rtime = have_time ? as<std::string>(as<CharacterVector>(rtime_name)[0]) : std::string();

  // effective counts (need >=2 to “exist” in this contrast)
  const int nm = have_met  ? (int)method_levels.size() : 0;
  const int nt = have_time ? (int)time_levels.size()   : 0;
  const int nm_eff = (nm >= 2 ? nm : 0);
  const int nt_eff = (nt >= 2 ? nt : 0);

  // p = number of beta columns
  const int p = colnames_X.size();
  if (p <= 0) stop("Empty design column names.");

  // map column name -> position
  std::unordered_map<std::string,int> pos;
  pos.reserve((size_t)p);
  for (int j=0; j<p; ++j) pos.emplace(std::string(colnames_X[j]), j);

  // no method contrasts -> NULL/NULL
  if (nm_eff == 0) {
    return Rcpp::List::create(
      _["L"]  = R_NilValue,
      _["Dm"] = R_NilValue,
      _["nm"] = nm_eff,
      _["nt"] = nt_eff
    );
  }

  // Helper to fetch Dsub (nt_eff==0 => 1x1)
  auto make_Dsub = [&](int nt_eff)->mat{
    if (nt_eff == 0) {
      mat D(1,1,arma::fill::ones);
      return D;
    }
    if (Dmat_global.isNotNull()) {
      NumericMatrix Dg(Dmat_global);
      if ((int)Dg.nrow() != nt_eff || (int)Dg.ncol() != nt_eff) {
        stop("Dmat_global must be %d x %d for the provided time levels.", nt_eff, nt_eff);
      }
      mat D(&Dg[0], Dg.nrow(), Dg.ncol(), /*copy_aux_mem*/ false);
      return D; // shallow view is fine (we immediately copy via return)
    } else {
      return arma::eye<mat>(nt_eff, nt_eff);
    }
  };

  // ---------- Case A: exactly two methods (fast path) ----------
    if (nm_eff == 2) {
      // Column for the 2nd method main effect: paste0(rmet, lev2)
      const std::string lev2 = as<std::string>(method_levels[1]); // 2nd level
      const std::string met2_name = rmet + lev2;
      const int met2_idx = find_col(pos, met2_name);
      if (met2_idx < 0) stop("Cannot find method column '%s' in design.", met2_name.c_str());

      mat L;
      if (nt_eff == 0) {
        // p x 1, picks beta_met2
        L.zeros(p, 1);
        L(met2_idx, 0) = 1.0;
      } else {
        // p x nt_eff, per-time diffs; baseline time uses only beta_met2
        L.zeros(p, nt_eff);
        for (int t = 0; t < nt_eff; ++t) {
          L(met2_idx, t) = 1.0; // all times include main effect
          if (has_interaction && t >= 1) {
            const std::string tlev = as<std::string>(time_levels[t]);
            // try "met2:rtimetlev" then "rtimetlev:met2"
            const std::string tname = rtime + tlev;
            int inter_idx = find_inter(pos, met2_name, tname);
            if (inter_idx < 0) {
              stop("Cannot find interaction column for '%s:%s' (or swapped) in design.",
                   met2_name.c_str(), tname.c_str());
            }
            L(inter_idx, t) = 1.0;
          }
        }
      }

      mat Dsub = make_Dsub(nt_eff);
      // For nm==2: Dm is just Dsub (nt_eff x nt_eff) or 1x1
      return Rcpp::List::create(
        _["L"]  = Rcpp::wrap(L),
        _["Dm"] = Rcpp::wrap(Dsub),
        _["nm"] = nm_eff,
        _["nt"] = nt_eff
      );
    }

  // ---------- Case B: >= 3 methods ----------
    // number of pairwise diffs per time
  const int nd = nm_eff * (nm_eff - 1) / 2;
  const int ntime_blocks = std::max(nt_eff, 1);
  const int q = nd * ntime_blocks;

  mat L(p, q, arma::fill::zeros);

  // Precompute method main-effect column indices (treatment contrasts):
    // level1 (baseline) has no column -> index = -1
  std::vector<int> met_col(nm_eff, -1);
  for (int mlev = 1; mlev < nm_eff; ++mlev) {
    const std::string mname = rmet + std::string(method_levels[mlev]);
    met_col[mlev] = find_col(pos, mname);
    if (met_col[mlev] < 0) {
      stop("Cannot find method column '%s' in design.", mname.c_str());
    }
  }
  // Precompute interaction indices for (method level >=2, time level >=2)
    // Missing entries stay -1 (will error if referenced).
  std::vector< std::vector<int> > inter_col(
    nm_eff, std::vector<int>(std::max(nt_eff,2), -1));
  if (has_interaction && nt_eff > 0) {
    for (int mlev = 1; mlev < nm_eff; ++mlev) {
      const std::string mname = rmet + std::string(method_levels[mlev]);
      for (int t = 1; t < nt_eff; ++t) {
        const std::string tname = rtime + std::string(time_levels[t]);
        int id = find_inter(pos, mname, tname);
        if (id < 0) {
          stop("Cannot find interaction column for '%s:%s' (or swapped).",
               mname.c_str(), tname.c_str());
        }
        inter_col[mlev][t] = id;
      }
    }
  }

  // Fill L: for each time block (or 1 if no time) and each pair (i<j)
  int col = 0;
  for (int tb = 0; tb < ntime_blocks; ++tb) {
    for (int i = 0; i < nm_eff - 1; ++i) {
      for (int j = i + 1; j < nm_eff; ++j, ++col) {
        // Main-effect part at baseline time (always contributes)
        // treatment coding: level1 has no column
        if (i == 0 && j > 0) {
          // difference (j - 1): +beta_mj
          L(met_col[j], col) += 1.0;
        } else if (i > 0 && j == 0) {
          // difference (0 - i): -beta_mi
          L(met_col[i], col) -= 1.0;
        } else { // i>0, j>0
          L(met_col[j], col) += 1.0;
          L(met_col[i], col) -= 1.0;
        }

        // Interaction part (only if time exists, tb>=1, and interactions on)
        if (has_interaction && nt_eff > 0 && tb >= 1) {
          if (i == 0 && j > 0) {
            L(inter_col[j][tb], col) += 1.0;
          } else if (i > 0 && j == 0) {
            L(inter_col[i][tb], col) -= 1.0;
          } else { // i>0, j>0
            L(inter_col[j][tb], col) += 1.0;
            L(inter_col[i][tb], col) -= 1.0;
          }
        }
      }
    }
  }

  // Dm = kron(Dsub, I_nd) when nm>=3 (nt_eff==0 => Dsub=1x1, so Dm = I_nd)
  mat Dsub = make_Dsub(nt_eff);
  mat Dm;
  if (nd == 1) {
    // rare case nm_eff==2 would have been caught above; still safe:
      Dm = Dsub;
  } else {
    Dm = arma::kron(Dsub, arma::eye<mat>(nd, nd));
  }

  return Rcpp::List::create(
    _["L"]  = Rcpp::wrap(L),
    _["Dm"] = Rcpp::wrap(Dm),
    _["nm"] = nm_eff,
    _["nt"] = nt_eff
  );
}
