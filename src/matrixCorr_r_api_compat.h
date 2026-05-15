// Compatibility declarations for R headers that may omit selected
// non-API extern symbols used transitively by older Rcpp headers.
//
// This header is force-included from Makevars/Makevars.win and is
// intentionally narrow: it only adds forward declarations and does not
// change runtime behavior for supported older R versions.
#ifndef MATRIXCORR_R_API_COMPAT_H
#define MATRIXCORR_R_API_COMPAT_H

#include <Rversion.h>

#if defined(R_VERSION) && (R_VERSION >= R_Version(4, 6, 0))

// This header is force-included before any Rcpp headers, so it must not pull
// in Rinternals.h directly. Forward-declaring SEXP keeps the shim compatible
// with Rcpp's own include order and avoids UB-prone header interleaving.
struct SEXPREC;
typedef struct SEXPREC* SEXP;

#ifdef __cplusplus
extern "C" {
#endif
extern SEXP R_NamespaceRegistry;
extern SEXP R_UnboundValue;
#ifdef __cplusplus
}
#endif

#endif

#endif
