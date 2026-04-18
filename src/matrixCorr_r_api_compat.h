// Compatibility declarations for R headers that may omit selected
// non-API extern symbols used transitively by older Rcpp headers.
//
// This header is force-included from Makevars/Makevars.win and is
// intentionally narrow: it only adds forward declarations and does not
// change runtime behavior for supported older R versions.
#ifndef MATRIXCORR_R_API_COMPAT_H
#define MATRIXCORR_R_API_COMPAT_H

#include <Rversion.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(R_VERSION) && (R_VERSION >= R_Version(4, 6, 0))
extern SEXP R_NamespaceRegistry;
extern SEXP R_UnboundValue;
#endif

#ifdef __cplusplus
}
#endif

#endif
