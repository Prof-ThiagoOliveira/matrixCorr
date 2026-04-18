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

// Prevent legacy R remap macros (e.g. length/isNull) from leaking into C++
// headers when this compatibility header is force-included.
#ifndef R_NO_REMAP
#define R_NO_REMAP
#define MATRIXCORR_UNDEF_R_NO_REMAP
#endif

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif
extern SEXP R_NamespaceRegistry;
extern SEXP R_UnboundValue;
#ifdef __cplusplus
}
#endif

#ifdef MATRIXCORR_UNDEF_R_NO_REMAP
#undef MATRIXCORR_UNDEF_R_NO_REMAP
#undef R_NO_REMAP
#endif

#endif

#endif
