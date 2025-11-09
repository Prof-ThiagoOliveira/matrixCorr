#!/usr/bin/env bash

# Custom rchk driver for matrixCorr.

set -euo pipefail

PKG=$(grep '^Package:' DESCRIPTION | awk '{ print $2 }')
LIB_ROOT=/opt/R/devel-rchk/packages/lib
PKG_LIB="${LIB_ROOT}/${PKG}"

rm -rf "${PKG_LIB}"
mkdir -p "${LIB_ROOT}"

R CMD INSTALL --preclean --no-byte-compile --install-tests --library="${LIB_ROOT}" .

if [ -x /opt/R/devel-rchk/bin/bcheck ]; then
  /opt/R/devel-rchk/bin/bcheck "${PKG}"
else
  echo "::warning::bcheck executable not found; skipping static analysis" >&2
  exit 0
fi

Rscript - <<'RSCRIPT'
pkg <- read.dcf("DESCRIPTION", "Package")[[1]]
bcheck <- file.path("/opt/R/devel-rchk/packages/lib", pkg, "libs", paste0(pkg, ".so.bcheck"))
if (!file.exists(bcheck)) stop("bcheck output file does not exist")
lines <- readLines(bcheck)
errs <- grep("^  \\[[A-Z][A-Z]\\] ", lines, value = TRUE)
ignored <- c(
  "unsupported form of unprotect",
  "has address taken, results will be incomplete",
  "possible protection stack imbalance",
  "negative depth",
  "attempt to unprotect more items",
  "protection/Armor.h",
  "protection/Shield.h",
  "Armor<SEXPREC*>",
  "Rcpp_protect(",
  "Rcpp_unprotect("
)
for (ign in ignored) {
  keep <- !grepl(ign, errs, fixed = TRUE)
  errs <- errs[keep]
}
if (length(errs)) {
  cat("::group::rchk issues\n")
  writeLines(errs)
  cat("::endgroup::\n")
  q(save = "no", status = 1)
}
cat("::group::rchk results\nNo actionable rchk findings after filtering known Rcpp false positives.\n::endgroup::\n")
RSCRIPT
