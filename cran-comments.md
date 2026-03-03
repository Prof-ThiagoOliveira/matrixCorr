## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* Using cli and rlang to control error and checks
* Expanded correlation test coverage to exceed 75% targets
* Fixed biweight heatmap labelling regressions in plot method
* Added partial correlation diagnostics and print/plot updates
* Normalised sparse bicor outputs to avoid Matrix method regressions
* Hardened partial correlation clustering against malformed distances
* Validated Kendall tau two-column path through shared input checks
* Locked correlation heatmap scales to [-1, 1] for consistent colour mapping
* Added interactive Shiny viewer with coverage across all correlation outputs
* Implemented a fast univariate O(n log n) dispatch for distance correlation
  with an exact unbiased O(n^2) fallback path for robustness
* Removed per-pair column-copy overhead in distance-correlation C++ kernels and
  tuned OpenMP scheduling for uniform upper-triangle workloads
* Verified distance-correlation before/after numerical equivalence on fixed-seed
  workloads; outputs match to floating-point tolerance
* Confirmed no new package dependency was introduced for the distance-correlation
  acceleration path (all changes are internal C++/R code)
