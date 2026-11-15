# matrixCorr 0.13.0

## New features

* `robust_ccc()` supports the package's standard matrix, sparse, and edge-list
  outputs; `error`, `complete`, and `pairwise` missing-data modes; and the
  shared `print()`, `summary()`, `plot()`, `estimate()`, `coef()`, `tidy()`,
  `ci()`, and `confint()` interfaces.
* Optional paired percentile-bootstrap confidence intervals are available via
  `ci = TRUE`, with bootstrap success counts retained in result diagnostics.

## Bug fixes and result-schema changes

* Native pair coordinates are authoritative. Duplicate-named edges carry a
  row-keyed `corr_pair_pos` map.
* Summary metadata now keeps the output filter threshold separate from latent
  correlation cutpoints. 
* Seed offsets use overflow-safe arithmetic that preserves the previously
  valid mapping and RNG streams.
* Only a Pearson summary without stored inference adds the normal
  approximation, and only for pairs with `n_complete > 3`. Rank, robust, and
  other estimators do not gain a Fisher column from that fallback.
* Summary and tidy append stored adjusted p-values, rejection decisions, and
  skipped-correlation counts and proportions when those matrices are present.
  Unadjusted p-values stay as stored.
* Unnamed matrix edge accessors use source coordinates rather than the edge
  frame's column names. Stored edge values and the `corr_dimnames` attribute
  are unchanged.
* Empty edge and dense accessor tables no longer fail when a scalar bootstrap
  count is attached. The count is repeated once per row, including zero rows.
  `confint()` on an empty edge result reports that no confidence limits are
  stored. Thin edge lists avoid allocating the estimate matrix. Denser edge
  lists reuse the existing dense path when that cost is preferable. Stored
  inference payloads retain their representation.
