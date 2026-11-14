# Contributing to AbacusIndexTesting

Thank you for taking the time to contribute!  We welcome bug reports,
feature requests, documentation fixes, and pull requests.

## 1. Report a bug or request a feature

Use the *New issue* button and choose the appropriate template:

* **Bug report** – reproducible error or unexpected behaviour.
* **Feature request** – enhancement or new function.
* **Question** – general usage question.

## 2. Set up your development environment

```r
# clone your fork, then:
install.packages(c('devtools', 'roxygen2', 'testthat'))
devtools::load_all()       # load package code
devtools::test()           # run the full test suite
```

## 3. Pull-request checklist

- [ ] `R CMD check` passes (0 ERROR, 0 WARNING, ≤ 1 NOTE).
- [ ] New code is documented (Rd + examples).
- [ ] Unit tests cover the new behaviour.
- [ ] NEWS.md entry added.

Please keep pull requests focused; open a new issue for any large redesigns.

## 4. Code style

* Base R + {data.table}.
* Follow the tidyverse style guide for naming and spacing.
* No trailing whitespace; max 100 characters per line.

## 5. License

By contributing, you agree that your code will be released under the
package's GPL (>= 3) license.

## Internal result infrastructure

`R/corr_output.R` owns output validation, construction, conversion, and
pair extraction. `R/corr_shared.R` orchestrates shared estimator output
and inference attachment. `R/corr_result_s3.R` owns representation
methods and summaries. `R/display_utils.R` owns presentation.
`R/accessors.R` owns accessors.

`.mc_robust_corr_pairwise_payload` shares only PB and Winsor payload
assembly, with statistical inference kept separate and pair and RNG
order unchanged. Specialised estimators and S4 Matrix support stay
distinct, as do summary inference and tidy inference.
