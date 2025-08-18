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
package's MIT license.
