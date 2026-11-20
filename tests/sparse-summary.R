# Run in a clean R CMD check session so package attachment is observable.
library(matrixCorr)

fit <- pearson_corr(iris[, 1:4], output = "sparse", ci = TRUE)
statistical <- base::summary(fit)
stopifnot(inherits(statistical, "summary.corr_sparse"),
          all(c("estimate", "lwr", "upr") %in% names(statistical)),
          identical(summary(fit), statistical))

library(Matrix)
stopifnot(identical(base::summary(fit), statistical),
          identical(summary(fit), Matrix::summary(fit)),
          identical(names(summary(fit)), c("i", "j", "x")))
