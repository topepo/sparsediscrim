context("predictor pre-processing")

library(testthat)
library(sparsediscrim)
library(modeldata)

# ------------------------------------------------------------------------------

data(scat, package = "modeldata")

# ------------------------------------------------------------------------------

test_that("formula method", {
   
   expect_warning(
      mod <- lda_diag(Species ~ ., data = scat[1:90, ]),
      "had zero variance"
   )
   expect_equal(
      mod$N,
      sum(complete.cases(scat[1:90, ]))
   )
   missing_rows <- which(!complete.cases(scat[-(1:90), -1]))
   
   pred <- predict(mod, newdata = scat[-(1:90), -1])
   expect_true(
      all(!is.na(pred$posterior[-missing_rows,]))
   )
   expect_true(
      all(is.na(pred$posterior[missing_rows,]))
   )
})
