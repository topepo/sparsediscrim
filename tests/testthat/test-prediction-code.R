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
   
   pred_cls <- predict(mod, newdata = scat[-(1:90), -1])
   pred_prb <- predict(mod, newdata = scat[-(1:90), -1], type = "prob")
   expect_true(
      all(!is.na(pred_prb[-missing_rows,]))
   )
   expect_true(
      all(is.na(pred_cls[missing_rows]))
   )
   expect_true(
      all(!is.na(pred_cls[-missing_rows]))
   )
   expect_true(
      all(is.na(pred_cls[missing_rows]))
   )
})


# ------------------------------------------------------------------------------

test_that("x/y method", {
   
   mod <- lda_diag(x = scat[1:90, 6:12], y = scat$Species[1:90])
   expect_equal(
      mod$N,
      sum(complete.cases(scat[1:90, 6:12]))
   )
   missing_rows <- which(!complete.cases(scat[-(1:90), 6:12]))
   
   pred_cls <- predict(mod, newdata = scat[-(1:90), 6:12])
   pred_prb <- predict(mod, newdata = scat[-(1:90), 6:12], type = "prob")
   expect_true(
      all(!is.na(pred_prb[-missing_rows,]))
   )
   expect_true(
      all(is.na(pred_cls[missing_rows]))
   )
   expect_true(
      all(!is.na(pred_cls[-missing_rows]))
   )
   expect_true(
      all(is.na(pred_cls[missing_rows]))
   )
})
