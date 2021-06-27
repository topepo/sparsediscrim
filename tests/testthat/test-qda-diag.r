library(testthat)
library(sparsediscrim)

context("The DQDA Classifier from Dudoit et al. (2002)")

test_that("The DQDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  dqda_out <- qda_diag(Species ~ ., data = iris[train, ])
  predicted <- predict(dqda_out, iris[-train, -5], type = "prob")

  dqda_out2 <- qda_diag(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(dqda_out2, iris[-train, -5], type = "prob")

  # Tests that the same prob results from the matrix and formula versions of
  # the DQDA classifier
  expect_equal(predicted, predicted2)
})

# Related to issue #41
test_that("The DQDA classifier works properly when 1 feature used", {
  require('MASS')

  train <- seq(1, 150, by = 3)

  dqda_out <- qda_diag(x = iris[train, 1, drop = FALSE], y = iris[train, 5])
  predicted_cls   <- predict(dqda_out, iris[-train, 1, drop = FALSE])
  predicted_prob  <- predict(dqda_out, iris[-train, 1, drop = FALSE], type = "prob")
  predicted_score <- predict(dqda_out, iris[-train, 1, drop = FALSE], type = "score")

  expect_equal(length(predicted_cls), 150 - length(train))
  expect_is(predicted_prob,  "data.frame")
  expect_is(predicted_score, "data.frame")
  expect_equal(dim(predicted_prob),  c(150 - length(train), nlevels(iris$Species)))
  expect_equal(dim(predicted_score), c(150 - length(train), nlevels(iris$Species))) 
})
