library(testthat)
library(sparsediscrim)

context("LDA with the Moore-Penrose Pseudo-Inverse")

test_that("The lda_pseudo classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  lda_pseudo_out <- lda_pseudo(Species ~ ., data = iris[train, ])
  predicted <- predict(lda_pseudo_out, iris[-train, -5], type = "prob")

  lda_pseudo_out2 <- lda_pseudo(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(lda_pseudo_out2, iris[-train, -5], type = "prob")

  # Tests that the same prob results from the matrix and formula versions of
  # the lda_pseudo classifier
  expect_equal(predicted, predicted2)
})

# Related to issue #41
test_that("The lda_pseudo classifier works properly when 1 feature used", {
  require('MASS')

  train <- seq(1, 150, by = 3)

  lda_pseudo_out <- lda_pseudo(x = iris[train, 1, drop = FALSE], y = iris[train, 5])
  predicted_cls   <- predict(lda_pseudo_out, iris[-train, 1, drop = FALSE])
  predicted_prob  <- predict(lda_pseudo_out, iris[-train, 1, drop = FALSE], type = "prob")
  predicted_score <- predict(lda_pseudo_out, iris[-train, 1, drop = FALSE], type = "score")

  expect_equal(length(predicted_cls), 150 - length(train))
  expect_is(predicted_prob,  "data.frame")
  expect_is(predicted_score, "data.frame")
  expect_equal(dim(predicted_prob),  c(150 - length(train), nlevels(iris$Species)))
  expect_equal(dim(predicted_score), c(150 - length(train), nlevels(iris$Species))) 
})
