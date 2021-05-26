library(testthat)
library(sparsediscrim)

context("The MDMEB Classifier from Srivistava and Kubokawa (2007)")

test_that("The MDMEB classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  mdmeb_out <- lda_emp_bayes_eigen(Species ~ ., data = iris[train, ])
  predicted <- predict(mdmeb_out, iris[-train, -5], type = "prob")

  mdmeb_out2 <- lda_emp_bayes_eigen(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(mdmeb_out2, iris[-train, -5], type = "prob")

  # Tests that the same prob results from the matrix and formula versions of
  # the MDMEB classifier
  expect_equal(predicted, predicted2)
})

# Related to issue #41
test_that("The MDMEB classifier works properly when 1 feature used", {
  require('MASS')

  train <- seq(1, 150, by = 3)

  mdmeb_out <- lda_emp_bayes_eigen(x = iris[train, 1, drop = FALSE], y = iris[train, 5])
  predicted_cls   <- predict(mdmeb_out, iris[-train, 1, drop = FALSE])
  predicted_prob  <- predict(mdmeb_out, iris[-train, 1, drop = FALSE], type = "prob")
  predicted_score <- predict(mdmeb_out, iris[-train, 1, drop = FALSE], type = "score")

  expect_equal(length(predicted_cls), 150 - length(train))
  expect_is(predicted_prob,  "data.frame")
  expect_is(predicted_score, "data.frame")
  expect_equal(dim(predicted_prob),  c(150 - length(train), nlevels(iris$Species)))
  expect_equal(dim(predicted_score), c(150 - length(train), nlevels(iris$Species))) 
})
