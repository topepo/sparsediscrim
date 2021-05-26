library(testthat)
library(sparsediscrim)

context("The SmDLDA Classifier from Tong et al. (2012)")

test_that("The SmDLDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  smdlda_out <- lda_shrink_mean(Species ~ ., data = iris[train, ])
  predicted <- predict(smdlda_out, iris[-train, -5])

  smdlda_out2 <- lda_shrink_mean(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(smdlda_out2, iris[-train, -5])

  # Tests that the same labels result from the matrix and formula versions of
  # the SmDLDA classifier
  expect_equal(predicted$class, predicted2$class)

  expect_is(predicted$posterior, "matrix")
})

# Related to issue #41
test_that("The SmDLDA classifier works properly when 1 feature used", {
  require('MASS')

  train <- seq(1, 150, by = 3)

  smdlda_out <- lda_shrink_mean(x = iris[train, 1, drop = FALSE], y = iris[train, 5])
  predicted <- predict(smdlda_out, iris[-train, 1, drop = FALSE])

  expect_equal(length(predicted$class), 150 - length(train))
  expect_is(predicted$posterior, "matrix")
  expect_equal(dim(predicted$posterior), c(150 - length(train), nlevels(iris$Species)))
})
