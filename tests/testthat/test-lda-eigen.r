library(testthat)
library(sparsediscrim)

context("The MDMP Classifier from Srivistava and Kubokawa (2007)")

test_that("The MDMP classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  mdmp_out <- lda_eigen(Species ~ ., data = iris[train, ])
  predicted <- predict(mdmp_out, iris[-train, -5])

  mdmp_out2 <- lda_eigen(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(mdmp_out2, iris[-train, -5])

  # Tests that the same labels result from the matrix and formula versions of
  # the MDMP classifier
  expect_equal(predicted$class, predicted2$class)

  expect_is(predicted$posterior, "matrix")
})

# Related to issue #41
test_that("The MDMP classifier works properly when 1 feature used", {
  require('MASS')

  train <- seq(1, 150, by = 3)

  mdmp_out <- lda_eigen(x = iris[train, 1, drop = FALSE], y = iris[train, 5])
  predicted <- predict(mdmp_out, iris[-train, 1, drop = FALSE])

  expect_equal(length(predicted$class), 150 - length(train))
  expect_is(predicted$posterior, "matrix")
  expect_equal(dim(predicted$posterior), c(150 - length(train), nlevels(iris$Species)))
})
