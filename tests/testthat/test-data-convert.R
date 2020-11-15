library(testthat)
library(sparsediscrim)

context("Data conversion")

test_that("Bad predictor input", {
   
   expect_error(
     sparsediscrim:::pred_to_matrix(iris), 
     "the matrix was no longer numeric"
   )
   
   
})
test_that("Bad outcome input", {
   
   expect_error(
      sparsediscrim:::outcome_to_factor(mtcars$am), 
      "data should be a character or factor vector"
   )
   expect_error(
      sparsediscrim:::outcome_to_factor(mtcars[, 1:2]), 
      "data should be a character or factor vector"
   )
   
})
