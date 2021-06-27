#' The Minimum Distance Empirical Bayesian Estimator (MDEB) classifier
#'
#' Given a set of training data, this function builds the MDEB classifier from
#' Srivistava and Kubokawa (2007). The MDEB classifier is an adaptation of the
#' linear discriminant analysis (LDA) classifier that is designed for
#' small-sample, high-dimensional data. Rather than using the standard maximum
#' likelihood estimator of the pooled covariance matrix, Srivastava and Kubokawa
#' (2007) have proposed an Empirical Bayes estimator where the eigenvalues of
#' the pooled sample covariance matrix are shrunken towards the identity matrix:
#' the shrinkage constant has a closed form and is quick to calculate.
#'
#' The matrix of training observations are given in `x`. The rows of `x`
#' contain the sample observations, and the columns contain the features for each
#' training observation.
#'
#' The vector of class labels given in `y` are coerced to a `factor`.
#' The length of `y` should match the number of rows in `x`.
#'
#' An error is thrown if a given class has less than 2 observations because the
#' variance for each feature within a class cannot be estimated with less than 2
#' observations.
#'
#' The vector, `prior`, contains the _a priori_ class membership for
#' each class. If `prior` is NULL (default), the class membership
#' probabilities are estimated as the sample proportion of observations belonging
#' to each class. Otherwise, `prior` should be a vector with the same length
#' as the number of classes in `y`. The `prior` probabilities should be
#' nonnegative and sum to one.
#'
#' @export
#'
#' @inheritParams lda_diag
#' @return `lda_emp_bayes` object that contains the trained MDEB classifier
#' @examples
#' library(modeldata)
#' data(penguins)
#' pred_rows <- seq(1, 344, by = 20)
#' penguins <- penguins[, c("species", "body_mass_g", "flipper_length_mm")]
#' mdeb_out <- lda_emp_bayes(species ~ ., data = penguins[-pred_rows, ])
#' predicted <- predict(mdeb_out, penguins[pred_rows, -1], type = "class")
#'
#' mdeb_out2 <- lda_emp_bayes(x = penguins[-pred_rows, -1], y = penguins$species[-pred_rows])
#' predicted2 <- predict(mdeb_out2, penguins[pred_rows, -1], type = "class")
#' all.equal(predicted, predicted2)
#' @references Srivastava, M. and Kubokawa, T. (2007). "Comparison of
#' Discrimination Methods for High Dimensional Data," Journal of the Japanese
#' Statistical Association, 37, 1, 123-134.
lda_emp_bayes <- function(x, ...) {
  UseMethod("lda_emp_bayes")
}

#' @rdname lda_emp_bayes
#' @export
lda_emp_bayes.default <- function(x, y, prior = NULL, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)
  complete <- complete.cases(x) & complete.cases(y)
  x <- x[complete,,drop = FALSE]
  y <- y[complete]

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  # Creates an object of type 'lda_emp_bayes'
  obj$col_names <- colnames(x)
  obj <- new_discrim_object(obj, "lda_emp_bayes")

  obj
}

#' @inheritParams lda_diag
#' @rdname lda_emp_bayes
#' @importFrom stats model.frame model.matrix model.response
#' @export
lda_emp_bayes.formula <- function(formula, data, prior = NULL, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)

  mf <- model.frame(formula = formula, data = data)
  .terms <- attr(mf, "terms")
  x <- model.matrix(.terms, data = mf)
  y <- model.response(mf)

  est <- lda_emp_bayes.default(x = x, y = y, prior = prior)
  est$.terms <- .terms
  est <- new_discrim_object(est, class(est))
  est
}

#' Outputs the summary for a MDEB classifier object.
#'
#' Summarizes the trained lda_emp_bayes classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @keywords internal
#' @export
print.lda_emp_bayes <- function(x, ...) {
  cat("Minimum Distance Empirical Bayesian Estimator\n\n")
  print_basics(x, ...)
  invisible(x)
}

#' Predicts of class membership of a matrix of new observations using the
#' Minimum Distance Empirical Bayesian Estimator (MDEB) classifier
#'
#' The MDEB classifier from Srivistava and Kubokawa (2007) is an adaptation of the
#' linear discriminant analysis (LDA) classifier that is designed for
#' small-sample, high-dimensional data. Rather than using the standard maximum
#' likelihood estimator of the pooled covariance matrix, Srivastava and Kubokawa
#' (2007) have proposed an Empirical Bayes estimator where the eigenvalues of
#' the pooled sample covariance matrix are shrunken towards the identity matrix:
#' the shrinkage constant has a closed form and is quick to calculate
#'
#' @rdname lda_emp_bayes
#' @export
#' @inheritParams predict.lda_diag

predict.lda_emp_bayes <- function(object, newdata, type = c("class", "prob", "score"), ...) {
  type <- rlang::arg_match0(type, c("class", "prob", "score"), arg_nm = "type")
  newdata <- process_newdata(object, newdata)

  # Calculates the MDEB shrinkage constant and then computes the inverse of the
  # MDEB covariance matrix estimator
  shrink <- with(object, sum(diag(cov_pool)) / min(N, p))
  cov_pool <- with(object, cov_pool + shrink * diag(p))

  # Calculates the discriminant scores for each test observation
  scores <- apply(newdata, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, quadform_inv(cov_pool, obs - xbar) + log(prior))
    })
  })
  
  if (type == "prob") {
    # Posterior probabilities via Bayes Theorem
    means <- lapply(object$est, "[[", "xbar")
    covs <- replicate(n=object$num_groups, cov_pool, simplify=FALSE)
    priors <- lapply(object$est, "[[", "prior")
    res <- posterior_probs(x = newdata, means = means, covs = covs, priors = priors)
    res <- as.data.frame(res)
    
  } else if (type == "class") {
    res <- score_to_class(scores, object)
  } else {
    res <- t(scores)
    res <- as.data.frame(res)
  }
  res
}
