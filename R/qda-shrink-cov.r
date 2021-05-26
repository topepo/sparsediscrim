#' Shrinkage-based Diagonal Quadratic Discriminant Analysis (SDQDA)
#'
#' Given a set of training data, this function builds the Shrinkage-based
#' Diagonal Quadratic Discriminant Analysis (SDQDA) classifier, which is based
#' on the DQDA classifier, often attributed to Dudoit et al. (2002). The DQDA
#' classifier belongs to the family of Naive Bayes classifiers, where the
#' distributions of each class are assumed to be multivariate normal. To improve
#' the estimation of the class variances, Pang et al. (2009) proposed the SDQDA
#' classifier which uses a shrinkage-based estimators of each class covariance
#' matrix.
#'
#' The DQDA classifier is a modification to the well-known QDA classifier, where
#' the off-diagonal elements of the pooled covariance matrix are assumed to be
#' zero -- the features are assumed to be uncorrelated. Under multivariate
#' normality, the assumption uncorrelated features is equivalent to the
#' assumption of independent features. The feature-independence assumption is a
#' notable attribute of the Naive Bayes classifier family. The benefit of these
#' classifiers is that they are fast and have much fewer parameters to estimate,
#' especially when the number of features is quite large.
#'
#' The matrix of training observations are given in `x`. The rows of
#' `x` contain the sample observations, and the columns contain the
#' features for each training observation.
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
#' probabilities are estimated as the sample proportion of observations
#' belonging to each class. Otherwise, `prior` should be a vector with the
#' same length as the number of classes in `y`. The `prior`
#' probabilities should be nonnegative and sum to one.
#'
#' @export
#'
#' @inheritParams lda_diag
#' @param num_alphas the number of values used to find the optimal amount of
#' shrinkage
#' @return `qda_shrink_cov` object that contains the trained SDQDA classifier
#'
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal
#' Discriminant Analysis and Its Applications in High-Dimensional Data,"
#' Biometrics, 65, 4, 1021-1029.
#' @examples
#' set.seed(42)
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' sdqda_out <- qda_shrink_cov(Species ~ ., data = iris[train, ])
#' predicted <- predict(sdqda_out, iris[-train, -5])$class
#'
#' sdqda_out2 <- qda_shrink_cov(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(sdqda_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
qda_shrink_cov <- function(x, ...) {
  UseMethod("qda_shrink_cov")
}

#' @rdname qda_shrink_cov
#' @export
qda_shrink_cov.default <- function(x, y, prior = NULL, num_alphas = 101, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)
  complete <- complete.cases(x) & complete.cases(y)
  x <- x[complete,,drop = FALSE]
  y <- y[complete]

  obj <- diag_estimates(x, y, prior, pool = FALSE)

  # Calculates the shrinkage-based estimator for each diagonal sample class
  # covariance matrix. We add these to the corresponding obj$est$var_shrink
  for(k in seq_len(obj$num_groups)) {
    obj$est[[k]]$var_shrink <- var_shrinkage(
      N = obj$est[[k]]$n,
      K = 1,
      var_feature = obj$est[[k]]$var,
      num_alphas = num_alphas,
      t = -1
    )
  }

  # Creates an object of type 'qda_shrink_cov'
  obj$col_names <- colnames(x)
  obj <- new_discrim_object(obj, "qda_shrink_cov")

  obj
}

#' @inheritParams lda_diag.formula
#' @importFrom stats model.frame model.matrix model.response
#' @rdname qda_shrink_cov
#' @export
qda_shrink_cov.formula <- function(formula, data, prior = NULL, num_alphas = 101, ...) {
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

  est <- qda_shrink_cov.default(x = x, y = y, prior = prior, num_alphas = num_alphas)

  est$.terms <- .terms
  est <- new_discrim_object(est, class(est))
  est
}

#' Outputs the summary for a SDQDA classifier object.
#'
#' Summarizes the trained SDQDA classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @export
print.qda_shrink_cov <- function(x, ...) {
  cat("Shrinkage-Based Diagonal QDA\n\n")
  print_basics(x, ...)
  invisible(x)
}

#' SDQDA prediction of the class membership of a matrix of new observations.
#'
#' The SDQDA classifier is a modification to QDA, where the off-diagonal
#' elements of the pooled sample covariance matrix are set to zero. To improve
#' the estimation of the pooled variances, we use a shrinkage method from Pang
#' et al.  (2009).
#'
#' @rdname qda_shrink_cov
#' @export
#'
#' @param object trained SDQDA object
#' @param newdata matrix of observations to predict. Each row corresponds to a
#' new observation.
#' @param ... additional arguments
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal
#' Discriminant Analysis and Its Applications in High-Dimensional Data,"
#' Biometrics, 65, 4, 1021-1029.
#' @return list predicted class memberships of each row in newdata
predict.qda_shrink_cov <- function(object, newdata, ...) {
  newdata <- process_newdata(object, newdata)

  scores <- apply(newdata, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, sum((obs - xbar)^2 / var_shrink + log(var_shrink))
           + log(prior))
    })
  })

  if (is.vector(scores)) {
    min_scores <- which.min(scores)
  } else {
    min_scores <- apply(scores, 2, which.min)
  }

  # Posterior probabilities via Bayes Theorem
  means <- lapply(object$est, "[[", "xbar")
  covs <- lapply(object$est, "[[", "var_shrink")
  priors <- lapply(object$est, "[[", "prior")
  posterior <- posterior_probs(x=newdata,
                               means=means,
                               covs=covs,
                               priors=priors)

  class <- factor(object$groups[min_scores], levels = object$groups)

  list(class = class, scores = scores, posterior = posterior)
}
