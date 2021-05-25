#' Diagonal Linear Discriminant Analysis (DLDA)
#'
#' Given a set of training data, this function builds the Diagonal Linear
#' Discriminant Analysis (DLDA) classifier, which is often attributed to Dudoit
#' et al. (2002). The DLDA classifier belongs to the family of Naive Bayes
#' classifiers, where the distributions of each class are assumed to be
#' multivariate normal and to share a common covariance matrix.
#'
#' The DLDA classifier is a modification to the well-known LDA classifier, where
#' the off-diagonal elements of the pooled sample covariance matrix are assumed
#' to be zero -- the features are assumed to be uncorrelated. Under multivariate
#' normality, the assumption uncorrelated features is equivalent to the
#' assumption of independent features. The feature-independence assumption is a
#' notable attribute of the Naive Bayes classifier family. The benefit of these
#' classifiers is that they are fast and have much fewer parameters to estimate,
#' especially when the number of features is quite large.
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
#' @param x Matrix or data frame containing the training data. The rows are the 
#' sample observations, and the columns are the features. Only complete data are 
#' retained. 
#' @param y Vector of class labels for each training observation. Only complete 
#' data are retained. 
#' @param prior Vector with prior probabilities for each class. If NULL
#' (default), then equal probabilities are used. See details.
#' @return `lda_diag` Object that contains the trained DLDA classifier
#'
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' dlda_out <- lda_diag(Species ~ ., data = iris[train, ])
#' predicted <- predict(dlda_out, iris[-train, -5])$class
#'
#' dlda_out2 <- lda_diag(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(dlda_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
lda_diag <- function(x, ...) {
  UseMethod("lda_diag")
}

#' @rdname lda_diag
#' @export
lda_diag.default <- function(x, y, prior = NULL, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)
  complete <- complete.cases(x) & complete.cases(y)
  x <- x[complete,,drop = FALSE]
  y <- y[complete]

  obj <- diag_estimates(x = x, y = y, prior = prior, pool = TRUE)

  # Creates an object of type 'lda_diag' 
  # Use columns that pass filters
  obj$col_names <- names(obj$est[[1]]$xbar)
  obj <- new_discrim_object(obj, "lda_diag")

  obj
}

#' @param formula A formula of the form `groups ~ x1 + x2 + ...` That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in `formula` are
#' preferentially to be taken.
#' @rdname lda_diag
#' @importFrom stats model.frame model.matrix model.response
#' @export
lda_diag.formula <- function(formula, data, prior = NULL, ...) {
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

  est <- lda_diag.default(x = x, y = y, prior = prior)
  est$.terms <- .terms
  est <- new_discrim_object(est, class(est))
  est
}

#' Outputs the summary for a DLDA classifier object.
#'
#' Summarizes the trained DLDA classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @export
print.lda_diag <- function(x, ...) {
  cat("Diagonal LDA\n\n")
  
  print_basics(x, ...)
  invisible(x)
}

#' DLDA prediction of the class membership of a matrix of new observations.
#'
#' The DLDA classifier is a modification to LDA, where the off-diagonal elements
#' of the pooled sample covariance matrix are set to zero.
#'
#' @rdname lda_diag
#' @export
#'
#' @param object trained DLDA object
#' @param newdata matrix of observations to predict. Each row corresponds to a new observation.
#' @param ... additional arguments
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @return list predicted class memberships of each row in newdata
predict.lda_diag <- function(object, newdata, ...) {
  if (!inherits(object, "lda_diag"))  {
    rlang::abort("object not of class 'lda_diag'")
  }
  if (is.vector(newdata)) {
    newdata <- as.matrix(newdata)
  }
  newdata <- process_newdata(object, newdata)

  scores <- apply(newdata, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, sum((obs - xbar)^2 / object$var_pool) + log(prior))
    })
  })

  if (is.vector(scores)) {
    min_scores <- which.min(scores)
  } else {
    min_scores <- unlist(apply(scores, 2, which.min))
  }

  # Posterior probabilities via Bayes Theorem
  means <- lapply(object$est, "[[", "xbar")
  covs <- replicate(n = object$num_groups, object$var_pool, simplify = FALSE)
  priors <- lapply(object$est, "[[", "prior")
  posterior <- posterior_probs(x=newdata,
                               means=means,
                               covs=covs,
                               priors=priors)

  class <- factor(object$groups[min_scores], levels = object$groups)

  list(class = class, scores = scores, posterior = posterior)
}
