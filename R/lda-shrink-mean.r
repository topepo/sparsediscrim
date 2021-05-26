#' Shrinkage-mean-based Diagonal Linear Discriminant Analysis (SmDLDA) from
#' Tong, Chen, and Zhao (2012)
#'
#' Given a set of training data, this function builds the Shrinkage-mean-based
#' Diagonal Linear Discriminant Analysis (SmDLDA) classifier from Tong, Chen,
#' and Zhao (2012). The SmDLDA classifier incorporates a Lindley-type shrunken
#' mean estimator into the DLDA classifier from Dudoit et al. (2002). For more
#' about the DLDA classifier, see [lda_diag()].
#'
#' The DLDA classifier belongs to the family of Naive Bayes classifiers, where
#' the distributions of each class are assumed to be multivariate normal and to
#' share a common covariance matrix.
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
#' @inheritParams lda_diag
#' @return `lda_shrink_mean` object that contains the trained SmDLDA classifier
#'
#' @references Tong, T., Chen, L., and Zhao, H. (2012), "Improved Mean
#' Estimation and Its Application to Diagonal Discriminant Analysis,"
#' Bioinformatics, 28, 4, 531-537.
#' \url{http://bioinformatics.oxfordjournals.org/content/28/4/531.long}
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' smdlda_out <- lda_shrink_mean(Species ~ ., data = iris[train, ])
#' predicted <- predict(smdlda_out, iris[-train, -5])$class
#'
#' smdlda_out2 <- lda_shrink_mean(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(smdlda_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
lda_shrink_mean <- function(x, ...) {
  UseMethod("lda_shrink_mean")
}

#' @rdname lda_shrink_mean
#' @export
lda_shrink_mean.default <- function(x, y, prior = NULL, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)
  complete <- complete.cases(x) & complete.cases(y)
  x <- x[complete,,drop = FALSE]
  y <- y[complete]

  obj <- diag_estimates(x = x, y = y, prior = prior, pool = TRUE,
                        est_mean = "tong")

  # Creates an object of type 'lda_shrink_mean'
  obj$col_names <- colnames(x)
  obj <- new_discrim_object(obj, "lda_shrink_mean")

  obj
}

#' @inheritParams lda_diag
#' @rdname lda_shrink_mean
#' @importFrom stats model.frame model.matrix model.response
#' @export
lda_shrink_mean.formula <- function(formula, data, prior = NULL, ...) {
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

  est <- lda_shrink_mean.default(x = x, y = y, prior = prior)
  est$.terms <- .terms
  est <- new_discrim_object(est, class(est))
  est
}

#' Outputs the summary for a SmDLDA classifier object.
#'
#' Summarizes the trained SmDLDA classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @keywords internal
#' @export
print.lda_shrink_mean <- function(x, ...) {
  cat("Shrinkage-Mean-Based Diagonal LDA\n\n")
  print_basics(x, ...)
  invisible(x)
}

#' SmDLDA prediction of the class membership of a matrix of new observations.
#'
#' The SmDLDA classifier is a modification to LDA, where the off-diagonal
#' elements of the pooled sample covariance matrix are set to zero.
#'
#' @rdname lda_shrink_mean
#' @export
#' @inheritParams predict.lda_diag

predict.lda_shrink_mean <- function(object, newdata, ...) {
  newdata <- process_newdata(object, newdata)

  scores <- apply(newdata, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, sum((obs - xbar)^2 / object$var_pool) + log(prior))
    })
  })

  # Posterior probabilities via Bayes Theorem
  means <- lapply(object$est, "[[", "xbar")
  covs <- replicate(n=object$num_groups, object$var_pool, simplify=FALSE)
  priors <- lapply(object$est, "[[", "prior")
  posterior <- posterior_probs(x=newdata,
                               means=means,
                               covs=covs,
                               priors=priors)

  class <- score_to_class(scores, object)

  list(class = class, scores = scores, posterior = posterior)
}
