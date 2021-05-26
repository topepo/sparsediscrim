#' Linear Discriminant Analysis using the Thomaz-Kitani-Gillies Covariance
#' Matrix Estimator
#'
#' Given a set of training data, this function builds the Linear Discriminant
#' Analysis (LDA) classifier, where the distributions of each class are assumed
#' to be multivariate normal and share a common covariance matrix. When the
#' pooled sample covariance matrix is singular, the linear discriminant function
#' is incalculable. This function replaces the pooled sample covariance matrix
#' with a regularized estimator from Thomaz et al. (2006), where the smallest
#' eigenvalues are replaced with the average eigenvalue. Specifically, small
#' eigenvalues here means that the eigenvalues are less than the average
#' eigenvalue.
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
#' @return `lda_thomaz` object that contains the trained classifier
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' lda_thomaz_out <- lda_thomaz(Species ~ ., data = iris[train, ])
#' predicted <- predict(lda_thomaz_out, iris[-train, -5])$class
#'
#' lda_thomaz_out2 <- lda_thomaz(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(lda_thomaz_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
#' @references Thomaz, C. E., Kitani, E. C., and Gillies, D. F. (2006). "A
#' maximum uncertainty LDA-based approach for limited sample size problems with
#' application to face recognition," J. Braz. Comp. Soc., 12, 2, 7-18.
lda_thomaz <- function(x, ...) {
  UseMethod("lda_thomaz")
}

#' @rdname lda_thomaz
#' @export
lda_thomaz.default <- function(x, y, prior = NULL, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)
  complete <- complete.cases(x) & complete.cases(y)
  x <- x[complete,,drop = FALSE]
  y <- y[complete]

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  # Computes eigenvalue decomposition of pooled sample covariance matrix
  # Then regularizes the estimator based on Thomaz et al.'s (2006) method
  cov_eigen <- eigen(obj$cov_pool, symmetric = TRUE)
  evals <- cov_eigen$values
  mean_eval <- mean(evals)
  evals[evals < mean_eval] <- mean_eval

  if (obj$p > 1) {
    obj$cov_pool <- with(cov_eigen,
                         tcrossprod(vectors %*% diag(evals), vectors))

    # Computes the inverse of the resulting covariance matrix estimator
    obj$cov_inv <- with(cov_eigen,
                        tcrossprod(vectors %*% diag(1 / evals), vectors))
  } else {
    obj$cov_pool <- with(cov_eigen,
                         tcrossprod(vectors %*% as.matrix(evals), vectors))
    obj$cov_inv <- with(cov_eigen,
                        tcrossprod(vectors %*% as.matrix(1 / evals), vectors))
  }

  # Creates an object of type 'lda_thomaz'
  obj$col_names <- colnames(x)
  obj <- new_discrim_object(obj, "lda_thomaz")

  obj
}

#' @inheritParams lda_diag
#' @rdname lda_thomaz
#' @importFrom stats model.frame model.matrix model.response
#' @export
lda_thomaz.formula <- function(formula, data, prior = NULL, ...) {
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

  est <- lda_thomaz.default(x = x, y = y, prior = prior)
  est$.terms <- .terms
  est <- new_discrim_object(est, class(est))
  est
}

#' Outputs the summary for a lda_thomaz classifier object.
#'
#' Summarizes the trained lda_thomaz classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @keywords internal
#' @export
print.lda_thomaz <- function(x, ...) {
  cat("LDA using the Thomaz-Kitani-Gillies Covariance Matrix Estimator\n\n")
  print_basics(x, ...)
  invisible(x)
}

#' Predicts of class membership of a matrix of new observations using Linear
#' Discriminant Analysis (LDA) using the Schafer-Strimmer Covariance Matrix
#' Estimator
#'
#' Given a set of training data, this function builds the Linear Discriminant
#' Analysis (LDA) classifier, where the distributions of each class are assumed
#' to be multivariate normal and share a common covariance matrix. When the
#' pooled sample covariance matrix is singular, the linear discriminant function
#' is incalculable. This function replaces the pooled sample covariance matrix
#' with a regularized estimator from Thomaz et al. (2006), where the smallest
#' eigenvalues are replaced with the average eigenvalue. Specifically, small
#' eigenvalues here means that the eigenvalues are less than the average
#' eigenvalue.
#'
#' @rdname lda_thomaz
#' @export
#' @inheritParams predict.lda_diag

predict.lda_thomaz <- function(object, newdata, ...) {
  newdata <- process_newdata(object, newdata)

  # Calculates the discriminant scores for each test observation
  scores <- apply(newdata, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, quadform(object$cov_inv, obs - xbar) + log(prior))
    })
  })

  # Posterior probabilities via Bayes Theorem
  means <- lapply(object$est, "[[", "xbar")
  covs <- replicate(n=object$num_groups, object$cov_pool, simplify=FALSE)
  priors <- lapply(object$est, "[[", "prior")
  posterior <- posterior_probs(x=newdata,
                               means=means,
                               covs=covs,
                               priors=priors)

  class <- score_to_class(scores, object)

  list(class = class, scores = scores, posterior = posterior)
}
