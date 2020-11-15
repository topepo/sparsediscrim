#' Linear Discriminant Analysis using the Schafer-Strimmer Covariance Matrix
#' Estimator
#'
#' Given a set of training data, this function builds the Linear Discriminant
#' Analysis (LDA) classifier, where the distributions of each class are assumed
#' to be multivariate normal and share a common covariance matrix. When the
#' pooled sample covariance matrix is singular, the linear discriminant function
#' is incalculable. This function replaces the inverse of pooled sample
#' covariance matrix with an estimator proposed by Schafer and Strimmer
#' (2005). The estimator is calculated via \code{\link[corpcor]{invcov.shrink}}.
#'
#' The matrix of training observations are given in \code{x}. The rows of \code{x}
#' contain the sample observations, and the columns contain the features for each
#' training observation.
#'
#' The vector of class labels given in \code{y} are coerced to a \code{factor}.
#' The length of \code{y} should match the number of rows in \code{x}.
#'
#' An error is thrown if a given class has less than 2 observations because the
#' variance for each feature within a class cannot be estimated with less than 2
#' observations.
#'
#' The vector, \code{prior}, contains the \emph{a priori} class membership for
#' each class. If \code{prior} is NULL (default), the class membership
#' probabilities are estimated as the sample proportion of observations belonging
#' to each class. Otherwise, \code{prior} should be a vector with the same length
#' as the number of classes in \code{y}. The \code{prior} probabilities should be
#' nonnegative and sum to one.
#'
#' @importFrom corpcor cov.shrink invcov.shrink
#' @export
#'
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param prior vector with prior probabilities for each class. If NULL
#' (default), then equal probabilities are used. See details.
#' @param ... additional arguments passed to
#' \code{\link[corpcor]{invcov.shrink}}
#' @return \code{lda_schafer} object that contains the trained classifier
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' lda_schafer_out <- lda_schafer(Species ~ ., data = iris[train, ])
#' predicted <- predict(lda_schafer_out, iris[-train, -5])$class
#'
#' lda_schafer_out2 <- lda_schafer(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(lda_schafer_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
#' @references Schafer, J., and Strimmer, K. (2005). "A shrinkage approach to
#' large-scale covariance estimation and implications for functional genomics,"
#' Statist. Appl. Genet. Mol. Biol. 4, 32.
lda_schafer <- function(x, ...) {
  UseMethod("lda_schafer")
}

#' @rdname lda_schafer
#' @export
lda_schafer.default <- function(x, y, prior = NULL, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = FALSE)

  # Centers the data matrix from each class and then estimates the lda_schafer
  # precision matrix
  x_centered <- tapply(seq_along(y), y, function(i) {
    scale(x[i, ], center = TRUE, scale = FALSE)
  })
  x_centered <- do.call(rbind, x_centered)
  obj$cov_pool <- cov.shrink(x_centered, verbose = FALSE, ...)
  obj$cov_inv <- invcov.shrink(x_centered, verbose = FALSE, ...)

  # The `corpcor` package returns a `shrinkage` object, which is actually a
  # matrix with some attributes.
  # Coerces the classes to avoid conflicts downstream.
  class(obj$cov_pool) <- "matrix"
  class(obj$cov_inv) <- "matrix"

  # Creates an object of type 'lda_schafer' and adds the 'match.call' to the object
  obj$call <- match.call()
  class(obj) <- "lda_schafer"

  obj
}

#' @param formula A formula of the form \code{groups ~ x1 + x2 + ...} That is,
#' the response is the grouping factor and the right hand side specifies the
#' (non-factor) discriminators.
#' @param data data frame from which variables specified in \code{formula} are
#' preferentially to be taken.
#' @rdname lda_schafer
#' @importFrom stats model.frame model.matrix model.response
#' @export
lda_schafer.formula <- function(formula, data, prior = NULL, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)

  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- lda_schafer.default(x = x, y = y, prior = prior)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' Outputs the summary for a lda_schafer classifier object.
#'
#' Summarizes the trained lda_schafer classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @export
print.lda_schafer <- function(x, ...) {
  cat("LDA using the Schafer-Strimmer Covariance Matrix Estimator\n\n")
  print_basics(x, ...)
  invisible(x)
}

#' Predicts of class membership of a matrix of new observations using Linear
#' Discriminant Analysis (LDA) using the Schafer-Strimmer Covariance Matrix
#' Estimator
#'
#' The Linear Discriminant Analysis (LDA) classifier involves the assumption
#' that the distributions of each class are assumed to be multivariate normal
#' and share a common covariance matrix. When the pooled sample covariance
#' matrix is singular, the linear discriminant function is incalculable. Here,
#' the inverse of the pooled sample covariance matrix is replaced with an
#' estimator from Schafer and Strimmer (2005).
#'
#' @rdname lda_schafer
#' @export
#'
#' @references Schafer, J., and Strimmer, K. (2005). "A shrinkage approach to
#' large-scale covariance estimation and implications for functional genomics,"
#' Statist. Appl. Genet. Mol. Biol. 4, 32.
#' @param object trained lda_schafer object
#' @param new_data matrix of observations to predict. Each row corresponds to a
#' new observation.
#' @return list predicted class memberships of each row in new_data
predict.lda_schafer <- function(object, new_data, ...) {
  if (!inherits(object, "lda_schafer"))  {
    rlang::abort("object not of class 'lda_schafer'")
  }

  new_data <- as.matrix(new_data)

  # Calculates the discriminant scores for each test observation
  scores <- apply(new_data, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, quadform(object$cov_inv, obs - xbar) + log(prior))
    })
  })

  if (is.vector(scores)) {
    min_scores <- which.min(scores)
  } else {
    min_scores <- apply(scores, 2, which.min)
  }

  # Posterior probabilities via Bayes Theorem
  means <- lapply(object$est, "[[", "xbar")
  covs <- replicate(n=object$num_groups, object$cov_pool, simplify=FALSE)
  priors <- lapply(object$est, "[[", "prior")
  posterior <- posterior_probs(x=new_data,
                               means=means,
                               covs=covs,
                               priors=priors)

  class <- factor(object$groups[min_scores], levels = object$groups)

  list(class = class, scores = scores, posterior = posterior)
}
