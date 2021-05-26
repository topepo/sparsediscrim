#' Linear Discriminant Analysis (LDA) with the Moore-Penrose Pseudo-Inverse
#'
#' Given a set of training data, this function builds the Linear Discriminant
#' Analysis (LDA) classifier, where the distributions of each class are assumed
#' to be multivariate normal and share a common covariance matrix.
#' When the pooled sample covariance matrix is singular, the linear discriminant
#' function is incalculable. A common method to overcome this issue is to
#' replace the inverse of the pooled sample covariance matrix with the
#' Moore-Penrose pseudo-inverse, which is unique and always exists. Note that
#' when the pooled sample covariance matrix is nonsingular, it is equal to the
#' pseudo-inverse.
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
#' @param tol tolerance value below which eigenvalues are considered numerically
#' equal to 0
#' @return `lda_pseudo` object that contains the trained lda_pseudo
#' classifier
#' @examples
#' n <- nrow(iris)
#' train <- sample(seq_len(n), n / 2)
#' lda_pseudo_out <- lda_pseudo(Species ~ ., data = iris[train, ])
#' predicted <- predict(lda_pseudo_out, iris[-train, -5])$class
#'
#' lda_pseudo_out2 <- lda_pseudo(x = iris[train, -5], y = iris[train, 5])
#' predicted2 <- predict(lda_pseudo_out2, iris[-train, -5])$class
#' all.equal(predicted, predicted2)
lda_pseudo <- function(x, ...) {
  UseMethod("lda_pseudo")
}

#' @rdname lda_pseudo
#' @export
lda_pseudo.default <- function(x, y, prior = NULL, tol = 1e-8, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)
  complete <- complete.cases(x) & complete.cases(y)
  x <- x[complete,,drop = FALSE]
  y <- y[complete]

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  # Calculates the Moore-Penrose pseudo inverse based on the pooled sample
  # covariance matrix
  cov_eigen <- eigen(obj$cov_pool, symmetric = TRUE)
  evals <- cov_eigen$values
  evals_inv <- rep.int(0, times = length(evals))
  evals_inv[evals > tol] <- 1 / evals[evals > tol]

  if (obj$p > 1) {
    obj$cov_pool <- with(cov_eigen,
                         tcrossprod(vectors %*% diag(1 / evals_inv), vectors))
    obj$cov_inv <- with(cov_eigen,
                        tcrossprod(vectors %*% diag(evals_inv), vectors))
  } else {
    obj$cov_pool <- with(cov_eigen,
                         tcrossprod(vectors %*% as.matrix(1 / evals_inv), vectors))
    obj$cov_inv <- with(cov_eigen,
                        tcrossprod(vectors %*% as.matrix(evals_inv), vectors))
  }
  # Creates an object of type 'lda_pseudo'
  obj$col_names <- colnames(x)
  obj <- new_discrim_object(obj, "lda_pseudo")

  obj
}

#' @inheritParams lda_diag.formula
#' @rdname lda_pseudo
#' @importFrom stats model.frame model.matrix model.response
#' @export
lda_pseudo.formula <- function(formula, data, prior = NULL, tol = 1e-8, ...) {
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

  est <- lda_pseudo.default(x = x, y = y, prior = prior)
  est$.terms <- .terms
  est <- new_discrim_object(est, class(est))
  est
}

#' Outputs the summary for a lda_pseudo classifier object.
#'
#' Summarizes the trained lda_pseudo classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @export
print.lda_pseudo <- function(x, ...) {
  cat("LDA with the Moore-Penrose Pseudo-Inverse\n\n")
  print_basics(x, ...)
  invisible(x)
}

#' Predicts of class membership of a matrix of new observations using Linear
#' Discriminant Analysis (LDA) with the Moore-Penrose Pseudo-Inverse
#'
#' The Linear Discriminant Analysis (LDA) classifier involves the assumption
#' that the distributions of each class are assumed to be multivariate normal
#' and share a common covariance matrix. When the pooled sample covariance
#' matrix is singular, the linear discriminant function is incalculable. A
#' common method to overcome this issue is to replace the inverse of the pooled
#' sample covariance matrix with the Moore-Penrose pseudo-inverse, which is
#' unique and always exists. Note that when the pooled sample covariance matrix
#' is nonsingular, it is equal to the pseudo-inverse.
#'
#' @rdname lda_pseudo
#' @export
#'
#' @param object trained lda_pseudo object
#' @param newdata matrix of observations to predict. Each row corresponds to a
#' new observation.
#' @param ... additional arguments
#' @return list predicted class memberships of each row in newdata
predict.lda_pseudo <- function(object, newdata, ...) {
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
