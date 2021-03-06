#' The Minimum Distance Rule using Modified Empirical Bayes (MDMEB) classifier
#'
#' Given a set of training data, this function builds the MDMEB classifier from
#' Srivistava and Kubokawa (2007). The MDMEB classifier is an adaptation of the
#' linear discriminant analysis (LDA) classifier that is designed for
#' small-sample, high-dimensional data. Srivastava and Kubokawa (2007) have
#' proposed a modification of the standard maximum likelihood estimator of the
#' pooled covariance matrix, where only the largest 95% of the eigenvalues and
#' their corresponding eigenvectors are kept. The resulting covariance matrix is
#' then shrunken towards a scaled identity matrix. The value of 95% is the
#' default and can be changed via the `eigen_pct` argument.
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
#' @param eigen_pct the percentage of eigenvalues kept
#' @return `lda_emp_bayes_eigen` object that contains the trained MDMEB classifier
#' @examples
#' library(modeldata)
#' data(penguins)
#' pred_rows <- seq(1, 344, by = 20)
#' penguins <- penguins[, c("species", "body_mass_g", "flipper_length_mm")]
#' mdmeb_out <- lda_emp_bayes_eigen(species ~ ., data = penguins[-pred_rows, ])
#' predicted <- predict(mdmeb_out, penguins[pred_rows, -1], type = "class")
#'
#' mdmeb_out2 <- lda_emp_bayes_eigen(x = penguins[-pred_rows, -1], y = penguins$species[-pred_rows])
#' predicted2 <- predict(mdmeb_out2, penguins[pred_rows, -1], type = "class")
#' all.equal(predicted, predicted2)
#' @references Srivastava, M. and Kubokawa, T. (2007). "Comparison of
#' Discrimination Methods for High Dimensional Data," Journal of the Japanese
#' Statistical Association, 37, 1, 123-134.
lda_emp_bayes_eigen <- function(x, ...) {
  UseMethod("lda_emp_bayes_eigen")
}

#' @rdname lda_emp_bayes_eigen
#' @export
lda_emp_bayes_eigen.default <- function(x, y, prior = NULL, eigen_pct = 0.95, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)
  complete <- complete.cases(x) & complete.cases(y)
  x <- x[complete,,drop = FALSE]
  y <- y[complete]

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  cov_eigen <- eigen(obj$cov_pool, symmetric = TRUE)

  # trace(cov_kept) / trace(cov_pool) \approx eigen_pct
  # as described in the middle of page 125
  kept_evals <- with(cov_eigen,
                     which(cumsum(values) / sum(values) < eigen_pct))

  # The eigenvalues are then shrunken towards their mean.
  kept_evals <- kept_evals + mean(kept_evals)

  # Computes the pseudoinverse of the resulting covariance matrix estimator
  evals_inv <- 1 / cov_eigen$values[kept_evals]
  obj$cov_pool <- with(cov_eigen,
                       tcrossprod(vectors[, kept_evals] %*% diag(1 / evals_inv),
                                  vectors[, kept_evals]))
  obj$cov_inv <- with(cov_eigen,
                      tcrossprod(vectors[, kept_evals] %*% diag(evals_inv),
                                 vectors[, kept_evals]))

  # Creates an object of type 'lda_emp_bayes_eigen'
  obj$col_names <- colnames(x)
  obj <- new_discrim_object(obj, "lda_emp_bayes_eigen")

  obj
}

#' @inheritParams lda_diag
#' @rdname lda_emp_bayes_eigen
#' @importFrom stats model.frame model.matrix model.response
#' @export
lda_emp_bayes_eigen.formula <- function(formula, data, prior = NULL, ...) {
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

  est <- lda_emp_bayes_eigen.default(x = x, y = y, prior = prior)
  est$.terms <- .terms
  est <- new_discrim_object(est, class(est))
  est
}

#' Outputs the summary for a MDMEB classifier object.
#'
#' Summarizes the trained lda_emp_bayes_eigen classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @keywords internal
#' @export
print.lda_emp_bayes_eigen <- function(x, ...) {
  cat("Minimum Distance Rule using Modified Empirical Bayes\n\n")
  print_basics(x, ...)
  invisible(x)
}

#' Predicts of class membership of a matrix of new observations using the
#' Minimum Distance Rule using Modified Empirical Bayes (MDMEB) classifier
#'
#' The MDMEB classifier is an adaptation of the linear discriminant analysis
#' (LDA) classifier that is designed for small-sample, high-dimensional
#' data. Srivastava and Kubokawa (2007) have proposed a modification of the
#' standard maximum likelihood estimator of the pooled covariance matrix, where
#' only the largest 95% of the eigenvalues and their corresponding eigenvectors
#' are kept. The resulting covariance matrix is then shrunken towards a scaled
#' identity matrix.
#'
#' @rdname lda_emp_bayes_eigen
#' @export
#' @inheritParams predict.lda_diag

predict.lda_emp_bayes_eigen <- function(object, newdata, type = c("class", "prob", "score"), ...) {
  type <- rlang::arg_match0(type, c("class", "prob", "score"), arg_nm = "type")
  newdata <- process_newdata(object, newdata)

  # Calculates the discriminant scores for each test observation
  scores <- apply(newdata, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, quadform(object$cov_inv, obs - xbar) + log(prior))
    })
  })
  
  if (type == "prob") {
    # Posterior probabilities via Bayes Theorem
    means <- lapply(object$est, "[[", "xbar")
    covs <- replicate(n=object$num_groups, object$cov_pool, simplify=FALSE)
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
