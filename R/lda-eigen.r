#' The Minimum Distance Rule using Moore-Penrose Inverse (MDMP) classifier
#'
#' Given a set of training data, this function builds the MDMP classifier from
#' Srivistava and Kubokawa (2007). The MDMP classifier is an adaptation of the
#' linear discriminant analysis (LDA) classifier that is designed for
#' small-sample, high-dimensional data. Srivastava and Kubokawa (2007) have
#' proposed a modification of the standard maximum likelihood estimator of the
#' pooled covariance matrix, where only the largest 95% of the eigenvalues and
#' their corresponding eigenvectors are kept. The value of 95% is the default
#' and can be changed via the \code{eigen_pct} argument.
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
#' @export
#'
#' @inheritParams lda_diag
#' @param eigen_pct the percentage of eigenvalues kept
#' @examples
#' data(cells, package = "modeldata")
#' 
#' cells$case <- NULL
#' 
#' cell_train <- cells[-(1:500), ]
#' cell_test  <- cells[  1:500 , ]
#' 
#' mod <- lda_eigen(class ~ ., data = cell_train)
#' mod
#' 
#' # ------------------------------------------------------------------------------
#' 
#' table(
#'    predict(mod, cell_test)$.pred_class, 
#'    cell_test$class
#' )
#' 
#' predict(mod, head(cell_test), type = "prob") 
#' @references Srivastava, M. and Kubokawa, T. (2007). "Comparison of
#' Discrimination Methods for High Dimensional Data," Journal of the Japanese
#' Statistical Association, 37, 1, 123-134.
lda_eigen <- function(x, ...) {
  UseMethod("lda_eigen")
}

#' @rdname lda_eigen
#' @export
lda_eigen.default <- function(x, y, prior = NULL, eigen_pct = 0.95, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  cov_eigen <- eigen(obj$cov_pool, symmetric = TRUE)

  # trace(cov_kept) / trace(cov_pool) \approx eigen_pct
  # as described in the middle of page 125
  kept_evals <- with(cov_eigen,
                     which(cumsum(values) / sum(values) < eigen_pct))

  # Computes the pseudoinverse of the resulting covariance matrix estimator
  evals_inv <- 1 / cov_eigen$values[kept_evals]
  obj$cov_pool <- with(cov_eigen,
                       tcrossprod(vectors[, kept_evals] %*% diag(1 / evals_inv),
                                  vectors[, kept_evals]))
  obj$cov_inv <- with(cov_eigen,
                      tcrossprod(vectors[, kept_evals] %*% diag(evals_inv),
                                 vectors[, kept_evals]))

  # Creates an object of type 'lda_eigen' and adds the 'match.call' to the object
  obj$predictors <- colnames(x)
  class(obj) <- "lda_eigen"

  obj
}

#' @rdname lda_eigen
#' @importFrom stats model.frame model.matrix model.response
#' @export
lda_eigen.formula <- function(formula, data, prior = NULL, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)

  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- lda_eigen.default(x = x, y = y, prior = prior)
  
  est$formula <- formula
  est
}

#' Outputs the summary for a MDMP classifier object.
#'
#' Summarizes the trained lda_eigen classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @export
print.lda_eigen <- function(x, ...) {
  cat("Minimum Distance Rule using Moore-Penrose Inverse\n\n")
  print_basics(x, ...)
  invisible(x)
}

#' Predicts of class membership of a matrix of new observations using the
#' Minimum Distance Rule using Moore-Penrose Inverse (MDMP) classifier
#'
#' The MDMP classifier from Srivistava and Kubokawa (2007) is an adaptation of the
#' linear discriminant analysis (LDA) classifier that is designed for
#' small-sample, high-dimensional data. Srivastava and Kubokawa (2007) have
#' proposed a modification of the standard maximum likelihood estimator of the
#' pooled covariance matrix, where only the largest 95% of the eigenvalues and
#' their corresponding eigenvectors are kept.
#'
#' @rdname lda_eigen
#' @export
#' @inheritParams predict.lda_diag
predict.lda_eigen <- function(object, new_data, type = "class", ...) {
  type <- match.arg(type, c("class", "score", "prob"))
  
  new_data <- process_new_data(new_data, object)
  
  if (type %in% c("score", "class")) {
    res <- apply(new_data, 1, function(obs) {
      sapply(object$est, function(class_est) {
        with(class_est, sum((obs - xbar)^2 / object$var_pool) + log(prior))
      })
    })
  }
  
  if (type == "prob") {
    # Posterior probabilities via Bayes Theorem
    means <- lapply(object$est, "[[", "xbar")
    covs <- replicate(n=object$num_groups, object$cov_pool, simplify=FALSE)
    priors <- lapply(object$est, "[[", "prior")
    res <- posterior_probs(x = new_data, means = means, covs = covs, priors = priors)
  } 
  
  if (type == "class") {
    if (is.vector(res)) {
      min_scores <- which.min(res)
    } else {
      min_scores <- apply(res, 2, which.min)
    }
    res <- factor(object$groups[min_scores], levels = object$groups)
  }
  
  format_predictions(res, type)
}
