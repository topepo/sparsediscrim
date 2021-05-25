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
#' @examples
#' data(cells, package = "modeldata")
#' 
#' cells$case <- NULL
#' 
#' cell_train <- cells[-(1:500), ]
#' cell_test  <- cells[  1:500 , ]
#' 
#' mod <- lda_emp_bayes(class ~ ., data = cell_train)
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
lda_emp_bayes <- function(x, ...) {
  UseMethod("lda_emp_bayes")
}

#' @rdname lda_emp_bayes
#' @export
lda_emp_bayes.default <- function(x, y, prior = NULL, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)

  obj <- regdiscrim_estimates(x = x, y = y, prior = prior, cov = TRUE)

  # Creates an object of type 'lda_emp_bayes' and adds the 'match.call' to the object
  obj$predictors <- colnames(x)
  class(obj) <- "lda_emp_bayes"

  obj
}

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
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- lda_emp_bayes.default(x = x, y = y, prior = prior)
  
  est$formula <- formula
  est
}

#' Outputs the summary for a MDEB classifier object.
#'
#' Summarizes the trained lda_emp_bayes classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
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
predict.lda_emp_bayes <- function(object, new_data, type = "class", ...) {
  type <- match.arg(type, c("class", "score", "prob"))
  
  new_data <- process_new_data(new_data, object)
  
  # Calculates the MDEB shrinkage constant and then computes the inverse of the
  # MDEB covariance matrix estimator
  shrink <- with(object, sum(diag(cov_pool)) / min(N, p))
  cov_pool <- with(object, cov_pool + shrink * diag(p))
  
  if (type %in% c("score", "class")) {
    res <- apply(new_data, 1, function(obs) {
      sapply(object$est, function(class_est) {
        with(class_est, quadform_inv(cov_pool, obs - xbar) + log(prior))
      })
    })
  }
  
  if (type == "prob") {
    # Posterior probabilities via Bayes Theorem
    means <- lapply(object$est, "[[", "xbar")
    covs <- replicate(n = object$num_groups, cov_pool, simplify = FALSE)
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
