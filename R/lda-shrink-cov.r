#' Shrinkage-based Diagonal Linear Discriminant Analysis (SDLDA)
#'
#' Given a set of training data, this function builds the Shrinkage-based
#' Diagonal Linear Discriminant Analysis (SDLDA) classifier, which is based on
#' the DLDA classifier, often attributed to Dudoit et al. (2002). The DLDA
#' classifier belongs to the family of Naive Bayes classifiers, where the
#' distributions of each class are assumed to be multivariate normal and to
#' share a common covariance matrix. To improve the estimation of the pooled
#' variances, Pang et al. (2009) proposed the SDLDA classifier which uses a
#' shrinkage-based estimators of the pooled covariance matrix.
#'
#' The DLDA classifier is a modification to the well-known LDA classifier, where
#' the off-diagonal elements of the pooled covariance matrix are assumed to be
#' zero -- the features are assumed to be uncorrelated. Under multivariate
#' normality, the assumption uncorrelated features is equivalent to the
#' assumption of independent features. The feature-independence assumption is a
#' notable attribute of the Naive Bayes classifier family. The benefit of these
#' classifiers is that they are fast and have much fewer parameters to estimate,
#' especially when the number of features is quite large.
#'
#' The matrix of training observations are given in \code{x}. The rows of
#' \code{x} contain the sample observations, and the columns contain the
#' features for each training observation.
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
#' probabilities are estimated as the sample proportion of observations
#' belonging to each class. Otherwise, \code{prior} should be a vector with the
#' same length as the number of classes in \code{y}. The \code{prior}
#' probabilities should be nonnegative and sum to one.
#'
#' @export
#'
#' @inheritParams lda_diag
#' @param num_alphas the number of values used to find the optimal amount of
#' shrinkage
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal
#' Discriminant Analysis and Its Applications in High-Dimensional Data,"
#' Biometrics, 65, 4, 1021-1029.
#' @examples
#' data(cells, package = "modeldata")
#' 
#' cells$case <- NULL
#' 
#' cell_train <- cells[-(1:500), ]
#' cell_test  <- cells[  1:500 , ]
#' 
#' mod <- lda_shrink_cov(class ~ ., data = cell_train)
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
lda_shrink_cov <- function(x, ...) {
  UseMethod("lda_shrink_cov")
}

#' @rdname lda_shrink_cov
#' @export
lda_shrink_cov.default <- function(x, y, prior = NULL, num_alphas = 101, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)

  obj <- diag_estimates(x, y, prior, pool = TRUE)

  # Calculates the shrinkage-based estimator of the pooled covariance matrix.
  obj$var_shrink <- var_shrinkage(
    N = obj$N,
    K = obj$num_groups,
    var_feature = obj$var_pool,
    num_alphas = num_alphas,
    t = -1
  )

  # Creates an object of type 'lda_shrink_cov' and adds the 'match.call' to the object
  obj$predictors <- colnames(x)
  class(obj) <- "lda_shrink_cov"

  obj
}

#' @rdname lda_shrink_cov
#' @importFrom stats model.frame model.matrix model.response
#' @export
lda_shrink_cov.formula <- function(formula, data, prior = NULL, num_alphas = 101, ...) {
  # The formula interface includes an intercept. If the user includes the
  # intercept in the model, it should be removed. Otherwise, errors and doom
  # happen.
  # To remove the intercept, we update the formula, like so:
  # (NOTE: The terms must be collected in case the dot (.) notation is used)
  formula <- no_intercept(formula, data)
  
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- lda_shrink_cov.default(x = x, y = y, prior = prior, num_alphas = num_alphas)

  
  est$formula <- formula
  est
}

#' Outputs the summary for a SDLDA classifier object.
#'
#' Summarizes the trained SDLDA classifier in a nice manner.
#'
#' @param x object to print
#' @param ... unused
#' @export
print.lda_shrink_cov <- function(x, ...) {
  cat("Shrinkage-based Diagonal LDA\n\n")
  print_basics(x, ...)
  invisible(x)
}

#' SDLDA prediction of the class membership of a matrix of new observations.
#'
#' The SDLDA classifier is a modification to LDA, where the off-diagonal
#' elements of the pooled sample covariance matrix are set to zero. To improve
#' the estimation of the pooled variances, we use a shrinkage method from Pang
#' et al.  (2009).
#' 
#' @rdname lda_shrink_cov
#' @export
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @references Pang, H., Tong, T., & Zhao, H. (2009). "Shrinkage-based Diagonal
#' Discriminant Analysis and Its Applications in High-Dimensional Data,"
#' Biometrics, 65, 4, 1021-1029.
predict.lda_shrink_cov <- function(object, new_data, type = "class", ...) {
  type <- match.arg(type, c("class", "score", "prob"))
  
  new_data <- process_new_data(new_data, object)
  
  if (type %in% c("score", "class")) {
    res <- apply(new_data, 1, function(obs) {
      sapply(object$est, function(class_est) {
        with(class_est, sum((obs - xbar)^2 / object$var_shrink) + log(prior))
      })
    })
  }
  
  if (type == "prob") {
    # Posterior probabilities via Bayes Theorem
    means <- lapply(object$est, "[[", "xbar")
    covs <- replicate(n = object$num_groups, object$var_shrink, simplify = FALSE)
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
