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
#' @param x matrix or data frame containing the training data. The rows are
#'  the sample observations, and the columns are the features. If a data frame
#'  is used, all columns should be numeric.
#' @param y vector of class labels for each training observation
#' @param prior vector with prior probabilities for each class. If NULL
#' (default), then equal probabilities are used. See details.
#' @return The modeling function returns object that contains the trained
#'  classifier. 
#'  
#'  The predict method always returns a tibble. When `type = "class"`, the 
#'  tibble has a single column names `.pred_class`. For the other
#'  types, there are columns in the tibble with names `.pred_{class-level}`.
#'
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). "Comparison of
#' Discrimination Methods for the Classification of Tumors Using Gene Expression
#' Data," Journal of the American Statistical Association, 97, 457, 77-87.
#' @examples
#' data(cells, package = "modeldata")
#' 
#' cells$case <- NULL
#' 
#' cell_train <- cells[-(1:500), ]
#' cell_test  <- cells[  1:500 , ]
#' 
#' mod <- lda_diag(class ~ ., data = cell_train)
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
lda_diag <- function(x, ...) {
  UseMethod("lda_diag")
}

#' @rdname lda_diag
#' @export
lda_diag.default <- function(x, y, prior = NULL, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)

  obj <- diag_estimates(x = x, y = y, prior = prior, pool = TRUE)

  # Creates an object of type 'lda_diag' and adds the 'match.call' to the object
  obj$predictors <- colnames(x)
  class(obj) <- "lda_diag"

  obj
}

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
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- lda_diag.default(x = x, y = y, prior = prior)
  
  est$formula <- formula
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
#' @param object trained model object
#' @param new_data matrix of observations to predict. Each row corresponds to a new observation.
#' @param ... additional arguments
predict.lda_diag <- function(object, new_data, type = "class", ...) {
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
    covs <- replicate(n=object$num_groups, object$var_pool, simplify=FALSE)
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
