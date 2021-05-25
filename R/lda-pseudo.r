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
#' @param tol tolerance value below which eigenvalues are considered numerically
#' equal to 0
#' @examples
#' data(cells, package = "modeldata")
#' 
#' cells$case <- NULL
#' 
#' cell_train <- cells[-(1:500), ]
#' cell_test  <- cells[  1:500 , ]
#' 
#' mod <- lda_pseudo(class ~ ., data = cell_train)
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

lda_pseudo <- function(x, ...) {
  UseMethod("lda_pseudo")
}

#' @rdname lda_pseudo
#' @export
lda_pseudo.default <- function(x, y, prior = NULL, tol = 1e-8, ...) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)

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
  # Creates an object of type 'lda_pseudo' and adds the 'match.call' to the object
  obj$predictors <- colnames(x)
  class(obj) <- "lda_pseudo"

  obj
}

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
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  est <- lda_pseudo.default(x = x, y = y, prior = prior)
  
  est$formula <- formula
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
predict.lda_pseudo <- function(object, new_data, type = "class", ...) {
  type <- match.arg(type, c("class", "score", "prob"))
  
  new_data <- process_new_data(new_data, object)
  
  if (type %in% c("score", "class")) {
    res <- apply(new_data, 1, function(obs) {
      sapply(object$est, function(class_est) {
        with(class_est, quadform(object$cov_inv, obs - xbar) + log(prior))
      })
    })
  }
  
  if (type == "prob") {
    # Posterior probabilities via Bayes Theorem
    means <- lapply(object$est, "[[", "xbar")
    covs <- replicate(n = object$num_groups, object$cov_pool, simplify = FALSE)
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
