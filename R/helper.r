#' Quadratic form of a matrix and a vector
#'
#' We compute the quadratic form of a vector and a matrix in an efficient
#' manner. Let `x` be a real vector of length `p`, and let `A` be
#' a p x p real matrix. Then, we compute the quadratic form \eqn{q = x' A x}.
#'
#' A naive way to compute the quadratic form is to explicitly write
#' `t(x) \%*\% A \%*\% x`, but for large `p`, this operation is
#' inefficient. We provide a more efficient method below.
#'
#' Note that we have adapted the code from:
#' \url{https://stat.ethz.ch/pipermail/r-help/2005-November/081940.html}
#'
#' @param A matrix of dimension p x p
#' @param x vector of length p
#' @return scalar value
quadform <- function(A, x) {
  drop(crossprod(x, A %*% x))
}

#' Quadratic Form of the inverse of a matrix and a vector
#'
#' We compute the quadratic form of a vector and the inverse of a matrix in an
#' efficient manner. Let `x` be a real vector of length `p`, and let
#' `A` be a p x p nonsingular matrix. Then, we compute the quadratic form
#' \eqn{q = x' A^{-1} x}.
#'
#' A naive way to compute the quadratic form is to explicitly write
#' `t(x) \%*\% solve(A) \%*\% x`, but for large `p`, this operation is
#' inefficient. We provide a more efficient method below.
#'
#' Note that we have adapted the code from:
#' \url{https://stat.ethz.ch/pipermail/r-help/2005-November/081940.html}
#'
#' @param A matrix that is p x p and nonsingular
#' @param x vector of length p
#' @return scalar value
quadform_inv <- function(A, x) {
  drop(crossprod(x, solve(A, x)))
}

#' Centers the observations in a matrix by their respective class sample means
#'
#' @inheritParams lda_diag
#' @param y vector of class labels for each training observation
#' @return matrix with observations centered by its corresponding class sample
#' mean
center_data <- function(x, y) {
  x <- pred_to_matrix(x)
  y <- outcome_to_factor(y)
  complete <- complete.cases(x) & complete.cases(y)
  x <- x[complete,,drop = FALSE]
  y <- y[complete]

  # Notice that the resulting centered data are sorted by class and do not
  # preserve the original ordering of the data.
  x_centered <- tapply(seq_along(y), y, function(i) {
    scale(x[i, ], center = TRUE, scale = FALSE)
  })
  x_centered <- do.call(rbind, x_centered)

  # Sorts the centered data to preserve the original ordering.
  orig_ordering <- do.call(c, tapply(seq_along(y), y, identity))
  x_centered[orig_ordering, ] <- x_centered
  x_centered
}

#' Computes the inverse of a symmetric, positive-definite matrix using the
#' Cholesky decomposition
#'
#' This often faster than [solve()] for larger matrices.
#' See, for example:
#' \url{http://blog.phytools.org/2012/12/faster-inversion-of-square-symmetric.html}
#' and
#' \url{https://stats.stackexchange.com/questions/14951/efficient-calculation-of-matrix-inverse-in-r}.
#'
#' @export
#' @param x symmetric, positive-definite matrix
#' @return the inverse of `x`
solve_chol <- function(x) {
  chol2inv(chol(x))
}

#' Computes the log determinant of a matrix.
#'
#' @export
#' @param x matrix
#' @return log determinant of `x`
log_determinant <- function(x) {
  # The call to 'as.vector' removes the attributes returned by 'determinant'
  as.vector(determinant(x, logarithm=TRUE)$modulus)
}

#' Computes multivariate normal density with a diagonal covariance matrix
#'
#' Alternative to `mvtnorm::dmvnorm`
#'
#' @importFrom stats dnorm
#' @param x matrix
#' @param mean vector of means
#' @param sigma vector containing diagonal covariance matrix
#' @return multivariate normal density
dmvnorm_diag <- function(x, mean, sigma) {
  exp(sum(dnorm(x, mean=mean, sd=sqrt(sigma), log=TRUE)))
}

#' Computes posterior probabilities via Bayes Theorem under normality
#'
#' @importFrom mvtnorm dmvnorm
#'
#' @param x matrix of observations
#' @param means list of means for each class
#' @param covs list of covariance matrices for each class
#' @param priors list of prior probabilities for each class
#' @return matrix of posterior probabilities for each observation
posterior_probs <- function(x, means, covs, priors) {
  if (is.vector(x)) {
    x <- matrix(x, nrow=1)
  }
  x <- pred_to_matrix(x)

  posterior <- mapply(function(xbar_k, cov_k, prior_k) {
    if (is.vector(cov_k)) {
      post_k <- apply(x, 1, function(obs) {
        dmvnorm_diag(x=obs, mean=xbar_k, sigma=cov_k)
      })
    } else {
      post_k <- dmvnorm(x=x, mean=xbar_k, sigma=cov_k)
    }
    prior_k * post_k
  }, means, covs, priors)

  if (is.vector(posterior)) {
    posterior <- posterior / sum(posterior)
  } else {
    posterior <- posterior / rowSums(posterior)
  }

  posterior
}



pred_to_matrix <- function(x) {
  if (is.vector(x)) {
    rlang::abort("'x' should be a matrix or data frame.")
  }
  
  if (is.null(colnames(x))) {
    rlang::abort("'x' should have column names.")
  }
  if (!is.matrix(x)) {
   x <- as.matrix(x)
   if (is.character(x)) {
     rlang::abort(
       paste(
         "When the predictors data were converted to a matrix, the matrix was",
         "no longer numeric. Were there non-numeric columns in the original",
         "data?"
       )
     )
   }
  }
  x
}

outcome_to_factor <- function(y) {
  if (is.numeric(y) | is.matrix(y) | is.data.frame(y)) {
    rlang::abort(
      paste(
        "The outcome data should be a character or factor vector."
      )
    )
  }
  if (!is.factor(y)) {
    y <- as.factor(y)
  }
  y
}

format_priors <- function(x) {
  priors <- vapply(x$est, function(x) x$prior, numeric(1))
  paste0(names(priors), " (", round(priors * 100, 2), "%)", collapse = ", ")
}

print_basics <- function(x, ...) {
  cat("Sample Size:", x$N, "\n")
  cat("Number of Features:", x$p, "\n\n")
  cat("Classes and Prior Probabilities:\n  ")
  cat(format_priors(x), "\n")
}

no_form_env <- function(x) {
  attr(x, ".Environment") <- rlang::base_env()
  x
}


new_discrim_object <- function(x, cls) {
  class(x) <- cls
  has_terms <- any(names(x) == ".terms")
  
  if (has_terms) {
    attr(x$.terms, ".Environment") <- rlang::base_env()
    class(x) <- c(paste0(cls, "_formula"), cls)
  }
  x
}

process_newdata <- function(object, x) {
  if (is.null(colnames(x))) {
    rlang::abort("'newdata' should have column names.")
  }
  has_terms <- any(names(object) == ".terms")
  if (has_terms) {
    .terms <- object$.terms
    .terms <- stats::delete.response(.terms)
    x <- stats::model.frame(.terms, x, na.action = stats::na.pass) #, xlev = object$xlevels)
    x <- model.matrix(.terms, x)
    attr(x, "contrasts") <- NULL
    attr(x, "assign") <- NULL
  } 
  x <- x[, object$col_names, drop = FALSE]
  as.matrix(x)
  
}

min_index <- function(x) {
  if (any(is.na(x))) {
    NA_integer_
  } else {
    which.min(x)
  }
}

score_to_class <- function(x, object) {
  if (is.vector(x)) {
    min_scores <- min_index(x)
  } else {
    min_scores <- apply(x, 2, min_index)
  }
  factor(object$groups[min_scores], levels = object$groups)
}

