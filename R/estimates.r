#' Computes estimates and ancillary information for diagonal classifiers
#'
#' Computes the maximum likelihood estimators (MLEs) for each class under the
#' assumption of multivariate normality for each class. Also, computes ancillary
#' information necessary for classifier summary, such as sample size, the number
#' of features, etc.
#'
#' This function computes the common estimates and ancillary information used in
#' all of the diagonal classifiers in the `sparsediscrim` package.
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
#' @inheritParams lda_diag
#' @param pool logical value. If TRUE, calculates the pooled sample variances
#' for each class.
#' @param est_mean the estimator for the class means. By default, we use the
#' maximum likelihood estimator (MLE). To improve the estimation, we provide the
#' option to use a shrunken mean estimator proposed by Tong et al. (2012).
#' @return named list with estimators for each class and necessary ancillary
#' information
#'
#' @references Tong, T., Chen, L., and Zhao, H. (2012), "Improved Mean
#' Estimation and Its Application to Diagonal Discriminant Analysis,"
#' Bioinformatics, 28, 4, 531-537.
#' \url{http://bioinformatics.oxfordjournals.org/content/28/4/531.long}
diag_estimates <- function(x, y, prior = NULL, pool = FALSE,
                           est_mean = c("mle", "tong")) {
  obj <- list()
  obj$labels <- y
  obj$N <- length(y)
  obj$p <- ncol(x)
  obj$groups <- levels(y)
  obj$num_groups <- nlevels(y)

  est_mean <- match.arg(est_mean)

  # Error Checking
  if (!is.null(prior)) {
    if (length(prior) != obj$num_groups) {
      rlang::abort("The number of 'prior' probabilities must match the number of classes in 'y'.")
    }
    if (any(prior <= 0)) {
      rlang::abort("The 'prior' probabilities must be nonnegative.")
    }
    if (sum(prior) != 1) {
      rlang::abort("The 'prior' probabilities must sum to one.")
    }
  }
  if (any(table(y) < 2)) {
    rlang::abort("There must be at least 2 observations in each class.")
  }

  # By default, we estimate the 'a priori' probabilities of class membership with
  # the MLEs (the sample proportions).
  if (is.null(prior)) {
    prior <- as.vector(table(y) / length(y))
  }

  # For each class, we calculate the MLEs (or specified alternative estimators)
  # for each parameter used in the DLDA classifier. The 'est' list contains the
  # estimators for each class.
  obj$est <- tapply(seq_along(y), y, function(i) {
    stats <- list()
    stats$n <- length(i)
    if (est_mean == "mle") {
      stats$xbar <- colMeans(x[i, , drop = FALSE])
    } else if (est_mean == "tong") {
      stats$xbar <- tong_mean_shrinkage(x[i, , drop = FALSE])
    }
    stats$var <- with(stats, (n - 1) / n * apply(x[i, , drop = FALSE], 2, var))
    stats
  })

  # Calculates the pooled variance across all classes.
  if (pool) {
    obj$var_pool <- Reduce('+', lapply(obj$est, function(x) x$n * x$var)) / obj$N
  }

  # Add each element in 'prior' to the corresponding obj$est$prior
  for(k in seq_len(obj$num_groups)) {
    obj$est[[k]]$prior <- prior[k]
  }
  obj
}

#' Computes estimates and ancillary information for regularized discriminant
#' classifiers
#'
#' Computes the maximum likelihood estimators (MLEs) for each class under the
#' assumption of multivariate normality for each class. Also, computes ancillary
#' information necessary for classifier summary, such as sample size, the number
#' of features, etc.
#'
#' This function computes the common estimates and ancillary information used in
#' all of the regularized discriminant classifiers in the `sparsediscrim`
#' package.
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
#' @inheritParams lda_diag
#' @param y vector of class labels for each training observation
#' @param cov logical. Should the sample covariance matrices be computed?
#' (Default: yes)
#' @param prior vector with prior probabilities for each class. If NULL
#' (default), then the sample proportions are used. See details.
#' @return named list with estimators for each class and necessary ancillary
#' information
regdiscrim_estimates <- function(x, y, cov = TRUE, prior = NULL) {
  obj <- list()
  obj$labels <- y
  obj$N <- length(y)
  obj$p <- ncol(x)
  obj$groups <- levels(y)
  obj$num_groups <- nlevels(y)

  # Error Checking
  if (!is.null(prior)) {
    if (length(prior) != obj$num_groups) {
      rlang::abort("The number of 'prior' probabilities must match the number of
            classes in 'y'.")
    }
    if (any(prior <= 0)) {
      rlang::abort("The 'prior' probabilities must be nonnegative.")
    }
    if (sum(prior) != 1) {
      rlang::abort("The 'prior' probabilities must sum to one.")
    }
  }
  if (any(table(y) < 2)) {
    rlang::abort("There must be at least 2 observations in each class.")
  }

  # By default, we estimate the 'a priori' probabilities of class membership with
  # the MLEs (the sample proportions).
  if (is.null(prior)) {
    prior <- as.vector(table(y) / length(y))
  }

  # For each class, we calculate the MLEs for each parameter used in the
  # 'regdiscrim' classifiers. The 'est' list contains the estimators for each
  # class.
  obj$est <- tapply(seq_along(y), y, function(i) {
    stats <- list()
    stats$n <- length(i)
    stats$xbar <- as.vector(colMeans(x[i, , drop = FALSE]))
    if (cov) {
      stats$cov <- with(stats, cov_mle(x[i, , drop = FALSE]))
    }
    stats
  })

  # Calculates the pooled sample covariance matrix.
  if (cov) {
    obj$cov_pool <- Reduce('+', lapply(obj$est, function(x) x$n * x$cov)) / obj$N
  }

  # Add each element in 'prior' to the corresponding obj$est$prior
  for(k in seq_len(obj$num_groups)) {
    obj$est[[k]]$prior <- prior[k]
  }
  obj
}
