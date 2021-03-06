% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lda-emp-bayes-eigen.r
\name{lda_emp_bayes_eigen}
\alias{lda_emp_bayes_eigen}
\alias{lda_emp_bayes_eigen.default}
\alias{lda_emp_bayes_eigen.formula}
\alias{predict.lda_emp_bayes_eigen}
\title{The Minimum Distance Rule using Modified Empirical Bayes (MDMEB) classifier}
\usage{
lda_emp_bayes_eigen(x, ...)

\method{lda_emp_bayes_eigen}{default}(x, y, prior = NULL, eigen_pct = 0.95, ...)

\method{lda_emp_bayes_eigen}{formula}(formula, data, prior = NULL, ...)

\method{predict}{lda_emp_bayes_eigen}(object, newdata, type = c("class", "prob", "score"), ...)
}
\arguments{
\item{x}{Matrix or data frame containing the training data. The rows are the
sample observations, and the columns are the features. Only complete data are
retained.}

\item{...}{additional arguments (not currently used).}

\item{y}{Vector of class labels for each training observation. Only complete
data are retained.}

\item{prior}{Vector with prior probabilities for each class. If NULL
(default), then equal probabilities are used. See details.}

\item{eigen_pct}{the percentage of eigenvalues kept}

\item{formula}{A formula of the form \code{groups ~ x1 + x2 + ...} That is,
the response is the grouping factor and the right hand side specifies the
(non-factor) discriminators.}

\item{data}{data frame from which variables specified in \code{formula} are
preferentially to be taken.}

\item{object}{Fitted model object}

\item{newdata}{Matrix or data frame of observations to predict. Each row
corresponds to a new observation.}

\item{type}{Prediction type: either `"class"`, `"prob"`, or `"score"`.}
}
\value{
\code{lda_emp_bayes_eigen} object that contains the trained MDMEB classifier
}
\description{
Given a set of training data, this function builds the MDMEB classifier from
Srivistava and Kubokawa (2007). The MDMEB classifier is an adaptation of the
linear discriminant analysis (LDA) classifier that is designed for
small-sample, high-dimensional data. Srivastava and Kubokawa (2007) have
proposed a modification of the standard maximum likelihood estimator of the
pooled covariance matrix, where only the largest 95\% of the eigenvalues and
their corresponding eigenvectors are kept. The resulting covariance matrix is
then shrunken towards a scaled identity matrix. The value of 95\% is the
default and can be changed via the \code{eigen_pct} argument.

The MDMEB classifier is an adaptation of the linear discriminant analysis
(LDA) classifier that is designed for small-sample, high-dimensional
data. Srivastava and Kubokawa (2007) have proposed a modification of the
standard maximum likelihood estimator of the pooled covariance matrix, where
only the largest 95\% of the eigenvalues and their corresponding eigenvectors
are kept. The resulting covariance matrix is then shrunken towards a scaled
identity matrix.
}
\details{
The matrix of training observations are given in \code{x}. The rows of \code{x}
contain the sample observations, and the columns contain the features for each
training observation.

The vector of class labels given in \code{y} are coerced to a \code{factor}.
The length of \code{y} should match the number of rows in \code{x}.

An error is thrown if a given class has less than 2 observations because the
variance for each feature within a class cannot be estimated with less than 2
observations.

The vector, \code{prior}, contains the \emph{a priori} class membership for
each class. If \code{prior} is NULL (default), the class membership
probabilities are estimated as the sample proportion of observations belonging
to each class. Otherwise, \code{prior} should be a vector with the same length
as the number of classes in \code{y}. The \code{prior} probabilities should be
nonnegative and sum to one.
}
\examples{
library(modeldata)
data(penguins)
pred_rows <- seq(1, 344, by = 20)
penguins <- penguins[, c("species", "body_mass_g", "flipper_length_mm")]
mdmeb_out <- lda_emp_bayes_eigen(species ~ ., data = penguins[-pred_rows, ])
predicted <- predict(mdmeb_out, penguins[pred_rows, -1], type = "class")

mdmeb_out2 <- lda_emp_bayes_eigen(x = penguins[-pred_rows, -1], y = penguins$species[-pred_rows])
predicted2 <- predict(mdmeb_out2, penguins[pred_rows, -1], type = "class")
all.equal(predicted, predicted2)
}
\references{
Srivastava, M. and Kubokawa, T. (2007). "Comparison of
Discrimination Methods for High Dimensional Data," Journal of the Japanese
Statistical Association, 37, 1, 123-134.
}
