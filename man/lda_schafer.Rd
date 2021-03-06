% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lda-schafer.r
\name{lda_schafer}
\alias{lda_schafer}
\alias{lda_schafer.default}
\alias{lda_schafer.formula}
\alias{predict.lda_schafer}
\title{Linear Discriminant Analysis using the Schafer-Strimmer Covariance Matrix
Estimator}
\usage{
lda_schafer(x, ...)

\method{lda_schafer}{default}(x, y, prior = NULL, ...)

\method{lda_schafer}{formula}(formula, data, prior = NULL, ...)

\method{predict}{lda_schafer}(object, newdata, type = c("class", "prob", "score"), ...)
}
\arguments{
\item{x}{Matrix or data frame containing the training data. The rows are the
sample observations, and the columns are the features. Only complete data are
retained.}

\item{...}{Options passed to \code{\link[corpcor:invcov.shrink]{corpcor::invcov.shrink()}}}

\item{y}{Vector of class labels for each training observation. Only complete
data are retained.}

\item{prior}{Vector with prior probabilities for each class. If NULL
(default), then equal probabilities are used. See details.}

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
\code{lda_schafer} object that contains the trained classifier
}
\description{
Given a set of training data, this function builds the Linear Discriminant
Analysis (LDA) classifier, where the distributions of each class are assumed
to be multivariate normal and share a common covariance matrix. When the
pooled sample covariance matrix is singular, the linear discriminant function
is incalculable. This function replaces the inverse of pooled sample
covariance matrix with an estimator proposed by Schafer and Strimmer
(2005). The estimator is calculated via \code{\link[corpcor:invcov.shrink]{corpcor::invcov.shrink()}}.

The Linear Discriminant Analysis (LDA) classifier involves the assumption
that the distributions of each class are assumed to be multivariate normal
and share a common covariance matrix. When the pooled sample covariance
matrix is singular, the linear discriminant function is incalculable. Here,
the inverse of the pooled sample covariance matrix is replaced with an
estimator from Schafer and Strimmer (2005).
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
lda_schafer_out <- lda_schafer(species ~ ., data = penguins[-pred_rows, ])
predicted <- predict(lda_schafer_out, penguins[pred_rows, -1], type = "class")

lda_schafer_out2 <- lda_schafer(x = penguins[-pred_rows, -1], y = penguins$species[-pred_rows])
predicted2 <- predict(lda_schafer_out2, penguins[pred_rows, -1], type = "class")
all.equal(predicted, predicted2)
}
\references{
Schafer, J., and Strimmer, K. (2005). "A shrinkage approach to
large-scale covariance estimation and implications for functional genomics,"
Statist. Appl. Genet. Mol. Biol. 4, 32.
}
