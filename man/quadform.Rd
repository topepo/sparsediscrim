% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.r
\name{quadform}
\alias{quadform}
\title{Quadratic form of a matrix and a vector}
\usage{
quadform(A, x)
}
\arguments{
\item{A}{matrix of dimension p x p}

\item{x}{vector of length p}
}
\value{
scalar value
}
\description{
We compute the quadratic form of a vector and a matrix in an efficient
manner. Let \code{x} be a real vector of length \code{p}, and let \code{A} be
a p x p real matrix. Then, we compute the quadratic form \eqn{q = x' A x}.
}
\details{
A naive way to compute the quadratic form is to explicitly write
\verb{t(x) \\\%*\\\% A \\\%*\\\% x}, but for large \code{p}, this operation is
inefficient. We provide a more efficient method below.

Note that we have adapted the code from:
\url{http://tolstoy.newcastle.edu.au/R/help/05/11/14989.html}
}
