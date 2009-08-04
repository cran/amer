\name{tp}
\alias{tp}
\title{Generate a truncated power basis for penalized spline smoothing.}
\usage{tp(x, degree=1, k=15, by, allPen=FALSE, varying, diag=FALSE,
    knots=seq(min(x), max(x), l = k - degree + 2)[-c(1, k -
    degree + 2)], centerscale=NULL, scaledknots=FALSE)}
\description{Generate a truncated power basis for penalized spline smoothing.}
\note{Truncated power bases have \code{degree} unpenalized basis functions, namely \eqn{x^1,\dots, x^{degree}} and \eqn{k-}\code{degree} penalized basis functions that contain the positive part \eqn{(x-\kappa_j)^degree} for knots \eqn{\kappa_j, j=1,dots,k-}\code{degree}.
This function can be used as a reference when implementing other \code{basisGenerators} that can be used for \code{\link{amer}}-fits. 
All such functions need to return a list of at least X (unpenalized basis functions, a matrix with zero columns if there are none of those), and Z (penalized basis functions)
that has a \code{call}-attribute with the expanded call returned by \code{\link{expand.call}()}. All such functions need to have at least arguments \code{x, by, allPen, diag} and \code{varying}.
See also section 4.4 in the vignette for an example on how to write your own basis-generating functions.}
\value{list with entries:
\code{"X"}: \code{n x degree} design matrix for unpenalized part (without intercept) (or a list of those for every level of by if allPen=F), 
\code{"Z"}: \code{n x (k-degree)} design matrix for penalized part (or a list of those for every level of by if allPen=F),}
\author{Fabian Scheipl}
\seealso{\code{\link{tpU}}}
\arguments{
\item{x}{covariate for the smooth function}
\item{degree}{integer: degree of truncated polynomials (0: piecewise constant, 1: piecewise linear etc..)}
\item{k}{integer: dimensionality of the basis (i.e.: number of knots + degree)}
\item{by}{factor variable: estimate separate functions for each level - this assumes standard treatment contrasts for the supplied factor.}
\item{allPen}{boolean: if TRUE, make design for group-specific curves with common smoothing parameter: all parameters (including the normally unpenalized basis functions in X) are penalized, every level of "by" has the same amount of smoothing 
if FALSE, make design for separate curves for each by-level: separate smoothing parameters for every level of "by", unpenalized estimates for the coefficients associated with X}
\item{varying}{numeric: if not NULL, a varying coefficient model is fit: f(x,varying) = f(x)*varying}
\item{diag}{logical: force a diagonal covariance-matrix for the random effects for X if \code{allPen=TRUE}?}
\item{knots}{vector of knot locations (optional). Defaults to equidistant knots. For quantile-based knots at the \eqn{(i+1)/(k+2)}-quantiles for \eqn{i=1,\dots,k} use \code{knots=quantile(x, probs = (2:(k - degree + 1))/(k - degree + 2))}}
\item{centerscale}{numeric(2), center&scale x by these values if not NULL}
\item{scaledknots}{boolean,	are knots given for the rescaled x-values?}
}
