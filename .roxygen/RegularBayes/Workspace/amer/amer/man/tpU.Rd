\name{tpU}
\alias{tpU}
\title{Generate a modified truncated power basis for penalized spline smoothing.}
\usage{tpU(x, degree=2, k=15, unpen=1, by, allPen=FALSE, varying, diag=FALSE,
    knots=quantile(x, probs = (2:(k - degree + 1))/(k - degree + 2)),
    centerscale, scaledknots=FALSE)}
\description{Generate a modified truncated power basis for penalized spline smoothing.}
\note{This is a more detailed implementation of the example on how to define additional basis generating functions in the vignette.}
\value{list with entries:
\code{"X"}: \code{n x unpen} design matrix for unpenalized part (without intercept)  
\code{"Z"}: \code{n x (k-unpen)} design matrix for penalized part}
\author{Fabian Scheipl}
\seealso{\code{\link{tp}}}
\arguments{\item{x}{covariate for the smooth function}
\item{degree}{integer: degree of truncated polynomials (0: piecewise constant, 1: piecewise linear etc..)}
\item{k}{integer: dimensionality of the basis (i.e.: number of knots + degree)}
\item{unpen}{integer: degree of the unpenalized nullspace, must be lower than degree: 1 for pen. deviations from linearity, 2 for pen. deviations from quadratic from etc.}
\item{by}{factor variable: estimate separate functions for each level - this assumes standard treatment contrasts for the supplied factor.}
\item{allPen}{boolean: if TRUE, make design for group-specific curves with common smoothing parameter: all parameters (including the normally unpenalized basis functions in X) are penalized, every level of "by" has the same amount of smoothing. If FALSE, make design for separate curves for each by-level: separate smoothing parameters for every level of "by", unpenalized estimates for the coefficients associated with X}
\item{varying}{numeric: if not NULL, a varying coefficient model is fit: f(x,varying) = f(x)*varying}
\item{diag}{logical: force a diagonal covariance-matrix for the random effects for X if \code{allPen=TRUE}?}
\item{knots}{vector of knot locations (optional). Defaults to quantile-based knot placement at the \eqn{(i+1)/(k+2)}-quantiles for \eqn{i=1,\dots,k}.}
\item{centerscale}{numeric(2): center&scale x by these values if not NULL}
\item{scaledknots}{boolean:	are knots given for the rescaled x-values?}}
