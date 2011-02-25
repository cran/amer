\name{bsp}
\alias{bsp}
\title{Generate a reparameterized b-spline basis...}
\usage{bsp(x, k=15, spline.degree=3, diff.ord=2, knots, by,
    allPen=FALSE, varying, diag=FALSE)}
\description{Generate a reparameterized b-spline basis}
\details{Generate a b-spline basis with equidistant knots in mixed model reparameterization}
\value{list with entries:
\code{"X"}: \code{n x (diff.ord - 1)} design matrix for unpenalized part (without intercept), 
\code{"Z"}: \code{n x (k - diff.ord+1)} design matrix for penalized part,}
\author{Fabian Scheipl}
\seealso{\code{\link{tp}}}
\arguments{\item{x}{covariate for the smooth function}
\item{k}{integer: dimensionality of the basis}
\item{spline.degree}{integer: degree of B-splines (defaults to cubic)}
\item{diff.ord}{integer: order of the difference penalty on the un-reparamerized spline coefficients. Defaults to 2, that is, penalized deviations from linearity.}
\item{knots}{double: vector of (inner & outer) knot locations}
\item{by}{factor variable: estimate separate functions for each level - this assumes standard treatment contrasts for the supplied factor.}
\item{allPen}{boolean: if TRUE, make design for group-specific curves with common smoothing parameter: all parameters (including the normally unpenalized basis functions in X) are penalized, every level of "by" has the same amount of smoothing 
if FALSE, make design for separate curves for each by-level: separate smoothing parameters for every level of "by", unpenalized estimates for the coefficients associated with X}
\item{varying}{numeric: if not NULL, a varying coefficient model is fit: f(x,varying) = f(x)*varying}
\item{diag}{logical: force a diagonal covariance-matrix for the random effects for X if \code{allPen=TRUE}?}}
