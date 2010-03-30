\name{tpU2d}
\alias{tpU2d}
\title{Generate a 2-D truncated power basis...}
\usage{tpU2d(x, z, kx=7, kz=7, degree=2, dimU=1, by, allPen=FALSE,
    varying, diag=FALSE, knotsx=quantile(x, probs = (2:(kx -
    dimU + 1))/(kx - dimU + 2)), knotsz=quantile(z, probs =
    (2:(kz - dimU + 1))/(kz - dimU + 2)), centerscale,
    scaledknots=FALSE)}
\description{Generate a 2-D truncated power basis}
\details{Generates a two-dimensional tensor product spline basis of marginal truncated power bases with a single smoothing parameter. 
The penalty matrix is the sum of the two kronecker products of the penalties for the marginal bases with the identity matrix.}
\value{list with entries:
\code{"X"}: \code{n x (2*degree + degree^2) )} design matrix for unpenalized part (without intercept): marginal polynomials plus their interactions, 
\code{"Z"}: \code{n x ((kx+1)*(kz+1) - (1 + 2*degree + degree^2))} design matrix for penalized part,}
\author{Fabian Scheipl}
\seealso{\code{\link{tp}}, \code{\link{tpU}}}
\arguments{\item{x}{covariate for the smooth function}
\item{z}{covariate for the smooth function}
\item{kx,kz}{integer: dimensionality of the marginal bases}
\item{dimU}{integer: dimension of marginal nullspaces: 1 means unpenalized terms up to linear-linear-interaction, 
2 means unpenalized terms up to quadratic-quadratic-interactions etc.}
\item{degree}{integer: degree of (truncated) polynomial}
\item{by}{factor variable: estimate separate functions for each level - this assumes standard treatment contrasts for the supplied factor.}
\item{allPen}{boolean: if TRUE, make design for group-specific curves with common smoothing parameter: all parameters (including the normally unpenalized basis functions in X) are penalized, every level of "by" has the same amount of smoothing 
if FALSE, make design for separate curves for each by-level: separate smoothing parameters for every level of "by", unpenalized estimates for the coefficients associated with X}
\item{varying}{numeric: if not NULL, a varying coefficient model is fit: f(x,varying) = f(x)*varying}
\item{diag}{logical: force a diagonal covariance-matrix for the random effects for X if \code{allPen=TRUE}?}
\item{knotsx,knotsz}{double: vector of knot locations (defaults to marginal-quantile-based knot locations)}
\item{centerscale}{numeric(2): center&scale x by these values if not NULL}
\item{scaledknots}{boolean:	are knot locations given for the rescaled x-values?}}
