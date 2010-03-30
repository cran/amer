\name{reparameterizeDesign}
\alias{reparameterizeDesign}
\title{Reparameterize a spline basis...}
\usage{reparameterizeDesign(B, K, tol=1e-10)}
\description{Reparameterize a spline basis}
\details{Reparameterize a spline basis via the spectral decomposition of its penalty matrix into a design for the unpenalized part of the function and 
a design for the penalized part with i.i.d. spline coefficients. A centering constraint is applied to the design of the unpenalized part,
which only makes sense if the intercept column is actually in the column space of X, i.e. constant functions are in the nullspace of the penalty.}
\value{list with entries:
\code{"X"}: \code{n x (p-1)} design matrix for unpenalized part (without intercept) 
\code{"Z"}: \code{n x (k-p-1)} design matrix for penalized part
p is the dimension+1 of the penalty nullspace (e.g. 2 for a linear TP-basis or a B-spline with 2nd order difference penalty)}
\author{Fabian Scheipl}
\seealso{source code of \code{tpU2D} for a usage example}
\arguments{\item{B}{matrix containing the basis functions}
\item{K}{penalty matrix (i.e. the penalty term for the fitted spline is the quadratic form of K and the spline coefficients)}
\item{tol}{eigenvalues smaller than tol are considered zero}}
